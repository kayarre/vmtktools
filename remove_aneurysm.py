from vmtkClassifyBifucation import *
from common import *
from argparse import ArgumentParser
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from os import path, listdir
from subprocess import STDOUT, check_output
import sys
import math
from time import time


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--d', '--dir_path', type=str, default=".", 
                        help="Path to the folder with all the cases")
    parser.add_argument('--case', type=str, default=None, help="Choose case")
    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram (surface) should be" + \
                        "smoothed before it is manipulated", metavar="smooth") 
    parser.add_argument('--a', '--angle', type=float, default=10,
                        help="Each daughter branch is rotated an angle a in the" + \
                        " bifurcation plane. a should be expressed in radians as" + \
                        " any math expression from the" + \
                        " math module in python", metavar="rotation_angle")
    parser.add_argument('--smooth_factor', type=float, default=0.25,
                         help="If smooth option is true then each voronoi point" + \
                         " that has a radius less then MISR*(1-smooth_factor) at" + \
                         " the closest centerline point is removes",
                         metavar="smoothening_factor")
    parser.add_argument("--leave1", type=bool, default=False,
                        help="Leave one branch untuched")
    parser.add_argument("--leave2", type=bool, default=False,
                        help="Leave one branch untuched")
    parser.add_argument("--bif", type=bool, default=False, 
                        help="interpolate bif as well")
    parser.add_argument("--addPoint", type=bool, default=False,
                        help="Add additional point for integration")
    parser.add_argument("--lower", type=bool, default=False,
                        help="Make a fourth line to interpolate along that" + \
                             " is lower than the other bif line.")
    parser.add_argument("--cylinder_factor", type=float, default=7.0,
                        help="Factor for choosing the smaller cylinder")

    args = parser.parse_args()
    ang_ = 0 if args.a == 0 else math.pi/args.a

    return args.s, ang_, args.smooth_factor, args.leave1, args.leave2, \
           args.bif, args.addPoint, args.d, args.case, args.lower, \
           args.cylinder_factor

def getData(centerline_par, centerline_dau1, centerline_dau2, tol):
    # Declear variables before loop incase values are not found
    data = {"par":{}, "dau1":{}, "dau2":{}}

    # List of points conected to ID
    points_ids_0 = vtk.vtkIdList()
    points_ids_1 = vtk.vtkIdList()

    data_list = [("dau1", centerline_dau1), ("dau2", centerline_dau2), ("par", centerline_par)]
    for key, centerline in data_list:
        # One is the branch to the left and the other is the one to the right
        centerline.GetCellPoints(0, points_ids_0)
        centerline.GetCellPoints(1, points_ids_1)

        # Find clipping points
        N = min(points_ids_0.GetNumberOfIds(), points_ids_1.GetNumberOfIds())
        for i in range(0, N):
            cell_point_0 = centerline.GetPoint(points_ids_0.GetId(i))
            cell_point_1 = centerline.GetPoint(points_ids_1.GetId(i))
            
            distance_between_points = math.sqrt(distance(cell_point_0, cell_point_1))
            if distance_between_points > tol:
                tmpI = i
                point_ID_0 = points_ids_0.GetId(i)
                point_ID_1 = points_ids_1.GetId(i)
                center = centerline.GetPoint(point_ID_0)
                r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID_0)
                break

        end, r_end = move_past_sphere(centerline, center, r, point_ID_0,
                                        stop=point_ID_0*100, step=1, X=1)

        data[key]["end_point"] = end
        data[key]["r_end"] = r_end
        data[key]["div_point"] = center
        data[key]["ID_div"] = point_ID_0
        data[key]["i_div"] = tmpI
        data[key]["r_div"] = r

    return data


def get_points(data, key, bif=False):
    div_points = np.asarray([data["par"][key], data["dau1"][key], data["dau2"][key]])

    points = vtk.vtkPoints()
    for point in div_points[bif:]:
        points.InsertNextPoint(point)
    
    return points, div_points[bif:]


def get_startpoint(centerline):
    line = ExtractSingleLine(centerline, 0)
    return line.GetPoints().GetPoint(0)


def merge_cl(centerline, end_point, div_point):
    merge = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    N_lines = centerline.GetNumberOfLines()

    arrays = []
    N_, names = get_number_of_arrays(centerline)
    for i in range(N_):
        tmp = centerline.GetPointData().GetArray(names[i])
        tmp_comp = tmp.GetNumberOfComponents()
        array = get_vtk_array(names[i], tmp_comp, centerline.GetNumberOfPoints())
        arrays.append(array)

    # Find lines to merge
    lines = [ExtractSingleLine(centerline, i) for i in range(N_lines)]
    locators = [get_locator(lines[i]) for i in range(N_lines)]
    div_ID = [locators[i].FindClosestPoint(div_point[0]) for i in range(N_lines)]
    end_ID = [locators[i].FindClosestPoint(end_point[0]) for i in range(N_lines)]
    dist = [np.sum(lines[i].GetPoint(end_ID[i]) - end_point[0]) for i in range(N_lines)]
    change = [j for j in range(N_lines) if dist[j] != 0 or dist[j] > 5]

    # Find the direction of each line
    map_other = {0: 1, 1: 0}
    ID0 = locators[0].FindClosestPoint(end_point[1])
    ID1 = locators[1].FindClosestPoint(end_point[1])
    dist0 = math.sqrt(np.sum((np.asarray(lines[0].GetPoint(ID0)) - end_point[1])**2))
    dist1 = math.sqrt(np.sum((np.asarray(lines[1].GetPoint(ID1)) - end_point[1])**2))
    end1 = 0 if dist0 < dist1 else 1
    end2 = int(not end1)
    for i in range(2, N_lines):
        ID1 = locators[i].FindClosestPoint(end_point[1])
        ID2 = locators[i].FindClosestPoint(end_point[2])
        dist1 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID1)) - end_point[1])**2))
        dist2 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID2)) - end_point[2])**2))
        map_other[i] = end1 if dist1 > dist2 else end2

    counter = 0
    for i in range(centerline.GetNumberOfLines()):
        line = lines[i]
        other = lines[map_other[i]]
        N = line.GetNumberOfPoints()
        cellArray.InsertNextCell(N)
        
        for j in range(N):
            # Add point
            if div_ID[i] < j < end_ID[i] and i in change:
                new = (np.asarray(other.GetPoint(j)) +
                    np.asarray(line.GetPoint(j))) / 2.
                points.InsertNextPoint(new)
            else:
                points.InsertNextPoint(line.GetPoint(j))

            cellArray.InsertCellPoint(counter)

            # Add array
            for k in range(N_):
                num = arrays[k].GetNumberOfComponents()
                if num == 1:
                    tmp = line.GetPointData().GetArray(names[k]).GetTuple1(j)
                    arrays[k].SetTuple1(counter, tmp)
                elif num == 3:
                    tmp = line.GetPointData().GetArray(names[k]).GetTuple3(j)
                    arrays[k].SetTuple3(counter, tmp[0], tmp[1], tmp[2])
                else:
                    print "Add more options"
                    sys.exit(0)

            counter += 1

    # Insert points, lines and arrays
    merge.SetPoints(points)
    merge.SetLines(cellArray)
    for i in range(N_):
        merge.GetPointData().AddArray(arrays[i])
    
    return merge


def main(dirpath, smooth, smooth_factor, angle, l1, l2, bif, addPoint, lower,
         cylinder_factor):
    # Input filenames
    model_path = path.join(dirpath, "surface", "model.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, "surface", "model_smooth.vtp")
    
    centerline_par_path = path.join(dirpath, "surface", "centerline_par.vtp")
    centerline_dau1_path = path.join(dirpath, "surface", "centerline_dau1.vtp")
    centerline_dau2_path = path.join(dirpath, "surface", "centerline_dau2.vtp")
    centerline_bif_path = path.join(dirpath, "surface", "centerline_anu_bif.vtp")
    centerline_clipped_path = path.join(dirpath, "surface", "centerline_clipped_anu.vtp")
    centerline_clipped_bif_path = path.join(dirpath, "surface", "centerline_clipped_bif_anu.vtp")
    centerline_new_path = path.join(dirpath, "surface", "centerline_interpolated_anu.vtp")
    centerline_new_bif_path = path.join(dirpath, "surface", "centerline_interpolated_bif_anu.vtp")
    centerline_new_bif_lower_path = path.join(dirpath, "surface", \
                                        "centerline_interpolated_bif_lower_anu.vtp")

    voronoi_path = path.join(dirpath, "surface", "voronoi.vtp")
    voronoi_smoothed_path = path.join(dirpath, "surface", "voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, "surface", "voronoi_clipped_anu.vtp")
    voronoi_anu_path = path.join(dirpath, "surface", "voronoi_anu.vtp")

    s = "_pi%s" % angle if angle == 0 else "_pi%s" % (1/(angle/math.pi))
    s += "" if not l1 else "_l1"
    s += "" if not l2 else "_l2"
    s += "" if not bif else "_bif"
    s += "" if not smooth else "_smooth"
    s += "" if not addPoint else "_extraPoint"
    s += "" if not lower else "_lower"
    s += "" if cylinder_factor == 7.0 else "_cyl%s" % cylinder_factor
    model_new_surface = path.join(dirpath, "surface", "model_anu"+s+".vtp")

    # Model
    if not path.exists(model_path):
        print "The given directory: %s did not contain surface/model.vtp" % dirpath
        sys.exit(0)

    # Centerline
    parameters = getParameters(dirpath)
    anu = max([int(k.split("outlet")[-1]) for k in parameters.keys() if "outlet" in k])
    rest = range(anu)
    centerline_par = makeCenterline(model_path, centerline_par_path,
                   #                 recompute=True,
                                    in_out=[-1] + rest, length=0.2, smooth=False)
    centerline_dau1 = makeCenterline(model_path, centerline_dau1_path,
                   #                 recompute=True,
                                    in_out=[0, anu, 1], length=0.2, smooth=False)
    centerline_dau2 = makeCenterline(model_path, centerline_dau2_path,
                   #                 recompute=True,
                                    in_out=[1, anu, 0], length=0.2, smooth=False)
    centerline_bif = makeCenterline(model_path, centerline_bif_path,
                   #                 recompute=True,
                                    in_out=[0, 1], length=0.2, smooth=False)


    # Voronoi
    print "Compute voronoi diagram"
    voronoi = makeVoronoi(model_path, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerline_par, 0.25)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)

        surface_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(surface_smoothed, model_smoothed_path)

    voronoi = voronoi if not smooth else ReadPolyData(voronoi_smoothed_path)

    # Create a tolerance for diverging
    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                    centerline_par.GetPoint(10), \
                                    centerline_par.GetPoint(11)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance * 3

    # Get data from centerlines
    data = getData(centerline_par, centerline_dau1, centerline_dau2, divergingTolerance)

    def write_spheres(points, radius=None, name="sphere%s.vtp", base=0.0005):
        radius = [base]*len(points) if radius is None else radius
        for counter, point in enumerate(points):
            sphere = vtk.vtkSphereSource()
            sphere.SetCenter(point)
            sphere.SetPhiResolution(500)
            sphere.SetThetaResolution(500)
            sphere.SetRadius(radius[counter])
            sphere_ = sphere.GetOutput()
    
            WritePolyData(sphere_, path.join(dirpath, "surface", name % counter))

    points = [data[k]["end_point"] for k in data.keys()]
    write_spheres(points, name="sphere_anu_end%s.vtp", base=0.15)

    points = [data[k]["div_point"] for k in data.keys()]
    radius = [data[k]["r_div"] for k in data.keys()]
    write_spheres(points, radius=radius, name="sphere_anu%s.vtp")
    write_spheres(points, name="sphere_anu_small%s.vtp", base=0.15)

    # divpoints or endpoints, rotated or not, and for bif or not
    key = "div_point"
    div_points = get_points(data, key)
    div_points_bif = get_points(data, key, bif=True)

    center = [((1/9.)*div_points[1][0] + (4/9.)*div_points[1][1] + \
                            (4/9.)*div_points[1][2]).tolist()]
    #print center
    #center = [[(points[0][i] + points[1][i] + points[2][i]) / 3 for i in range(3)]]
    write_spheres(center, name="sphere_center%s.vtp", base=0.15)
    #sys.exit(0)

    key = "end_point"
    end_points = get_points(data, key)
    end_points_bif = get_points(data, key, bif=True)

    print "Clipping centerlines and voronoi diagram."
    patch_cl = CreateParentArteryPatches(centerline_par, end_points[0])
    WritePolyData(patch_cl, centerline_clipped_path)

    patch_bif_cl = CreateParentArteryPatches(centerline_bif, end_points_bif[0])
    WritePolyData(patch_bif_cl, centerline_clipped_bif_path)

    # Clipp the voronoi diagram
    masked_voronoi = MaskVoronoiDiagram(voronoi, patch_cl)
    voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
    WritePolyData(voronoi_clipped, voronoi_clipped_path)

    # Interpolate the centerline
    print "Interpolate centerlines and voronoi diagram."
    version = True #bif * lower
    interpolated_cl = InterpolatePatchCenterlines(patch_cl, centerline_par,
                                                  div_points[0].GetPoint(0),
                                                  None, version)
    WritePolyData(interpolated_cl, centerline_new_path.replace(".vtp", "1.vtp"))

    if bif:
        interpolated_bif = InterpolatePatchCenterlines(patch_bif_cl, centerline_bif,
                                                        None, "bif", version)
        WritePolyData(interpolated_bif, centerline_new_bif_path)

    if lower:
        #center = np.sum(div_points[1], axis=0)/3.
        # Skewed center
        center = (3/9.)*div_points[1][0] + (3/9.)*div_points[1][1] + \
                 (3/9.)*div_points[1][2]
        div_points_bif[0].SetPoint(0, center[0], center[1], center[2])
        interpolated_bif_lower = InterpolatePatchCenterlines(patch_bif_cl, centerline_bif,
                                                             div_points_bif[0].GetPoint(0),
                                                             "lower", version)
        WritePolyData(interpolated_bif_lower, centerline_new_bif_lower_path)

    interpolated_cl = merge_cl(interpolated_cl, div_points[1],
                               end_points[1])
    WritePolyData(interpolated_cl, centerline_new_path)
    
    # Interpolate voronoi diagram
    if lower and bif:
        bif = [interpolated_bif, interpolated_bif_lower, patch_bif_cl]
    elif bif:
        bif = [interpolated_bif, patch_bif_cl]
    elif lower:
        bif = [interpolated_bif_lower, patch_bif_cl]
    else:
        bif = []
    
    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, patch_cl, 
                                                       voronoi_clipped,
                                                       end_points[0],
                                                       bif, lower, cylinder_factor)
    interpolated_voronoi = remove_distant_points(interpolated_voronoi, interpolated_cl)

    WritePolyData(interpolated_voronoi, voronoi_anu_path) 

    # Write a new surface from the new voronoi diagram
    print "Write new surface"
    new_surface = create_new_surface(interpolated_voronoi)
    WritePolyData(new_surface, model_new_surface)


if  __name__ == "__main__":
    smooth, angle, smooth_factor, l1, l2, bif, addPoint, basedir, case, lower, \
    cylinder_factor = read_command_line()
    folders = listdir(basedir) if case is None else [case]
    for folder in folders:
        if folder[:2] in ["P0", "C0"]:
            main(path.join(basedir, folder), smooth, smooth_factor, angle, l1,
                  l2, bif, addPoint, lower, cylinder_factor)
