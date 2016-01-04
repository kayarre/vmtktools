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
    parser.add_argument("--w", type=bool, default=False, help="Use the weighted" + \
                        "center")
    parser.add_argument("--version", type=bool, default=True, help="Type of" + \
                        "interpolation")
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
           args.cylinder_factor, args.w, args.version


def rotate_voronoi(clipped_voronoi, patch_cl, div_points, m, R):
    numberOfPoints = clipped_voronoi.GetNumberOfPoints()
    distance = vtk.vtkMath.Distance2BetweenPoints
    I = np.eye(3)
    R_inv = np.linalg.inv(R)

    locator = []
    cellLine = []
    not_rotate = [0]
    for i in range(patch_cl.GetNumberOfCells()):
        cellLine.append(ExtractSingleLine(patch_cl, i))
        tmp_locator = get_locator(cellLine[-1])
        locator.append(tmp_locator)

    for i in range(1, patch_cl.GetNumberOfCells()):
        pnt = cellLine[i].GetPoints().GetPoint(0)
        new = cellLine[0].GetPoints().GetPoint(locator[0].FindClosestPoint(pnt))
        dist = math.sqrt(distance(pnt, new)) < divergingRatioToSpacingTolerance
        if dist:
            not_rotate.append(i)

    def check_rotate(point):
        dist = []
        for i in range(len(locator)):
            tmp = locator[i].FindClosestPoint(point)
            tmp = cellLine[i].GetPoints().GetPoint(tmp)
            dist.append(math.sqrt(distance(tmp, point)))

        if dist.index(min(dist)) not in not_rotate:
            pnt = cellLine[dist.index(min(dist))].GetPoints().GetPoint(0)
            if math.sqrt(distance(pnt, div_points[1])) >  \
                    math.sqrt(distance(pnt, div_points[2])):
                m_ = m[2]
                div = div_points[2]
            else:
                m_ = m[1]
                div = div_points[1]
            return m_, div
        else:
            return I, np.array([0, 0, 0])

    maskedVoronoi = vtk.vtkPolyData()
    maskedPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPoints) 

    count = 0
    for i in range(numberOfPoints):
        point = [0.0, 0.0, 0.0]
        clipped_voronoi.GetPoint(i, point)

        pointRadius = clipped_voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        M, O = check_rotate(point)
        tmp = np.dot(np.dot(np.dot(np.asarray(point) - O, R), M), R_inv) + O
        maskedPoints.InsertNextPoint(tmp)
        radiusArray.SetTuple1(i, pointRadius)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

    maskedVoronoi.SetPoints(maskedPoints)
    maskedVoronoi.SetVerts(cellArray)
    maskedVoronoi.GetPointData().AddArray(radiusArray)

    return maskedVoronoi


def rotate_cl(patch_cl, div_points, rotation_matrix, R):
    distance = vtk.vtkMath.Distance2BetweenPoints
    I = np.eye(3)
    R_inv = np.linalg.inv(R)
    
    numberOfPoints = patch_cl.GetNumberOfPoints()

    centerline = vtk.vtkPolyData()
    centerlinePoints = vtk.vtkPoints()
    centerlineCellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPoints)

    line0 = ExtractSingleLine(patch_cl, 0)
    locator0 = get_locator(line0)

    count = 0
    for i in range(patch_cl.GetNumberOfCells()):
        cell = ExtractSingleLine(patch_cl, i)
        centerlineCellArray.InsertNextCell(cell.GetNumberOfPoints())

        start = cell.GetPoint(0)
        dist = line0.GetPoint(locator0.FindClosestPoint(start))
        test = math.sqrt(distance(start, dist)) > divergingRatioToSpacingTolerance

        if test or len(div_points) == 2:
            locator = get_locator(cell)

            pnt1 = cell.GetPoint(locator.FindClosestPoint(div_points[-2]))
            pnt2 = cell.GetPoint(locator.FindClosestPoint(div_points[-1]))
            dist1 = math.sqrt(distance(pnt1, div_points[-2]))
            dist2 = math.sqrt(distance(pnt2, div_points[-1]))
            k = -2 if dist1 < dist2 else -1
            O = div_points[k]
            m = rotation_matrix[k+3]

        else:
            m = I
            O = np.array([0, 0, 0])
        
        getData = cell.GetPointData().GetArray(radiusArrayName).GetTuple1
        for j in range(cell.GetNumberOfPoints()):
            point = np.asarray(cell.GetPoints().GetPoint(j))
            tmp = np.dot(np.dot(np.dot(point - O, R), m), R_inv) + O
            centerlinePoints.InsertNextPoint(tmp)
            radiusArray.SetTuple1(count, getData(j))
            centerlineCellArray.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerlinePoints)
    centerline.SetLines(centerlineCellArray)
    centerline.GetPointData().AddArray(radiusArray)

    return centerline


def rotationMatrix(data, angle, leave1, leave2):
    d = (np.asarray(data[0]["div_point"]) + \
        np.asarray(data[1]["div_point"]) + \
        np.asarray(data["bif"]["div_point"])) / 3.
    vec = np.eye(3)
    for i in range(2):
        e = np.asarray(data[i]["end_point"])
        tmp = e - d
        len = math.sqrt(np.dot(tmp, tmp))
        vec[:,i] = tmp / len

    R = GramSchmidt(vec)

    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    m1 = np.asarray([[cos_a, -sin_a, 0],
                     [sin_a,  cos_a, 0],
                     [    0,      0, 1]])
    m2 = np.asarray([[ cos_a, sin_a, 0],
                     [-sin_a, cos_a, 0],
                     [     0,     0, 1]])

    m = {1: m1, 2: m2}
    tmp1 = data[0]["div_point"] - d
    tmp2 = data[1]["div_point"] - d

    I = np.eye(3)

    if tmp1[0] < tmp2[0]:
        m = {1: m2, 2: m1}

    if leave1:
        k = 1 if data[0]["r_end"] > data[1]["r_end"] else 2
        m[k] = I
    if leave2:
        k = 1 if data[0]["r_end"] < data[1]["r_end"] else 2
        m[k] = I

    return R, m


def get_points(data, key, R, m, rotated=True, bif=False):
    div_points = np.asarray([data["bif"][key], data[0][key], data[1][key]])

    # Origo of the bifurcation
    O_key = "div_point"
    O = np.asarray([data["bif"][O_key], data[0][O_key], data[1][O_key]])
    O = np.sum(np.asarray(O),axis=0)/3.

    if rotated:
        R_inv = np.linalg.inv(R)
        for i in range(len(div_points)):
            m_ = m[i] if i > 0 else np.eye(3)
            div_points[i] = np.dot(np.dot(np.dot(div_points[i] - O, R), m_), R_inv) + O

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
    change = [j for j in range(N_lines) if dist[j] == 0]

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
         cylinder_factor, w, version):
    # Input filenames
    model_path = path.join(dirpath, "surface", "model.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, "surface", "model_smooth.vtp")
    centerline_path = path.join(dirpath, "surface", "centerline_for_rotating.vtp")
    centerline_bif_path = path.join(dirpath, "surface", "centerline_bifurcation.vtp")
    centerline_new_path = path.join(dirpath, "surface", "centerline_interpolated.vtp")
    centerline_new_bif_path = path.join(dirpath, "surface", "centerline_interpolated_bif.vtp")
    centerline_new_bif_lower_path = path.join(dirpath, "surface", \
                                      "centerline_interpolated_bif_lower.vtp")
    centerline_clipped_path = path.join(dirpath, "surface", "centerline_clipped.vtp")
    centerline_rotated_path = path.join(dirpath, "surface", "centerline_rotated.vtp")
    centerline_rotated_bif_path = path.join(dirpath, "surface", "centerline_rotated_bif.vtp")
    centerline_bif_clipped_path = path.join(dirpath, "surface", "centerline_clipped_bif.vtp")
    voronoi_path = path.join(dirpath, "surface", "voronoi.vtp")
    voronoi_smoothed_path = path.join(dirpath, "surface", "voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, "surface", "voronoi_clipped.vtp")
    voronoi_rotated_path = path.join(dirpath, "surface", "voronoi_rotated.vtp")
    voronoi_angle_path = path.join(dirpath, "surface", "voronoi_angle.vtp")

    s = "_pi%s" % angle if angle == 0 else "_pi%s" % (1/(angle/math.pi))
    s += "" if not l1 else "_l1"
    s += "" if not l2 else "_l2"
    s += "" if not bif else "_bif"
    s += "" if not smooth else "_smooth"
    s += "" if not addPoint else "_extraPoint"
    s += "" if not lower else "_lower"
    s += "" if cylinder_factor == 7.0 else "_cyl%s" % cylinder_factor
    model_new_surface = path.join(dirpath, "surface", "model_angle"+s+".vtp")

    # Model
    if not path.exists(model_path):
        print "The given directory: %s did not contain surface/model.vtp" % dirpath
        sys.exit(0)

    # Centerline
    centerline = makeCenterline(model_path, centerline_path, length=0.1, smooth=False)
    centerline_bif = makeCenterline(model_path, centerline_bif_path,
                                    in_out=[0,1], length=0.1, smooth=False)

    # Voronoi
    print "Compute voronoi diagram"
    voronoi = makeVoronoi(model_path, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerline, 0.25)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)

        surface_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(surface_smoothed, model_smoothed_path)

    voronoi = voronoi if not smooth else ReadPolyData(voronoi_smoothed_path)

    # Create a tolerance for diverging
    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                    centerline.GetPoint(10), \
                                    centerline.GetPoint(11)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance

    # Get data from centerlines
    data = getData(centerline, centerline_bif, divergingTolerance)
    R, m = rotationMatrix(data, angle, l1, l2) 

    # divpoints or endpoints, rotated or not, and for bif or not
    key = "div_point"
    div_points_rotated_bif = get_points(data, key, R, m, rotated=True, bif=True)
    div_points_rotated = get_points(data, key, R, m, rotated=True, bif=False)
    div_points = get_points(data, key, R, m, rotated=False, bif=False)

    key = "end_point"
    end_points_bif = get_points(data, key, R, m, rotated=False, bif=True)
    end_points = get_points(data, key, R, m, rotated=False, bif=False)
    end_points_rotated = get_points(data, key, R, m, rotated=True, bif=False)
    end_points_rotated_bif = get_points(data, key, R, m, rotated=True, bif=True)

    print "Clipping centerlines and voronoi diagram."
    patch_cl = CreateParentArteryPatches(centerline, end_points[0])
    patch_bif = CreateParentArteryPatches(centerline_bif, end_points_bif[0])
    WritePolyData(patch_cl, centerline_clipped_path)
    WritePolyData(patch_bif, centerline_bif_clipped_path)
    
    # Clipp the voronoi diagram
    masked_voronoi = MaskVoronoiDiagram(voronoi, patch_cl)
    voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
    WritePolyData(voronoi_clipped, voronoi_clipped_path)
    
    # Rotate branches (Both centerline and voronoi)
    print "Rotate centerlines and voronoi diagram."
    rotated_cl = rotate_cl(patch_cl, end_points[1], m, R)
    rotated_bif = rotate_cl(patch_bif, end_points_bif[1], m, R)
    WritePolyData(rotated_cl, centerline_rotated_path)
    WritePolyData(rotated_bif, centerline_rotated_bif_path)

    rotated_voronoi = rotate_voronoi(voronoi_clipped, patch_cl, end_points[1], m, R)
    WritePolyData(rotated_voronoi, voronoi_rotated_path)
    
    # Interpolate the centerline
    print "Interpolate centerlines and voronoi diagram."
    interpolated_cl = InterpolatePatchCenterlines(rotated_cl, centerline,
                                                  div_points_rotated[0].GetPoint(0),
                                                  None, version)

    if bif:
        interpolated_bif = InterpolatePatchCenterlines(rotated_bif, centerline_bif,
                                                        None, "bif", True)
        WritePolyData(interpolated_bif, centerline_new_bif_path)

    if lower:
        if not w:
            center = np.sum(div_points[1], axis=0)/3.
        else:
            # Skewed center
            center = (1/9.)*div_points[1][0] + (4/9.)*div_points[1][1] + \
                     (4/9.)*div_points[1][2]
        div_points_rotated_bif[0].SetPoint(0, center[0], center[1], center[2])
        interpolated_bif_lower = InterpolatePatchCenterlines(rotated_bif, centerline_bif,
                                                             div_points_rotated_bif[0].GetPoint(0),
                                                             "lower", True)
        WritePolyData(interpolated_bif_lower, centerline_new_bif_lower_path)

    interpolated_cl = merge_cl(interpolated_cl, div_points_rotated[1],
                               end_points_rotated[1])
    WritePolyData(interpolated_cl, centerline_new_path)
    
    # Interpolate voronoi diagram
    if lower and bif:
        bif = [interpolated_bif, interpolated_bif_lower, rotated_bif]
    elif bif:
        bif = [interpolated_bif, rotated_bif]
    elif lower:
        bif = [interpolated_bif_lower, rotated_bif]
    else:
        bif = []

    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, rotated_cl, 
                                                       rotated_voronoi,
                                                       end_points_rotated[0],
                                                       bif, lower, cylinder_factor)
    interpolated_voronoi = remove_distant_points(interpolated_voronoi, interpolated_cl)

    WritePolyData(interpolated_voronoi, voronoi_angle_path) 

    # Write a new surface from the new voronoi diagram
    print "Write new surface"
    new_surface = create_new_surface(interpolated_voronoi)
    WritePolyData(new_surface, model_new_surface)


if  __name__ == "__main__":
    smooth, angle, smooth_factor, l1, l2, bif, addPoint, basedir, case, lower, \
    cylinder_factor, w, version = read_command_line()
    folders = listdir(basedir) if case is None else [case]
    for folder in folders:
        if folder[:2] in ["P0", "C0"]:
            main(path.join(basedir, folder), smooth, smooth_factor, angle, l1,
                  l2, bif, addPoint, lower, cylinder_factor, w, version)
