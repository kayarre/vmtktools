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

    args = parser.parse_args()
    ang_ = math.pi/args.a

    return args.s, ang_, args.smooth_factor, args.leave1, args.leave2, \
           args.bif, args.addPoint, args.d, args.case, args.lower


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
        k = 1 if data[0][r_end] > data[1][r_end] else 2
        m[2] = I
    if leave2:
        k = 1 if data[0][r_end] < data[1][r_end] else 2
        m[2] = I

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


def main(dirpath, smooth, smooth_factor, angle, l1, l2, bif, addPoint, lower):
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
    s = "_pi%s" % (1/(angle/math.pi))
    s += "" if not l1 else "_l1"
    s += "" if not l2 else "_l2"
    s += "" if not bif else "_bif"
    s += "" if not smooth else "_smooth"
    s += "" if not addPoint else "_extraPoint"
    s += "" if not lower else "_lower"
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
    addPoint_ = False
    lower_ = [False]
    interpolated_cl = InterpolatePatchCenterlines(rotated_cl, centerline,
                                                  end_points_rotated[0],
                                                  div_points_rotated[0],
                                                  addPoint)
    interpolated_bif = InterpolatePatchCenterlines(rotated_bif, centerline_bif,
                                                   end_points_rotated_bif[0],
                                                   div_points_rotated_bif[0],
                                                   False,)
    if lower:
        print "Interpolate the lower thing"
        center = np.sum(div_points[1], axis=0)/3.
        div_points_rotated_bif[0].SetPoint(0, center[0], center[1], center[2])
        interpolated_bif_lower = InterpolatePatchCenterlines(rotated_bif, centerline_bif,
                                                       end_points_rotated_bif[0],
                                                       div_points_rotated_bif[0],
                                                       addPoint)
        WritePolyData(interpolated_bif_lower, centerline_new_bif_lower_path)
    WritePolyData(interpolated_cl, centerline_new_path)
    WritePolyData(interpolated_bif, centerline_new_bif_path)

    # Interpolate voronoi diagram
    if lower and bif:
        bif = [interpolated_bif, interpolated_bif_lower, rotated_bif]
    elif bif:
        bif = [interpolated_bif, rotated_bif] if bif else None
    elif lower:
        bif = [interpolated_bif_lower, rotated_bif] if lower else None

    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, rotated_cl, 
                                                       rotated_voronoi,
                                                       end_points_rotated[0],
                                                       bif, lower)
    new_surface = create_new_surface(interpolated_voronoi)
    interpolated_voronoi = remove_distant_points(interpolated_voronoi, interpolated_cl)

    # Write a new surface from the new voronoi diagram
    print "Write new surface"
    new_surface = create_new_surface(interpolated_voronoi)
    WritePolyData(new_surface, model_new_surface)


if  __name__ == "__main__":
    smooth, angle, smooth_factor, l1, l2, bif, addPoint, basedir, case, lower = read_command_line()
    folders = listdir(basedir) if case is None else [case]
    for folder in folders:
        if folder[:2] in ["P0", "C0"]:
            main(path.join(basedir, folder), smooth, smooth_factor, angle, l1,
                  l2, bif, addPoint, lower)
