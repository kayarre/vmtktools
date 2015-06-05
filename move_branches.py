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

    args = parser.parse_args()
    ang_ = math.pi/args.a

    return args.s, ang_, args.smooth_factor, args.leave1, args.leave2, args.bif


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
        tmp = np.dot(np.dot(np.dot(np.asarray(point) - O, M), R), R_inv) + O
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
        cell = vtk.vtkGenericCell()
        patch_cl.GetCell(i, cell)
        
        centerlineCellArray.InsertNextCell(cell.GetNumberOfPoints())
        
        start = cell.GetPoints().GetPoint(0)
        dist = line0.GetPoints().GetPoint(locator0.FindClosestPoint(start))
        test = math.sqrt(distance(start, dist)) > divergingRatioToSpacingTolerance
        
        if test:
            cellLine = ExtractSingleLine(patch_cl, i)

            locator = get_locator(cellLine)

            pnt1 = cellLine.GetPoints().GetPoint(locator.FindClosestPoint(div_points[1]))
            pnt2 = cellLine.GetPoints().GetPoint(locator.FindClosestPoint(div_points[2]))
            dist1 = math.sqrt(distance(pnt1, div_points[1]))
            dist2 = math.sqrt(distance(pnt2, div_points[2]))
            k = 1 if dist1 < dist2 else 2
            O = div_points[k]
            m = rotation_matrix[k]

        else:
            m = I
            O = np.array([0, 0, 0])
        
        for j in range(cell.GetNumberOfPoints()):
            point = np.asarray(cell.GetPoints().GetPoint(j))
            tmp = np.dot(np.dot(np.dot(point - O, m), R), R_inv) + O
            centerlinePoints.InsertNextPoint(tmp)
            radiusArray.SetTuple1(count, patch_cl.GetPointData().\
                            GetArray(radiusArrayName).GetTuple1(cell.GetPointId(j)))
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

    if tmp1[0] > tmp2[0]:
        m = {1: m2, 2: m1}

    if leave1:
        m[1] = I
    if leave2:
        m[2] = I

    return R, m


def get_startpoint(centerline):
    line = ExtractSingleLine(centerline, 0)
    return line.GetPoints().GetPoint(0)


def main(dirpath, smooth, smooth_factor, angle, l1, l2, bif):
    # Input filenames
    model_path = path.join(dirpath, "surface", "model.vtp")

    # Output names
    centerline_path = path.join(dirpath, "surface", "centerline_for_rotating.vtp")
    centerline_bif_path = path.join(dirpath, "surface", "centerline_bifurcation.vtp")
    centerline_new_path = path.join(dirpath, "surface", "centerline_interpolated.vtp")
    centerline_new_bif_path = path.join(dirpath, "surface", "centerline_interpolated.vtp")
    centerline_clipped_path = path.join(dirpath, "surface", "centerline_patch.vtp")
    centerline_rotated_path = path.join(dirpath, "surface", "centerline_rotated.vtp")
    centerline_rotated_path = path.join(dirpath, "surface", "centerline_rotated_bif.vtp")
    centerline_bif_clipped_path = path.join(dirpath, "surface", "centerline_bif_clipped.vtp")
    if smooth:
        voronoi_path = path.join(dirpath, "surface", "voronoi.vtp")
    else:
        voronoi_path = path.join(dirpath, "surface", "voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, "surface", "voronoi_clipped.vtp")
    voronoi_rotated_path = path.join(dirpath, "surface", "voronoi_rotated.vtp")
    if l1 or l2:
        s = "l1" if l1 else "l2"
        model_new_surface = path.join(dirpath, "surface", "model_angle_%.04f_%s.vtp" % (angle, s))
    else:
        model_new_surface = path.join(dirpath, "surface", "model_angle_%s.vtp" % angle)

    # Check if data already exsists, if not it is created
    # Model
    if not path.exists(model_path):
        print "The given directory: %s did not contain surface/model.vtp" % dirpath
        sys.exit(0)

    # Voronoi
    print "Compute voronoi diagram"
    voronoi = makeVoronoi(model_path, voronoi_path)

    # Centerline
    centerline = makeCenterline(model_path, centerline_path, length=0.1, smooth=False)
    centerline_bif = makeCenterline(model_path, centerline_bif_path,
                                    in_out=[0,1], length=0.1, smooth=False)

    # TODO: Control the curvature in the bifurcation
    
    # Create a tolerance for diverging
    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                    centerline.GetPoint(10), \
                                    centerline.GetPoint(11)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance

    # Get data from centerlines
    data = getData(centerline, centerline_bif, divergingTolerance)

    key = "end_point"
    div_points = [data["bif"][key], data[0][key], data[1][key]]
    clippingPoints = vtk.vtkPoints()
    for point in div_points:
        clippingPoints.InsertNextPoint(point)
    #div_points.append(get_startpoint(centerline))

    print "Clipping centerlines and voronoi diagram."
    patch_cl = CreateParentArteryPatches(centerline, clippingPoints)
    patch_bif = CreateParentArteryPatches(centerline_bif_clipped_path, clippingPoints)
    WritePolyData(patch_cl, centerline_clipped_path)
    
    # Clipp the voronoi diagram
    masked_voronoi = MaskVoronoiDiagram(voronoi, patch_cl)
    voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
    WritePolyData(voronoi_clipped, voronoi_clipped_path)
    
    # Rotate branches (Both centerline and voronoi)
    print "Rotate centerlines and voronoi diagram."
    R, m = rotationMatrix(data, angle, l1, l2)
    
    rotated_cl = rotate_cl(patch_cl, div_points, m, R)
    rotated_bif = rotate_cl(patch_bif, div_points, m, R)
    WritePolyData(rotated_cl, centerline_rotated_path)
    WritePolyData(rotated_bif, centerline_rotated_bif_path)

    rotated_voronoi = rotate_voronoi(voronoi_clipped, patch_cl, div_points, m, R)
    WritePolyData(rotated_voronoi, voronoi_rotated_path)

    # Interpolate the centerline
    print "Interpolate centerlines and voronoi diagram."
    interpolated_cl = InterpolatePatchCenterlines(rotated_cl, centerline, clippingPoints)
    interpolated_bif = InterpolatePatchCenterlines(rotated_bif, centerline, clippingPoints)
    WritePolyData(interpolated_cl, centerline_new_path)
    WritePolyData(interpolated_bif, centerline_new_bif_path)
    
    # Interpolate voronoi diagram
    bif = [interpolated_bif, rotated_bif] if bif else None
    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, patch_cl, 
                                                       rotated_voronoi,
                                                       clippingPoints, bif)
    WritePolyData(interpolated_voronoi, "voronoi_test.vtp")

    # Write a new surface from the new voronoi diagram
    print "Write new surface"
    new_surface = create_new_surface(interpolated_voronoi)
    WritePolyData(new_surface, model_new_surface)


if  __name__ == "__main__":
    smooth, angle, smooth_factor, l1, l2, bif = read_command_line()
    #basedir = "."
    #for folder in listdir(basedir):
    #    if folder[:2] in ["P0", "C0"]:
    #        main(path.join(basedir, folder), smooth, smooth_factor, angle)
    main("C0001", smooth, smooth_factor, angle, l1, l2, bif)
