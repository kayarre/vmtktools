from vmtkClassifyBifucation import *
from common import *
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from os import path, listdir
from subprocess import STDOUT, check_output
import sys


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

    radiusArray = vtk.vtkDoubleArray()
    radiusArray.SetNumberOfComponents(1)
    radiusArray.SetNumberOfTuples(numberOfPoints)
    radiusArray.SetName(radiusArrayName)
    radiusArray.FillComponent(0, 0.0)

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

    radiusArray = vtk.vtkDoubleArray()
    radiusArray.SetName(radiusArrayName)
    radiusArray.SetNumberOfComponents(1)
    radiusArray.SetNumberOfTuples(numberOfPoints)
    radiusArray.FillComponent(0,0.0)

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


def rotationMatrix(data):
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
    m1 = np.asarray([[np.cos(np.pi / 8), -np.sin(np.pi / 8), 0],
                     [np.sin(np.pi / 8),  np.cos(np.pi / 8), 0],
                     [0, 0, 1]])
    m2 = np.asarray([[ np.cos(np.pi / 8), np.sin(np.pi / 8), 0],
                     [-np.sin(np.pi / 8), np.cos(np.pi / 8), 0],
                     [0, 0, 1]])

    m = {1: m1, 2: m2}
    tmp1 = data[0]["div_point"] - d
    tmp2 = data[1]["div_point"] - d

    if tmp1[0] > tmp[2]:
        m = {1: m2, 2: m1}

    return R, m


def get_startpoint(centerline):
    line = ExtractSingleLine(centerline, 0)
    return line.GetPoints().GetPoint(0)


def main(dirpath, smooth=1, dir_out_path=None):
    # Input filenames
    model_path = path.join(dirpath, "surface", "model.vtp")

    # Output names
    dirpath = dirpath if dir_out_path is None else dir_out_path
    centerline_bif_path = path.join(dirpath, "surface", "cl_bif.vtp")
    centerline_path = path.join(dirpath, "surface", "cl.vtp")
    centerline_new_path = path.join(dirpath, "surface", "cl_new.vtp")
    centerline_clipped_path = path.join(dirpath, "surface", "cl_patch.vtp")
    centerline_rotated_path = path.join(dirpath, "surface", "cl_rotated.vtp")
    voronoi_path = path.join(dirpath, "surface", "voronoi.vtp")
    voronoi_clipped_path = path.join(dirpath, "surface", "voronoi_clipped.vtp")
    voronoi_rotated_path = path.join(dirpath, "surface", "voronoi_rotated.vtp")
    voronoi_new_surface = path.join(dirpath, "surface", "model_new.vtp")

    # Check if data already exsists, if not it is created
    # Model
    if not path.exists(model_path):
        print "The given directory: %s did not contain surface/model.vtp" % dirpath
        sys.exit(0)

    # Voronoi
    print "Compute voronoi diagram"
    voronoi = makeVoronoi(model_path, voronoi_path)

    # Centerline
    if not path.exists(centerline_path):
        print "Compute centerlines"
        print "The program need a centerline that goes from the inlet to one \nthe outlet at" + \
              "inter part and one outlet at the anterior part. \nPleace select" + \
              "the inlet as source point and the two oters as outlets. \nPress any key when" + \
              " you want to continue."
        raw_input()
        a = check_output(("vmtk vmtkcenterlines -ifile %s --pipe " + \
                          "vmtkcenterlineresampling -length 0.1 -ofile %s") %
                          (model_path, centerline_path), stderr=STDOUT, shell=True)
        if not success(a):
            print "Something went wrong when computing the centerline"
            print a
            sys.exit(0)

    # Centerline bif
    if not path.exists(centerline_bif_path):
        points = get_endpoints(centerline_path)
        a = check_output(("vmtk vmtkcenterlines -ifile %s -seedselector "
                          "pointlist -sourcepoints %s %s %s -targetpoints %s %s %s" + \
                          " --pipe vmtkcenterlineresampling -length 0.1 -ofile %s") % 
                          (model_path, points[0][0], points[0][1], points[0][2], 
                           points[1][0], points[1][1], points[1][2], centerline_bif_path), 
                          stderr=STDOUT, shell=True)
        if not success(a):
            print "Something went wrong when computing the centerline"
            print a
            sys.exit(0)

    # Read data
    voronoi = ReadPolyData(voronoi_path)
    centerline = ReadPolyData(centerline_path)
    centerline_bif = ReadPolyData(centerline_bif_path)

    # TODO: Read options from commandline
    # TODO: Smooth option
    # TODO: Optional dirpath
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
    div_points.append(get_startpoint(centerline))

    # Clipp the centerline
    print "Clipping centerlines and voronoi diagram."
    patch_cl = CreateParentArteryPatches(centerline, clippingPoints)
    WritePolyData(patch_cl, centerline_clipped_path)
    
    # Clipp the voronoi diagram
    masked_voronoi = MaskVoronoiDiagram(voronoi, patch_cl)
    voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
    WritePolyData(voronoi_clipped, voronoi_clipped_path)

    # Rotate branches (Both centerline and voronoi)
    print "Rotate centerlines and voronoi diagram."
    R, m = rotationMatrix(data)
    
    rotated_cl = rotate_cl(patch_cl, div_points, m, R)
    WritePolyData(rotated_cl, centerline_rotated_path)

    rotated_voronoi = rotate_voronoi(voronoi_clipped, patch_cl, div_points, m, R)
    WritePolyData(rotated_voronoi, voronoi_rotated_path)

    # Interpolate the centerline
    print "Interpolate centerlines and voronoi diagram."
    interpolated_cl = InterpolatePatchCenterlines(rotated_cl, centerline, clippingPoints)
    WritePolyData(interpolated_cl, centerline_new_path)
    
    # Interpolate voronoi diagram
    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, patch_cl, 
                                                       rotated_voronoi, clippingPoints)

    # Write a new surface from the new voronoi diagram
    print "Write new surface"
    new_surface = create_new_surface(interpolated_voronoi)
    WritePolyData(new_surface, voronoi_new_surface)


if  __name__ == "__main__": 
    main("C0001")
