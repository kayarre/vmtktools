from common import *
import time
from argparse import ArgumentParser
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from scipy.interpolate import splrep, splev
from os import path, listdir
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelextrema

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()
   
    parser.add_argument('--case', type=str, default=None, help="Case")
    parser.add_argument('--dir_path', type=str, default=".", help="Path")
    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram", metavar="smooth")
    parser.add_argument('--beta', type=float, default=1.5, 
            help="The new voronoi diagram is computed as (A/mean)**beta*r_old," + \
            " over the respective area. If beta < -1 the geometry will be more even, and" + \
            " if beta > 1, the differences in the geometry will be larger")
    parser.add_argument("--stats", type=bool, default=False,
                        help="Collect stats")

    args = parser.parse_args()

    return args.s, args.beta, args.stats, args.dir_path, args.case


def move_past_sphere(centerline, i):
    center = centerline.GetPoints().GetPoint(i)
    r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
    MISphere = vtk.vtkSphere()
    MISphere.SetCenter(center)
    MISphere.SetRadius(r)

    for j in range(i, 0, -1):
        value = MISphere.EvaluateFunction(centerline.GetPoint(j))
        if value >= 0:
            break

    return j


def splineCenterline(line, knots):
    # Allocate data structure to store centerline points
    data = np.zeros((line.GetNumberOfPoints(), 3))
    MISR = get_array(radiusArrayName, line)

    # Collect data from centerline
    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i,:3] = line.GetPoints().GetPoint(i)

    t = np.linspace(curv_coor[0], curv_coor[-1], knots+2)[1:-1]
    fx = splrep(curv_coor, data[:,0], k=4, t=t)
    fy = splrep(curv_coor, data[:,1], k=4, t=t)
    fz = splrep(curv_coor, data[:,2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    # Store data for converting to vtkPolyLine
    data = np.zeros((len(curv_coor), 4))
    data[:,0] = fx_
    data[:,1] = fy_
    data[:,2] = fz_
    data[:,3] = MISR[:,0]

    header = ["X", "Y", "Z", radiusArrayName]
    line = data_to_vtkPolyData(data, header)

    return line


def get_stats(centerline_area, folder, centerline):
    WritePolyData(centerline_area, "test_area.vtp")

    area = get_array("CenterlineSectionArea", centerline_area)
    MISR_ = get_array(radiusArrayName, centerline)**2*math.pi
    MISR = np.array([MISR_[i] for i in range(MISR_.shape[0] - 1, -1, -1)])[:area.shape[0]]
    length = get_curvilinear_coordinate(centerline_area)

    for i in range(2):
        area = gaussian_filter(area, 5)
        MISR = gaussian_filter(MISR, 5)
    
    # Area shape - "circleness"
    circleness = MISR / area
    max_circleness = circleness.max()
    min_circleness = circleness.min()
    mean_circleness = circleness.mean()
        
    # Local extremas, ignore min or max on boundary
    local_min_MISR_ID = argrelextrema(MISR, np.less)[0]
    local_max_MISR_ID = argrelextrema(MISR, np.greater)[0]
    local_min_area_ID = argrelextrema(area, np.less)[0]
    local_max_area_ID = argrelextrema(area, np.greater)[0]

    local_min_MISR = MISR[local_min_MISR_ID]
    local_max_MISR = MISR[local_max_MISR_ID]
    local_min_area = area[local_min_area_ID]
    local_max_area = area[local_max_area_ID]

    global_min_MISR = local_min_MISR.min()
    global_max_MISR = local_max_MISR.max()
    global_min_area = local_min_area.min()
    global_max_area = local_max_area.max()

    mean_area = area.mean()
    mean_MISR = MISR.mean()

    number_of_max = local_max_area.shape[0]

    # Min max derived parameters
    length_min_max = abs(length[(area == global_min_area).nonzero()[0]] - \
                         length[(area == global_max_area).nonzero()[0]])[0]

    max_mean_ratio_area = global_max_area / mean_area
    min_mean_ratio_area = global_min_area / mean_area
    max_mean_ratio_MSIR = global_max_MISR / mean_MISR
    min_mean_ratio_MISR = global_min_MISR / mean_MISR

    # Global and local disent
    global_min_max_disent = abs(math.sqrt(global_max_area)/math.pi -
                                math.sqrt(global_min_area) / math.pi) / \
                            length_min_max
    local_max_stepest_disent = 0

    if length[(area == local_min_area[0]).nonzero()[0]] > \
       length[(area == local_max_area[0]).nonzero()[0]]:
        start = 1
    else:
        start = 0

    N = min(number_of_max, local_min_area.shape[0] - start)
    for i in range(N):
        min_ = local_min_area[start + i]
        max_ = local_max_area[i]

        h = math.sqrt(max_)/math.pi - math.sqrt(min_) / math.pi
        l = abs(length[(area == max_).nonzero()[0]] - length[(area == min_).nonzero()[0]])
        if h/l > local_max_stepest_disent:
            local_max_stepest_disent = h / l
    
    
    # Max point disent (max |derivative|)
    knots = 40
    t = np.linspace(length[0], length[-1], knots+2)[1:-1]
    spline_area = splrep(length, area, k=4, t=t)
    spline_area_ = splev(length, spline_area)
    darea_dx = splev(length, spline_area, der=1)
    max_derivative = abs(darea_dx).max()

    stats = {"max_derivative": max_derivative,
             "local_max_stepest_disent": local_max_stepest_disent,
             "max_mean_ratio_area": max_mean_ratio_area,
             "min_mean_ratio_area" : min_mean_ratio_area,
             "mean_area": mean_area,
             "max_min_ratio_area": global_max_area / global_min_area,
             "length_min_max": length_min_max,
             "global_min_area": global_min_area,
             "global_max_area": global_max_area,
             "max_circleness": max_circleness,
             "min_circleness": min_circleness,
             "mean_circleness": mean_circleness,
             "number_of_max": number_of_max}


    writeParameters(stats, folder)


def get_lineToChange(centerline, tol, type):
    # NEW: Choose the lowest bif, this could either be something in the siphon,
    # the aneurism or ICA terminus.
    # OLD: Check type of case and choose end point of siphon. If it is a normal case
    # or has a terminal aneurism, use ICA terminus as end point. If it is a
    # leteral case, use the aneurism as endpoint.

    #if type is None or type == "terminal":
    #    split = 1
    #elif type == "LAT":
    #    split = centerline.GetNumberOfLines() - 1

    line2 = ExtractSingleLine(centerline, 0)
    numberOfPoints2 = line2.GetNumberOfPoints()

    pointIDs = []
    for j in range(1, centerline.GetNumberOfCells()):
        line1 = ExtractSingleLine(centerline, j)
        numberOfPoints1 = line1.GetNumberOfPoints()

        N = min(numberOfPoints1, numberOfPoints2)
        for i in range(N):
            point1 = line1.GetPoints().GetPoint(i)
            point2 = line2.GetPoints().GetPoint(i)
            if math.sqrt(distance(point1, point2)) > tol:
                pointID = i
                break

        pointIDs.append(move_past_sphere(line2, pointID))

    pointID = min(pointIDs)
    lineToChange = ExtractSingleLine(centerline, 0, endID=pointID)

    return lineToChange


def change_area(voronoi, lineToChange, tol, beta, type):
    # NEW: Check if voronoi point is within 2*MISR of centerline.
    #      So, this means change MISR of lineToChange, due to extreme
    #      circleness
    # OLD: Tube function to evaluate if the voronoi point should be changed
    arrayForTube = get_vtk_array("TubeRadius", 1, lineToChange.GetNumberOfPoints())
    MISR = get_array(radiusArrayName, lineToChange)*1.5
    for i in range(MISR.shape[0]):
        arrayForTube.SetTuple1(i, MISR[i])
    lineToChange.GetPointData().AddArray(arrayForTube)

    tubeFunction = vtkvmtk.vtkvmtkPolyBallLine()
    tubeFunction.SetInput(lineToChange)
    tubeFunction.SetPolyBallRadiusArrayName("TubeRadius")

    # Make a sphere at the end of the line
    pointID = lineToChange.GetNumberOfPoints()
    c1 = lineToChange.GetPoints().GetPoint(pointID)
    c2 = lineToChange.GetPoints().GetPoint(pointID - 1)
    r = get_array(radiusArrayName, lineToChange)[-1]
    t = [c1[i] - c2[i] for i in range(len(c1))]
    lastSphere = vtk.vtkSphere()
    lastSphere.SetRadius(r * 1.5)
    lastSphere.SetCenter(c1)

    # Array to change the radius
    area = get_array("CenterlineSectionArea", lineToChange)
    for i in range(10):
        area = gaussian_filter(area, 6)
    mean = np.mean(area)
    factor_ = (area / mean)**beta

    # A linear transition of the old and new geometry
    k = round(factor_.shape[0] * 0.10, 0)
    l = factor_.shape[0] - k
    trans = np.asarray(np.zeros(l).tolist() + np.linspace(0, 1, k).tolist())
    factor = factor_[:,0]*(1-trans) + trans
    one = np.zeros(factor_.shape[0]) + 1

    # Locator to find closest point on centerline
    locator = get_locator(lineToChange)

    # Voronoi diagram
    N = voronoi.GetNumberOfPoints()
    newVoronoi = vtk.vtkPolyData()
    voronoiPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, N)

    # If inside MISR tube and inside plane, change r and move point.
    point = [0., 0., 0.]
    for i in range(N):
        voronoi.GetPoint(i, point)
        pointRadius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        
        tubeValue = tubeFunction.EvaluateFunction(point)
        sphereValue = lastSphere.EvaluateFunction(point)
        voronoiVector = [point[j] - c1[j] for j in range(3)]
        vectorValue = vtk.vtkMath.Dot(voronoiVector, t)

        if (tubeValue <= 0.0) and not ((sphereValue < 0.0) and (vectorValue < 0.0)):
            # FIXME: Move point or project vector down in to prependicular
            # plane, or find vector from spline for continous evaluation
            tmp_ID = locator.FindClosestPoint(point)
            v1 = np.asarray(lineToChange.GetPoint(tmp_ID)) - np.asarray(point)
            v2 = v1 * (1 - factor[tmp_ID])
            point = (np.asarray(point) + v2).tolist()

            # Change radius
            pointRadius = pointRadius*factor[tmp_ID]

        voronoiPoints.InsertNextPoint(point)
        radiusArray.SetTuple1(i, pointRadius)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

    newVoronoi.SetPoints(voronoiPoints)
    newVoronoi.SetVerts(cellArray)
    newVoronoi.GetPointData().AddArray(radiusArray)

    return newVoronoi


def main(folder, beta, smooth, stats):
    # Naming convention
    model_path = path.join(folder, "surface", "model.vtp")
    model_smoothed_path = path.join(folder, "surface", "model_smoothed.vtp")
    s = "%s" % beta if not smooth else "%s_smooth" % beta
    model_area_path = path.join(folder, "surface", "model_area_%s.vtp" % s)
    voronoi_path = path.join(folder, "surface", "voronoi.vtp")
    voronoi_smoothed_path = path.join(folder, "surface", "voronoi_smoothed.vtp")
    voronoi_new_path = path.join(folder, "surface", "voronoi_area.vtp")
    centerlines_path = path.join(folder, "surface", "centerline_for_area.vtp")
    centerline_area_path = path.join(folder, "surface", "centerline_area.vtp")
    centerline_area_spline_path = path.join(folder, "surface", "centerline_area_spline.vtp")

    # Smooth voronoi diagram
    centerlines = makeCenterline(model_path, centerlines_path, length=0.1,
            smooth=False)
    voronoi = makeVoronoi(model_path, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerlines, 0.25)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)

        surface_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(surface_smoothed, model_smoothed_path)

    voronoi_smoothed = ReadPolyData(voronoi_smoothed_path)

    # Use smoothed voronoi or not
    voronoi = voronoi_smoothed if smooth else voronoi
    model_path = model_smoothed_path if smooth else model_path

    # Tolerance for finding diverging point
    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                centerlines.GetPoint(10), \
                                centerlines.GetPoint(11)))
    tol = centerlineSpacing / divergingRatioToSpacingTolerance

    # Get line to change
    type_ = getParameters(folder)["aneurysmType"]

    if not path.exists(centerline_area_spline_path):
        centerline_to_change = get_lineToChange(centerlines, tol, type_)
        WritePolyData(centerline_to_change, centerline_area_spline_path)
        
        centerline_splined = splineCenterline(centerline_to_change, 20)
        centerline_area = makeCenterlineSections(model_path,
                                        centerline_area_spline_path,
                                        centerline_area_path, recompute=True)
        WritePolyData(centerline_area, centerline_area_spline_path)
    else:
        centerline_area = ReadPolyData(centerline_area_spline_path)

    # Compute stats
    if stats:
        get_stats(centerline_area, folder, centerlines)
    else:
        # Change and compute the new voronoi diagram
        newvoronoi = change_area(voronoi, centerline_area, tol, beta, type_)
        WritePolyData(newvoronoi, voronoi_new_path) 

        # Make new surface
        surface_smoothed = create_new_surface(newvoronoi)
        WritePolyData(surface_smoothed, model_area_path)


if __name__ == '__main__':
    smooth, beta, stats, basefolder, case = read_command_line()
    folders = listdir(basefolder) if case is None else [case]
    for folder in folders:
        if (folder.startswith("C0") or folder.startswith("P0")) and not ".png" in folder:
            print folder
            case = path.join(basefolder, folder)
            main(case, beta, smooth, stats)
