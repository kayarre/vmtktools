from common import *
import time
from argparse import ArgumentParser
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from os import path, listdir
import numpy as np
from scipy.ndimage.filters import gaussian_filter

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram", metavar="smooth")
    parser.add_argument('--beta', type=float, default=1.5, 
            help="The new voronoi diagram is computed as (A/mean)**beta*r_old," + \
            " over the respective area. If beta < -1 the geometry will be more even, and" + \
            " if beta > 1, the differences in the geometry will be larger")
    parser.add_argument("--stats", type=bool, default=False,
                        help="Collect stats")

    args = parser.parse_args()

    return args.s, args.beta, args.stats


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


def get_stats(centerline_area, folder):
    WritePolyData(centerline_area, "test_area.vtp")

    area = get_array("CenterlineSectionArea", centerline_area)
    length = get_curvilinear_coordinate(centerline_area)
    MISR = get_array(radiusArrayName, centerline_area)

    # Basic parameters
    max = area.max()
    min = area.min()
    mean = area.mean()

    # Min max derived parameters
    length_min_max = abs(length[(area == min).nonzero()[0]] - \
                            length[(area == max).nonzero()[0]])
    max_mean_ratio = max / mean
    min_mean_ratio = min / mean

    # Area shape - "circleness"
    circleness = MISR**2*math.pi / area
    max_circleness = circleness.max()
    min_circleness = circleness.min()
    mean_circleness = circleness.mean()
        
    # Local extremas
    area_smooth = gaussian_filter(area, 0.1)

    from matplotlib.pyplot import plot, show, hold

    plot(length, area)
    hold("on")
    #plot(length, area_smooth)
    plot(length, MISR**2*math.pi)
    show()
    #sys.exit(0)

    # TODO: Need smoothening

    number_of_local_max = 2
    mean_dist_local_extrema = 2


def get_lineToChange(centerline, tol, type):
    # Check type of case and choose end point of siphon. If it is a normal case
    # or has a terminal aneurism, use ICA terminus as end point. If it is a
    # leteral case, use the aneurism as endpoint.
    if type in ["None", "terminal"]:
        split = 1
    elif type == "LAT":
        split = centerline.GetNumberOfLines() - 1

    line1 = ExtractSingleLine(centerline, split) 
    numberOfPoints1 = line1.GetNumberOfPoints()

    line2 = ExtractSingleLine(centerline, 0)
    numberOfPoints2 = line2.GetNumberOfPoints()

    N = min(numberOfPoints1, numberOfPoints2)
    for i in range(N):
        point1 = line1.GetPoints().GetPoint(i)
        point2 = line2.GetPoints().GetPoint(i)
        if math.sqrt(distance(point1, point2)) > tol:
            pointID = i
            break

    pointID = move_past_sphere(line1, i)
    lineToChange = ExtractSingleLine(centerline, 0, endID=pointID)

    return lineToChange


def change_area(voronoi, centerline, tol, beta, type):
    # Extract line to change 
    lineToChange = get_lineToChange(centerline, tol, type)
   
    # Tube function to evaluate if the voronoi point should be changed
    tubeFunction = vtkvmtk.vtkvmtkPolyBallLine()
    tubeFunction.SetInput(lineToChange)
    tubeFunction.SetPolyBallRadiusArrayName(radiusArrayName)

    # Make a sphere at the end of the line
    c1 = lineToChange.GetPoints().GetPoint(pointID)
    c2 = lineToChange.GetPoints().GetPoint(pointID - 1)
    r = get_array(radiusArrayName, lineToChange)[-1]
    t = [c1[i] - c2[i] for i in range(len(c1))]
    lastSphere = vtk.vtkSphere()
    lastSphere.SetRadius(r * 1.5)
    lastSphere.SetCenter(c1)

    # Array to change the radius
    #area_ = get_array("CenterlineSectionArea", lineToChange)
    #cylyndrenes = area_ - area
    #area_mean = (area_ + area) / 2.
    #length = get_curvilinear_coordinate(lineToChange)

    area = get_array(radiusArrayName, lineToChange)**2 * np.pi
    mean = np.mean(area)
    factor = (area / mean)**beta

    # A linear transition of the old and new geometry
    k = round(factor.shape[0] * 0.10, 0)
    l = factor.shape[0] - k
    trans = np.asarray(np.zeros(l).tolist() + np.linspace(0, 1, k).tolist())
    factor_ = factor[:,0]*(1-trans) + trans
    one = np.zeros(factor.shape[0]) + 1

    # Locator to find closest point on centerline
    locator = get_locator(lineToChange)

    # Voronoi diagram
    N = voronoi.GetNumberOfPoints()
    newVoronoi = vtk.vtkPolyData()
    voronoiPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, N)

    # If inside MISR tube and inside plane, change r.
    point = [0., 0., 0.]
    for i in range(N):
        voronoi.GetPoint(i, point)
        pointRadius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        
        tubeValue = tubeFunction.EvaluateFunction(point)
        sphereValue = lastSphere.EvaluateFunction(point)
        voronoiVector = [point[j] - c1[j] for j in range(3)]
        vectorValue = vtk.vtkMath.Dot(voronoiVector, t)

        if (sphereValue < 0.0) & (vectorValue < 0.0):
            pointRadius = pointRadius
        elif (tubeValue <= 0.0):
            pointRadius = pointRadius*factor[locator.FindClosestPoint(point)]

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

    # Smooth voronoi diagram
    centerlines = makeCenterline(model_path, centerlines_path, length=0.1,
                                 smooth=True, factor=0.1)
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

    centerline_area = makeCenterlineSections(model_path,
                                             centerlines_path,
                                             centerline_area_path)
    
    # Tolerance for finding diverging point
    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                centerline_area.GetPoint(10), \
                                centerline_area.GetPoint(11)))
    tol = centerlineSpacing / divergingRatioToSpacingTolerance

    # Get line to change
    type_ = getParameters(folder)["aneurysmType"]
    centerline_to_change = get_lineToChange(centerline_area, tol, type_)

    # Compute stats
    if stats:
        get_stats(centerline_to_change, folder)

    # Change and compute the new voronoi diagram
    #newvoronoi = change_area(voronoi, centerline_to_change, tol, beta, type_)
    #WritePolyData(newvoronoi, voronoi_new_path) 

    # Make new surface
    #surface_smoothed = create_new_surface(newvoronoi)
    #WritePolyData(surface_smoothed, model_area_path)


if __name__ == '__main__':
    beta, smooth, stats = read_command_line()
    basefolder = "/home/aslak/master/src/aneurysms"
    for folder in listdir(basefolder):
        if folder.startswith("C0"):
            print folder
            case = path.join(basefolder, folder)
            main(case, beta, smooth, stats)
