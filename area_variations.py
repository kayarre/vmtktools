from common import *
import time
from argparse import ArgumentParser
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from scipy.interpolate import splrep, splev
from os import path, listdir
import numpy as np
import math
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelextrema

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()
   
    parser.add_argument('--case', type=str, default=None, help="Case")
    parser.add_argument('--dir_path', type=str, default=".", help="Path")
    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram", metavar="smooth")
    parser.add_argument('--beta', type=float, default=0.5, 
            help="The new voronoi diagram is computed as (A/mean)**beta*r_old," + \
            " over the respective area. If beta < -1 the geometry will be more even, and" + \
            " if beta > 1, the differences in the geometry will be larger")
    parser.add_argument("--ratio", type=float, default=None, help="Wanted ratio" + \
                       " A_max/A_min, when this is given beta will be ignored" + \
                       " and beta computed such that this will (approxematly)" + \
                       " be the result")
    parser.add_argument("--stats", type=bool, default=False,
                        help="Collect stats")
    parser.add_argument("--make_plot", type=bool, default=False, help="Make plots")

    args = parser.parse_args()

    if args.ratio is not None and args.beta != 0.5:
        print "WARNING: The beta value you provided will be ignored."

    return args.s, args.beta, args.stats, args.dir_path, args.case, args.ratio, args.make_plot


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

    #print "max_min_ratio_area:", global_max_area / global_min_area
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
    
    return length, area


def get_lineToChange(centerline, tol):
    line2 = ExtractSingleLine(centerline, 0)
    numberOfPoints2 = line2.GetNumberOfPoints()

    n = 2
    pointIDs = []
    for j in range(1, n):
        line1 = ExtractSingleLine(centerline, j)
        numberOfPoints1 = line1.GetNumberOfPoints()

        N = min(numberOfPoints1, numberOfPoints2)
        for i in range(N):
            point1 = line1.GetPoints().GetPoint(i)
            point2 = line2.GetPoints().GetPoint(i)
            if math.sqrt(distance(point1, point2)) > tol:
                pointID = i
                break

        pointIDs.append(pointID) #move_past_sphere(line2, pointID))

    pointID = min(pointIDs)
    lineToChange = ExtractSingleLine(centerline, 0, endID=pointID)
    lineToChange = splineCenterline(lineToChange, nknots=25)
    
    curvature = get_array("Curvature", line)
    curvature = get_array("Torsion", line)
    length = get_curvilinear_coordinate(line)
    
    from matplotlib.pylab import *
    length = get_curvilinear_coordinate(lineToChange)
    plot(length, curvature)
    hold("on")
    plot(length, torsion)
    show()
    
    sys.exit(0)

    return lineToChange


def get_factor(lineToChange, beta, ratio):
    # Array to change the radius
    area = get_array("CenterlineSectionArea", lineToChange)
    for i in range(2):
        area = gaussian_filter(area, 5)
    mean = np.mean(area)

    if ratio is not None:
        # Inital guess
        R_old = area.max() / area.min()
        beta = 0.5 * math.log(ratio / R_old) / math.log(ratio) + 1

        # Parameters for algorithm
        R = 1e10
        a = 0
        b = 2
        sign = 0 if ratio > R_old else 1
        max_iter = 30
        iter = 0

        # Exclude first and last 10 %
        area_ = area[int(area.shape[0]*0.02):-int(area.shape[0]*0.02)]

        #print "Area 80% max:", area_.max(), "area min", area_.min()
        #print "Area full max:", area.max(), "area min", area.min()

        while abs(R - ratio) >= 0.001 and iter < max_iter:
            #print beta - 1, "Diff:", abs(R - ratio)
            factor_ = (area / mean)**(beta-1)
            k = round(factor_.shape[0] * 0.10, 0)
            l = factor_.shape[0] - k*2
            trans = np.asarray(np.linspace(1, 0, k).tolist() +
            np.zeros(l).tolist() + np.linspace(0, 1, k).tolist())
            factor = factor_[:,0]*(1-trans) + trans

            area_new = (np.sqrt(area[:,0]/math.pi)*factor)**2 * math.pi
            R = area_new.max() / area_new.min()
            
            print "R now:", R, "Want:", ratio, "R_old:", R_old
            if R < ratio:
                a = beta
                beta = a + (b - a) / 2.
            else:
                b = beta
                beta = a + (b - a) / 2.

            iter += 1

        print beta - 1
        beta = beta - 1

    else:
        factor_ = (area / mean)**beta

        # A linear transition of the old and new geometry
        k = round(factor_.shape[0] * 0.10, 0)
        l = factor_.shape[0] - k*2
        trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(l).tolist() + np.linspace(0, 1, k).tolist())
        factor = factor_[:,0]*(1-trans) + trans

    return factor


def change_area(voronoi, lineToChange, tol, beta, ratio):
    # NEW: Check if voronoi point is within 1.5*MISR of centerline.
    #      So, this means change MISR of lineToChange, due to extreme
    #      circleness
    # OLD: Tube function to evaluate if the voronoi point should be changed

    arrayForTube = get_vtk_array("TubeRadius", 1, lineToChange.GetNumberOfPoints())
    MISR = get_array(radiusArrayName, lineToChange)*1.7
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

    # Get factor    
    factor = get_factor(lineToChange, beta, ratio)

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
            # TODO: Move point or project vector down in to prependicular
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


def main(folder, beta, smooth, stats, r_change):
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

    if not path.exists(centerline_area_spline_path):
        centerline_to_change = get_lineToChange(centerlines, tol)
        WritePolyData(centerline_to_change, centerline_area_spline_path)
        
        centerline_splined = splineCenterline(centerline_to_change, 20)
        centerline_area = makeCenterlineSections(model_path,
                                        centerline_area_spline_path,
                                        centerline_area_path, recompute=True)
        WritePolyData(centerline_area, centerline_area_spline_path)
    else:
        centerline_area = ReadPolyData(centerline_area_spline_path)

    # Compute stats
    area = None
    length = None
    if stats:
        length, area = get_stats(centerline_area, folder, centerlines)
    else:
        # Change and compute the new voronoi diagram
        print "Change Voronoi diagram"
        newvoronoi = change_area(voronoi, centerline_area, tol, beta, ratio)
        print "Write Voronoi diagram"
        WritePolyData(newvoronoi, voronoi_new_path)

        # Make new surface
        print "Create surface"
        surface_smoothed = create_new_surface(newvoronoi)
        print "Write surface", model_area_path
        WritePolyData(surface_smoothed, model_area_path)

    return length, area


def make_plots(smooth, beta, stats, basefolder, case, ratio):
    folders = []
    folders = folders + listdir(path.join(basefolder, "new_cases"))
    folders = ["N0124", "N0123", "N0125"]
    folders = ["A0026", "A0027"]
    color = ["c", "k", "r"]
    label = {"N0123": "Orginal", "N0124": "-1 SD", "N0125": "+1 SD", "N0132": "+2 SD"}
    from matplotlib.pylab import figure, legend, plot, savefig, show, hold, \
                                 xlabel, ylabel, tight_layout, axis, subplot, \
                                 tick_params, xlim

    #f = figure(figsize=(10, 3))
    for i, folder in enumerate(folders):
        if folder.startswith("B0"): #or folder.startswith("P0")) and \
          #not ".png" in folder or folder.startswith("N0"):
            #if folder in ["C0023", "C0099", "C0057b", "C0093", "C0087"]: continue
            if folder in ["B0010", "C0087", "C0093"]:
                continue
            print "Working on:", folder
            case = path.join(basefolder, folder)
            length, a = main(case, beta, smooth, stats, ratio)

            # Plot
            f = plot(length, a, label=label[folder], linewidth=2, color=color[i])
            hold("on")

        ylabel("Area", fontsize="large")

    #### FOR AREA VARIATION PLOT ####
    """
    tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        left='off',        # ticks along the bottom edge are off
        right='off',       # ticks along the top edge are off
        labelleft='off')   # labels along the bottom edge are off

    tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        top='off')         # ticks
    """

    #xlabel("Length", fontsize="large")
    #xlim([0, 140])
    #legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
    #                borderaxespad=0., fontsize="large")
    #tight_layout()
    #savefig("area_variations_%s.eps" % folder) # pad_size=0
    #show()
