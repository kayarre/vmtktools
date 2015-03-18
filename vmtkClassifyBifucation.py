import vmtk
import vtk
import math
import sys
import time
import subprocess 
import numpy as np
from os import path, listdir


# Global variables
divergingRatioToSpacingTolerance = 2.0
radiusArrayName = 'MaximumInscribedSphereRadius'

def ReadPolyData(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def WritePolyData(input,filename):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInput(input)
    writer.Write()


def getData(centerline, centerline_bif, toll):
    # For speed up, create local function
    distance = vtk.vtkMath.Distance2BetweenPoints

    # Declear variables before loop incase values are not found
    diverging_point_ID = -1
    diverging_point = [0.0, 0.0, 0.0]
    diverging_point_MISR = -1

    clipping_point_ID = -1
    clipping_point = [0.0, 0.0, 0.0]
    data = {"bif":{}, 0:{}, 1:{}}

    # List of points conected to ID
    points_ids_0 = vtk.vtkIdList()
    points_ids_1 = vtk.vtkIdList()

    # One is the branch to the left and the other is the one to the right
    centerline.GetCellPoints(0, points_ids_0)
    centerline.GetCellPoints(1, points_ids_1)

    # Find lower clipping point
    N = min(points_ids_0.GetNumberOfIds(), points_ids_1.GetNumberOfIds())
    for i in range(0, N):
        cell_point_0 = centerline.GetPoint(points_ids_0.GetId(i))
        cell_point_1 = centerline.GetPoint(points_ids_1.GetId(i))
        
        distance_between_points = math.sqrt(distance(cell_point_0, cell_point_1))
        if distance_between_points > toll:
            tmpI = i
            point_ID_0 = points_ids_0.GetId(i)
            point_ID_1 = points_ids_1.GetId(i)
            center = centerline.GetPoint(point_ID_0)
            r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID_0)
            break

    end, r_end = move_past_sphere(centerline, center, r, point_ID_0,
                                  stop=point_ID_0*100, step=1)
    data["bif"]["end_point"] = end
    data["bif"]["r_end"] = r_end
    data["bif"]["div_point"] = center
    data["bif"]["ID_div"] = point_ID_0
    data["bif"]["i_div"] = tmpI
    data["bif"]["r_div"] = r

    # Find the diverging point for anterior and midt bifurcation
    # continue further downstream in each direction and stop when
    # a point is closer than toll, than move point MISR * X
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(centerline_bif)
    locator.BuildLocator()

    counter = 0
    for point_ids in [points_ids_0, points_ids_1]:
        for i in range(tmpI, point_ids.GetNumberOfIds(), 1):
            tmp_point = centerline.GetPoint(point_ids.GetId(i))
            closest_point_ID = locator.FindClosestPoint(tmp_point)
            closest_point = centerline_bif.GetPoint(closest_point_ID)
            distance_between_points = distance(tmp_point, closest_point)
            if distance_between_points < toll:
                point_ID = point_ids.GetId(i)
                center = centerline.GetPoint(point_ID)
                r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID)
                break
        
        end, r_end = move_past_sphere(centerline, center, r, point_ID)
        data[counter]["end_point"] = end
        data[counter]["r_end"] = r_end
        data[counter]["r_div"] = r
        data[counter]["ID_end"] = locator.FindClosestPoint(data[counter]["end_point"])
        data[counter]["ID_div"] = locator.FindClosestPoint(center)
        data[counter]["div_point"] = center
        
        counter += 1
        
    return data


def move_past_sphere(centerline, center, r, start, step=-1, stop=0, X=2.5):
    """Moves a point along the centerline until it as outside MIS"""
    # Create the minimal inscribed sphere
    MISphere = vtk.vtkSphere()
    MISphere.SetCenter(center)
    MISphere.SetRadius(r * X)
    tempPoint = [0.0, 0.0, 0.0]

    # Go the length of one MISR backwards
    for i in range(start, stop, step):
        value = MISphere.EvaluateFunction(centerline.GetPoint(i))
        if (value>=0.0):
            tempPoint = centerline.GetPoint(i)
            break

    r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
    return tempPoint, r


def get_endpoints(centerline):
    """Returns the endpoints of the two streamlines"""
    end = []
    N_ = min(centerline.GetCell(0).GetNumberOfPoints(),
             centerline.GetCell(1).GetNumberOfPoints())
    
    for i in range(2):
        cells = centerline.GetCell(i)
        points = cells.GetPoints()
        N = points.GetNumberOfPoints() - 1
        end.append(points.GetPoint(N))

    return end


def viz(centerline, centerline_bif, points):
    """Help method during development to view the results"""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for i in range(3):
        if i == 2:
            i = 0
            centerline = centerline_bif
        point_ids = vtk.vtkIdList()
        centerline.GetCellPoints(i, point_ids)
        points0 = []
        for k in range(point_ids.GetNumberOfIds()):
            points0.append(centerline.GetPoint(point_ids.GetId(k)))
        arr = np.asarray(points0)
        x = arr[:,0]
        y = arr[:,1]
        z = arr[:,2]
        ax.plot(x, y, z, label=i)
        ax.legend()
        plt.hold("on")
    
    for p in points:
        ax.plot([p[0]], [p[1]], [p[2]],"o", label=3)
        ax.legend()
        plt.hold("on")

    plt.show()


def find_angle(data, u_vec, method, dirname):
    # Following the syntacs of Ingebrigtsen
    phi_0 = math.acos(np.dot(u_vec["bif"], u_vec[0])) * 180 / math.pi
    phi_1 = math.acos(np.dot(u_vec["bif"], u_vec[1])) * 180 / math.pi

    if data[0]["r_div"] > data[1]["r_div"]:
        phi1 = phi_0
        phi2 = phi_1
    else:
        phi1 = phi_1
        phi2 = phi_0

    text = "\nphi1_" + method + ": " + str(phi1) + "\n"
    text += "phi2_" + method + ": " + str(phi2)
    file_path = path.join(dirname, "manifest.txt")
    if not text_exists(text, file_path):
        f = open(path.join(dirname, "manifest.txt"), "a")
        f.write(text)
        f.close()


def text_exists(text, file_path):
    file = open(file_path, "r")
    t = file.read()
    file.close()
    return text in t


def curvature_stats(centerline, centerline_bif, data, dirname, tol):
    #curvature = centerline_geo.GetPointData().GetArray("Curvature")
    opt = data.keys()
    cases = [(opt[0], opt[1]), (opt[0], opt[2]), (opt[1], opt[2])]
    file_path = path.join(dirname, "manifest.txt")

    # List of points conected to ID
    points_ids = vtk.vtkIdList()

    for case in cases:
        stats = []
        #print case

        if "bif" not in case:
            name = "bif"
            start = data[case[0]]["ID_div"]
            end = data[case[1]]["ID_div"]
            curvature = centerline_bif.GetPointData().GetArray("Curvature")

            for i in range(min(start, end), max(start, end) + 1, 1):
                stats.append(curvature.GetTuple1(i))

        else:
            key = case[0] if case[0] != "bif" else case[1]
            key_comp = 0 if key == 1 else 1
            name = "1" if data[key]["r_end"] > data[key_comp]["r_end"] else "2"   
            #print key
            centerline.GetCellPoints(key, points_ids)
            #print centerline.GetNumberOfLines()
            curvature = centerline.GetPointData().GetArray("Curvature")
             
            end_point = data[key]["div_point"]
            i = data["bif"]["i_div"]
            point = (0,0,0)
            dist = 1e10
            dist_prev = 1e10
            points = []

            # Collect curvature from desired area
            while dist_prev >= dist:
                #print i
                #print dist, tol
                stats.append(curvature.GetTuple1(points_ids.GetId(i)))
                #points.append(centerline.GetPoint(points_ids.GetId(i)))
                point = centerline.GetPoint(points_ids.GetId(i))

                #dist_prev_2 = dist_prev_1
                dist_prev = dist
                dist = math.sqrt(np.sum(abs(np.asarray(end_point) - 
                                            np.asarray(point))))
                i += 1
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        arr = np.asarray(points)
        x = arr[:,0]
        y = arr[:,1]
        z = arr[:,2]
        ax.plot(x, y, z, label=i)
        plt.hold("on")
        ax.plot([end_point[0]], [end_point[1]], [end_point[2]], "o", label="point")
        p = data[key_comp]["div_point"]
        ax.plot([p[0]], [p[1]], [p[2]], "o", label="point_other")
        ax.legend()
        plt.show()

        sys.exit(0) 
        """

        # Get stats
        stats = np.asarray(stats)
        mean = np.mean(stats)
        max_ = np.max(stats)

        # Write to file
        text = "\ncurvature_max_%s: %s\ncurvature_mean_%s: %s" % \
                (name, max_, name, mean)
        #print text
        if not text_exists(text, file_path):
            f = open(path.join(dirname, "manifest.txt"), "a")
            f.write(text)
            f.close()
    

def angle_stats(data):
    # Create unit vectors
    for method in ["mean", "div", "bif_half"]:
        u_vec = {}

        # Define 'center' of bifurcation
        if method == "mean":
            d = (np.asarray(data[0]["div_point"]) + \
                 np.asarray(data[1]["div_point"]) + \
                 np.asarray(data["bif"]["div_point"])) / 3

        elif method == "bif_half":
            id = int(round(abs(data[0]["ID_div"] - data[1]["ID_div"]) / 2))
            id = min(data[0]["ID_div"] + id, data[1]["ID_div"] + id)
            d = np.asarray(centerline_bif.GetPoint(id))

        for key in data.keys():
            if method == "div":
                d = np.asarray(data[key]["div_point"])
            e = np.asarray(data[key]["end_point"])
            if key == "bif":
                vec = d - e
            else:
                vec = e - d
            length = math.sqrt(np.dot(e - d, e - d))
            u_vec[key] = (np.asarray(vec)/length)

        find_angle(data, u_vec, method, dirpath)


def main(dirpath):
    """Takes filepath to a unziped 'case' folder from ANURISK database"""
    
    # Naming convention on files
    model_path = path.join(dirpath, "surface", "model.vtp")
    centerline_path = path.join(dirpath, "surface", "model_usr_centerline.vtp")
    centerline_path_bif = path.join(dirpath, "surface", "model_usr_centerline_bif.vtp")
    centerline_path_geo = path.join(dirpath, "surface", "model_usr_centerline_geo.vtp")
    centerline_path_bif_geo = path.join(dirpath, "surface", "model_usr_centerline_bif_geo.vtp")

    # Make the user create the bifurcation centerline
    if not path.exists(centerline_path_bif):
        print "Pick the inflow as a source and then pick one outflow on each side of" + \
            " the ICA(?) bifraction.\nHit enter to continue."
        #raw_input()
        subprocess.check_output("vmtk vmtkcenterlines -ifile %s -ofile %s" \
                                 % (model_path, centerline_path),
                                stderr=subprocess.STDOUT,
                                shell=True)

        # The user defined senterline to see the correct bifucation
        centerline = ReadPolyData(centerline_path)

        # Get the senterline from the two endpoints
        points = get_endpoints(centerline)
        subprocess.check_output(("vmtk vmtkcenterlines -ifile %s -ofile %s" + \
                                 " -seedselector pointlist -sourcepoints %s %s %s -targetpoints " + \
                                 "%s %s %s") % (model_path, centerline_path_bif, points[0][0], 
                                                points[0][1], points[0][2], points[1][0],
                                                points[1][1], points[1][2]),
                                 stderr=subprocess.STDOUT,
                                 shell=True)
                    
    # If the bifurcation allready have been defined
    else:
        centerline = ReadPolyData(centerline_path)
    
    centerline_bif = ReadPolyData(centerline_path_bif)

    # Creat a tolerance for diverging
    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                  centerline.GetPoint(10), \
                                  centerline.GetPoint(11)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance

    # Find diverging points
    data = getData(centerline, centerline_bif, divergingTolerance)

    # Compute centerline properties
    if not path.exists(centerline_path_geo):
        subprocess.check_output(("vmtk vmtkcenterlinegeometry -ifile %s -smoothing 1 " + \
                                 "-iterations 300 -factor 0.1 -outputsmoothed 1 -ofile %s") \
                                 % (centerline_path, centerline_path_geo),
                                    stderr=subprocess.STDOUT, shell=True)

    if not path.exists(centerline_path_bif_geo):
        subprocess.check_output(("vmtk vmtkcenterlinegeometry -ifile %s -smoothing 1 " + \
                                 "-iterations 300 -factor 0.1 -outputsmoothed 1 -ofile %s") \
                                 % (centerline_path_bif, centerline_path_bif_geo),
                                    stderr=subprocess.STDOUT, shell=True)

    centerline_geo = ReadPolyData(centerline_path_geo)
    centerline_bif_geo = ReadPolyData(centerline_path_bif_geo)
    curvature_stats(centerline_geo, centerline_bif_geo, data, dirpath,
                    divergingTolerance)


def csv_to_txt(folder):
    """Make it easier to access data with a normal txt file"""
    csv = path.join(folder, "manifest.csv")
    txt = path.join(folder, "manifest.txt")
    reader = open(csv, "r")
    header = reader.readline().split(",")
    row = reader.readline().split(",")
    for i in range(len(header)):
        header[i] = ": ".join([header[i].replace("\n",""),
                               row[i].replace("\n", "")])
    text = "\n".join(header)
    reader.close()
    writer = open(csv, "w")
    writer.write(text)
    writer.close()
    subprocess.check_output("mv " + csv + " " + txt, stderr=subprocess.STDOUT, shell=True)
    

if __name__ == "__main__":
    for folder in listdir("."):
        if path.isdir(folder) and folder != "backup":
            print "Looking at case", folder
            main(folder)
    #main("C0039")
