import vtk
import math
import subprocess
from common import *
import numpy as np
from os import path, listdir


def getData(centerline, centerline_bif, tol):
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
        if distance_between_points > tol:
            tmpI = i
            point_ID_0 = points_ids_0.GetId(i)
            point_ID_1 = points_ids_1.GetId(i)
            center = centerline.GetPoint(point_ID_0)
            r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID_0)
            break

    end, r_end = move_past_sphere(centerline, center, r, point_ID_0, stop=point_ID_0*100, step=1)
    data["bif"]["end_point"] = end
    data["bif"]["r_end"] = r_end
    data["bif"]["div_point"] = center
    data["bif"]["ID_div"] = point_ID_0
    data["bif"]["i_div"] = tmpI
    data["bif"]["r_div"] = r

    # Find the diverging point for anterior and midt bifurcation
    # continue further downstream in each direction and stop when
    # a point is closer than tol, than move point MISR * X
    locator = get_locator(centerline_bif)

    counter = 0
    for point_ids in [points_ids_0, points_ids_1]:
        for i in range(tmpI, point_ids.GetNumberOfIds(), 1):
            tmp_point = centerline.GetPoint(point_ids.GetId(i))
            closest_point_ID = locator.FindClosestPoint(tmp_point)
            closest_point = centerline_bif.GetPoint(closest_point_ID)
            distance_between_points = distance(tmp_point, closest_point)
            if distance_between_points < tol:
                point_ID = point_ids.GetId(i)
                center = centerline.GetPoint(point_ID)
                r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID)
                break
        
        end, r_end = move_past_sphere(centerline, center, r, point_ID, X=1.5)
        data[counter]["end_point"] = end
        data[counter]["r_end"] = r_end
        data[counter]["r_div"] = r
        data[counter]["ID_end"] = locator.FindClosestPoint(data[counter]["end_point"])
        data[counter]["ID_div"] = locator.FindClosestPoint(center)
        data[counter]["div_point"] = center
        
        counter += 1
        
    return data


def move_past_sphere(centerline, center, r, start, step=-1, stop=0, X=0.8):
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

    new_data = {"phi1_%s" % method: phi1, "phi2_%s" % method: phi2}
    writeParameters(dirname)


def curvature_stats(centerline, centerline_bif, data, dirname, tol):
    opt = data.keys()
    cases = [(opt[0], opt[1]), (opt[0], opt[2]), (opt[1], opt[2])]

    curvature_bif = get_array("Curvature", centerline_bif)
    curvature = get_array("Curvature", centerline)

    for case in cases:
        if "bif" not in case:
            name = "bif"
            start = data[case[0]]["ID_end"]
            end = data[case[1]]["ID_end"]
            stats = curvature_bif[start:end+1]

        else:
            key = case[0] if case[0] != "bif" else case[1]
            key_comp = 0 if key == 1 else 1
            name = "1" if data[key]["r_end"] > data[key_comp]["r_end"] else "2"   
            
            locator = get_locator(centerline)
            end = locator.FindClosestPoint(data[key]["end_point"])
            stats = curvature[data["bif"]["ID_end"]:end+1]

        # Get stats
        mean = np.mean(stats)
        max_ = np.max(stats)

        # Write to file
        new_data = {"curvature_max_%s" % name: max_, "curvature_mean_%s" % name: mean}
        writeParameters(new_data, dirname)
    

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
    centerline_path = path.join(dirpath, "surface", "centerline_for_angles.vtp")
    centerline_path_bif = path.join(dirpath, "surface", "centerline_bifurcation.vtp")
    centerline_path_geo = path.join(dirpath, "surface", "centerline_geo.vtp")
    centerline_path_bif_geo = path.join(dirpath, "surface", "centerline_geo_bif.vtp")

    # Make centerlines, will only read if they allready exists
    centerline = makeCenterline(model_path, centerline_path, smooth=False, resampling=True)
    centerline_bif = makeCenterline(model_path, cenerline_path_bif, length=0.1
                                    smooth=False, in_out=[0,1])
                    
    # Tolerance for diverging
    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                  centerline.GetPoint(10), \
                                  centerline.GetPoint(11)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance

    # Diverging points
    data = getData(centerline, centerline_bif, divergingTolerance)

    # TODO: Extract lines and spline

    # Compute centerline properties
    centerline_geo = CenterlineAttribiute(centerline, remove=False,
                                          filename=centerline_path_geo,
                                          smooth=True)
    centerline_bif_geo = CenterlineAttribiutecenterline_path_bif_geo(centerline_bif, 
                                          remove=False, filename=centerline_path_geo,
                                          smooth=True)

    curvature_stats(centerline_geo, centerline_bif_geo, data, dirpath,
                    divergingTolerance)
    

if __name__ == "__main__":
    for folder in listdir("."):
        if path.isdir(folder) and folder != "backup":
            print "Looking at case", folder
            main(folder)
    #main("C0039")
