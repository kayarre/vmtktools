import vtk
import math
from argparse import ArgumentParser
import subprocess
from common import *
import numpy as np
from os import path, listdir
from landmark import splineCenterline


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--dir_path', '--d', type=str, default=".",
                        help="Path to the folder with all the cases",
                        metavar="PATH")

    args = parser.parse_args()

    return args.dir_path


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


def curvature_stats(centerline, centerline_bif, data, folder):
    key = "end_point"
    for i in range(3):
        if i != 2:
            startPoint = data["bif"][key]
            endPoint = data[i][key]
            tmp = ExtractSingleLine(centerline, i)

            # Find largest and smallest daugther branch
            other = 0 if i == 1 else 0
            key_line = "1" if data[i]["r_end"] > data[other]["r_end"] else "2"

        else:
            startPoint = data[0][key]
            endPoint = data[1][key]
            tmp = ExtractSingleLine(centerline_bif, 0)
            key_line = "bif"

        locator = get_locator(tmp)
        startID_ = locator.FindClosestPoint(startPoint)
        endID_ = locator.FindClosestPoint(endPoint)
        startID = min(startID_, endID_)
        endID = max(startID_, endID_)
        
        line = ExtractSingleLine(tmp, 0, startID=startID, endID=endID)
        line = splineCenterline(line, nknots=4)
        WritePolyData(line[0], path.join(folder, "surface", "spline_%s.vtp" % key_line))

        curvature = get_array("Curvature", line[0])
        mean = np.mean(curvature)
        max_ = np.max(curvature)

        new = {"curvature_max_new_%s" % key_line: max_, "curvature_mean_new_%s" % key_line: mean}
        writeParameters(new, folder)


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
    centerline = makeCenterline(model_path, centerline_path, smooth=False, resampling=False)
    centerline_bif = makeCenterline(model_path, centerline_path_bif, resampling=False,
                                    smooth=False, in_out=[0,1])
                    
    # Tolerance for diverging
    centerlineSpacing = math.sqrt(distance(centerline.GetPoint(10), \
                                           centerline.GetPoint(11)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance

    # Diverging points
    data = getData(centerline, centerline_bif, divergingTolerance)
    angle_stats(data)
    curvature_stats(centerline, centerline_bif, data, dirpath)
    

if __name__ == "__main__":
    basedir = read_command_line()
    for folder in listdir(basedir):
        folder_path = path.join(basedir, folder)
        if path.isdir(folder_path) and "P0" in folder or "C0" in folder:
            print "Looking at case", folder
            main(folder_path)
    #main("C0001")
