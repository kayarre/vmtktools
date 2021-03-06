from common import *
from subprocess import check_output, STDOUT
from paralleltransportvoronoidiagram import create_new_surface
from clipvoronoidiagram import SmoothClippedVoronoiDiagram
from os import listdir, path
from argparse import ArgumentParser


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--dir_path', '--d', type=str, default=".",
                        help="Path to the folder with all the cases",
                        metavar="PATH")

    args = parser.parse_args()

    return args.dir_path


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
    check_output("mv " + csv + " " + txt, stderr=STDOUT, shell=True)


def start_case(casepath, folder):
    if folder.endswith(".tar.gz"):
        folder_ = path.join(casepath, folder[:-14])
    else:
        folder_ = path.join(casepath, folder)

    model_path = path.join(folder_, "surface", "model.vtp")
    model_smoothed_path = path.join(folder_, "surface", "model_smoothed.vtp")
    voronoi_path = path.join(folder_, "surface", "voronoi.vtp")
    voronoi_smoothed_path = path.join(folder_, "surface", "voronoi_smoothed.vtp")
    centerlines_path = path.join(folder_, "surface", "centerline_complete.vtp")
   
    # Untar
    if not path.isdir(folder_):
        check_output("bash untar.sh %s %s" % (casepath, folder), stderr=STDOUT, shell=True)
        check_output("cd -", stderr=STDOUT, shell=True)

    # Convert case info to plain text
    if not path.exists(path.join(folder_, "manifest.txt")):
        csv_to_txt(folder_)

    # Compute centerlines and voronoi diagram
    if not path.exists(centerlines_path):
        centerlines = makeCenterline(model_path, centerlines_path, smooth=False,
                                    resampling=False)
    else:
        centerlines = ReadPolyData(centerlines_path)

    if not path.exists(voronoi_path):
        voronoi = makeVoronoi(model_path, voronoi_path)
    else:
        voronoi = ReadPolyData(voronoi_path)

    if not path.exists(voronoi_smoothed_path):
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerlines, 0.25)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)
    else:
        voronoi_smoothed = ReadPolyData(voronoi_smoothed_path)

    if not path.exists(model_smoothed_path):
        model_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(model_smoothed, model_smoothed_path)


if __name__ == "__main__":
    dir = read_command_line()
    cases = ["P0086"]
    #for folder in listdir(dir):
    for folder in cases:
        print folder
        start_case(dir, folder)
