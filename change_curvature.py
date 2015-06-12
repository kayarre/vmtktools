from common import *
from argparse import ArgumentParser
import math


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram (surface) should be" + \
                             "smoothed before it is manipulated", metavar="smooth")
    parser.add_argument('--D2', type=bool, default=True,
                        help="Find the best fit plan of the part that should be" + \
                        "smoothed and change the distance to that particular" + \
                        "part", metavar="best_fit_plane")
    parser.add_argument('--global_smooth', type=bool, default=False,
                        help="Takes the entire siphon and do a global
                        "smoothing, this does not preserve torsion and the local" + \
                        "geometry.", metavar="global_smooth")
    parser.add_argument('--D2_factor', type=float, default=0.75)
    parser.add_argument('--global_smooth_factor', type=float, default=0.5)
    parser.add_argument('global_smooth_iterations', type=int, default=100)
    parser.add_argument('--rotate', type=bool, default=False)
    parser.add_argument('--rotate_angle', type=int, default=8)
    parser.add_argument('--D3', type=bool, defualt=False)
    parser.add_argument('--D3_factor', type=float, defualt=0.75)

    args = parser.parse_args()
    ang_ = math.pi / args.rotate_angle

    return args.s, args.D2, args.global_smooth, args.D2_factor,
            args.global_smooth_factor, args.global_smooth_iterations, args.rotate,
            args.rotate_angle, args.D3, args.D3_factor


def main(folder, smooth, D2, global_smooth, D2_factor, global_smooth_factor,
        global_smooth_iterations, rotate_angle, D3, D3_factor):
    # Input filenames
    model_path = path.join(dirpath, "surface", "model.vtp")

    # Outpate filenames
    centerline_path = path.join(dirpath, "surface", "centerline_for_curv.vtp")
    centerline_curv_path = path.join(dirpath, "surface", "centerline_curv_change.vtp")
    centerline_changed_path = path.join(dirpath, "surface", "centerline_curv_new.vtp")
    if smooth:
        voronoi_path = path.join(dirpath, "surface", "voronoi.vtp")
    else:
        voronoi_path = path.join(dirpath, "surface", "voronoi_smoothed.vtp")
    voronoi_curvature_path = path.join(dirpath, "surface", "voronoi_curvature.vtp")
    model_curvature_path = path.join(dirpath, "surface", "model_curvature.vtp")

    # TODO: Start, stop.
    centerline = makeCenterline(model_path, centerline_path, length=0.1, smooth=False)
    centerline 
    centerline_curve = ExtractSingleLine(
    




if __name__ == "__main__":
    smooth, D2, global_smooth, D2_factor, global_smooth_factor,
    global_smooth_iterations, rotate, rotate_angle, D3, D3_factor = read_command_line()
    basedir = "."
    for folder in listdir(basedir):
        if folder[:2] in ["P0", "C0"]:
            main(path.join(basedir, folder), smooth, D2, global_smooth,
                 D2_factor, global_smooth_factor, global_smooth_iterations,
                 rotate, rotate_angle, D3, D3_factor)

