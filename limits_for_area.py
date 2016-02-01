from stats import find_SD, get_info
from os import path, listdir
from argparse import ArgumentParser


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    args = parser.parse_args()

    return args.d


RED   = "\033[1;37;31m%s\033[0m"
BLUE  = "\033[1;37;34m%s\033[0m"
GREEN = "\033[1;37;32m%s\033[0m"

def info_blue(s, check=True):
    print BLUE % s

def info_green(s, check=True):
    print GREEN % s

def info_red(s, check=True):
    print RED % s

def main(dir_path):
    var = []
    fd = []
    #key = "max_mean_ratio_area"
    #key = "max_circleness"
    #key = "global_max_area"
    #key = "length_min_max"
    #key = "min_circleness"
    #key = "mean_area"
    #key = "max_derivative"
    #key = "min_mean_ratio_area"
    key = "max_min_ratio_area"
    #key = "global_min_area"
    #key = "local_max_stepest_disent"
    #key = "mean_circleness"
    folders = listdir(dir_path)
    folders.sort()
    for folder in folders:
        if folder.startswith("C0"):
            manifest_path = path.join(dir_path, folder, "manifest.txt")
            info = get_info(manifest_path)
            if info["aneurysmType"] == "TER" and info["aneurysmLocation"] == "ICA" and info.has_key(key):
                var.append(float(info[key]))
                fd.append(folder)
    
    mean, SD = find_SD(var)
    upper = mean + SD
    lower = mean - SD
    for i in range(len(fd)):
        if lower <= var[i] <= upper:
            info_green(fd[i] + ": " + str(var[i]))
        else:
            info_red(fd[i] + ": " + str(var[i]))

    print "Mean: ", mean, "SD:", SD


if __name__ == "__main__":
    dir_path = read_command_line()
    main(dir_path)
