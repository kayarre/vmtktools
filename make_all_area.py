from subprocess import *
from os import path, listdir, system
from numpy import linspace
from common import *

dirpath = "/home/aslak/master/src/aneurysms/"
a = 0
cases = listdir(dirpath)

cases = [c for c in cases if c.startswith("P0")]
cases.sort()

base_command = "python area_variations.py --dir_path ../master/src/aneurysms" + \
               " --smooth True --ratio %s --case %s"

for case in cases:
    # Get area ratio
    data = getParameters(path.join(dirpath, case))
    ratio = data["max_min_ratio_area"]
    sd = 0.92
    leg = ["-SD", "+SD"] #, "+2SD"]

    for j, i in enumerate([-1, 1]):
        ratio_tmp = ratio + i*sd
        print "="*50, "\n", case, leg[j], "\n", "="*50
        cmd = base_command % (ratio_tmp, case)
        print cmd
        system(cmd)
        print ""

    print ""
