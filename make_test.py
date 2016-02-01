from subprocess import *
from os import path, listdir
from numpy import linspace
import sys

dirpath = "/home/aslak/master/src/aneurysms/"
a = 0
cases = listdir(dirpath)
#cases = {"N0120": [0.2625, -0.4358] , "N0122": [0.1743, -0.2554], "P0121": [0.2470, -0.4036]}
#cases = {"N0122": [0.2652, -0.4315]}
#cases = ["N01%s" % str(i) for i in range(12, 18, 1)]
cases = [c for c in cases if c.startswith("P0")]
cases.sort()
print cases
cases = ["P0220"]
for case in cases:
    for bif in [True, False]:
        for lower in [True, False]:
            print "="*50
            print "      Running:   Bif =", bif, "lower = ", lower
            print "="*50

            
            cmd = ("python move_branches.py --dir_path %s --smooth True" + \
                " --case %s --addPoint True --a 0") % (dirpath, case)
            
            if bif:
                cmd += " --bif True"
            if lower:
                cmd += " --lower True"

            print cmd
            output = check_output(cmd, stderr=STDOUT, shell=True)
            print output
            sys.exit(0)
            #print "="*50
            #print "    Old version"
            #print "="*50

            #cmd = ("python move_branches.py --dir_path %s --smooth True" + \
            #    " --case %s --addPoint True --a 0") % (dirpath, case)
            #print cmd
            #output = check_output(cmd, stderr=STDOUT, shell=True)
            #print output
