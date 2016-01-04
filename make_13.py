from subprocess import *
from os import path, listdir
from numpy import linspace

dirpath = "/home/aslak/master/src/aneurysms/"
a = 0
cases = listdir(dirpath)
#cases = {"N0120": [0.2625, -0.4358] , "N0122": [0.1743, -0.2554], "P0121": [0.2470, -0.4036]}
#cases = {"N0122": [0.2652, -0.4315]}
#cases = ["N01%s" % str(i) for i in range(12, 18, 1)]

cases = [c for c in cases if c.startswith("P0")]
cases.sort()
print cases
for case in cases:
    if case == cases[0]:
        continue
    print "="*50
    print "     Rotate up", case
    print "="*50

    cmd = ("python move_branches.py --dir_path %s --smooth True" + \
           " --case %s --lower %s --bif %s --addPoint True --a 10 --w True") \
           % (dirpath, case, True, True)
    print cmd
    output = check_output(cmd, stderr=STDOUT, shell=True)
    print output

    print "="*50
    print "    Rotate down", case
    print "="*50

    cmd = ("python move_branches.py --dir_path %s --smooth True" + \
           " --case %s --lower %s --bif %s --addPoint True --a -10 --w True") \
           % (dirpath, case, True, True)
    print cmd
    output = check_output(cmd, stderr=STDOUT, shell=True)
    print output
