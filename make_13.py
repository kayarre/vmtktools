from subprocess import *
from os import path, listdir
from numpy import linspace
from glob import glob

dirpath = "/home/aslak/master/src/aneurysms/"
a = 0
cases = glob(path.join(dirpath, "P0*")
cases.sort()

print cases
for case in cases:
    #if case == cases[0]:
    #    continue

    # Rotate up
    print "="*50
    print "     Rotate up", case
    print "="*50

    cmd = ("python move_branches.py --dir_path %s --smooth True" + \
           " --case %s --lower %s --bif %s --addPoint True --a 10 --w True") \
           % (dirpath, case, True, True)
    print cmd
    
    output = check_output(cmd, stderr=STDOUT, shell=True)
    print output

    # Rotate down
    print "="*50
    print "    Rotate down", case
    print "="*50

    cmd = ("python move_branches.py --dir_path %s --smooth True" + \
           " --case %s --lower %s --bif %s --addPoint True --a -10 --w True") \
           % (dirpath, case, True, True)
    print cmd
    
    output = check_output(cmd, stderr=STDOUT, shell=True)
    print output
