from subprocess import *
from os import path, listdir

dirpath = "/home/aslak/master/src/aneurysms/"
a = 0
cases = listdir(dirpath)
cases = [c for c in cases if c.startswith("P0") and path.isdir(path.join(dirpath, c))]

for case in cases:
    for b in [True, False]:
        for l in [True, False]:
            bif = " --bif True" if b else ""
            l = " --lower True" if l else ""
            cmd = ("python move_branches.py --dir_path %s --smooth True%s" + \
                   " --case %s --angle %s%s --addPoint True") % \
                    (dirpath, bif, case, a, l)

            print cmd
            output = check_output(cmd, stderr=STDOUT, shell=True)
            print output

            if (not l) and (not b):
                cmd = cmd.replace(" --addPoint True", "")
                print cmd
                output = check_output(cmd, stderr=STDOUT, shell=True)
                print output
