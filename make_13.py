from subprocess import *

path = "../master/src/aneurysms/"
case = "C0001"
for a in [15, 10, 7, 5, 3]:
    for b in [True, False]:
        for l in [True, False]:
            bif = "--bif True" if b else ""
            cmd = ("python move_branches.py --dir_path %s --smooth True %s" + \
                   " --case %s --angle %s --lower %s --addPoint True") % \
                    (path, bif, case, a, l)
            print cmd
            output = check_output(cmd, stderr=STDOUT, shell=True)
            print output


