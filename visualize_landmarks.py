from common import *
from os import path, listdir
from subprocess import *

def dir(f, m=0):
    s = "surface" if m == 0 else "morphology"
    return path.join(dirpath, s, f)

mainfolder = "/home/aslak/master/src/aneurysms"
for folder in listdir(mainfolder):
    dirpath = path.join(mainfolder, folder)
    if path.isdir(dirpath) and folder != "backup" and "test" not in folder and folder != ".git":
        print "Dooing case", folder
        if not path.exists(dir("model_split.vtp")):
            print "Extracting the siphon"
            cmd = ("vmtkbranchextractor -ifile %s -radiusarray" + \
                    " MaximumInscribedSphereRadius --pipe vmtkbranchclipper " + \
                    "-radiusarray MaximumInscribedSphereRadius -ifile %s" + \
                    " -ofile %s -groupids 0") % \
                    (dir("model_usr_centerline.vtp"), dir("model.vtp"), dir("model_split.vtp"))

            a = check_output(cmd, stderr=STDOUT, shell=True)

        line = ReadPolyData(dir("centerlines.vtp", 1))

        locator = get_locator(line)

        parameters = getParameters(dirpath)
        key = ["inf_end", "sup_ant", "ant_post", "post_inf"]
        key_name = ["lower", "inferior", "posterior", "anterior", "superior"]

        points = [parameters[k] for k in key]
        new_points = []
        for point in points:
            new_points.append(line.GetPoints().GetPoint(locator.FindClosestPoint(point)))
        
        for i in range(4):
            ifile = dir("model_split.vtp")
            ofile = dir("landmark_" + key_name[i] + ".vtp")
            cmd = ("vmtkpointsplitextractor -ifile %s -splitpoint %s %s %s -tolerance 1E-4" + \
                    " -gaplength 1 -radiusarray MaximumInscribedSphereRadius" + \
                    " --pipe vmtkbranchclipper -ifile %s -groupids 0" + \
                    " -ofile %s -radiusarray MaximumInscribedSphereRadius") % \
                    ((dir("centerlines.vtp", 1),) + new_points[i] + (ifile, ofile))
            
            a = check_output(cmd, stderr=STDOUT, shell=True)
