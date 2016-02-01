from os import path, listdir, makedirs, system
import sys

if len(sys.argv) != 2:
    print "Usage %s <model_to_test>" % sys.argv[0]
    sys.exit(0)

dirpath = "../master/src/aneurysms/"

model_path= sys.argv[1]
centerline_path = path.join(path.dirname(model_path), "centerline_complete.vtp")

newfolder = listdir(dirpath)
newfolder = [path.join(dirpath, f) for f in newfolder if path.isdir(path.join(dirpath, f))]
newfolder = [f for f in newfolder if "N0" in f]
newfolder.sort()
number = int(newfolder[-1][-3:]) + 1
newfolder = newfolder[-1][:-3] + str(number)
makedirs(path.join(newfolder, "surface"))

system("touch %s" % path.join(newfolder, "manifest.txt"))
system("cp %s %s" % (model_path, path.join(newfolder, "surface", "model.vtp")))
system("cp %s %s" % (centerline_path, path.join(newfolder, "surface")))
