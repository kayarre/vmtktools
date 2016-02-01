from os import listdir, path, system
from common import *

base = "/home/aslak/master/src/aneurysms"
folders = listdir(base)
folders = [path.join(base, f, "surface") for f in folders if f.startswith("P0")]
j = 1

for folder in folders:
    # Get three models
    files = listdir(folder)
    files = [f for f in files if f.startswith("model_area") and len(f) > 33]
    #print "Length of file:", len(files)

    for file in files:
        # Destination folder
        if j > 9:
            i = "00%s" % j
        elif j < 10:
            i = "000%s" % j
        elif j > 99:
            i = "0%s" % j

        folder_name = path.join(base, "A%s" % i, "surface")
        file_path = path.join(folder, file)
        cl_path = path.join(folder, "centerline_complete.vtp")
        system("cp %s %s" % (file_path, path.join(folder_name, "model.vtp")))
        system("cp %s %s" % (cl_path, folder_name))
        new_text = open(path.join(base, "A%s" % i, "manifest.txt"), "w")
        new_text.write("origin: %s" % file_path)

        # Get inlet and outlet from centerline
        cl = ReadPolyData(cl_path)
        numberOfLines = cl.GetNumberOfCells()
        new_text.write("\ninlet: %s" % str(ExtractSingleLine(cl, 0).GetPoints().GetPoint(0)))
        for k in range(numberOfLines):
            tmp_cl =  ExtractSingleLine(cl, k)
            N = tmp_cl.GetNumberOfPoints()
            new_text.write("\noutlet%s: %s" % (k, str(tmp_cl.GetPoints().GetPoint(N - 1))))
        new_text.close()
        j += 1


