import sys
import os

if len(sys.argv) == 3:
    g = str(sys.argv[1])  # Team number 1-30(ish)
    c = str(sys.argv[2])  # Case number 1-5
else:
    print "Usage %s [team number] [case number]" % sys.argv[0]
    sys.exit(0)

# Move files to folder
os.system("cp *.py " + g)

# Change directory
os.chdir(g)
os.system("echo Now in $PWD")

# Naming convention
case_path = path.join(c, "Case" + c)
domain_path = case_path + "domain.vtp"
domain_name = path.basename(domain_path)

# Convert stl file
os.system("vmtksurfacereader -ifile " + case_path + "domain.stl -ofile " + domain_path) 

# Make voronoi diagram
os.system("vmtkdelaunayvoronoi -ifile "+c+"/Case"+c+"domain.vtp -removesubresolution 1 -voronoidiagramfile "+c+"/Case"+c+"domain_voronoi.vtp")

# Remove aneurysm and reconstruct the deometry
os.system("python make_centerlines.py --domain_path " + domain_path)

os.system("python patchandinterpolatecenterlines.py "+c+" Case"+c+"domain terminal")
os.system("python clipvoronoidiagram.py "+c+"/ Case"+c+"domain")
os.system("python paralleltransportvoronoidiagram.py "+c+"/Case"+c+"domain 1 ")

# Centerlines and voronoi for the new geometry
os.system("python make_centerlines.py --domain_path " + case_path + "domain_reconstructedmodel.vtp --n True")
os.system("vmtkdelaunayvoronoi -ifile "+c+"/Case"+c+"domain_reconstructedmodel.vtp -removesubresolution 1 -voronoidiagramfile "+c+"/Case"+c+"domain_reconstructedmodel_vor.vtp")

# Extract parameters from geometry
os.system("python extractaneurysmneckplanesection.py %s %s terminal" % (c, domain_name))
os.system("python computeaneurysmsaccenterline.py %s %s domain" % (c, domain_name))
os.system("python computesurfacesandvolumes.py %s %s domain" % (c, domain_name))
os.system("python extractparentvesselparameters.py %s %s terminal" % (c, domain_name))
os.system("rm *.py ")
