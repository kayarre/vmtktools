from common import *
from patchandinterpolatecenterlines import *
from vmtkClassifyBifucation import getData
import numpy as np
from os import path, listdir
from scipy.interpolate import splrep, splev
from scipy.signal import argrelextrema, gaussian
import math
import subprocess


def CenterlineAttribiute(line):
    tmpFileName = "tmp_cl.vtp"
    WritePolyData(line, "tmp_cl.vtp")
    command = ('vmtkcenterlineattributes -ifile %s --pipe vmtkcenterlinegeometry ' + \
    '-ofile %s -smoothing 0') % (tmpFileName, tmpFileName)
    subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
    line = ReadPolyData(tmpFileName)
    subprocess.check_output('rm ' + tmpFileName, stderr=subprocess.STDOUT, shell=True)
    return line


def splineCenterline(line):
    # Allocate data structure to store centerline points
    data = np.zeros((line.GetNumberOfPoints(), 3))

    # Collect data from centerline
    for i in range(data.shape[0]):
        curv_coor = get_curvilinear_coordinate(line)
        data[i,:] = line.GetPoints().GetPoint(i)

    t = np.linspace(curv_coor[0], curv_coor[-1], 25)[1:-1]
    fx = splrep(curv_coor, data[:,0], k=4, t=t)
    fy = splrep(curv_coor, data[:,1], k=4, t=t)
    fz = splrep(curv_coor, data[:,2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    data = np.zeros((len(curv_coor), 3))
    data[:,0] = fx_
    data[:,1] = fy_
    data[:,2] = fz_

    header = ["X", "Y", "Z"]
    line = data_to_vtkPolyData(data, header)

    # Let vmtk compute curve attributes
    line = CenterlineAttribiute(line)

    # Compute curvature from the 'exact' spline to get a robust way of
    # finding max / min points on the centerline
    dlsfx = splev(curv_coor, fx, der=1)
    dlsfy = splev(curv_coor, fy, der=1)
    dlsfz = splev(curv_coor, fz, der=1)

    ddlsfx = splev(curv_coor, fx, der=2)
    ddlsfy = splev(curv_coor, fy, der=2)
    ddlsfz = splev(curv_coor, fz, der=2)

    C1xC2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
    C1xC2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
    C1xC2_3 = ddlsfy*dlsfx - ddlsfx * dlsfy

    curvature_ = np.sqrt(C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2) / \
                        (dlsfx**2 + dlsfy**2 + dlsfz**2)**1.5

    max_point_ids = list(argrelextrema(curvature_, np.greater)[0])
    min_point_ids = list(argrelextrema(curvature_, np.less)[0])

    # TODO: Replace the locator with curv_coor = length
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(line)
    locator.BuildLocator()

    min_points = [[fx_[i], fy_[i], fz_[i]] for i in min_point_ids]
    max_points = [[fx_[i], fy_[i], fz_[i]] for i in max_point_ids]
    min_point_ids = []
    max_point_ids = []

    for point_min, point_max in zip(min_points, max_points):
        min_point_ids.append(locator.FindClosestPoint(point_min))
        max_point_ids.append(locator.FindClosestPoint(point_max))

    # The ParallelTransportNormals and the FrenetTangent is not orthonormal
    # (but close) from vmtk. Using the GramSchmidt proses gives E2 and fixes
    # the non-orthogonality
    E1 = get_array("ParallelTransportNormals", line, k=3)
    T = get_array("FrenetTangent", line, k=3)
    E2 = np.zeros((E1.shape[0], 3))

    for i in range(E1.shape[0]):
        V = np.eye(3)
        V[:, 0] = T[i,:]
        V[:, 1] = E1[i,:]
        V = GramSchmidt(V)

        E1[i,:] = V[:,1]
        E2[i,:] = V[:,2]

    # Compute k_1, k_2 furfilling T' = curv(s)*N(s) = k_1(s)*E_1(s) + k_2(s)*E_2(s).
    # This is simply a change of basis for the curvature vector N. The room
    # of k_1 and k_2 can be used to express both the curvature and the
    # torsion.
    N = get_array("FrenetNormal", line, k=3)
    curvature = get_array("Curvature", line)

    k2 = (curvature.T * (E1[:,1]*N[:,0] - N[:,1]*E1[:,0]) / \
                        (E2[:,1]*E1[:,0] - E2[:,0]*E1[:,1]))[0]
    k1 = (-(curvature.T * N[:,0] + k2*E2[:,0]) / E1[:,0])[0]

    for k in [(k1, "k1"), (k2, "k2")]:
        k_array = create_vtk_array(k[0], k[1])
        line.GetPointData().AddArray(k_array)

    # Compute a improved torsion that does not include all the noice from
    # vmtk. Here the analytical representation of the spline is accessebole.
    length = get_curvilinear_coordinate(line)
    dddlsfx = splev(length, fx, der=3)
    dddlsfy = splev(length, fy, der=3)
    dddlsfz = splev(length, fz, der=3)

    torsion_spline = (dddlsfx*C1xC2_1 + dddlsfy*C1xC2_2 + dddlsfz * C1xC2_3) / \
                        (C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2)
    torsion_array = create_vtk_array(torsion_spline, "Torsion")
    line.GetPointData().AddArray(torsion_array)

    return line, max_point_ids, min_point_ids


def landmarking(centerline, folder):
    """Takes the file path to a ceterline patch created by
    patchandinterpolatecenterlines.CreateParentArteryPatches, and spline the
    centerline and uses Bogunevic et al. (2012) to do automated landmarking"""

    line, max_point_ids, min_point_ids = splineCenterline(centerline)
    WritePolyData(line, path.join(folder, "surface", "centerline_splined.vtp"))
    curvature = get_array("Curvature", line)
    length = get_curvilinear_coordinate(line)
    k1 = get_array("k1", line)
    k2 = get_array("k2", line)

    # Remove a min / max point that is in reality a saddle point
    for i in min_point_ids:
        for j in max_point_ids:
            if abs(i-j) < 5 and abs(curvature[i] - curvature[j]) < 0.01:
                min_point_ids.remove(i)
                max_point_ids.remove(j)

    k1_points = k1[max_point_ids]
    k2_points = k2[max_point_ids]
    k_points = np.zeros((k1_points.shape[0], 2))
    k_points[:,0] = k1_points[:,0]
    k_points[:,1] = k2_points[:,0]
    tetha = np.zeros(k1_points.shape[0] - 1)
    for i in range(tetha.shape[0]):
        a = k_points[i,:] / np.sqrt(np.sum(k_points[i,:]**2))
        b = k_points[i+1,:] / np.sqrt(np.sum(k_points[i+1,:]**2))
        tetha[i] = math.acos(np.dot(a, b))
        tetha[i] = tetha[i] * 180 / math.pi
    
    Z = np.zeros(length.shape[0])
    Y = np.zeros(length.shape[0])
    X = np.zeros(length.shape[0])
    for i in range(Z.shape[0]):
        X[i] = line.GetPoints().GetPoint(i)[0]
        Y[i] = line.GetPoints().GetPoint(i)[1]
        Z[i] = line.GetPoints().GetPoint(i)[2]


    from matplotlib.pyplot import plot, hold, show, legend
    plot(length, X, label="X")
    hold("on")
    plot(length, Y, label="Y")
    plot(length, Z, label="Z")
    legend()

    # Tolerance parameters from Bogunevic et al. (2012)
    tol_ant_post = 60
    tol_sup_ant = 45
    tol_post_inf = 45
    tol_inf_end = 110

    # Find max coronal coordinate
    if Z.min() in Z[argrelextrema(Z, np.less)[0]]:
        value_index = Z[argrelextrema(Z, np.less)[0]].min()
        max_coronal_ids = np.array(Z.tolist().index(value_index))
    else:
        value_index = Z[argrelextrema(Z, np.greater)[0]].max()
        max_coronal_ids = np.array(Z.tolist().index(value_index))

    plot(length[max_coronal_ids], Z[max_coronal_ids], "g^")
    show()

    m = max_coronal_ids
    viz(line, [[X[m], Y[m], Z[m]]])
    plot(length, curvature, [length[m], length[m]], [0,1])
    hold("on")
    plot(length[max_point_ids], curvature[max_point_ids], 'g^')
    show()

    sys.exit(0)
    # Find all interfaces
    def find_interface(start, dir, tol, part):
        stop = dir if dir == -1 else tetha.shape[0]
        sucess = False
        for i in range(start-1, stop, dir):
            if tetha[i] > tol:
                sucess = True
                break

        if sucess:
            start = max_point_ids[i]
            stop = max_point_ids[i + 1]
            index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
            min_point = min_point_ids[index]
            interfaces[part] = min_point

        elif not sucess and part == "sup_ant":
            print "Where not able to identify the interface between the" + \
                  "anterior and superior bend. Chekc the coronal coordinates"
            return None

        elif not sucess and part != "inf_end":
            print "The geometry is to short to be classified with superior" + \
                    ", anterior, posterior and inferior."
            return None

        elif not sucess and part == "inf_end":
            interfaces["inf_end"] = 0
            i = 0
            print "End of inferior is at the end of the geometry, this might" + \
                  "affect the geometry stats"
        else:
            print "Something happend, idea: some bend ended at the last point"
            return None

        return i

    interfaces = {}
    min_point_ids = np.array(min_point_ids)
    index = np.array((max_coronal_ids > max_point_ids).nonzero()[0]).max()
    start = find_interface(index, -1, tol_ant_post, "ant_post")
    if start is None:
        return None
    start = find_interface(start, -1, tol_post_inf, "post_inf")
    if start is None:
        return None
    start = find_interface(start, -1, tol_inf_end, "inf_end")
    start = find_interface(index + 1, 1, tol_sup_ant, "sup_ant")
    if start is None:
        return None

    # Find a "center" of each bend
    bends = ["inferior", "posterior", "anterior", "superior"]
    values = [interfaces["inf_end"], interfaces["post_inf"],
              interfaces["ant_post"], interfaces["sup_ant"], 
              curvature.shape[0]]

    max_landmarks = {}
    for i in range(len(bends)):
        curv_part = curvature[values[i]: values[i+1] + 1][:, 0]
        max_new = []

        for j in range(len(max_point_ids)):
            if values[i] < max_point_ids[j] < values[i+1] + 1:
                max_new.append(max_point_ids[j])
        max_new = np.array(max_new)

        while max_new.shape[0] > 1:
            Gauss = gaussian(curv_part.shape[0], std=curv_part.shape[0]//2)    
            new = np.convolve(curv_part, Gauss, 'same')
            max_new = argrelextrema(new, np.greater)[0]

        max_landmarks[bends[i]] = max_new + values[i]

    # Following Bogunovic
    max_landmarks.pop("superior")
    landmarks = {}
    for k, v in max_landmarks.iteritems():
        landmarks[k] = line.GetPoints().GetPoint(float(v))

    for k, v in interfaces.iteritems():
        landmarks[k] = line.GetPoints().GetPoint(float(v))

    return landmarks


def main(folder):
    centerline_path = path.join(folder, "surface", "model_usr_centerline.vtp")
    centerline_bif_path = path.join(folder, "surface", "model_usr_centerline_bif.vtp")
    centerline_bif = ReadPolyData(centerline_bif_path)
    centerline = ReadPolyData(centerline_path)

    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                  centerline.GetPoint(10), \
                                  centerline.GetPoint(11)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance
    data = getData(centerline, centerline_bif, centerlineSpacing)
    line = ExtractSingleLine(centerline, 0, startID=0, endID=data["bif"]["ID_div"])
    WritePolyData(line, path.join(folder, "surface", "carotid_siphon.vtp"))
    
    landmarks = landmarking(line, folder)
    if landmarks is not None:
        writeParameters(landmarks, folder)
    

def check_landmarked(dirpath):
    parameters = getParameters(dirpath)
    return parameters.has_key("sup_ant")


if __name__ == '__main__':
    #maindir = "/home/aslak/master/src/aneurysms/"
    maindir = "."
    for folder in listdir(maindir):
        dirpath = path.join(maindir, folder)
        if path.isdir(dirpath) and folder != "backup" and folder != ".git" and "test" not in folder:
            print folder,
            #if check_landmarked(dirpath):
            #    print "is allready marked, moving on!"
            #elif folder in ["C0075", "C0063", "P0134", "C0035", "C0086",
            #        "C0088a", "C0067", "C0085", "P0252", "C0005", "C0022", "C0074a"]:
            #    print "special case, moving on!"
            #else:
            #    print "starting to landmark..."
            main(dirpath)
