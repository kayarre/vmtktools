from common import *
from patchandinterpolatecenterlines import *
import numpy as np
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


def splineCenterline(centerline):
    # TODO: Remove the top part
    # Allocate data structure to store centerline points
    data = []
    curv_coor = []
    for i in range(centerline.GetNumberOfCells()):
        data.append(np.zeros((centerline.GetCell(i).GetNumberOfPoints(), 3)))

    # Collect data from centerline
    for i in range(centerline.GetNumberOfCells()):
        line = ExtractSingleLine(centerline, i)
        if i == 0:
            WritePolyData(line, "original_line.vtp")
        curv_coor.append(get_curvilinear_coordinate(line))
        cell = vtk.vtkGenericCell()
        centerline.GetCell(i, cell)
        for j in range(cell.GetNumberOfPoints()):
            data[i][j,:] = cell.GetPoints().GetPoint(j)

    for i in range(len(data)):
        t = np.linspace(curv_coor[i][0], curv_coor[i][-1], 25)[1:-1]
        fx = splrep(curv_coor[i], data[i][:,0], k=4, t=t)
        fy = splrep(curv_coor[i], data[i][:,1], k=4, t=t)
        fz = splrep(curv_coor[i], data[i][:,2], k=4, t=t)

        fx_ = splev(curv_coor[i], fx)
        fy_ = splev(curv_coor[i], fy)
        fz_ = splev(curv_coor[i], fz)

        data = np.zeros((len(curv_coor[i]), 3))
        data[:,0] = fx_
        data[:,1] = fy_
        data[:,2] = fz_

        header = ["X", "Y", "Z"]
        line = data_to_vtkPolyData(data, header)

        # Let vmtk compute curve attributes
        line = CenterlineAttribiute(line)

        # Compute curvature from the 'exact' spline to get a robust way of
        # finding max / min points on the centerline
        dlsfx = splev(curv_coor[i], fx, der=1)
        dlsfy = splev(curv_coor[i], fy, der=1)
        dlsfz = splev(curv_coor[i], fz, der=1)

        ddlsfx = splev(curv_coor[i], fx, der=2)
        ddlsfy = splev(curv_coor[i], fy, der=2)
        ddlsfz = splev(curv_coor[i], fz, der=2)

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
        # vmtk due to an analytical representation of the spline
        length = get_curvilinear_coordinate(line)
        dddlsfx = splev(length, fx, der=3)
        dddlsfy = splev(length, fy, der=3)
        dddlsfz = splev(length, fz, der=3)

        torsion_spline = (dddlsfx*C1xC2_1 + dddlsfy*C1xC2_2 + dddlsfz * C1xC2_3) / \
                         (C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2)
        torsion_array = create_vtk_array(torsion_spline, "Torsion")
        line.GetPointData().AddArray(torsion_array)

        return line, max_point_ids, min_point_ids


def landmarking(centerline_patch_path):
    """Takes the file path to a ceterline patch created by
    patchandinterpolatecenterlines.CreateParentArteryPatches, and spline the
    centerline and uses Bogunevic et al. (2012) to do automated landmarking"""

    centerlines = ReadPolyData(centerline_patch_path)
    centerline = ExtractSingleLine(centerlines, 0)
    line, max_point_ids, min_point_ids = splineCenterline(centerline)
   
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
    

    X = np.zeros(length.shape[0])
    Y = np.zeros(length.shape[0])
    Z = np.zeros(length.shape[0])
    for i in range(X.shape[0]):
        X[i] = line.GetPoints().GetPoint(i)[0]
        Y[i] = line.GetPoints().GetPoint(i)[1]
        Z[i] = line.GetPoints().GetPoint(i)[2]

    # Tolerance parameters from Bogunevic et al. (2012)
    tol_ant_post = 60
    tol_sup_ant = 45
    tol_post_inf = 45
    tol_inf_end = 110

    # TODO: Test not flipped
    if Z.min() not in Z[argrelextrema(X, np.less)[0]]:
        value_index = Z[argrelextrema(Z, np.less)[0]].min()
        max_coronal_ids = np.array(Z.tolist().index(value_index))
    else:
        value_index = Z[argrelextrema(Z, np.greater)[0]].max()
        max_coronal_ids = np.array(Z.tolist().index(value_index))

    #from matplotlib.pyplot import *
    #viz(line, [[X[max_coronal_ids], Y[max_coronal_ids], Z[max_coronal_ids]]])
    #plot(length, Z, length[max_coronal_ids], Z[max_coronal_ids], 'g^')
    #show()

    # TODO: Make some recursive function of this algorithm
    index = np.array((max_coronal_ids > max_point_ids).nonzero()[0]).max()
    for i in range(index-1, -1, -1):
        print "i", i, tetha[i]
        if tetha[i] > tol_ant_post:
            break

    start = max_point_ids[i]
    stop = max_point_ids[i + 1]
    min_point_ids = np.array(min_point_ids)
    index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
    min_point = min_point_ids[index]
    interfaces = {'ant_post': min_point}

    for j in range(i-1, -1, -1):
        print "j", j, tetha[j]
        if tetha[j] > tol_post_inf:
            break

    start = max_point_ids[j]
    stop = max_point_ids[j + 1]
    min_point_ids = np.array(min_point_ids)
    index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
    min_point = min_point_ids[index]
    interfaces["post_inf"] = min_point

    stopped = False
    for k in range(j-1, -1, -1):
        print "k", k, tetha[k]
        if tetha[k] > tol_inf_end:
            stopped = True
            break

    # TODO: If the geometry is not long enough let inferior bend be at end
    #       point. The same test needs to be done for each landmark.
    if stopped:
        start = max_point_ids[k]
        stop = max_point_ids[k + 1]
        min_point_ids = np.array(min_point_ids)
        index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
        min_point = min_point_ids[index]
        interfaces["inf_end"] = min_point
    else:
        print "End is end"

    index = np.array((max_coronal_ids > max_point_ids).nonzero()[0]).max()
    for l in range(index, tetha.shape[0]):
        print "l", l, tetha[l]
        if tetha[l] > tol_sup_ant:
            break

    start = max_point_ids[l]
    stop = max_point_ids[l + 1]
    min_point_ids = np.array(min_point_ids)
    index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
    min_point = min_point_ids[index]
    interfaces['sup_ant'] = min_point

    bends = ["inferior", "posterior", "anterior", "superior"]
    values = [interfaces["inf_end"], interfaces["post_inf"],
              interfaces["ant_post"], interfaces["sup_ant"], 
              curvature.shape[0]]

    from matplotlib.pyplot import *
    max_landmarks = {}
    for i in range(len(bends)):
        curv_part = curvature[values[i]: values[i+1] + 1]
        max_new = argrelextrema(curv_part, np.greater)[0]
        print max_new
        plot(length[values[i]: values[i+1] +1], curv_part) 
        while max_new.shape[0] > 1:
            Gauss = gaussian(curv_part.shape[0], std=5)
            new = np.convolve(curv_part, Gauss, 'same')
            max_new = argrelextrema(new, np.greater)[0]

        max_landmarks[bends[i]] = max_new

    print max_landmarks

    if 0:
        #print max_coronal_ids
        #plot(length, X)
        #plot(length[max_coronal_ids], X[max_coronal_ids], 'bs')
        #show()
    
        #point = [[float(X[max_coronal_ids]), float(Y[max_coronal_ids]),
        #    float(Z[max_coronal_ids])]]
        #print point
        point = [line.GetPoints().GetPoint(k) for k in interfaces.values()]# +
        #        [max_coronal_ids]]
        #viz(line, point)
        fig = figure()
        ax = fig.add_subplot(111)
    
        #plot(k1[:-10], k2[:-10])
        plot(length, curvature)
        hold("on")
        #for i in min_point_ids:
            #plot(length[i], curvature[i], 'g^')
        plot([length[max_coronal_ids], length[max_coronal_ids]], [0,1]) 
            #plot(k1[i], k2[i], 'g^')
        c = 0
        for key, value in interfaces.iteritems():
            plot(length[value], curvature[value], 'bs')
            #plot(k1[i], k2[i], 's', label=c)
            #print (float(k1[i, k2[i])
            ax.annotate('%s' % c, xy=(float(length[value]),
                float(curvature[value])), textcoords='offset points')
            c += 1
        for p in max_point_ids:
            plot(length[p], curvature[p], "b^")
    
        torsion = get_array("Torsion", line)
        #plot(length, torsion)
        legend()
        grid()
        show()
    
        viz(line, point)

if __name__ == '__main__':
    landmarking("C0013/surface/cl_patch.vtp")
