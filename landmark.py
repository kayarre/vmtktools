from common import *
from patchandinterpolatecenterlines import *
import numpy as np
from scipy.interpolate import splrep, splev
import math

def splineCenterline(centerline):
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
        fx = splrep(curv_coor[i], data[i][:,0], k=5, t=t)
        fy = splrep(curv_coor[i], data[i][:,1], k=5, t=t)
        fz = splrep(curv_coor[i], data[i][:,2], k=5, t=t)

        fx_ = splev(curv_coor[i], fx)
        fy_ = splev(curv_coor[i], fy)
        fz_ = splev(curv_coor[i], fz)

        dlsfx = splev(curv_coor[i], fx, der=1)
        dlsfy = splev(curv_coor[i], fy, der=1)
        dlsfz = splev(curv_coor[i], fz, der=1)

        ddlsfx = splev(curv_coor[i], fx, der=2)
        ddlsfy = splev(curv_coor[i], fy, der=2)
        ddlsfz = splev(curv_coor[i], fz, der=2)

        dddlsfx = splev(curv_coor[i], fx, der=3)
        dddlsfy = splev(curv_coor[i], fy, der=3)
        dddlsfz = splev(curv_coor[i], fz, der=3)

        C1xC2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
        C1xC2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
        C1xC2_3 = ddlsfy*dlsfx - ddlsfx * dlsfy

        # Compute curvature and torsion
        lscurvature = np.sqrt(C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2) / \
                             (dlsfx**2 + dlsfy**2 + dlsfz**2)**1.5
        lstorsion = (dddlsfx*C1xC2_1 + dddlsfy*C1xC2_2 + dddlsfz * C1xC2_3) / \
                      (C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2)

        # Compute Frenet normals
        T = np.zeros((len(curv_coor[i]), 3))
        N = np.zeros((len(curv_coor[i]), 3))
        B = np.zeros((len(curv_coor[i]), 3))

        dx_length = np.sqrt(dlsfx**2 + dlsfy**2 + dlsfz**2)
        T[:,0] = dlsfx / dx_length
        T[:,1] = dlsfy / dx_length
        T[:,2] = dlsfz / dx_length

        ddx_length = np.sqrt(ddlsfx**2 + ddlsfy**2 + ddlsfz**2)
        N[:,0] = ddlsfx / ddx_length
        N[:,1] = ddlsfy / ddx_length
        N[:,2] = ddlsfz / ddx_length
        print np.max(N)
        print np.min(abs(N))
        print N

        T_cross_N = np.cross(N, T)
        B_length = np.sqrt(np.sum(T_cross_N**2, axis=1))
        B[:,0] = T_cross_N[:, 0] / B_length
        B[:,1] = T_cross_N[:, 1] / B_length
        B[:,2] = T_cross_N[:, 2] / B_length

        # Compute parallel transport frame
        E = np.eye(3)
        E[:,0] = T[0,:]
        E = GramSchmidt(E)
        E1_ = E[:, 1]
        E2_ = E[:, 2]

        E1 = np.zeros((T.shape[0], 3))
        E2 = np.zeros((T.shape[0], 3))
        E1[0,:] = E1_.T
        E2[0,:] = E2_.T
        
        H = np.cross(T[:-1], T[1:])
        length = np.sqrt(np.sum(H**2, axis=1))
        H[:,0] = H[:,0] / length
        H[:,1] = H[:,1] / length
        H[:,2] = H[:,2] / length

        for k in range(T.shape[0] - 1):
            tetha = np.arccos(np.dot(T[k, :], T[k+1, :]))
            R_ = R(H[k, :], tetha)
            E1[k+1, :] = np.dot(R_, E1[k, :].T).T
            E2[k+1, :] = np.dot(R_, E2[k, :].T).T

        # Compute k_1, k_2 furfilling T' = curv(s)*N(s) = k_1(s)*E_1(s) + k_2(s)*E_2(s)
        k2 = lscurvature * (E1[:,1]*N[:,0] - N[:,1]*E1[:,0]) / \
                           (E2[:,1]*E1[:,0] - E2[:,0]*E1[:,1])
        k1 = -(lscurvature * N[:,0] + k2*E2[:,2]) / E1[:,0]

        data_ = np.zeros((len(fx_), 8))
        data_[:, 0] = fx_
        data_[:, 1] = fy_
        data_[:, 2] = fz_
        data_[:, 3] = lscurvature
        data_[:, 4] = lstorsion
        data_[:, 5] = k1
        data_[:, 6] = k2
        data_[:, 7] = curv_coor[i]

        header = ["X", "Y", "Z", "Curvature", "Torsion", "k1", "k2", "Length", "FrenetTanget",
                  "FrenetNormal", "FrenetBinormal", "ParallelTransportFrame1",
                  "ParallelTransportFrame2"]
        line = data_to_vtkPolyData(data_, header, TNB=[T, N, B], PT=[E1, E2])

        return line
        #WritePolyData(line, "spline_analytical.vtp")


def landmarking(centerline_patch_path):
    """Takes the file path to a ceterline patch created by
    patchandinterpolatecenterlines.CreateParentArteryPatches, and spline the
    centerline and uses Bogunevic et al. (2012) to do automated landmarking"""

    centerlines = ReadPolyData(centerline_patch_path)
    centerline = ExtractSingleLine(centerlines, 0)
    line = splineCenterline(centerline)
   
    curvature  = get_array("Curvature", line)
    length = get_array("Length", line)
    k1 = get_array("k1", line)
    k2 = get_array("k2", line)
    max_curvature_id = []
    min_curvature_id = []

    for i in range(1, curvature.shape[0] - 1):
        if curvature[i-1] < curvature[i] > curvature[i+1]:
            max_curvature_id.append(i)
        if curvature[i-1] > curvature[i] < curvature[i+1]:
            min_curvature_id.append(i)

    # TODO: Compute tetha for all k
    # TODO: Compute coronal coordinate, distance to x-y plane?
    # TODO: Implement smoothing to find "max" curvature point of a bend.

    # Remove "sadle point" that are very close in space and value
    for i in min_curvature_id:
        for j in max_curvature_id:
            if abs(i-j) < 10 and abs(curvature[i] - curvature[j]) < 0.01:
                min_curvature_id.remove(i)
                max_curvature_id.remove(j)
     
    k1_points = k1[max_curvature_id]
    k2_points = k2[max_curvature_id]
    k_points = np.zeros((k1_points.shape[0], 2))
    k_points[:,0] = k1_points
    k_points[:,1] = k2_points
    tetha = np.zeros(len(k1_points) - 1)
    for i in range(tetha.shape[0]):
        a = k_points[i]
        b = k_points[i+1]
        tetha[i] = math.cosh(np.dot(a, b) / math.sqrt(np.sum(b**2)) / math.sqrt(np.sum(a**2)))
        tetha[i] = tetha[i] * 180 / math.pi

    print tetha
    

    from matplotlib.pyplot import *
    plot(k1, k2)

    #plot(length, curvature)
    hold("on")
    for i in min_curvature_id:
        #plot(length[i], curvature[i], 'g^')
        plot(k1[i], k2[i], 'g^')
    for i in max_curvature_id:
        #plot(length[i], curvature[i], 'bs')
        plot(k1[i], k2[i], 'bs')
    show()

if __name__ == '__main__':
    landmarking("cl_patch.vtp")
