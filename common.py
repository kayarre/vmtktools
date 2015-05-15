import vtk
import numpy as np
import sys
from sympy.interpolate import splrep, splev

# Global names
radiusArrayName = 'MaximumInscribedSphereRadius'
parallelTransportNormalsArrayName = 'ParallelTransportNormals'
AbscissasArrayName = 'Abscissas'
divergingRatioToSpacingTolerance = 2.0


def ReadPolyData(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def WritePolyData(input, filename):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInput(input)
    writer.Write()


def get_curvilinear_coordinate(line):
    curv_coor = np.zeros(line.GetNumberOfPoints())
    for i in range(line.GetNumberOfPoints() - 1):
        pnt1 = np.asarray(line.GetPoints().GetPoint(i))
        pnt2 = np.asarray(line.GetPoints().GetPoint(i+1))
        curv_coor[i+1] = np.sum(np.sqrt((pnt1 - pnt2)**2)) + curv_coor[i]

    return curv_coor


def get_array(arrayName, line):
    array = np.zeros(line.GetNumberOfPoints())
    for i in range(line.GetNumberOfPoints()):
        array[i] = line.GetPointData().GetArray(arrayName).GetTuple1(i)

    return array


def GramSchmidt(V):
    V = 1.0 * V
    U = np.copy(V)

    def proj(u, v):
        return u * np.dot(v,u) / np.dot(u,u)

    for i in xrange(1, V.shape[1]):
        for j in xrange(i):
            U[:,i] -= proj(U[:,j], V[:,i])

    # normalize column
    den=(U**2).sum(axis=0)**0.5
    E = U/den
    return E


def read_dat(filename):
    f = open(filename, "r")
    text = f.readline()
    
    header = text.split(" ")
    header[-1] = header[-1][:-2]

    lines = f.readlines()
    f.close()

    data = np.zeros((len(lines), len(header)))
    col_len = len(lines[0].split(" "))

    counter = 0
    for line in lines:
        values = line.split(" ")
        for i in range(col_len):
            data[counter, i] = float(values[i])
        counter += 1

    return data, header
        
        
def data_to_vtkPolyData(data, header):
    line = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    cellArray.InsertNextCell(data.shape[0])
    linePoints = vtk.vtkPoints()

    info_array = []
    for i in range(3, data.shape[1]):
        radiusArray = vtk.vtkDoubleArray()
        radiusArray.SetNumberOfComponents(1)
        radiusArray.SetNumberOfTuples(data.shape[0])
        radiusArray.SetName(header[i])
        radiusArray.FillComponent(0, 0.0)
        info_array.append(radiusArray)

    for i in range(data.shape[0]):
        cellArray.InsertCellPoint(i)
        linePoints.InsertNextPoint(data[i,:3])
        for j in range(3, data.shape[1]):
            info_array[j-3].SetTuple1(i, data[i, j])
    
    line.SetPoints(linePoints)
    line.SetLines(cellArray)
    for i in range(data.shape[1] - 3):
        line.GetPointData().AddArray(info_array[i])

    return line


def ExtractSingleLine(centerlines, id):
   cell = vtk.vtkGenericCell()
   centerlines.GetCell(id, cell)

   line = vtk.vtkPolyData()
   cellArray = vtk.vtkCellArray()
   cellArray.InsertNextCell(cell.GetNumberOfPoints())
   linePoints = cell.GetPoints()

   radiusArray = vtk.vtkDoubleArray()
   radiusArray.SetNumberOfComponents(1)
   radiusArray.SetNumberOfTuples(cell.GetNumberOfPoints())
   radiusArray.SetName(radiusArrayName)
   radiusArray.FillComponent(0, 0.0)

   for i in range(cell.GetNumberOfPoints()):
      cellArray.InsertCellPoint(i)   
      radius = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(cell.GetPointId(i))
      radiusArray.SetTuple1(i, radius)
       
   line.SetPoints(linePoints)
   line.SetLines(cellArray)
   line.GetPointData().AddArray(radiusArray)

   return line


def splineCenterline(centerline):
    # Alloimport matplotlib.pyplot as pltcate data structure to store centerline points
    data = []
    curv_coor = []
    for i in range(centerline.GetNumberOfCells()):
        data.append(zeros((centerline.GetCell(i).GetNumberOfPoints(), 3)))

    # Collect data from centerline
    for i in range(centerline.GetNumberOfCells()):
        line = ExtractSingleLine(centerline, i)
        curv_coor.append(get_curvilinear_coordinate(line))
        cell = vtk.vtkGenericCell()
        centerline.GetCell(i, cell)
        for j in range(cell.GetGetNumberOfPoints()):
            data[i][i,:] = cell.GetPoints().GetPoint(j)


    # least square fit of the centerline
    interpolated_data = []
    interpolated_data_1 = []
    interpolated_data_2 = []
    for i in len(data):
        fx = splrep(curv_coor[i], data[i][:,0], k=7)
        fy = splrep(curv_coor[i], data[i][:,1], k=7)
        fz = splrep(curv_coor[i], data[i][:,2], k=7)

        interpolated_data.append([fx, fy, fz])

        dlsfx = splev(curv_coor[i], fx, der=1)
        dlsfy = splev(curv_coor[i], fy, der=1)
        dlsfz = splev(curv_coor[i], fz, der=1)

        ddlsfx = splev(curv_coor[i], fx, der=2)
        ddlsfy = splev(curv_coor[i], fy, der=2)
        ddlsfz = splev(curv_coor[i], fz, der=2)

        C1xC2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
        C1xC2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
        C1xC2_3 = ddlsfy*dlsfx - ddlsfx * dlsfy

        lscurvature = np.sqrt(C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2) / \
                             (dlsfx**2 + dlsfy**2 + dlsfz**2)**1.5

        from matplotlib.pyplot import *
        plot(curv_coor[i], lscurvature)
        show()
        sys.exit(0)


def viz(centerline, points):
    """Help method during development to view the results"""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    N = centerline.GetNumberOfCells()
    for i in range(N):
        point_ids = vtk.vtkIdList()
        centerline.GetCellPoints(i, point_ids)
        points0 = []
        for k in range(point_ids.GetNumberOfIds()):
            points0.append(centerline.GetPoint(point_ids.GetId(k)))
        arr = np.asarray(points0)
        x = arr[:,0]
        y = arr[:,1]
        z = arr[:,2]
        ax.plot(x, y, z, label=i)
        ax.legend()
        plt.hold("on")
    
    counter = 0
    for p in points:
        ax.plot([p[0]], [p[1]], [p[2]],"o", label=N+counter)
        ax.legend()
        plt.hold("on")
        counter += 1

    plt.show()
