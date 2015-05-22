import vtk
import numpy as np
import sys
import math

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


def get_array(arrayName, line, k=1):
    array = np.zeros((line.GetNumberOfPoints(), k))
    vtkArray = line.GetPointData().GetArray(arrayName)
    if k == 1:
        getData = line.GetPointData().GetArray(arrayName).GetTuple1
    elif k == 2:
        getData = line.GetPointData().GetArray(arrayName).GetTuple2
    elif k ==3:
        getData = line.GetPointData().GetArray(arrayName).GetTuple3

    for i in range(line.GetNumberOfPoints()):
        array[i,:] = getData(i)

    return array


def create_vtk_array(values, name, k=1):
    vtkArray = vtk.vtkDoubleArray()
    vtkArray.SetNumberOfComponents(k)
    vtkArray.SetNumberOfTuples(values.shape[0])
    vtkArray.SetName(name)
    for i in range(k):
        vtkArray.FillComponent(i, 0.0)

    if k == 1:
        for i in range(values.shape[0]):
            vtkArray.SetTuple1(i, values[i])
    elif k == 2:
        for i in range(values.shape[0]):
            vtkArray.SetTuple2(i, values[i,0], values[i,1])
    elif k == 3:
        for i in range(values.shape[0]):
            vtkArray.SetTuple3(i, values[i,0], values[i,1], values[i,2])

    return vtkArray


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
        
        
def data_to_vtkPolyData(data, header, TNB=None, PT=None):
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

    if TNB is not None:
        for i in range(3):
            radiusArray = vtk.vtkDoubleArray()
            radiusArray.SetNumberOfComponents(3)
            radiusArray.SetNumberOfTuples(data.shape[0])
            radiusArray.SetName(header[i+data.shape[1]])
            radiusArray.FillComponent(0, 0.0)
            radiusArray.FillComponent(1, 0.0)
            radiusArray.FillComponent(2, 0.0)
            info_array.append(radiusArray)

    if PT is not None:
        start = data.shape[1] if TNB is None else data.shape[1] + 3
        for i in range(2):
            radiusArray = vtk.vtkDoubleArray()
            radiusArray.SetNumberOfComponents(3)
            radiusArray.SetNumberOfTuples(PT[0].shape[0])
            radiusArray.SetName(header[i+start])
            radiusArray.FillComponent(0, 0.0)
            radiusArray.FillComponent(1, 0.0)
            radiusArray.FillComponent(2, 0.0)
            info_array.append(radiusArray)

    for i in range(data.shape[0]):
        cellArray.InsertCellPoint(i)
        linePoints.InsertNextPoint(data[i,:3])
        for j in range(3, data.shape[1]):
            info_array[j-3].SetTuple1(i, data[i, j])

    if TNB is not None:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]-3, data.shape[1], 1):
                tnb_ = TNB[j - data.shape[1]][i,:]
                info_array[j].SetTuple3(i, tnb_[0], tnb_[1], tnb_[2])

    if PT is not None:
        start = data.shape[1]-3 if TNB is None else data.shape[1]
        for i in range(PT[-1].shape[0]):
            for j in range(start, start + 2, 1):
                pt_ = PT[j - start][i, :]
                info_array[j].SetTuple3(i, pt_[0], pt_[1], pt_[2])

    line.SetPoints(linePoints)
    line.SetLines(cellArray)
    for i in range(len(header) - 3):
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


def R(n, t):
    cos = math.cos
    sin = math.sin
    n1 = n[0]; n2 = n[1]; n3 = n[2]
    r = np.array([[cos(t) + n1**2 * (1 - cos(t)),   \
                   n1*n2*(1 - cos(t)) - sin(t)*n3,  \
                   n3*n1*(1 - cos(t)) + sin(t)*n2], \
                  [n1*n2*(1 - cos(t)) + sin(t)*n3,  \
                   cos(t) + n2**2*(1 - cos(t)),     \
                   n3*n2*(1 - cos(t)) - sin(t)*n1], \
                  [n1*n3*(1 - cos(t)) - sin(t)*n2,  \
                   n2*n3*(1 - cos(t)) + sin(t)*n1,  \
                   cos(t) + n3**2*(1 - cos(t))]])
    return r


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
