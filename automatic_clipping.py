from argparse import ArgumentParser
from common import *
import numpy as np

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('--m', type=str, default="model_smoothed.vtp", help="Name of the model file")
    parser.add_argument('--c', type=str, default="centerline_complete.vtp", help="Name of the centerline file")
    parser.add_argument('--anu', type=bool, default=False)

    args = parser.parse_args()

    return args.d, args.m, args.c, args.anu


def move_past_sphere(cl, points, direction):
    i = 0 if not direction else cl.GetNumberOfPoints() - 1
    j = 0 if direction else cl.GetNumberOfPoints() - 1
    r = cl.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
    center = cl.GetPoints().GetPoint(i)

    MISphere = vtk.vtkSphere()
    MISphere.SetCenter(center)
    MISphere.SetRadius(r*(1./3))
    direction = -1 if direction else 1
    for k in range(i, j, direction):
        value = MISphere.EvaluateFunction(cl.GetPoint(k))
        if value >= 0:
            break

    return cl.GetPoint(k), cl.GetPointData().GetArray(radiusArrayName).GetTuple1(i)


def getBoundingBox(cl, inlet):
    endPoint = cl.GetPoint(cl.GetNumberOfPoints() - 1) if inlet else cl.GetPoint(0)
    bottom, bottom_r = move_past_sphere(cl, endPoint, inlet)

    line = CenterlineAttribiutes(cl)

    E1 = get_array("ParallelTransportNormals", line, k=3)
    E1 = E1[E1.shape[0]-1,:]
    T = get_array("FrenetTangent", line, k=3)
    T = T[T.shape[0]-1,:]
    E2 = np.zeros(3)

    V = np.eye(3)
    V[:, 0] = T
    V[:, 1] = E1
    V = GramSchmidt(V)

    E1 = V[:,1] * bottom_r * 1.5
    E2 = V[:,2] * bottom_r * 1.6
    T = T * bottom_r * 3 if not inlet else T * bottom_r * 3 * (-1)


    corners = []
    for O in [bottom, bottom + T]:
        for dir1 in [1, -1]:
            for dir2 in [1, -1]:
                corners.append(O + dir1*E1 + dir2*E2)

    viz(line, [bottom, endPoint] + corners)

    corners = np.array(corners)
    limits = []
    for i in range(3):
        for f in [np.min, np.max]:
            limits.append(f(corners[:,i]))

    return limits


def clipp(dir_path, model, centerline, anu):
    cl = ReadPolyData(path.join(dir_path, centerline))
    surface = ReadPolyData(path.join(dir_path, model))

    #clipper = vtk.vtkBoxClipDataSet()
    box = vtk.vtkBox()
    clipper = vtk.vtkClipPolyData()
    clipper.SetInput(surface)
    clipper.SetClipFunction(box)

    inlet = True
    for i in [0] + range(cl.GetNumberOfLines() - anu):
        limits = getBoundingBox(ExtractSingleLine(cl, i), inlet)
        inlet = False
        box.SetBonds(limits[0], limits[1], limits[2], limits[3], limits[4], limits[5])
        clipper.Update()
        #clipper.SetBoxClip(limits[0], limits[1], limits[2], limits[3], limits[4], limits[5])
        clipper.GenerateClippedOutputOn()

        clipper.Update()
        filter = vtk.vtkGeometryFilter()
        filter.SetInput(clipper.GetClippedOutput())
        filter.Update()
        surface.DeepCopy(filter.GetOutput())
        clipper.Update()
        
        #TODO: Deep copy of surface and update clipper

        WritePolyData(surface, "test_clipping.vtp")
        sys.exit(0)


if __name__ == "__main__":
    dir_path, model, centerline, anu = read_command_line()
    clipp(dir_path, model, centerline, anu)
