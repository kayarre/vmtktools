from argparse import ArgumentParser
from os import path, listdir, sep
from subprocess import STDOUT, check_output
import sys


def ReadPolyData(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('--n', '--new', type=bool, default=False,
                        help="Only do the new geometry")
    
    args = parser.parse_args()

    return args.d, args.n


def ExtractSingleLine(centerlines, id, startID=0, endID=None):
    cell = vtk.vtkGenericCell()
    centerlines.GetCell(id, cell)
    N = cell.GetNumberOfPoints() if endID is None else endID + 1

    line = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    cellArray.InsertNextCell(N - startID)
    linePoints = vtk.vtkPoints()

    arrays = []
    N_, names = get_number_of_arrays(centerlines)
    for i in range(N_):
        tmp = centerlines.GetPointData().GetArray(names[i])
        tmp_comp = tmp.GetNumberOfComponents()
        radiusArray = get_vtk_array(names[i], tmp_comp, N - startID)
        arrays.append(radiusArray)

    getArray = []
    for i in range(N_):
        getArray.append(centerlines.GetPointData().GetArray(names[i]))

    count = 0
    for i in range(startID, N):
        cellArray.InsertCellPoint(count)
        linePoints.InsertNextPoint(cell.GetPoints().GetPoint(i))

        for j in range(N_):
            num = getArray[j].GetNumberOfComponents()
            if num == 1:
                tmp = getArray[j].GetTuple1(i)
                arrays[j].SetTuple1(count, tmp)
            elif num == 2:
                tmp = getArray[j].GetTuple2(i)
                arrays[j].SetTuple2(count, tmp[0], tmp[1])
            elif num == 3:
                tmp = getArray[j].GetTuple3(i)
                arrays[j].SetTuple3(count, tmp[0], tmp[1], tmp[2])
            elif num == 9:
                tmp = getArray[j].GetTuple9(i)
                arrays[j].SetTuple9(count, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
                                       tmp[5], tmp[6], tmp[7], tmp[8])
        count += 1

    line.SetPoints(linePoints)
    line.SetLines(cellArray)
    for j in range(N_):
        line.GetPointData().AddArray(arrays[j])
    return line


def makeCenterline(ifile, ofile, in_out=None, length=0.1, resampling=True):
    if in_out is None:
        source = ""        
    else:
        inlet = in_out[0]
        source = " -seedselector pointlist -sourcepoints %s %s %s -targetpoints " \
                    % (inlet[0], inlet[1], inlet[2])

        points = []
        for p in in_out[1:]:
            points += [str(p[0]), str(p[1]), str(p[2])]

        text = " ".join(points)
        source += text

    # Add resampling
    resampling = " -resampling 1 -resamplingstep %s" % length

    # Execute command
    a = check_output(("vmtkcenterlines -ifile %s%s%s%s -ofile %s") % \
                    (ifile, source, resampling, smooth, ofile),
                    stderr=STDOUT, shell=True)

    return ReadPolyData(ofile)


def main(d, new):
    case_folder = path.dirname(d)
    ofile = path.join(case_folder, "Case" + case_folder)

    if not new:
        cl = makeCenterline(d, ofile+"full_centerline.vtp")

        # Get points from centerline
        end_points = []
        for i in range(cl.GetNumberOfLines()):
            tmp_line = ExtractSingleLine(centerline, i)
            tmp_N = tmp_line.GetNumberOfPoints()
            end_points.append(tmp_line.GetPoint(tmp_N - 1))
        start_point = tmp_line.GetPoint(0)

        makeCenterline(d, "parentvessel.vtp", in_out=[start_point, end_points[0], end_points[1]])
        makeCenterline(d, "dau1cl.vtp", in_out=[end_points[0], end_points[-1], end_points[1]])
        makeCenterline(d, "dau1cl.vtp", in_out=[end_points[1], end_points[-1], end_points[0]])
    else:
        cl = ReadPolyData(ofile+"full_centerline.vtp")

        # Get points from centerline
        end_points = []
        for i in range(cl.GetNumberOfLines()):
            tmp_line = ExtractSingleLine(centerline, i)
            tmp_N = tmp_line.GetNumberOfPoints()
            end_points.append(tmp_line.GetPoint(tmp_N - 1))
        start_point = tmp_line.GetPoint(0)
        in_out = [p for p in end_points]
        in_out.append(start_point)
        
        cl = makeCenterline(d, d[:-4]+"_cl.vtp", in_out=in_out)


if __name__ == "__main__":
    d, new = read_command_line()
    main(d)
