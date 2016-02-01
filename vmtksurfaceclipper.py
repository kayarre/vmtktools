#!/usr/bin/env python

import vtk
import sys
import vmtkrenderer
from argparse import ArgumentParser


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('--m', type=str, default="model.vtp", help="Name of the model file")
    parser.add_argument('--c', type=str, default="centerline.vtp",
                        help="Name of the centerline file")


Surface = None

def ClipCallback(ClipWidget, Clipper, ClipFunction, Surface, ClippedSurface,
                 Cutter, CutLines)
    ClipWidget.GetPlanes(ClipFunction)
    Clipper.Update()
    Surface.DeepCopy(Clipper.GetOutput())
    Surface.Update()
    ClippedSurface.DeepCopy(Clipper.GetClippedOutput())
    ClippedSurface.Update()
    Cutter.Update()
    CutLines.DeepCopy(Cutter.GetOutput())


def Execute(Surface, cl):
    if Surface == None:
        print 'Error: no Surface.'
        sys.exit(0)

    # Initialize
    Clipper = vtk.vtkClipPolyData()
    Clipper.SetInput(Surface)
    Clipper.GenerateClippedOutputOn()
    Clipper.SetInsideOut(0)
 
    ClipFunction = vtk.vtkPlanes()
       
    Clipper.SetClipFunction(ClipFunction)

    Cutter = vtk.vtkCutter()
    Cutter.SetInput(Surface)
    Cutter.SetCutFunction(ClipFunction)

    ClippedSurface = vtk.vtkPolyData()
    CutLines = vtk.vtkPolyData()

    ClipWidget = vtk.vtkBoxWidget()
    ClipWidget.GetFaceProperty().SetColor(0.6,0.6,0.2)
    ClipWidget.GetFaceProperty().SetOpacity(0.25)

    Transform = vtk.vtkTransform()
    ClipWidget.GetTransform(Transform)

    for i in range(cl.GetNumberOfLines()):
        # TODO: Implement thing here.

    # Clean surface
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInput(Surface)
    cleaner.Update()
    Surface = cleaner.GetOutput()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInput(ClippedSurface)
    cleaner.Update()
    ClippedSurface = cleaner.GetOutput()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInput(CutLines)
    cleaner.Update()
    stripper = vtk.vtkStripper()
    stripper.SetInput(cleaner.GetOutput())
    stripper.Update()
    CutLines = stripper.GetOutput()

    if Surface.GetSource():
        Surface.GetSource().UnRegisterAllOutputs()


if __name__=='__main__':
    surface = read_command_line()
    Execute()
