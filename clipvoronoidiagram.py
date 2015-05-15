#!/usr/bin/env python

import vtk
import sys
import math
from common import *
from vmtk import vtkvmtk
from os import path


def MaskVoronoiDiagram(voronoi, centerlines):
   numberOfCenterlinesPatches = centerlines.GetNumberOfCells()
   numberOfVoronoiPoints = voronoi.GetNumberOfPoints()

   maskArray = vtk.vtkIntArray()
   maskArray.SetNumberOfComponents(1)
   maskArray.SetNumberOfTuples(numberOfVoronoiPoints)
   maskArray.FillComponent(0, 0)

   for i in range(numberOfCenterlinesPatches):
     tangent, center, centerMISR = ComputePatchEndPointParameters(i, centerlines) 
     MaskWithPatch(i, tangent, center, centerMISR, maskArray, centerlines, voronoi)

   return maskArray


def ComputePatchEndPointParameters(id, centerlines):
   point0 = [0.0, 0.0, 0.0]
   point1 = [0.0, 0.0, 0.0]
   tan = [0.0, 0.0, 0.0]
   radius0 = -1

   cell = vtk.vtkGenericCell()
   centerlines.GetCell(id, cell)

   if (id == 0):
     point0 = cell.GetPoints().GetPoint(cell.GetNumberOfPoints() - 1)
     point1 = cell.GetPoints().GetPoint(cell.GetNumberOfPoints() - 2)
     radius0 = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(cell.GetPointId(cell.GetNumberOfPoints() - 1))
   
   else:
     point0 = cell.GetPoints().GetPoint(0)
     point1 = cell.GetPoints().GetPoint(1)
     radius0 = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(cell.GetPointId(0))
   
   tan[0] = point1[0] - point0[0]
   tan[1] = point1[1] - point0[1]
   tan[2] = point1[2] - point0[2]
   vtk.vtkMath.Normalize(tan)

   return tan, point0, radius0 


def MaskWithPatch(id, t, c, r, maskArray, centerlines, voronoi):
   patch = ExtractSingleLine(centerlines, id)
   
   tubeFunction = vtkvmtk.vtkvmtkPolyBallLine()
   tubeFunction.SetInput(patch)
   tubeFunction.SetPolyBallRadiusArrayName(radiusArrayName)

   lastSphere = vtk.vtkSphere()
   lastSphere.SetRadius(r * 1.5)
   lastSphere.SetCenter(c)

   for i in range(voronoi.GetNumberOfPoints()):
        point = [0.0, 0.0, 0.0]
        voronoiVector = [0.0, 0.0, 0.0]

        voronoi.GetPoint(i, point)
        voronoiVector = [point[j] - c[j] for j in range(3)]
        voronoiVectorDot = vtk.vtkMath.Dot(voronoiVector, t)

        tubevalue = tubeFunction.EvaluateFunction(point)
        spherevalue = lastSphere.EvaluateFunction(point)
        
        if (spherevalue < 0.0) & (voronoiVectorDot < 0.0): continue
        elif (tubevalue <= 0.0):
            maskArray.SetTuple1(i, 1)


def ComputeNumberOfMaskedPoints(dataArray):
   numberOfPoints = 0
   for i  in range(dataArray.GetNumberOfTuples()):
      value = dataArray.GetTuple1(i)
      if (value == 1): numberOfPoints += 1
   return numberOfPoints   


def ExtractMaskedVoronoiPoints(voronoi,maskArray):
   numberOfPoints = ComputeNumberOfMaskedPoints(maskArray)

   maskedVoronoi = vtk.vtkPolyData()
   maskedPoints = vtk.vtkPoints()
   cellArray = vtk.vtkCellArray()

   radiusArray = vtk.vtkDoubleArray()
   radiusArray.SetNumberOfComponents(1)
   radiusArray.SetNumberOfTuples(numberOfPoints)
   radiusArray.SetName(radiusArrayName)
   radiusArray.FillComponent(0, 0.0)

   count = 0
   for i in range(voronoi.GetNumberOfPoints()):
      point = [0.0, 0.0, 0.0]
      voronoi.GetPoint(i, point)
      pointRadius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i) 
      value = maskArray.GetTuple1(i)
      if (value == 1):
        maskedPoints.InsertNextPoint(point)
        radiusArray.SetTuple1(count, pointRadius)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(count)
        count += 1

   maskedVoronoi.SetPoints(maskedPoints)
   maskedVoronoi.SetVerts(cellArray)
   maskedVoronoi.GetPointData().AddArray(radiusArray)
   return maskedVoronoi


def SmoothClippedVoronoiDiagram(voronoi, centerlines, smoothingFactor):
   numberOfPoints = voronoi.GetNumberOfPoints()
   numberOfCenterlinesPoints = centerlines.GetNumberOfPoints()

   maskArray = vtk.vtkIntArray()
   maskArray.SetNumberOfComponents(1)
   maskArray.SetNumberOfTuples(numberOfPoints)
   maskArray.FillComponent(0,0)

   for i in range(numberOfCenterlinesPoints):
     localRadius = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

     threshold = localRadius * (1.0 - smoothingFactor) 

     sphere = vtk.vtkSphere()
     sphere.SetRadius(localRadius)
     sphere.SetCenter(centerlines.GetPoint(i))

     localMaskArray = vtk.vtkIntArray()
     localMaskArray.SetNumberOfComponents(1)
     localMaskArray.SetNumberOfTuples(numberOfPoints)
     localMaskArray.FillComponent(0,0)

     for j in range(numberOfPoints):
        value = sphere.EvaluateFunction(voronoi.GetPoint(j))    
        if (value <= 0.0):
           localMaskArray.SetTuple1(j, 1) 

     for j in range(numberOfPoints):
       value = localMaskArray.GetTuple1(j)
       if (value == 1):
          r = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(j) 
          if (r > threshold):
             maskArray.SetTuple1(j, 1)

   finalNumberOfMaskedPoints = ComputeNumberOfMaskedPoints(maskArray)
   print 'from original number of points ', numberOfPoints,'to ',finalNumberOfMaskedPoints

   smoothedDiagram = vtk.vtkPolyData()
   points = vtk.vtkPoints()
   cellArray = vtk.vtkCellArray()

   radiusArray = vtk.vtkDoubleArray()
   radiusArray.SetNumberOfComponents(1)
   radiusArray.SetNumberOfTuples(finalNumberOfMaskedPoints)
   radiusArray.FillComponent(0, 0.0)
   radiusArray.SetName(radiusArrayName)

   count = 0
   for i in range(numberOfPoints):
      value = maskArray.GetTuple1(i)
      if (value == 1):
         radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
      
         points.InsertNextPoint(voronoi.GetPoint(i))
         cellArray.InsertNextCell(1)
         cellArray.InsertCellPoint(count) 
         radiusArray.SetTuple1(count, radius)
         count += 1

   smoothedDiagram.SetPoints(points)
   smoothedDiagram.SetVerts(cellArray)
   smoothedDiagram.GetPointData().AddArray(radiusArray)

   return smoothedDiagram


#SOME COMMON VMTK DATA ARRAY NAMES
smoothVoronoiDiagram = 1	# recommended; perform smoothing of the Voronoi Diagram to 

if __name__ == "__main__":
    if (smoothVoronoiDiagram==1):
        smoothVoronoiDiagram = SmoothClippedVoronoiDiagram(clippedVoronoiDiagram,patchCenterlines, 0.25)
