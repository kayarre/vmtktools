#!/usr/bin/env python

from common import *
import vtk
import sys
import math
from vmtk import vtkvmtk


def CreateParentArteryPatches(parentCenterlines, clipPoints):
   numberOfDaughterPatches = parentCenterlines.GetNumberOfCells()

   patchedCenterlines = vtk.vtkPolyData()
   patchedCenterlinesPoints = vtk.vtkPoints()
   patchedCenterlinesCellArray = vtk.vtkCellArray()
   patchedRadiusArray = vtk.vtkDoubleArray()

   clipIds, numberOfPatchedCenterlinesPoints = ExtractPatchesIds(parentCenterlines, clipPoints)
   pnt = []

   radiusArray = vtk.vtkDoubleArray()
   radiusArray.SetNumberOfComponents(1)
   radiusArray.SetName(radiusArrayName)
   radiusArray.SetNumberOfTuples(numberOfPatchedCenterlinesPoints) 
   radiusArray.FillComponent(0,0.0)

   numberOfCommonPatch = clipIds[0]+1
   patchedCenterlinesCellArray.InsertNextCell(numberOfCommonPatch)

   count = 0
   for i in range(0, numberOfCommonPatch):
      patchedCenterlinesPoints.InsertNextPoint(parentCenterlines.GetPoint(i))
      patchedCenterlinesCellArray.InsertCellPoint(i)
      radiusArray.SetTuple1(i, parentCenterlines.GetPointData().GetArray(radiusArrayName).GetTuple1(i))
      count+=1 
 
   for j in range(numberOfDaughterPatches):
      cell = vtk.vtkGenericCell()
      parentCenterlines.GetCell(j,cell)

      numberOfCellPoints = cell.GetNumberOfPoints()
      startId = clipIds[j+1]
      patchNumberOfPoints = numberOfCellPoints-startId
      patchedCenterlinesCellArray.InsertNextCell(patchNumberOfPoints)

      for i in range(startId, cell.GetNumberOfPoints()):
         point = cell.GetPoints().GetPoint(i)
         patchedCenterlinesPoints.InsertNextPoint(point)
         patchedCenterlinesCellArray.InsertCellPoint(count)
         radiusArray.SetTuple1(count,parentCenterlines.GetPointData().GetArray(radiusArrayName).GetTuple1(cell.GetPointId(i)))
         count+=1

   patchedCenterlines.SetPoints(patchedCenterlinesPoints)
   patchedCenterlines.SetLines(patchedCenterlinesCellArray)
   patchedCenterlines.GetPointData().AddArray(radiusArray)

   return patchedCenterlines


def ExtractPatchesIds(parentCl, clipPts):
    distance = vtk.vtkMath.Distance2BetweenPoints
    clipIds = []
    numberOfPoints = 0

    commonPoint = clipPts.GetPoint(0)
    pnt_1 = clipPts.GetPoint(1)
    pnt_2 = clipPts.GetPoint(2)
    for j in range(parentCl.GetNumberOfCells()):
        cellLine = ExtractSingleLine(parentCl, j)

        locator = vtk.vtkPointLocator()
        locator.SetDataSet(cellLine)
        locator.BuildLocator()

        if j==0:
            upstreamId = locator.FindClosestPoint(commonPoint)
            clipIds.append(upstreamId)
            numberOfPoints += upstreamId + 1

        ID1 = locator.FindClosestPoint(pnt_1) 
        ID2 = locator.FindClosestPoint(pnt_2)

        distance1 = math.sqrt(distance(pnt_1, cellLine.GetPoints().GetPoint(ID1)))
        distance2 = math.sqrt(distance(pnt_2, cellLine.GetPoints().GetPoint(ID2)))

        ID = ID1 if distance1 < distance2 else ID2

        clipIds.append(ID)
        numberOfPoints += cellLine.GetNumberOfPoints() - ID

    return clipIds, numberOfPoints


def InterpolatePatchCenterlines(patchCenterlines, parentCenterlines, clippingPoints):
   additionalPoint = [-1.0,-1.0,-1.0]
   additionalPointIds = [] 

   if (useAdditionalInterpolationPoint == 1):
      additionalPoint = divergingPoints.GetPoint(0)
      line1 = ExtractSingleLine(parentCenterlines,0) 
      line2 = ExtractSingleLine(parentCenterlines,1) 
      additionalPointIds.append(line1.FindPoint(additionalPoint))
      additionalPointIds.append(line2.FindPoint(additionalPoint))
   else:
      for i in range(parentCenterlines.GetNumberOfCells()):	   
          additionalPoint = clippingPoints.GetPoint(0)
          line1 = ExtractSingleLine(parentCenterlines,0)
          additionalPointIds.append(line1.FindPoint(additionalPoint)) 
  
   interpolatedLines = vtk.vtkPolyData()
   interpolatedPoints = vtk.vtkPoints()
   interpolatedCellArray = vtk.vtkCellArray()

   pointsInserted = 0
   interpolatedCellArray.Initialize()

   for i in range(parentCenterlines.GetNumberOfCells()):
      startingCell = vtk.vtkGenericCell()
      endingCell = vtk.vtkGenericCell()

      numberOfInterpolationPoints = parentCenterlines.GetCell(i).GetNumberOfPoints()

      patchCenterlines.GetCell(0,startingCell)
      patchCenterlines.GetCell(i+1,endingCell)
     
      splinePoints = InterpolateTwoCells(startingCell,endingCell,numberOfInterpolationPoints,additionalPointIds[i],additionalPoint)

      interpolatedCellArray.InsertNextCell(splinePoints.GetNumberOfPoints())
      for j in range(splinePoints.GetNumberOfPoints()):
         interpolatedPoints.InsertNextPoint(splinePoints.GetPoint(j))
         interpolatedCellArray.InsertCellPoint(pointsInserted + j)
      pointsInserted += splinePoints.GetNumberOfPoints()

   interpolatedLines.SetPoints(interpolatedPoints)
   interpolatedLines.SetLines(interpolatedCellArray)
 
   attributeFilter = vtkvmtk.vtkvmtkCenterlineAttributesFilter()
   attributeFilter.SetInput(interpolatedLines)
   attributeFilter.SetAbscissasArrayName(AbscissasArrayName)
   attributeFilter.SetParallelTransportNormalsArrayName(parallelTransportNormalsArrayName)
   attributeFilter.Update()

   attributeInterpolatedLines = attributeFilter.GetOutput()
  
   return attributeInterpolatedLines


def InterpolateTwoCells(startCell,endCell,numberOfSplinePoints,additionalPointId,additionalPoint):
   points = vtk.vtkPoints()
   xspline = vtk.vtkCardinalSpline()
   yspline = vtk.vtkCardinalSpline()
   zspline = vtk.vtkCardinalSpline()
   
   numberOfStartCellPoints = startCell.GetNumberOfPoints()
   numberOfEndCellPoints = endCell.GetNumberOfPoints()
   endCellFirstId = numberOfSplinePoints - numberOfEndCellPoints
   numberOfClippedCenterlinePoints = numberOfStartCellPoints + numberOfEndCellPoints

   for i in range(numberOfStartCellPoints):
      point = startCell.GetPoints().GetPoint(i)
      xspline.AddPoint(float(i),point[0])
      yspline.AddPoint(float(i),point[1])
      zspline.AddPoint(float(i),point[2])

   if (useAdditionalInterpolationPoint == 1):
      xspline.AddPoint(float(additionalPointId),additionalPoint[0])
      yspline.AddPoint(float(additionalPointId),additionalPoint[1])
      zspline.AddPoint(float(additionalPointId),additionalPoint[2])

   for i in range(numberOfEndCellPoints):
      point = endCell.GetPoints().GetPoint(i)
      index = endCellFirstId+i
      xspline.AddPoint(float(endCellFirstId + i),point[0])
      yspline.AddPoint(float(endCellFirstId + i),point[1])
      zspline.AddPoint(float(endCellFirstId + i),point[2])
   xspline.Compute()
   yspline.Compute()  
   zspline.Compute()
   
   points.SetNumberOfPoints(numberOfSplinePoints)
   for i in range(numberOfSplinePoints):
      points.SetPoint(i,xspline.Evaluate(float(i)),yspline.Evaluate(float(i)),zspline.Evaluate(float(i)))
   return points


useAdditionalInterpolationPoint = 0      # automatically set to 1 for terminal aneurysms  