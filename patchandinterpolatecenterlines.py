#!/usr/bin/env python

from common import *
import vtk
import sys
import math
from vmtk import vtkvmtk


def CreateParentArteryPatches(parentCenterlines, clipPoints):
    numberOfDaughterPatches = parentCenterlines.GetNumberOfCells()
    clipIds, numberOfPatchedCenterlinesPoints = ExtractPatchesIds(parentCenterlines, clipPoints)
    pnt = []

    patchedCenterlines = vtk.vtkPolyData()
    patchedCenterlinesPoints = vtk.vtkPoints()
    patchedCenterlinesCellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPatchedCenterlinesPoints)

    numberOfCommonPatch = clipIds[0]+1
    patchedCenterlinesCellArray.InsertNextCell(numberOfCommonPatch)

    count = 0
    line = ExtractSingleLine(parentCenterlines, 0)
    getData = line.GetPointData().GetArray(radiusArrayName).GetTuple1
    for i in range(0, numberOfCommonPatch):
        patchedCenterlinesPoints.InsertNextPoint(line.GetPoint(i))
        patchedCenterlinesCellArray.InsertCellPoint(i)
        radiusArray.SetTuple1(i, getData(i))
        count+=1

    for j in range(numberOfDaughterPatches):
        cell = ExtractSingleLine(parentCenterlines, j)

        getData = cell.GetPointData().GetArray(radiusArrayName).GetTuple1
        numberOfCellPoints = cell.GetNumberOfPoints()
        startId = clipIds[j+1]

        patchNumberOfPoints = numberOfCellPoints-startId
        patchedCenterlinesCellArray.InsertNextCell(patchNumberOfPoints)

        for i in range(startId, cell.GetNumberOfPoints()):
            point = cell.GetPoint(i)
            patchedCenterlinesPoints.InsertNextPoint(point)
            patchedCenterlinesCellArray.InsertCellPoint(count)
            radiusArray.SetTuple1(count, getData(i))
            count+=1

    patchedCenterlines.SetPoints(patchedCenterlinesPoints)
    patchedCenterlines.SetLines(patchedCenterlinesCellArray)
    patchedCenterlines.GetPointData().AddArray(radiusArray)
    
    return patchedCenterlines


def ExtractPatchesIds(parentCl, clipPts):
    distance = vtk.vtkMath.Distance2BetweenPoints
    clipIds = []
    numberOfPoints = 0
    N = clipPts.GetNumberOfPoints()
    
    if N == 3:
        commonPoint = clipPts.GetPoint(0)
        pnt_1 = clipPts.GetPoint(1)
        pnt_2 = clipPts.GetPoint(2)
    else:
        pnt_1 = clipPts.GetPoint(0)
        pnt_2 = clipPts.GetPoint(1)

    for j in range(parentCl.GetNumberOfCells()):
        cellLine = ExtractSingleLine(parentCl, j)
        locator = get_locator(cellLine)

        if j == 0 and N == 3:
            upstreamId = locator.FindClosestPoint(commonPoint)
            clipIds.append(upstreamId)
            numberOfPoints += upstreamId + 1

        ID1 = locator.FindClosestPoint(pnt_1)
        ID2 = locator.FindClosestPoint(pnt_2)
        
        distance1 = math.sqrt(distance(pnt_1, cellLine.GetPoints().GetPoint(ID1)))
        distance2 = math.sqrt(distance(pnt_2, cellLine.GetPoints().GetPoint(ID2)))

        if distance1 > 1 and distance2 > 1:
            ID = 0
        else:
            ID = ID1 if distance1 < distance2 else ID2
        
        if N == 2:
            clipIds = [ID1, ID2]
            numberOfPoints = cellLine.GetNumberOfPoints()
        else:
            clipIds.append(ID)
            numberOfPoints += cellLine.GetNumberOfPoints() - ID

    return clipIds, numberOfPoints


def InterpolatePatchCenterlines(patchCenterlines, parentCenterlines,
                                clippingPoints, divergingPoints, addPoint):
    additionalPoint = [-1.0, -1.0, -1.0]
    additionalPointIds = [] 

    if addPoint:
        additionalPoint = divergingPoints.GetPoint(0)
        for i in range(parentCenterlines.GetNumberOfCells()):
            line = ExtractSingleLine(parentCenterlines, i)
            additionalPointIds.append(line.FindPoint(additionalPoint))
    else:
        additionalPoint = clippingPoints.GetPoint(0)
        for i in range(parentCenterlines.GetNumberOfCells()):	   
            line = ExtractSingleLine(parentCenterlines, 0)
            additionalPointIds.append(line.FindPoint(additionalPoint)) 
  
    interpolatedLines = vtk.vtkPolyData()
    interpolatedPoints = vtk.vtkPoints()
    interpolatedCellArray = vtk.vtkCellArray()

    pointsInserted = 0
    interpolatedCellArray.Initialize()

    for i in range(parentCenterlines.GetNumberOfCells()):
        startingCell = vtk.vtkGenericCell()
        endingCell = vtk.vtkGenericCell()

        numberOfInterpolationPoints = parentCenterlines.GetCell(i).GetNumberOfPoints()

        patchCenterlines.GetCell(0, startingCell)
        patchCenterlines.GetCell(i+1, endingCell)
    
        splinePoints = InterpolateTwoCells(startingCell, endingCell, \
                                            numberOfInterpolationPoints, \
                                            additionalPointIds[i],
                                            additionalPoint, addPoint)

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


def InterpolateTwoCells(startCell, endCell, numberOfSplinePoints, additionalPointId,
                        additionalPoint, addPoint):
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

    if addPoint:
        xspline.AddPoint(float(additionalPointId), additionalPoint[0])
        yspline.AddPoint(float(additionalPointId), additionalPoint[1])
        zspline.AddPoint(float(additionalPointId), additionalPoint[2])

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
