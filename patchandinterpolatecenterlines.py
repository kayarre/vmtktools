#!/usr/bin/env python

from common import *
import vtk
from scipy.interpolate import splrep, splev
import sys
import math
from vmtk import vtkvmtk
import numpy as np

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
        count += 1

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
                                additionalPoint, lower, version):
   
    if additionalPoint is not None:
        additionalPointIds = []
        for i in range(parentCenterlines.GetNumberOfCells()):
            line = ExtractSingleLine(parentCenterlines, i)
            additionalPointIds.append(line.FindPoint(additionalPoint))
    else:
        additionalPointIds = ["" for i in range(parentCenterlines.GetNumberOfCells())]
  
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
    
        if version:
            splinePoints = InterpolateSpline(startingCell, endingCell, additionalPoint)
        else:
            splinePoints = InterpolateTwoCells(startingCell, endingCell, \
                                               numberOfInterpolationPoints, \
                                               additionalPointIds[i],
                                               additionalPoint, lower)

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


def InterpolateSpline(startCell, endCell, additionalPoint):

    # If the centerline does not pass the bifurcation, return the centerline
    if startCell.GetPoints().GetPoint(0) == endCell.GetPoints().GetPoint(0):
        return endCell.GetPoints()

    # Get number of cells
    points = []
    num_start = startCell.GetNumberOfPoints()
    num_end = endCell.GetNumberOfPoints()
    get_startCell = startCell.GetPoints()
    get_endCell = endCell.GetPoints()

    points = []
    n = 5
    N = 100
    num_centerline_points = 3
    
    for i in range(num_centerline_points -1, -1, -1):
        points.append(get_startCell.GetPoint(num_start - n*i - 1))
    
    if additionalPoint is not None:
        points.append(additionalPoint)

    for i in range(num_centerline_points):
        points.append(get_endCell.GetPoint(i*n))

    curv_coor = np.zeros(len(points))
    for i in range(len(points)-1):
        curv_coor[i+1] = curv_coor[i] + math.sqrt(distance(points[i], points[i+1]))

    # FIXME: Test the effect of adding an extra distance
    # Increse the longest distances
    #base_line = curv_coor[0]
    #for i, c in enumerate(curv_coor):
    #    if c > 5*base_line:
    #        print "Making distance longer, prev:", c, "  new:", 1.2*c
    #        curv_coor[i] = c*1.2

    points = np.asarray(points)
    
    fx = splrep(curv_coor, points[:,0], k=3)
    fy = splrep(curv_coor, points[:,1], k=3)
    fz = splrep(curv_coor, points[:,2], k=3)

    curv_coor = np.linspace(curv_coor[0], curv_coor[-1], N)
    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    tmp = []
    for i in range(num_start - n*num_centerline_points):
        tmp.append(get_startCell.GetPoint(i))

    for j in range(N):
        tmp.append([fx_[j], fy_[j], fz_[j]])

    for k in range(n*num_centerline_points, num_end):
        tmp.append(get_endCell.GetPoint(k))

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(tmp))
    for l in range(len(tmp)):
        points.SetPoint(l, tmp[l])

    return points
    

def InterpolateTwoCells(startCell, endCell, numberOfSplinePoints, additionalPointId,
                        additionalPoint, type):
    points = vtk.vtkPoints()
    xspline = vtk.vtkCardinalSpline()
    yspline = vtk.vtkCardinalSpline()
    zspline = vtk.vtkCardinalSpline()
   
    numberOfStartCellPoints = startCell.GetNumberOfPoints()
    numberOfEndCellPoints = endCell.GetNumberOfPoints()
    endCellFirstId = numberOfSplinePoints - numberOfEndCellPoints

    for i in range(numberOfStartCellPoints):
        point = startCell.GetPoints().GetPoint(i)
        xspline.AddPoint(float(i), point[0])
        yspline.AddPoint(float(i), point[1])
        zspline.AddPoint(float(i), point[2])

    if additionalPoint is not None:
        xspline.AddPoint(float(additionalPointId), additionalPoint[0])
        yspline.AddPoint(float(additionalPointId), additionalPoint[1])
        zspline.AddPoint(float(additionalPointId), additionalPoint[2])

    for i in range(numberOfEndCellPoints):
        point = endCell.GetPoints().GetPoint(i)
        index = float(endCellFirstId + i)
        xspline.AddPoint(index, point[0])
        yspline.AddPoint(index, point[1])
        zspline.AddPoint(index, point[2])

    xspline.Compute()
    yspline.Compute()  
    zspline.Compute()
   
    points.SetNumberOfPoints(numberOfSplinePoints)
    for i in range(numberOfSplinePoints):
        points.SetPoint(i,xspline.Evaluate(float(i)),yspline.Evaluate(float(i)),zspline.Evaluate(float(i)))
    return points
