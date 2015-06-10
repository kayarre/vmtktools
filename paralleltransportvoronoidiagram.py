#!/usr/bin/env python

import vtk
import sys
import math
from common import *
from vmtk import vtkvmtk


def ExtractCylindricInterpolationVoronoiDiagram(cellId,pointId,cylinderRadius,voronoi,centerlines):
    isInside = 0

    if (cellId == 0):
        cylinderTop = centerlines.GetPoint(pointId) 
        cylinderCenter = centerlines.GetPoint(pointId - interpolationHalfSize)
        cylinderBottom = centerlines.GetPoint(pointId - 2*interpolationHalfSize) 
    else:
        cylinderTop = centerlines.GetPoint(pointId) 
        cylinderCenter = centerlines.GetPoint(pointId + interpolationHalfSize)
        cylinderBottom = centerlines.GetPoint(pointId + 2*interpolationHalfSize) 

    interpolationDataset = vtk.vtkPolyData()
    interpolationDatasetPoints = vtk.vtkPoints()
    interpolationDatasetCellArray = vtk.vtkCellArray()

    maskArray = vtk.vtkIntArray()
    maskArray.SetNumberOfComponents(1)
    maskArray.SetNumberOfTuples(voronoi.GetNumberOfPoints())
    maskArray.FillComponent(0,0)

    for i in range(voronoi.GetNumberOfPoints()):
        point = voronoi.GetPoint(i)
        isInside = IsPointInsideInterpolationCylinder(point, cylinderTop, cylinderCenter, cylinderBottom, cylinderRadius)

        if (isInside == 1):
            maskArray.SetTuple1(i, 1)

    numberOfInterpolationPoints = ComputeNumberOfMaskedPoints(maskArray)

    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfInterpolationPoints)

    count = 0
    for i in range(voronoi.GetNumberOfPoints()):
        value = maskArray.GetTuple1(i)
        if (value == 1):
            interpolationDatasetPoints.InsertNextPoint(voronoi.GetPoint(i))
            interpolationDatasetCellArray.InsertNextCell(1)
            interpolationDatasetCellArray.InsertCellPoint(count)
            radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
            radiusArray.SetTuple1(count,radius)
            count +=1

    interpolationDataset.SetPoints(interpolationDatasetPoints)
    interpolationDataset.SetVerts(interpolationDatasetCellArray)
    interpolationDataset.GetPointData().AddArray(radiusArray)
    
    return interpolationDataset


def IsPointInsideInterpolationCylinder(x, t, c, b, r):
    halfheigth = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(b, t)) / 2.0

    xc = [x[i] - c[i] for i in range(len(x))]
    tb = [t[i] - b[i] for i in range(len(t))]

    xcnorm = vtk.vtkMath.Norm(xc)
    vtk.vtkMath.Normalize(xc)  
    vtk.vtkMath.Normalize(tb)  

    alpha = math.acos(vtk.vtkMath.Dot(xc, tb))

    parallelxc = xcnorm * math.cos(alpha) 
    perpendicularxc = xcnorm * math.sin(alpha)

    thetamin = math.atan(r / halfheigth)
    thetamax = thetamin + (math.pi - 2*thetamin)

    inside = 0
    if (thetamin <= alpha <= thetamax):
        if (abs(perpendicularxc)<=r):
            inside = 1
    else:
        if (abs(parallelxc)<=halfheigth):
            inside = 1

    return inside   


def ComputeNumberOfMaskedPoints(dataArray):
   numberOfPoints = 0
   for i  in range(dataArray.GetNumberOfTuples()):
      value = dataArray.GetTuple1(i)
      if (value ==1): numberOfPoints +=1
   return numberOfPoints


def VoronoiDiagramInterpolation(interpolationcellid, id0, id1, voronoiDataset0, 
                                voronoiDataset1, centerlines, step,
                                clippingPoints):
    cellLine = ExtractSingleLine(centerlines, interpolationcellid)

    startPoint = clippingPoints.GetPoint(id0)
    endPoint = clippingPoints.GetPoint(id1)

    startId = cellLine.FindPoint(startPoint)
    endId = cellLine.FindPoint(endPoint)
   
    gapStartId = startId + 1*step
    gapEndId = endId - 1*step
    arrivalId = gapEndId + 1*step
    endSavingInterval = gapEndId + 1*step

    numberOfGapPoints = int(math.fabs(gapEndId - gapStartId)) + 1
    numberOfInterpolationPoints = voronoiDataset0.GetNumberOfPoints()
    numberOfCenterlinesPoints = cellLine.GetNumberOfPoints()
    numberOfAddedPoints = numberOfGapPoints*numberOfInterpolationPoints

    finalNewVoronoiPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    finalRadiusArray = get_vtk_array(radiusArrayName, 1, numberOfAddedPoints)
    count = 0

    for i in range(numberOfInterpolationPoints):
       voronoiPoint = voronoiDataset0.GetPoint(i)
       voronoiPointRadius = voronoiDataset0.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

       centerlinePointLocator = get_locator(cellLine)

       closestPointId = centerlinePointLocator.FindClosestPoint(voronoiPoint)
       closestPoint = cellLine.GetPoint(closestPointId)
    
       voronoiVector = [0.0,0.0,0.0]
       voronoiVector[0] = voronoiPoint[0] - closestPoint[0]
       voronoiVector[1] = voronoiPoint[1] - closestPoint[1]
       voronoiVector[2] = voronoiPoint[2] - closestPoint[2]
       voronoiVectorNorm = vtk.vtkMath.Norm(voronoiVector)
       rotationAngle = ComputeVoronoiVectorToCenterlineAngle(closestPointId,voronoiVector,cellLine)

       PTPoints = vtk.vtkPoints()
       for j in range(closestPointId, arrivalId, step):
          localtangent = [0.0,0.0,0.0]
          localnormal = [0.0,0.0,0.0]
          newVoronoiVector = [0.0,0.0,0.0]
          newVoronoiPoint = [0.0,0.0,0.0]
          
          transform = vtk.vtkTransform()
          point0 = [0.0,0.0,0.0]
          point0 = cellLine.GetPoint(j)

          if (j < numberOfCenterlinesPoints - 1):
             point1 = [0.0,0.0,0.0]
             cellLine.GetPoint(j+1,point1)   
             localtangent[0] += point1[0] - point0[0]  
             localtangent[1] += point1[1] - point0[1]  
             localtangent[2] += point1[2] - point0[2]  
          if (j > 0):
             point2 = [0.0,0.0,0.0]
             cellLine.GetPoint(j-1,point2) 
             localtangent[0] += point0[0] - point2[0]
             localtangent[1] += point0[1] - point2[1]
             localtangent[2] += point0[2] - point2[2]

          localnormal = cellLine.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(j)
          localnormaldot = vtk.vtkMath.Dot(localtangent,localnormal)

          localtangent[0] -= localnormaldot*localnormal[0]
          localtangent[1] -= localnormaldot*localnormal[1]
          localtangent[2] -= localnormaldot*localnormal[2]
          vtk.vtkMath.Normalize(localtangent)

          transform.RotateWXYZ(rotationAngle,localtangent)
          transform.TransformNormal(localnormal,newVoronoiVector)
          vtk.vtkMath.Normalize(newVoronoiVector)

          newVoronoiPoint[0] = point0[0] + voronoiVectorNorm*newVoronoiVector[0]
          newVoronoiPoint[1] = point0[1] + voronoiVectorNorm*newVoronoiVector[1]
          newVoronoiPoint[2] = point0[2] + voronoiVectorNorm*newVoronoiVector[2]

          PTPoints.InsertNextPoint(newVoronoiPoint)

       numberOfPTPoints = PTPoints.GetNumberOfPoints() 
       lastPTPoint = PTPoints.GetPoint(PTPoints.GetNumberOfPoints()-1)
      
       voronoiPointLocator = get_locator(voronoiDataset1)
       
       arrivalVoronoiPointId = voronoiPointLocator.FindClosestPoint(lastPTPoint)
       arrivalVoronoiPoint = voronoiDataset1.GetPoint(arrivalVoronoiPointId)
       arrivalVoronoiPointRadius = voronoiDataset1.GetPointData().GetArray(radiusArrayName).GetTuple1(arrivalVoronoiPointId)
       
       arrivalCenterlinePointLocator = get_locator(cellLine)
      
       arrivalCenterlineClosestPointId = arrivalCenterlinePointLocator.FindClosestPoint(arrivalVoronoiPoint) 
       arrivalCenterlineClosestPoint = cellLine.GetPoint(arrivalCenterlineClosestPointId)
      
       arrivalVoronoiVector = [0.0,0.0,0.0]
       arrivalVoronoiVector[0] = arrivalVoronoiPoint[0] - arrivalCenterlineClosestPoint[0]  
       arrivalVoronoiVector[1] = arrivalVoronoiPoint[1] - arrivalCenterlineClosestPoint[1]  
       arrivalVoronoiVector[2] = arrivalVoronoiPoint[2] - arrivalCenterlineClosestPoint[2]  
       arrivalVoronoiVectorNorm = vtk.vtkMath.Norm(arrivalVoronoiVector)

       radiusArray = ComputeSpline(voronoiPointRadius,arrivalVoronoiPointRadius,numberOfPTPoints)
       vectorNormArray = ComputeSpline(voronoiVectorNorm,arrivalVoronoiVectorNorm, numberOfPTPoints)

       pointsToGap = (gapStartId - closestPointId)*step
         
       pointId = pointsToGap
       for k in range(gapStartId,endSavingInterval,step):
          ptpoint = PTPoints.GetPoint(pointsToGap)
          clpoint = cellLine.GetPoint(k)
    
          vector = [0.0,0.0,0.0]
          vector[0] = ptpoint[0] - clpoint[0]
          vector[1] = ptpoint[1] - clpoint[1]
          vector[2] = ptpoint[2] - clpoint[2]
          vtk.vtkMath.Normalize(vector)
          
          norm = vectorNormArray.GetTuple1(pointsToGap)

          newvector = [0.0,0.0,0.0]
          newvector[0] = norm*vector[0]
          newvector[1] = norm*vector[1]
          newvector[2] = norm*vector[2]
          
          newpoint = [0.0,0.0,0.0]
          newpoint[0] = clpoint[0] + newvector[0] 
          newpoint[1] = clpoint[1] + newvector[1] 
          newpoint[2] = clpoint[2] + newvector[2] 

          finalNewVoronoiPoints.InsertNextPoint(newpoint)
          cellArray.InsertNextCell(1)
          cellArray.InsertCellPoint(count)
          finalRadiusArray.SetTuple1(count,radiusArray.GetTuple1(pointsToGap))
          pointsToGap +=1
          count +=1
       
    return finalNewVoronoiPoints, finalRadiusArray


def ComputeVoronoiVectorToCenterlineAngle(pointId, vector, centerline):
    point0 = centerline.GetPoint(pointId)
    point1 = centerline.GetPoint(pointId + 1)
    point2 = centerline.GetPoint(pointId - 1)

    tangent = [0.0,0.0,0.0]
    for i in range(3):
        tangent[i] += point1[i]-point0[i]
        tangent[i] += point0[i]-point2[i]

    ptnnormal = centerline.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(pointId)
    alpha = ComputeAngleBetweenVectors(ptnnormal,tangent,vector)

    return alpha


def ComputeAngleBetweenVectors(normal, tangent, vector):
    #compute the tangent component orthogonal to normal 
    otangent = [0.0,0.0,0.0]
    normalDot = vtk.vtkMath.Dot(tangent, normal)
    otangent[0] = tangent[0] - normalDot*normal[0]
    otangent[1] = tangent[1] - normalDot*normal[1]
    otangent[2] = tangent[2] - normalDot*normal[2]
    vtk.vtkMath.Normalize(otangent)

    #compute the vector component orthogonal to otangent, i.e. parallel to normal
    vtk.vtkMath.Normalize(vector)
    ovector = [0.0,0.0,0.0]
    vectorDot = vtk.vtkMath.Dot(vector,otangent)
    ovector[0] = vector[0] - vectorDot*otangent[0]
    ovector[1] = vector[1] - vectorDot*otangent[1]
    ovector[2] = vector[2] - vectorDot*otangent[2]
    vtk.vtkMath.Normalize(ovector)
 
    theta = vtkvmtk.vtkvmtkMath.AngleBetweenNormals(normal, ovector)
    theta = vtk.vtkMath.DegreesFromRadians(theta)
   
    cross = [0.0,0.0,0.0]
    vtk.vtkMath.Cross(ovector, normal, cross)
    tangentDot = vtk.vtkMath.Dot(otangent, cross)

    if (tangentDot < 0.0): 
        theta = -1.0*theta
    
    angle = -theta
    
    return angle 


def ComputeSpline(startValue, endValue, numberOfPoints):
   averageValue = (startValue + endValue) / 2.0
   averageId = int(numberOfPoints / 2) 

   splineArray = vtk.vtkDoubleArray()
   splineArray.SetNumberOfComponents(1)
   splineArray.SetNumberOfTuples(numberOfPoints)
   splineArray.FillComponent(0,0.0)

   spline = vtk.vtkCardinalSpline()
   spline.AddPoint(float(0), startValue)
   spline.AddPoint(float(averageId), averageValue)
   spline.AddPoint(float(numberOfPoints), endValue)
   spline.Compute()   

   for i in range(numberOfPoints):
      splineArray.SetTuple1(i,spline.Evaluate(float(i)))

   return splineArray 


def InsertNewVoronoiPoints(oldDataset,newPoints,newArray):
   numberOfDatasetPoints = oldDataset.GetNumberOfPoints()
   numberOfNewPoints = newPoints.GetNumberOfPoints()
   numberOfNewVoronoiPoints = numberOfDatasetPoints + numberOfNewPoints

   newDataset = vtk.vtkPolyData()
   cellArray = vtk.vtkCellArray()  
   points = vtk.vtkPoints()
   radiusArray = get_vtk_array(radiusArrayName, 1, numberOfNewVoronoiPoints)

   for i in range(numberOfDatasetPoints):
      point = [0.0,0.0,0.0]
      oldDataset.GetPoint(i,point)
      points.InsertNextPoint(point)
      cellArray.InsertNextCell(1)
      cellArray.InsertCellPoint(i)

      value = oldDataset.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
      radiusArray.SetTuple1(i,value)

   for i in range(numberOfNewPoints):
      point = [0.0,0.0,0.0]
      newPoints.GetPoint(i,point)
      points.InsertNextPoint(point)
      cellArray.InsertNextCell(1)
      cellArray.InsertCellPoint(numberOfDatasetPoints+i)

      value = newArray.GetTuple1(i)
      radiusArray.SetTuple1(numberOfDatasetPoints+i,value)

   newDataset.SetPoints(points)
   newDataset.SetVerts(cellArray)
   newDataset.GetPointData().AddArray(radiusArray)

   return newDataset


def create_new_surface(completeVoronoiDiagram):
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    modeller.SetInput(completeVoronoiDiagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(polyBallImageSize)
    modeller.Update()

    # Write the new surface
    marchingCube = vtk.vtkMarchingCubes()
    marchingCube.SetInput(modeller.GetOutput())
    marchingCube.SetValue(0,0.0)
    marchingCube.Update()
    envelope = marchingCube.GetOutput()

    return envelope


def interpolate_voronoi_diagram(interpolatedCenterlines, patchCenterlines, 
                                clippedVoronoi, clippingPoints, bif):
    # Copy the voronoi diagram
    completeVoronoiDiagram = vtk.vtkPolyData()
    completeVoronoiDiagram.DeepCopy(clippedVoronoi)

    numberOfInterpolatedCenterlinesCells = interpolatedCenterlines.GetNumberOfCells()

    for j in range(1, 3): #numberOfInterpolatedCenterlinesCells + 1):
        interpolationCellId = j - 1
        startId = 0
        endId = j
        
        startCell = vtk.vtkGenericCell()
        patchCenterlines.GetCell(startId, startCell)

        startCellPointId = startCell.GetPointId(startCell.GetNumberOfPoints() - 1)
        startCellPoint = patchCenterlines.GetPoint(startCellPointId)
        startCellPointRadius = patchCenterlines.GetPointData().GetArray(radiusArrayName)\
                                                        .GetTuple1(startCellPointId)
        startCellPointHalfRadius = startCellPointRadius / 7.0
        
        startInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId, 
                                                    startCellPointId, startCellPointRadius, 
                                                    clippedVoronoi, patchCenterlines)
        startHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId, 
                                                    startCellPointId, startCellPointHalfRadius, 
                                                    clippedVoronoi, patchCenterlines)

        endCell = vtk.vtkGenericCell()
        patchCenterlines.GetCell(endId,endCell)

        endCellPointId = endCell.GetPointId(0)
        endCellPoint = patchCenterlines.GetPoint(endCellPointId)
        endCellPointRadius = patchCenterlines.GetPointData().GetArray(radiusArrayName)\
                                                        .GetTuple1(endCellPointId)
        endCellPointHalfRadius = endCellPointRadius / 7.0

        endInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId, endCellPointId, 
                                                    endCellPointRadius, clippedVoronoi,
                                                    patchCenterlines)
        endHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId, 
                                                    endCellPointId, endCellPointHalfRadius, 
                                                    clippedVoronoi, patchCenterlines)
   
        # Find and insert new points
        newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId,
                                                    startId, endId, startInterpolationDataset,
                                                    endHalfInterpolationDataset, 
                                                    interpolatedCenterlines, 1,
                                                    clippingPoints)
        completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints, 
                                                    newVoronoiPointsMISR)
        
        newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId, 
                                                    endId, startId, endInterpolationDataset, 
                                                    startHalfInterpolationDataset, 
                                                    interpolatedCenterlines, -1,
                                                    clippingPoints)
        completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints,
                                                    newVoronoiPointsMISR)

    if bif is not None:
        interpolationCellId = 0
        startId_ = 0
        endId_ = 1
        bif_clipped = bif[1]
        bif = bif[0]

        startCell = ExtractSingleLine(bif_clipped, startId_)
        locator = get_locator(startCell)
        startPoint = clippingPoints.GetPoint(1)
        startId = locator.FindClosestPoint(startPoint)
        startPoint = startCell.GetPoint(startId)
        startR = startCell.GetPointData().GetArray(radiusArrayName).GetTuple1(startId)
        startRHalf = startR / 7.0

        endCell = ExtractSingleLine(bif_clipped, endId_)
        locator = get_locator(endCell)
        endPoint = endCell.GetPoint(0)
        endId = locator.FindClosestPoint(endPoint)
        endPoint = endCell.GetPoint(endId)
        endR = endCell.GetPointData().GetArray(radiusArrayName).GetTuple1(endId)
        endRHalf = endR / 7.0

	startInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId_,
                                                    startId, startR,
                                                    clippedVoronoi, startCell)
        startHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId_, 
                                                    startId, startRHalf,
                                                    clippedVoronoi, startCell)

	endInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId_, endId, 
                                                    endR, clippedVoronoi,
                                                    endCell)
        endHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId_, endId,
                                                    endRHalf, clippedVoronoi,
                                                    endCell)

        newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId,
                                                    1, 2, startInterpolationDataset,
                                                    endHalfInterpolationDataset, 
                                                    bif, 1, clippingPoints)
        completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints, 
                                                    newVoronoiPointsMISR)

        newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId, 
                                                    2, 1, endInterpolationDataset, 
                                                    startHalfInterpolationDataset, 
                                                    bif, -1, clippingPoints)
        completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints,
                                                    newVoronoiPointsMISR)

    return completeVoronoiDiagram
