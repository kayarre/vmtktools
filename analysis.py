#!/usr/bin/env python 

import sys
import vtk
import scipy
import numpy
import math
import os
import matplotlib
from pylab import *

from vmtk import pypes, vtkvmtk
from numpy import zeros, array, linspace
from scipy.interpolate import splprep, splev 
from numpy import *


def Integration(line, arrayName):
   array = line.GetPointData().GetArray(arrayName)
   sum = 0.0
   for i in range(0, array.GetNumberOfTuples()):
      sum += pow(array.GetTuple1(i),2)


def Summation(line, arrayName):
   array = line.GetPointData().GetArray(arrayName)
   sum = 0.0
   for i in range(0, array.GetNumberOfTuples()):
      sum += abs(array.GetTuple1(i))
   return sum


def ClipSurfaceTopOrBottom(surface, origin, normal, flag):
    #if flag > 0:
    #    origin[2]=origin[2]-2.   # clip the bottom part
	
    plane = vtk.vtkPlane()
    plane.SetOrigin(origin)
    plane.SetNormal(normal)

    clipper = vtk.vtkClipPolyData()
    clipper.SetInput(surface)
    clipper.SetClipFunction(plane)

    if flag == 0:                      # clip the top part
        clipper.GenerateClippedOutputOn()

    clipper.GenerateClippedOutputOn()
    clipper.GetClippedOutput()

    if flag == 0:
        clipper.InsideOutOn()

    clipper.Update()
    cutSurface = clipper.GetClippedOutput()

    stripper = vtk.vtkStripper()
    stripper.SetInput(cutSurface)
    stripper.Update()

    connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
    connectivityFilter.SetInput(stripper.GetOutput())
    connectivityFilter.SetClosestPoint(origin)
    connectivityFilter.SetExtractionModeToClosestPointRegion()
    connectivityFilter.Update()

    contour = connectivityFilter.GetOutput()

    return cutSurface


def ReorderPoints(clippedsurf, origin):
    # this filter reorder the points
    stripper = vtk.vtkStripper()
    stripper.SetInput(clippedsurf.GetOutput())
    stripper.Update()

    connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
    connectivityFilter.SetInput(stripper.GetOutput())
    connectivityFilter.SetClosestPoint(origin)
    connectivityFilter.SetExtractionModeToClosestPointRegion()
    connectivityFilter.Update()

    contour = connectivityFilter.GetOutput()

    return contour


def CutWithPlane(input, origin, normal):
   plane = vtk.vtkPlane()
   plane.SetOrigin(origin)
   plane.SetNormal(normal)

   cutEdges = vtk.vtkCutter()
   cutEdges.SetInput(input)
   cutEdges.SetCutFunction(plane)
   cutEdges.GenerateCutScalarsOn()
   cutEdges.SetValue(0,0.0)

   cutStrips = vtk.vtkStripper()
   cutStrips.SetInputConnection(cutEdges.GetOutputPort())
   cutStrips.Update()
   cutPoly = vtk.vtkPolyData()
   cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
   cutPoly.SetPolys(cutStrips.GetOutput().GetLines())

   cutTriangles = vtk.vtkTriangleFilter()
   cutTriangles.SetInput(cutPoly)

   section2 = cutTriangles.GetOutput()
   return section2
    

def SplineCenterline(centerline, resampledCl, numberOfOutputPoints):
   inputPoints = resampledCl.GetPoints()
   aSplineX = vtk.vtkCardinalSpline()
   aSplineY = vtk.vtkCardinalSpline()
   aSplineZ = vtk.vtkCardinalSpline()   

   for i in range(0, inputPoints.GetNumberOfPoints()):
       point=[0.0, 0.0, 0.0]
       inputPoints.GetPoint(i, point)
       aSplineX.AddPoint(i, point[0])
       aSplineY.AddPoint(i, point[1])
       aSplineZ.AddPoint(i, point[2])

   # Generate the polyline for the spline.
   splinePoints = vtk.vtkPoints()
   splineVtp = vtk.vtkPolyData()
   splineSphereArray = vtk.vtkDoubleArray()
   
   # Interpolate x, y and z by using the three spline filters and
   # create new points
   for i in range(0, numberOfOutputPoints):
       t = (inputPoints.GetNumberOfPoints()-1.0)/(numberOfOutputPoints-1.0)*i
       splinePoints.InsertPoint(i, aSplineX.Evaluate(t), aSplineY.Evaluate(t),
                          aSplineZ.Evaluate(t))
    
   # Create the polyline.
   linesSpline = vtk.vtkCellArray()
   linesSpline.InsertNextCell(numberOfOutputPoints)
   print numberOfOutputPoints
   for i in range(0, numberOfOutputPoints):
       linesSpline.InsertCellPoint(i)
    
   # spline the radius array  
   radiusArray=centerline.GetPointData.GetArray('MaximumInscribedSphereRadius')
   radiusArrayResampled=vtk.vtkDoubleArray()
   sphereSpline=vtk.vtkCardinalSpline()

   j = 0
   stepsOfInputSphereSpline = 9
   for i in range(0, radiusArray.GetNumberOfTuples(), stepsOfInputSphereSpline):
      radius = radiusArray.GetTuple1(i)
      sphereSpline.AddPoint(j, radius)
      j += 1

   numberOfInputSpheres = j

   # Generate poly line for spline for the radius
   splineSphereArray.SetName('MaximumInscribedSphereRadius')
   splineSphereArray.SetNumberOfComponents(1)
   splineSphereArray.SetNumberOfTuples(numberOfOutputPoints)

   # Create new spheres
   print numberOfOutputPoints
   for i in range (0,numberOfOutputPoints):
       t = (numberOfInputSpheres-1.0)/(numberOfOutputPoints-1.0)*i
       splineSphereArray.SetValue(i,sphereSpline.Evaluate(t))

   splineVtp.SetPoints(splinePoints)  
   splineVtp.SetLines(linesSpline)          
   splineVtp.GetPointData().AddArray(splineSphereArray)

   return splineVtp


def AverageCenterline2(outdirectory,apexCl,icaCl,ecaCl):
    """This method of averaging centerine is the intilligent one, tfor finding the
       closest point on the centerline for each point of the other centerline , the
       cut plane cuts the first centerline first"""
   from numpy import zeros,array, linspace
   icaPointLocator=vtk.vtkPointLocator()
   ecaPointLocator=vtk.vtkPointLocator()
   
   averageCell=vtk.vtkCellArray()
   averageCell.InsertNextCell(apexCl.GetNumberOfPoints()-1)
   averageCl=vtk.vtkPolyData()
   averagePoints=vtk.vtkPoints()
   
   ecaPointLocator.SetDataSet(ecaCl)
    
   for i in range(1, apexCl.GetNumberOfPoints()):

      ecaPoint=[0.0,0.0,0.0]
      icaPoint=[0.0,0.0,0.0]
      averagePoint=[0.0,0.0,0.0]
      apexPoint=[0.0,0.0,0.0]
      
      apexCl.GetPoints().GetPoint(i,apexPoint)
      point1=numpy.zeros(3)
      point2=numpy.zeros(3)
     
      apexCl.GetPoints().GetPoint(i,point1)
      apexCl.GetPoints().GetPoint(i-1,point2)

      #better way fo tind icaPoin
      infofile = file(outdirectory+'_ref.dat',"r")
      infolines=infofile.readlines()
      infoLine=infolines[1]
      bifurcationCenter=numpy.zeros((3,1))
      bifurcationCenter[0]=float(infoLine.split()[0]);bifurcationCenter[1]=float(infoLine.split()[1]);bifurcationCenter[2]=float(infoLine.split()[2])
      bifurcationUpNormal=numpy.zeros((3,1))
      bifurcationUpNormal[0]=float(infoLine.split()[6]);bifurcationUpNormal[1]=float(infoLine.split()[7]);bifurcationUpNormal[2]=float(infoLine.split()[8])

      cutPlane=vtk.vtkPlane()
      cutPlane.SetOrigin(point1)       
      cutPlane.SetNormal(bifurcationUpNormal)       
      clipper=vtk.vtkClipPolyData()
      #clipper.SetInputConnection(apexCl)
      clipper.SetInput(icaCl)
      clipper.SetClipFunction(cutPlane)
      clipper.GenerateClipScalarsOff()
      clipper.GenerateClippedOutputOff();
      clipper.InsideOutOff()
      xxx = vtk.vtkPolyData()
      WritePolyData(clipper.GetOutput(),outdirectory+'_xx.vtp')
      xxx = ReadPolyData(outdirectory+'_xx.vtp')       
      icaPointLocator.SetDataSet(xxx)
      icaPointId=icaPointLocator.FindClosestPoint(apexPoint)
      ecaPointId=ecaPointLocator.FindClosestPoint(apexPoint)
      
      xxx.GetPoints().GetPoint(icaPointId,icaPoint)
      ecaCl.GetPoints().GetPoint(ecaPointId,ecaPoint)

      averagePoint[0]=(icaPoint[0]+ecaPoint[0])/2
      averagePoint[1]=(icaPoint[1]+ecaPoint[1])/2
      averagePoint[2]=(icaPoint[2]+ecaPoint[2])/2
      averagePoints.InsertNextPoint(averagePoint[0],averagePoint[1],averagePoint[2])
      averageCell.InsertCellPoint(i-1)

   averageCl.SetPoints(averagePoints)
   averageCl.SetLines(averageCell)
   return averageCl


def Triangle(point0, point1, point2):
   triangles=vtk.vtkCellArray()
   points = vtk.vtkPoints()
   points.InsertNextPoint(point0)
   points.InsertNextPoint(point1)
   points.InsertNextPoint(point2)
    
   triangle = vtk.vtkTriangle()
   triangle.GetPointIds().SetId(0,0)
   triangle.GetPointIds().SetId(1,1)
   triangle.GetPointIds().SetId(2,2)
   triangles.InsertNextCell(triangle)

#   point0=[0,0,0];point1=[1,0,0];point2=[0,1,0]

   area=triangle.TriangleArea(point0,point1,point2)
   return area 


def translate(circle,rotatAngle,meanPoint):
   trnasformFilter3=vtk.vtkTransformFilter()

   rotationF=vtk.vtkTransform()
   rotationF.RotateWXYZ(rotatAngle,[0,1,0])
   rotationF.Translate(meanPoint)

   trnasformFilter3.SetInput(circle)
   trnasformFilter3.SetTransform(rotationF)

   trnasformFilter3.Update()
   projClTrans3=trnasformFilter3.GetOutput()

   #trnasformFilter2=vtk.vtkTransformFilter()
   #trnasformFilter2.SetInput(projClTrans3)
   #trnasformFilter2.SetTransform(rotation)
   #trnasformFilter2.Update()
   #projClTrans2=trnasformFilter2.GetOutput()
   return projClTrans3
   #WritePolyData(projClTrans,'FOOO.vtp')


def ReadPolyData(filename):
   reader=vtk.vtkXMLPolyDataReader()
   reader.SetFileName(filename)
   reader.Update()                            
   output=reader.GetOutput()
   return output


def WritePolyData(input,filename):
   writer = vtk.vtkXMLPolyDataWriter()
   writer.SetFileName(filename)
   writer.SetInput(input)
   writer.Write()


def ExtractIcaCenterline(centerline,branch):
    centerlineComplete=vtk.vtkPolyData()
    icaCcaRadiusArray=vtk.vtkDoubleArray()
    icaCca=vtk.vtkGenericCell()
    icaCcaPoints=vtk.vtkPoints()
    
    centerlineComplete=centerline
    centerlineComplete.GetCell(branch,icaCca)
    icaCcaNumberOfPoints=icaCca.GetNumberOfPoints()
    icaCcaPointIds=vtk.vtkIdList()
    centerlineComplete.GetCellPoints(branch,icaCcaPointIds)
    icaCcaRadiusArray.SetName('MaximumInscribedSphereRadius')
    icaCcaRadiusArray.SetNumberOfComponents(1)
    icaCcaRadiusArray.SetNumberOfTuples(icaCcaNumberOfPoints)

    icaCcaRadiusArray.FillComponent(0,1.0)
    sphereArray=centerlineComplete.GetPointData().GetArray('MaximumInscribedSphereRadius')
    #extracting the sphere radius of icaCca + the point coordinates

    for i in range(0,icaCcaPointIds.GetNumberOfIds()):
       sphereId=icaCcaPointIds.GetId(i)
       sphere=sphereArray.GetTuple1(sphereId)
       icaCcaRadiusArray.SetTuple1(i,sphere)
          
       point=[0.0,0.0,0.0]
       icaCca.GetPoints().GetPoint(i,point)
       icaCcaPoints.InsertNextPoint(point[0],point[1],point[2])
 
    #creating the icaCca polyData
    icaCcaLine=vtk.vtkCellArray()
    icaCcaCenterline=vtk.vtkPolyData()
    icaCcaLine.InsertNextCell(icaCcaNumberOfPoints)

    for i in range(0,icaCcaNumberOfPoints):
       icaCcaLine.InsertCellPoint(i)
    
    icaCcaCenterline.SetPoints(icaCcaPoints)
    icaCcaCenterline.SetLines(icaCcaLine)
    icaCcaCenterline.GetPointData().AddArray(icaCcaRadiusArray)
    return icaCcaCenterline
      

def FindPoint(line,point):               
    pointId=line.FindPoint(point)
    return pointId


def ResamplingCenterline(inputLine, numberOfInputPoints, numberOfOutputPoints):
   #resampling
   cellArray=vtk.vtkCellArray()
   points=vtk.vtkPoints()
   resampledPoly=vtk.vtkPolyData()
   j=0
   samplingFactor=int(inputLine.GetNumberOfPoints()/numberOfInputPoints)
   cellArray.InsertNextCell(numberOfInputPoints)
   for i in range(0,inputLine.GetNumberOfPoints(),samplingFactor):
      point=[0.0,0.0,0.0]
      inputLine.GetPoints().GetPoint(i,point)
      points.InsertNextPoint(point[0],point[1],point[2]) 
      cellArray.InsertCellPoint(j)
      j+=1
   resampledPoly.SetPoints(points)
   resampledPoly.SetLines(cellArray)
   
   return resampledPoly


def CenterlineAttribiute(clFileName, clAttribiuteFileName, length, factor):
   command='vmtkcenterlineresampling -ifile '+clFileName+' -length '+length+' --pipe vmtkcenterlineattributes --pipe vmtkcenterlinegeometry -ofile '+clAttribiuteFileName+'  -smoothing 1 -factor '+factor+' -outputsmoothed 1 -iterations 100 '
   os.system(command)


def tortuosity(attribute):
    frenetNormal=attribute.GetPointData().GetArray('FrenetNormal')
    paralelTransportNormal=attribute.GetPointData().GetArray('ParallelTransportNormals')

    #creating new array for tortuosity
    tortuosity=vtk.vtkDoubleArray() 
    tortuosity.SetName('tortuosityArray')
    tortuosity.SetNumberOfComponents(1) 
    tortuosity.SetNumberOfTuples(attribute.GetNumberOfPoints())
    tortuosityAngle=0

    for i in range (0,attribute.GetNumberOfPoints()):
        xfrenetNormal=frenetNormal.GetComponent(i,0)
        yfrenetNormal=frenetNormal.GetComponent(i,1)
        zfrenetNormal=frenetNormal.GetComponent(i,2)
      
        xparallelTransportNormal=paralelTransportNormal.GetComponent(i,0)
        yparallelTransportNormal=paralelTransportNormal.GetComponent(i,1)
        zparallelTransportNormal=paralelTransportNormal.GetComponent(i,2)
      
        xTortuosity=xparallelTransportNormal-xfrenetNormal
        yTortuosity=yparallelTransportNormal-yfrenetNormal
        zTortuosity=zparallelTransportNormal-zfrenetNormal

        tortuosityAngle=numpy.arcsin(yfrenetNormal/numpy.sqrt(xfrenetNormal*xfrenetNormal+yfrenetNormal*yfrenetNormal))
        if xfrenetNormal>0 and yfrenetNormal >0:
            angle=numpy.degrees(tortuosityAngle)
                      
        if xfrenetNormal<0 and yfrenetNormal >0:
            angle=180-numpy.degrees(tortuosityAngle)

        if  xfrenetNormal<0 and yfrenetNormal <0:
            angle=180-numpy.degrees(tortuosityAngle)
       
        if  xfrenetNormal>0 and yfrenetNormal <0:
            angle=360+numpy.degrees(tortuosityAngle)

      #tortuosityAngle=numpy.arccos(xfrenetNormal/numpy.sqrt(xfrenetNormal*xfrenetNormal+yfrenetNormal*yfrenetNormal))
      #tortuosityAngle=(tortuosityAngle/math.pi)*180
      #tortuosity.SetTuple1(i,tortuosityAngle)
      tortuosity.SetTuple1(i,angle)
    
    attribute.GetPointData().AddArray(tortuosity)
    return attribute


def Teta(attribute):
    frenetNormal = attribute.GetPointData().GetArray('FrenetNormal')
    paralelTransportNormal = attribute.GetPointData().GetArray('ParallelTransportNormals')

    tortuosity = vtk.vtkDoubleArray()
    tortuosity.SetName('Teta')
    tortuosity.SetNumberOfComponents(1) 
    tortuosity.SetNumberOfTuples(attribute.GetNumberOfPoints())
    tortuosityAngle = 0

    for i in range (0,attribute.GetNumberOfPoints()):
      xfrenetNormal = frenetNormal.GetComponent(i, 0)
      yfrenetNormal = frenetNormal.GetComponent(i, 1)
      zfrenetNormal = frenetNormal.GetComponent(i, 2)
      
      xparallelTransportNormal = paralelTransportNormal.GetComponent(i, 0)
      yparallelTransportNormal = paralelTransportNormal.GetComponent(i, 1)
      zparallelTransportNormal = paralelTransportNormal.GetComponent(i, 2)
      
      xTortuosity = xparallelTransportNormal - xfrenetNormal
      yTortuosity = yparallelTransportNormal - yfrenetNormal
      zTortuosity = zparallelTransportNormal - zfrenetNormal
      tortuosityAngle = numpy.arctan(zparallelTransportNormal / \
                        numpy.sqrt((xparallelTransportNormal * xparallelTransportNormal + \
                                    yparallelTransportNormal * yparallelTransportNormal + \
                                    zparallelTransportNormal * zparallelTransportNormal)))
      tortuosityAngle=(tortuosityAngle/math.pi)*180
      tortuosity.SetTuple1(i,tortuosityAngle)
    
    attribute.GetPointData().AddArray(tortuosity)
    return attribute


def createPlane(normal,center):
   plane=vtk.vtkPlaneSource()
   plane.SetOrigin([0,0,0])
   #plane.SetPoint1(center[0]-10,center[1]-10,0)
   plane.SetPoint1(18,0,0)
   #plane.SetPoint2(center[0]-10,center[1]+10,0)
   plane.SetPoint2(0,35,0)
  
   plane.SetNormal(normal)
   plane.SetCenter(center)
  
   return plane                                                                   


def BestFitPlane(points):
   [eigValue,eigVec]=numpy.linalg.eig(points)
   normalVec=eigVec[:,2]
   #print normalVec      


def initialCenterFinder(cl):
   p1=[0.0,0.0,0.0]
   p2=[0.0,0.0,0.0]
   p3=[0.0,0.0,0.0]

   # first point and last point
   p1=cl.GetPoints().GetPoint(0)
   p3=cl.GetPoints().GetPoint(cl.GetNumberOfPoints()-1)

   # the average point is the middle point
   point=[0.0,0.0,0.0]
   pointArray=numpy.zeros((cl.GetNumberOfPoints(),3))
   for i in range(0,cl.GetNumberOfPoints()):
   
      cl.GetPoints().GetPoint(i,point)
      pointArray[i]=point
   pointArray=pointArray.T
   p2[0]=pointArray[0,:].mean()
   p2[1]=pointArray[1,:].mean()
   p2[2]=pointArray[2,:].mean()

   # becasue cl is projected on y-z plane
   p1=[p1[1],p1[2]]
   p2=[p2[1],p2[2]]
   p3=[p3[1],p3[2]]                                
   ma=(p2[1]-p1[1])/(p2[0]-p1[0])
   mb=(p3[1]-p2[1])/(p3[0]-p2[0])

   xc=(ma*mb*(p1[1]-p3[1])+mb*(p1[0]+p2[0])-ma*(p2[0]+p3[0]))/(2*(mb-ma))

   #yc=(xc-0.5*(p1[0]+p2[0]))/ma + 0.5*(p1[1]+p2[1])
   yc2=(-(xc-((p2[0]+p3[0])*0.5))/mb)+(0.5*(p2[1]+p3[1]))

   center=numpy.array((xc,yc2))
   return center


def initialCircle2D(cl):
   #giving a projected cl on the y-z plane this module extract three points and
   #measur the circle R,Center passing these three points. this is used as the
   #initial guess for best fit circle module.
   p1=[0.0,0.0,0.0]
   p2=[0.0,0.0,0.0]
   p3=[0.0,0.0,0.0]

   #first point and last point
   p1=cl.GetPoints().GetPoint(0)
   p3=cl.GetPoints().GetPoint(cl.GetNumberOfPoints()-1)

   #the average point is the middle point
   point=[0.0,0.0,0.0]
   pointArray=numpy.zeros((cl.GetNumberOfPoints(),3))
   for i in range(0,cl.GetNumberOfPoints()):
      cl.GetPoints().GetPoint(i,point)
      pointArray[i]=point

   pointArray=pointArray.T
   p2[0]=pointArray[0,:].mean()
   p2[1]=pointArray[1,:].mean()
   p2[2]=pointArray[2,:].mean()

   #becasue cl is projected on y-z plane
   p1=[p1[1], p1[2]]
   p2=[p2[1], p2[2]]
   p3=[p3[1], p3[2]]                                
   ma=(p2[1] - p1[1]) / (p2[0] - p1[0])
   mb=(p3[1] - p2[1]) / (p3[0] - p2[0])

   xc=(ma*mb*(p1[1]-p3[1])+mb*(p1[0]+p2[0])-ma*(p2[0]+p3[0]))/(2*(mb-ma))

   #yc=(xc-0.5*(p1[0]+p2[0]))/ma + 0.5*(p1[1]+p2[1])
   yc2=(-(xc-((p2[0]+p3[0])*0.5))/mb)+(0.5*(p2[1]+p3[1]))

   center = numpy.array((xc, yc2))
   b = numpy.array((p1[0], p1[1]))
   radius = numpy.linalg.norm(center - b)
   #print 'radius =',radius,' center= ', center

   return radius,center


#creating the 2d circle   
def Circle2D(radius,center):
   circle=vtk.vtkRegularPolygonSource()
   circle.SetNormal([1, 0, 0])
   circle.GeneratePolygonOff()
   circle.SetNumberOfSides(360) 
   circle.SetCenter(center)
   circle.SetRadius(radius)
   circle.Update()
   circlePoly=circle.GetOutput()
   return circlePoly


# this module search the centerline between a ccaId and icaId with the
# calculated tortuosity and return the start of the twisting.
def FindTwistPoint2(centerline,ccaId,icaId):

   ccaPoint=[0.0,0.0,0.0]
   icaPoint=[0.0,0.0,0.0]

   point=[0.0,0.0,0.0]
   points=numpy.zeros(3*numpy.abs((icaId-ccaId+1)))
   points=points.reshape(numpy.abs(icaId-ccaId+1),3)                                                                 

   tortuosityArray=numpy.zeros(numpy.abs(icaId-ccaId+1))
   tortuosityArray.reshape(icaId-ccaId+1)

   for i in range(ccaId,icaId+1):
      centerline.GetPoints().GetPoint(ccaId,ccaPoint)
      centerline.GetPoints().GetPoint(icaId,icaPoint)
      centerline.GetPoints().GetPoint(i,point)
      points[i-ccaId]=point
      tortuosityArray[i-ccaId]=centerline.GetPointData().GetArray('tortuosityArray').GetComponent(i,1)
   index=[]
   value=[]
   index.append(ccaId+1)
   value0=tortuosityArray[ccaId+1]
   value.append(value0)
   ignoreL=10

   for i in range (ccaId+1,icaId-1):
      if ((tortuosityArray[i-ccaId]-tortuosityArray[i-ccaId-1])*(-tortuosityArray[i-ccaId]+tortuosityArray[i-ccaId+1])<0) and (i-index[len(index)-1])>ignoreL:

      #if ((tortuosityArray[i-ccaId]-tortuosityArray[i-ccaId-1])*(-tortuosityArray[i-ccaId]+tortuosityArray[i-ccaId+1])<0):
         index.append(i)
         value.append(tortuosityArray[i-ccaId])
   print index
   print value

   length=[]
   ignoreL=10
   index2=[]
   for i in range (0,len(index)):
      if  value[i]<350 and value[i]>10:
         index2.append(index[i])      

   for i in range (0,len(index2)-1):
      length.append(int(index2[i+1]-index2[i]))
   max=0.0;start=0.0
   for i in range(0,len(length)):
      if length[i]>max:
         max=length[i]
         start=index2[i]
   print index                    
   print index2
   print length
   print 'start, llength', start, max
   return start,int(max)


def FindTwistPoint(centerline, ccaId, icaId):
   ccaPoint=[0.0, 0.0, 0.0]
   icaPoint=[0.0, 0.0, 0.0]

   point=[0.0, 0.0, 0.0]
   points=numpy.zeros(3 * numpy.abs((icaId - ccaId + 1)))
   points=points.reshape(numpy.abs(icaId - ccaId + 1), 3)

   tortuosityArray=numpy.zeros(numpy.abs(icaId-ccaId+1))
   tortuosityArray.reshape(icaId-ccaId+1)

   for i in range(ccaId,icaId+1):
      centerline.GetPoints().GetPoint(ccaId,ccaPoint)
      centerline.GetPoints().GetPoint(icaId,icaPoint)
      centerline.GetPoints().GetPoint(i,point)
      points[i-ccaId]=point
      tortuosityArray[i-ccaId]=centerline.GetPointData().GetArray('tortuosityArray').GetComponent(i,1)

   tortuosityDiff = numpy.ediff1d(tortuosityArray)
   x = tortuosityDiff

   endIndex = 0
   startIndex = ccaId
   lengthArray = numpy.zeros(numpy.abs(icaId-ccaId))
   endIndexArray = numpy.zeros(numpy.abs(icaId-ccaId))
   valueArray = numpy.zeros(numpy.abs(icaId-ccaId))
   infoArray = numpy.zeros(3*(numpy.abs(icaId-ccaId)))
   infoArray = infoArray.reshape(icaId-ccaId,3)
   realLength = 0
   ignorLength = 20     # any bump less than this is ognored
   ignoreAngle1 = 10    # if it goes less than this, we assume it reaches the -90  
   ignoreAngle2 = 345   # if it goes more than this, we assume it reaches the 90
   for i in range(ccaId+1,icaId):
            
      if  (x[i-ccaId]*x[i-ccaId-1] <= 0 ):
        
         endIndex=i
         length=(endIndex-startIndex)
         endValue=tortuosityArray[i-ccaId]
         startValue=tortuosityArray[startIndex-ccaId]

         if ((length) >  ignorLength and endValue < ignoreAngle1 ):
            realLength=length

         if  ((length) >  ignorLength and startValue <ignoreAngle1 ):
            length=length+realLength
            realLength=0             

         if length<=8.0:
            length=0
            
         lengthArray[i-ccaId]=length
         endIndexArray[i-ccaId]=endIndex
         valueArray[i-ccaId]=endValue
         infoArray[i-ccaId]=[endIndex,length,endValue]
         startIndex=endIndex

                                      
   maxIndexes=numpy.argmax(infoArray,0)
   maxLengthIndex=maxIndexes[1]
   maxLength=infoArray[maxLengthIndex,1]
   twistPointIndex=infoArray[maxLengthIndex,0]-maxLength
   print  'max length =',maxLength,'twist point =',int(twistPointIndex)
   return twistPointIndex,int(maxLength)


def BranchAngle2(surfaceFilename, centerlineFilename, outdirectory):
   command='/das/pbijari/mainCodes/lucalizer.py'+ ' '+surfaceFilename+' '+centerlineFilename+' '+outdirectory  
   os.system(command)
 

def BranchAngle(ICAECA0, ICAECA1, BifurcationFile, CCA0, CCA1):
   from numpy import * 
   from scipy import *
   from math import *                            

   ICAECA0Lines = ICAECA0.readlines()                                                   
   ICA0Line=ICAECA0Lines[2]
   
   ICA0x=float(ICA0Line.split()[5])
   ICA0y=float(ICA0Line.split()[6])
   ICA0z=float(ICA0Line.split()[7])
   
   ECA0Line=ICAECA0Lines[3]
   ECA0x=float(ECA0Line.split()[5])
   ECA0y=float(ECA0Line.split()[6])
   ECA0z=float(ECA0Line.split()[7])
   
   ICAECA1Lines = ICAECA1.readlines()
   ICA1Line=ICAECA1Lines[2]
   print ICA1Line
   ICA1x=float(ICA1Line.split()[5])
   ICA1y=float(ICA1Line.split()[6])
   ICA1z=float(ICA1Line.split()[7])
   
   ECA1Line=ICAECA1Lines[3]
   ECA1x=float(ECA1Line.split()[5])
   ECA1y=float(ECA1Line.split()[6])
   ECA1z=float(ECA1Line.split()[7])
   
   CCA0Lines = CCA0.readlines()
   CCA0Line=CCA0Lines[1]
   
   CCA0x=float(CCA0Line.split()[5])
   CCA0y=float(CCA0Line.split()[6])
   CCA0z=float(CCA0Line.split()[7])
   
   CCA1Lines = CCA1.readlines()
   CCA1Line=CCA1Lines[1]
   
   CCA1x=float(CCA1Line.split()[5])
   CCA1y=float(CCA1Line.split()[6])
   CCA1z=float(CCA1Line.split()[7])
   
   BifurcationLines=BifurcationFile.readlines()
   BifurcationLine=BifurcationLines[1]
   BifurcationOrigionx=float(BifurcationLine.split()[0])
   BifurcationOrigiony=float(BifurcationLine.split()[1])
   BifurcationOrigionz=float(BifurcationLine.split()[2])
   BifurcationNormalx=float(BifurcationLine.split()[3])
   BifurcationNormaly=float(BifurcationLine.split()[4])
   BifurcationNormalz=float(BifurcationLine.split()[5])
   BifurcationUpnormalx=float(BifurcationLine.split()[6])
   BifurcationUpnormaly=float(BifurcationLine.split()[7])
   BifurcationUpnormalz=float(BifurcationLine.split()[8])
   
   icaOrientationx=ICA1x-ICA0x
   icaOrientationy=ICA1y-ICA0y
   icaOrientationz=ICA1z-ICA0z
   
   ecaOrientationx=ECA1x-ECA0x
   ecaOrientationy=ECA1y-ECA0y
   ecaOrientationz=ECA1z-ECA0z
   
   ccaOrientationx=CCA0x-CCA1x # be careful this vector is from CCA1 to CCA0
   ccaOrientationy=CCA0y-CCA1y
   ccaOrientationz=CCA0z-CCA1z
   
   #ccaOrientationx=CCA1x-CCA0x # be careful this vector is from CCA1 to CCA0
   #ccaOrientationy=CCA1y-CCA0y
   #ccaOrientationz=CCA1z-CCA0z
   
   icaOrientationProjectionx=icaOrientationx
   icaOrientationProjectiony=icaOrientationy
   icaOrientationProjectionz=0
   
   icaOrientationProjectionOutx=0
   icaOrientationProjectionOuty=icaOrientationy
   icaOrientationProjectionOutz=icaOrientationz
   
   ecaOrientationProjectionx=ecaOrientationx
   ecaOrientationProjectiony=ecaOrientationy
   ecaOrientationProjectionz=0
   
   ccaOrientationProjectionx=ccaOrientationx
   ccaOrientationProjectiony=ccaOrientationy
   ccaOrientationProjectionz=0
      
   ccaOrientationProjectionOutx=0
   ccaOrientationProjectionOuty=ccaOrientationy
   ccaOrientationProjectionOutz=ccaOrientationz
   
   internalMultiply=(icaOrientationProjectionx*ecaOrientationProjectionx)+(icaOrientationProjectiony*ecaOrientationProjectiony)+(icaOrientationProjectionz*ecaOrientationProjectionz)
   
   icaOrientationProjectionMagnitude=sqrt(icaOrientationProjectionx**2+icaOrientationProjectiony**2+icaOrientationProjectionz**2)
   ecaOrientationProjectionMagnitude=sqrt(ecaOrientationProjectionx**2+ecaOrientationProjectiony**2+ecaOrientationProjectionz**2)
   
   cosBranchAngle=internalMultiply/(icaOrientationProjectionMagnitude*ecaOrientationProjectionMagnitude)
   branchAngleRadian=acos(cosBranchAngle)
   branchAngleDegree=branchAngleRadian*180/pi
   
   internalMultiplyCCAICA=(ccaOrientationProjectionx*icaOrientationProjectionx)+(ccaOrientationProjectiony*icaOrientationProjectiony)+(ccaOrientationProjectionz*icaOrientationProjectionz)
   ccaOrientationProjectionMagnitude=sqrt(ccaOrientationProjectionx**2+ccaOrientationProjectiony**2+ccaOrientationProjectionz**2)
   
   cosCCAICAangle=internalMultiplyCCAICA/(ccaOrientationProjectionMagnitude*icaOrientationProjectionMagnitude)
   CCAICAAngleRadian=acos(cosCCAICAangle)
   CCAICAAngleDegree=CCAICAAngleRadian*180/pi
   
   internalMultiplyCCAICAOutPlane = (ccaOrientationProjectionOutx * icaOrientationProjectionOutx) + \
                                    (ccaOrientationProjectionOuty * icaOrientationProjectionOuty) + \
                                    (ccaOrientationProjectionOutz * icaOrientationProjectionOutz)
   ccaOrientationProjectionOutMagnitude = sqrt(ccaOrientationProjectionOutx**2 + ccaOrientationProjectionOuty**2 + ccaOrientationProjectionOutz**2)
   icaOrientationProjectionOutMagnitude = sqrt(icaOrientationProjectionOutx**2 + icaOrientationProjectionOuty**2 + icaOrientationProjectionOutz**2)
   
   cosCCAICOutPlaneAangle=internalMultiplyCCAICAOutPlane/(ccaOrientationProjectionOutMagnitude*icaOrientationProjectionOutMagnitude)
   CCAICAOutPlaneAngleRadian=acos(cosCCAICOutPlaneAangle)
   CCAICAOutPlaneAngleDegree=CCAICAOutPlaneAngleRadian*180/pi
   return    branchAngleDegree,CCAICAAngleDegree,CCAICAOutPlaneAngleDegree
   outputfile=open(sys.argv[2],'a')
   outputfile.write('%s %f %f %f '%(modelname, branchAngleDegree,CCAICAAngleDegree,CCAICAOutPlaneAngleDegree))


def AverageCenterline(centerline,apexCl,icaCl,ecaCl):

   from numpy import zeros,array, linspace
   icaPointLocator=vtk.vtkPointLocator()
   ecaPointLocator=vtk.vtkPointLocator()

   #averageCl=vtk.vtkPolyData()
   #averagePoints=vtk.vtkPoints()
   #averageSpheres=vtk.vtkDoubleArray()

   averageCell=vtk.vtkCellArray()
   averageCell.InsertNextCell(apexCl.GetNumberOfPoints()-1)
   averageCl=vtk.vtkPolyData()
   averagePoints=vtk.vtkPoints()

   cl=centerline
   apexCl=apexCl

   icaPointLocator.SetDataSet(icaCl)
   ecaPointLocator.SetDataSet(ecaCl)
   
   for i in range(0,apexCl.GetNumberOfPoints()):

      ecaPoint=[0.0,0.0,0.0]
      icaPoint=[0.0,0.0,0.0]
      averagePoint=[0.0,0.0,0.0]
      apexPoint=[0.0,0.0,0.0]
      
      apexCl.GetPoints().GetPoint(i,apexPoint)

      icaPointId=icaPointLocator.FindClosestPoint(apexPoint)
      ecaPointId=ecaPointLocator.FindClosestPoint(apexPoint)

      icaCl.GetPoints().GetPoint(icaPointId,icaPoint)
      ecaCl.GetPoints().GetPoint(ecaPointId,ecaPoint)

      averagePoint[0]=(icaPoint[0]+ecaPoint[0])/2
      averagePoint[1]=(icaPoint[1]+ecaPoint[1])/2
      averagePoint[2]=(icaPoint[2]+ecaPoint[2])/2
      averagePoints.InsertNextPoint(averagePoint[0],averagePoint[1],averagePoint[2])
      print i,averagePoint
      averageCell.InsertCellPoint(i)
   averageCl.SetPoints(averagePoints)
   averageCl.SetLines(averageCell)


   # replcing the average point 
   #apexCl.GetPoints().SetPoint(i,averagePoint)
   return averageCl


def OffProjection(line):
   from numpy import zeros, array, linspace

   points=(line.GetPoints()) 
   point =[0.0,0.0,0.0]
   pointArray=numpy.zeros((points.GetNumberOfPoints(),3))
   
   for i in range(0,points.GetNumberOfPoints()):
   
      points.GetPoint(i,point)
      pointArray[i]=point
   
   pointArray=pointArray.T
   covPoints=numpy.cov(pointArray)   
   [eigValue,eigVec]=numpy.linalg.eig(covPoints)
   
   # double check to see if the right vector is chosen
   # this is OFF normal vector
   # this is OFF NORMAL SO the Name is Not Correct
   normalVec=eigVec[0,:]
   
   meanPoint=numpy.zeros((3,1))
   meanPoint[0]=pointArray[0,:].mean()
   meanPoint[1]=pointArray[1,:].mean()
   meanPoint[2]=pointArray[2,:].mean()

   plane1 = createPlane(eigVec[2,:],meanPoint)
   WritePolyData(plane1.GetOutput(),'bfPlane.vtp')     

   plane2 = createPlane(eigVec[0,:],meanPoint)
   WritePolyData(plane2.GetOutput(),'bfPlaneOff.vtp')     

   plane3 = createPlane(eigVec[1,:],meanPoint)
   WritePolyData(plane3.GetOutput(),outdirectory+'_XXPlaneOff.vtp')     

   # projecting the points on the plane normal to the best fit plane   
   projPoints=vtk.vtkPoints()
   projPointIds=vtk.vtkIdList()
   line.GetCellPoints(0,projPointIds)
   
   plane=vtk.vtkPlane()
   for i in range (0,points.GetNumberOfPoints()):
      projPoint=[0.0,0.0,0.0]
      plane.ProjectPoint(points.GetPoint(i),meanPoint,normalVec,projPoint)
      projPoints.InsertNextPoint(projPoint[0],projPoint[1],projPoint[2])
   
   projCl=vtk.vtkCellArray()
   projCl.InsertNextCell(projPoints.GetNumberOfPoints())
   for i in range(0,projPoints.GetNumberOfPoints()):
      projCl.InsertCellPoint(i)
   
   projCenterline=vtk.vtkPolyData()
   projCenterline.SetLines(projCl)
   projCenterline.SetPoints(projPoints)
   
   # rotation of the projected line over the mean point to be pralllel to Y-Z
   # and then projecting on Y-Z --> x=0.
   rotatAngle=math.degrees(numpy.arccos(numpy.dot([1,0,0],normalVec)))
   rotation=vtk.vtkTransform()
   rotation.RotateWXYZ(rotatAngle,meanPoint)
   trnasformFilter=vtk.vtkTransformFilter()
   trnasformFilter.SetInput(projCenterline)
   trnasformFilter.SetTransform(rotation)
   trnasformFilter.Update()
   projClTrans=trnasformFilter.GetOutput()
   
   #projecting on Y-Z 
   projPoints=vtk.vtkPoints()
   
   projPointIds=vtk.vtkIdList()
   projClTrans.GetCellPoints(0,projPointIds)
   meanPointProj=meanPoint
   meanPointProj[0]=0
   plane=vtk.vtkPlane()
   for i in range (0,projClTrans.GetNumberOfPoints()):
      projPoint=[0.0,0.0,0.0]
      plane.ProjectPoint(projClTrans.GetPoints().GetPoint(i),meanPointProj,[1,0,0],projPoint)
      projPoints.InsertNextPoint(projPoint[0],projPoint[1],projPoint[2])
   
   projCl=vtk.vtkCellArray()
   projCl.InsertNextCell(projPoints.GetNumberOfPoints())
   for i in range(0,projPoints.GetNumberOfPoints()):
      projCl.InsertCellPoint(i)
    
   projCenterlineYZ=vtk.vtkPolyData()
   projCenterlineYZ.SetLines(projCl)
   projCenterlineYZ.SetPoints(projPoints)
   
   return projCenterline,projClTrans,projCenterlineYZ,plane1,plane2



def Projection(outdirectory,line,partOfCl):
   from numpy import zeros,array, linspace 
   points=(line.GetPoints()) 
   point =[0.0,0.0,0.0]
   pointArray=numpy.zeros((points.GetNumberOfPoints(),3))
   
   for i in range(0,points.GetNumberOfPoints()):
   
      points.GetPoint(i,point)
      pointArray[i]=point

   pointArray=pointArray.T
   covPoints=numpy.cov(pointArray)

   [eigValue,eigVec]=numpy.linalg.eig(covPoints)
   indexNormVec=0
   indexUpNormVec=0

   for i in range(0,3):
      if eigValue[i]==min(eigValue):
         indexNormVec=i
      if eigValue[i]==max(eigValue):
         indexUpNormVec=i

   xnormalVec=eigVec[:,indexUpNormVec]
   normalVec=eigVec[:,indexNormVec]
   print 'normal to the BFP of the '+partOfCl+'is: ',normalVec[0],normalVec[1],normalVec[2]
#   outputfile = open('./normalToBFP_'+partOfCl+'.dat','a')
#   outputfile.write('%s %s %f %f %f \n' % (outdirectory,partOfCl,normalVec[0],normalVec[1],normalVec[2]))
   #creating the normal vector and writing as .vtp
   icaInPlanePoints=vtk.vtkPoints()
   icaInPlanePoints.InsertPoint(0,[0.0,0.0,0.0])
   icaInPlanePoints.InsertPoint(1,normalVec)
   icaInPlaneLine=vtk.vtkCellArray()
   icaInPlaneLine.InsertNextCell(2)
   icaInPlaneLine.InsertCellPoint(0)
   icaInPlaneLine.InsertCellPoint(1)
   icaInPlaneV=vtk.vtkPolyData()
   icaInPlaneV.SetPoints(icaInPlanePoints)
   icaInPlaneV.SetLines(icaInPlaneLine)
   WritePolyData(icaInPlaneV,outdirectory+'_'+partOfCl+'Normal.vtp')

   xxnormalVec=numpy.cross(xnormalVec,normalVec)
   meanPoint=numpy.zeros((3,1))
# defining the center of the BFP (method 1)
   meanPoint[0]=pointArray[0,:].mean()
   meanPoint[1]=pointArray[1,:].mean()
   meanPoint[2]=pointArray[2,:].mean()
# defining the center of the BFP (method 2)
   midIndex=int(points.GetNumberOfPoints()/2)
   firstIndex=int(points.GetNumberOfPoints()/3)
#   meanPoint[0]=pointArray[0,midIndex]
#   meanPoint[1]=pointArray[1,midIndex]
#   meanPoint[2]=pointArray[2,midIndex]


   #BFP
   plane1 = createPlane(normalVec,meanPoint)
#   plane1 = createPlane(normalVec,[0.0,0.0,0.0])
   WritePolyData(plane1.GetOutput(),outdirectory+partOfCl+'_bfPlaneBigger1.vtp')     
   #sys.exit()
   print 'the normal to the '+partOfCl+'is: ',normalVec[0], normalVec[1],normalVec[2]

   #BFC test, using R=81.87, Center (x,53,-61)
   correctionV=[-81.8783578685*xxnormalVec[0]-0.5,81.8783578685*xxnormalVec[1]+5,-81.8783578685*xxnormalVec[2]]
   polygon=vtk.vtkRegularPolygonSource()
   polygon.SetNumberOfSides(100)
   polygon.SetRadius(81.8783578685);
   newCenter=[meanPoint[0]+correctionV[0],meanPoint[1]+correctionV[1],meanPoint[2]+correctionV[2]+0.05]
   polygon.SetCenter(newCenter);
#   polygon.SetCenter(meanPoint[0]+correctionV[0],meanPoint[1]+correctionV[1],meanPoint[2]+correctionV[2]);
#   polygon.SetCenter(meanPoint);
   polygon.SetNormal(normalVec);
   polygon.Update();
   circleTest=polygon.GetOutput()
   #WritePolyData(circleTest,'circle3D.vtp')

   plane2 = createPlane(xnormalVec,meanPoint)
   #WritePolyData(plane2.GetOutput(),outdirectory+partOfCl+'_bfPlaneOff1.vtp')     
   
#off BFP
   plane3 = createPlane(xxnormalVec,meanPoint)
   #WritePolyData(plane3.GetOutput(),outdirectory+partOfCl+'_bfPlaneOff2.vtp')     

   # bifurcation plane (the center and the normal)
   infofile = file(outdirectory+'_ref.dat',"r")
   infolines=infofile.readlines()
   infoLine=infolines[1]
   bifurcationCenter=numpy.zeros((3,1))
   bifurcationCenter[0]=float(infoLine.split()[0]);bifurcationCenter[1]=float(infoLine.split()[1]);bifurcationCenter[2]=float(infoLine.split()[2])
   bifurcationNormal=numpy.zeros((3,1))
   bifurcationNormal[0]=float(infoLine.split()[3]);bifurcationNormal[1]=float(infoLine.split()[4]);bifurcationNormal[2]=float(infoLine.split()[5])
   #writing the bifurcation plane
   plane4 = createPlane(bifurcationNormal,bifurcationCenter)
   #WritePolyData(plane4.GetOutput(),outdirectory+'_bifurcationPlane.vtp')     

   #we will project the centerlineCut on the following planes
   vectors=[normalVec,xnormalVec,xxnormalVec,bifurcationNormal]
   meanPointList=[meanPoint,meanPoint,meanPoint,bifurcationCenter]
   projCenterlineList=[0,0,0,0];projClTrans3List=[0,0,0,0];projCenterlineYZList=[0,0,0,0]

   for j in range (0,4):
      #projecting the points on the best fit plane
      projPoints=vtk.vtkPoints()
      projPointIds=vtk.vtkIdList()
      line.GetCellPoints(0,projPointIds)
      plane=vtk.vtkPlane()

      for i in range (0,points.GetNumberOfPoints()):
         projPoint=[0.0,0.0,0.0]
         plane.ProjectPoint(points.GetPoint(i),meanPointList[j],vectors[j],projPoint)
         projPoints.InsertNextPoint(projPoint[0],projPoint[1],projPoint[2])
      
      projCl=vtk.vtkCellArray()
      projCl.InsertNextCell(projPoints.GetNumberOfPoints())
      for i in range(0,projPoints.GetNumberOfPoints()):
         projCl.InsertCellPoint(i)
      
      projCenterline=vtk.vtkPolyData()
      projCenterline.SetLines(projCl)
      projCenterline.SetPoints(projPoints)

      nProj=[0.0,0.0,0.0]
#      nProj[0]=normalVec[0];nProj[1]=normalVec[1]
      nProj[0]=vectors[j][0];nProj[1]=vectors[j][1]
      nProjPre=[0.0,0.0,0.0]
      nProjPre[0]=nProj[1];nProjPre[1]=-nProj[0]
#      tiltAngle=math.degrees(numpy.arcsin(normalVec[2]))
      tiltAngle=math.degrees(numpy.arcsin(vectors[j][2]))
   
      rotation=vtk.vtkTransform()
      rotation.RotateWXYZ(-tiltAngle,nProjPre)
      rotation.Update()
   
      trans= vtk.vtkTransform()
      trans.Translate(-meanPointList[j])
   
      trnasformFilter=vtk.vtkTransformFilter()
      trnasformFilter.SetInput(projCenterline)
      trnasformFilter.SetTransform(trans)
      trnasformFilter.Update()
      projClTrans=trnasformFilter.GetOutput()

      trnasformFilter2=vtk.vtkTransformFilter()
      trnasformFilter2.SetInput(projClTrans)
      trnasformFilter2.SetTransform(rotation)
      trnasformFilter2.Update()
      projClTrans2=trnasformFilter2.GetOutput()
      
      rotation2=vtk.vtkTransform()
      rotatAngle=-(math.degrees(numpy.arctan(nProj[1]/nProj[0])))
      rotation2.RotateWXYZ(rotatAngle,[0,0,1])
      rotation.Update()
   
      trnasformFilter3=vtk.vtkTransformFilter()
      trnasformFilter3.SetInput(projClTrans2)
      trnasformFilter3.SetTransform(rotation2)
      trnasformFilter3.Update()
      projClTrans3=trnasformFilter3.GetOutput()
   
      #projecting on Y-Z 
      projPoints=vtk.vtkPoints()
      projPointIds=vtk.vtkIdList()
      projClTrans.GetCellPoints(0,projPointIds)
      meanPointProj=[0.0,0.0,0.0]
      meanPointProj[1]=meanPointList[j][1];meanPointProj[2]=meanPointList[j][2]
      plane=vtk.vtkPlane()
   
      for i in range (0,projClTrans.GetNumberOfPoints()):
         projPoint=[0.0,0.0,0.0]
         plane.ProjectPoint(projClTrans.GetPoints().GetPoint(i),meanPointProj,[1,0,0],projPoint)
         projPoints.InsertNextPoint(projPoint[0],projPoint[1],projPoint[2])
      
      projCl=vtk.vtkCellArray()
      projCl.InsertNextCell(projPoints.GetNumberOfPoints())
   
      for i in range(0,projPoints.GetNumberOfPoints()):
         projCl.InsertCellPoint(i)
       
      projCenterlineYZ=vtk.vtkPolyData()
      projCenterlineYZ.SetLines(projCl)
      projCenterlineYZ.SetPoints(projPoints)
     
      projCenterlineList[j]=projCenterline
      projClTrans3List[j]=projClTrans3
      projCenterlineYZList[j]=projCenterlineYZ
   return projCenterlineList,projClTrans3List,projCenterlineYZList


def BestFitCircle(lineProjectedOnY_Z):
    #this module finds the best arc fit for the 2D Cl which MUST BE projected on y-z plane  
    line =lineProjectedOnY_Z
    points=(line.GetPoints())
    point =[0.0,0.0,0.0]
    pointArray=numpy.zeros((points.GetNumberOfPoints(),3))
    
    for i in range(0,points.GetNumberOfPoints()):
    
       points.GetPoint(i,point)
       pointArray[i]=point
    
    xCl=pointArray[:,0]#which is zero since the line is projetced on y-z plane
    yCl=pointArray[:,1]
    zCl=pointArray[:,2]
    
    #defining the initial values for Radis and Center (the line MUST BE projected on y-z plane)
    [radius,center] = initialCircle2D(line)
    
    print 'initial radius = ',radius,'initial center =',center
    p=[center[0],center[1],radius]
    
    #points 
    x=yCl
    y=zCl
    X=zeros((2,len(xCl)))
    X[0,:]=x
    X[1,:]=y
    
    #the residual function
    def residuals(p,y,x):
    
       radius=p[2]
       err=sqrt((x-p[0])**2+(y-p[1])**2)-radius
       
       return err 

    p0=p
    #y_meas=X
    from scipy.optimize import leastsq
    plsq = leastsq(residuals, p0, args=(y,x))
    a=plsq[0]

    circle = Circle2D(a[2],[0,a[0],a[1]])
    #WritePolyData(circle,'circleFinal.vtp')
    #circle = Circle2D(p0[2],[0,p0[0],p0[1]])
    rotatedCircle = translate(circle,70,[10.0,10.0,0.0])
    WritePolyData(rotatedCircle,'rotatedCircl2.vtp')
    #WritePolyData(circle,'circle.vtp')
    return a[2],a[0],a[1]
    

# we want to the maximum curvature point closest to origin
def FindMaxClosestCurvePoint(modelName,centerline,ccaId,icaId,origin):
   ccaPoint=[0.0,0.0,0.0]
   icaPoint=[0.0,0.0,0.0]

   point=[0.0,0.0,0.0]
   points=numpy.zeros(3*(icaId-ccaId+1))
   points=points.reshape(icaId-ccaId+1,3)                                                                 

   tortuosityArray=numpy.zeros(icaId-ccaId+1)
   tortuosityArray.reshape(icaId-ccaId+1)
   maximum=[0,0]
   for i in range(ccaId,icaId+1):
      centerline.GetPoints().GetPoint(ccaId,ccaPoint)
      centerline.GetPoints().GetPoint(icaId,icaPoint)
      centerline.GetPoints().GetPoint(i,point)
      points[i-ccaId]=point
      tortuosityArray[i-ccaId]=abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i,1))
   tortuosityDiff=numpy.ediff1d(1000*tortuosityArray)
   x=tortuosityDiff
   peaks=zeros(len(x))
   localMin=zeros(len(x))
   localMinId=[]
   localMaxId=[]
   distance=[]
   distanceForPeaks=[]
   localMinPoint=[0.0,0.0,0.0]
   localMaxPoint=[0.0,0.0,0.0]

   #find the peaks which are maximum and local minimums
   for i in range (1,len(tortuosityDiff)):
      if x[i]*x[i-1]<0 and x[i]<0:
         peaks[i]=tortuosityArray[i]
         localMaxId.append(i) 
         centerline.GetPoints().GetPoint(ccaId+i,localMaxPoint)
         distanceForPeaks.append(math.sqrt(vtk.vtkMath.Distance2BetweenPoints(localMaxPoint,origin)))

   for i in range (1,len(tortuosityDiff)):
      if x[i]*x[i-1]<0 and x[i+1]>0:
         localMin[i]=tortuosityArray[i]
         localMinId.append(i) 
         centerline.GetPoints().GetPoint(ccaId+i,localMinPoint)
         distance.append(math.sqrt(vtk.vtkMath.Distance2BetweenPoints(localMinPoint,origin)))
   print 'local Maxs',localMaxId
   return localMaxId[distanceForPeaks.index(min(distanceForPeaks))]


# we want to find all the local maximum easy way 
def FindMinClosestCurvePoint(modelName, centerline, ccaId,icaId, origin):
   ccaPoint=[0.0,0.0,0.0]
   icaPoint=[0.0,0.0,0.0]

   point=[0.0,0.0,0.0]
   points=numpy.zeros(3*(icaId-ccaId+1))
   points=points.reshape(icaId-ccaId+1,3)                                                                 

   tortuosityArray=numpy.zeros(icaId-ccaId+1)
   tortuosityArray.reshape(icaId-ccaId+1)
   maximum=[0,0]
   for i in range(ccaId,icaId+1):
      centerline.GetPoints().GetPoint(ccaId,ccaPoint)
      centerline.GetPoints().GetPoint(icaId,icaPoint)
      centerline.GetPoints().GetPoint(i,point)
      points[i-ccaId]=point
      tortuosityArray[i-ccaId]=abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i,1))
#      if (abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i,1)))>abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i-1,1)) and (abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i,1)))<(abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i+1,1))):
#         print 'local max',point
#         raw_input("local max")
   tortuosityDiff=numpy.ediff1d(1000*tortuosityArray)
   x=tortuosityDiff
   peaks=zeros(len(x))
   localMin=zeros(len(x))
   localMinId=[]
   distance=[]
   localMinPoint=[0.0,0.0,0.0]
#find the peaks which are maximum 
#   for i in range (1,len(tortuosityDiff)):
#      if x[i]*x[i-1]<0 and x[i]<0:
#         peaks[i]=tortuosityArray[i]

#find the peaks which are local minimums
   for i in range (1,len(tortuosityDiff)):

      if x[i]*x[i-1]<0 and x[i+1]>0:
         localMin[i]=tortuosityArray[i]
         localMinId.append(i) 
         centerline.GetPoints().GetPoint(ccaId+i,localMinPoint)
         distance.append(math.sqrt(vtk.vtkMath.Distance2BetweenPoints(localMinPoint,origin)))
#   print 'local Mins',localMinId
   return localMinId[distance.index(min(distance))]


#this module finds the maximum torsion point and return its Idm nem
def FindCurvaturePoint(modelName,centerline,ccaId,icaId,origin):
   ccaPoint=[0.0,0.0,0.0]
   icaPoint=[0.0,0.0,0.0]

   point=[0.0,0.0,0.0]
   points=numpy.zeros(3*(icaId-ccaId+1))
   points=points.reshape(icaId-ccaId+1,3)                                                                 

   tortuosityArray=numpy.zeros(icaId-ccaId+1)
   tortuosityArray.reshape(icaId-ccaId+1)
   maximum=[0,0]
   for i in range(ccaId,icaId+1):
      centerline.GetPoints().GetPoint(ccaId,ccaPoint)
      centerline.GetPoints().GetPoint(icaId,icaPoint)
      centerline.GetPoints().GetPoint(i,point)
      points[i-ccaId]=point
      tortuosityArray[i-ccaId]=abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i,1))
      
   tortuosityDiff=numpy.ediff1d(1000*tortuosityArray)
   x=tortuosityDiff
   peaks=zeros(len(x))
   localMax=zeros(len(x))
   localMaxId=[]
   distance=[]
   localMaxPoint=[0.0,0.0,0.0]

   #find the peaks
   for i in range (1,len(tortuosityDiff)):
      if x[i]*x[i-1]<0 and x[i+1]>0:
         peaks[i]=tortuosityArray[i]
         localMax[i]=tortuosityArray[i]
         localMaxId.append(i) 
         centerline.GetPoints().GetPoint(ccaId+i,localMaxPoint)
         distance.append(math.sqrt(vtk.vtkMath.Distance2BetweenPoints(localMaxPoint,origin)))
   maxCurve1Id=localMaxId[distance.index(min(distance))]                                                                                                                     
   distance[distance.index(min(distance))]=100000
   maxCurve2Id=localMaxId[distance.index(min(distance))]                                                                                                                     
   return maxCurve1Id,maxCurve2Id

'''
#find the two biggest peaks and retuen the Ids
   maximum1=[0.0,0.0]
   for i in range (0,len(peaks)):
      if peaks[i]>maximum1[0]:
         maximum1[0]=peaks[i]
         maximum1[1]=i                             
   peaks[int(maximum1[1])]=0.0

   maximum2=[0.0,0.0]
   for i in range (0,len(peaks)):
      if peaks[i]>maximum2[0]:
         maximum2[0]=peaks[i]
         maximum2[1]=i                             
   # transfer these two Id to the main id based on the ica centelrine
   maximum1[1]=int(ccaId+maximum1[1])
   maximum2[1]=int(ccaId+maximum2[1])                                                  
   xlabel("pointId")
   ylabel("Curvature")
   title('Max:'+str(maximum1[1])+','+str(maximum2[1])+','+str(ccaId))
   plot(peaks)
   savefig('/das/pbijari/geometricFactors/figures/curvature/'+modelName+'.png')
   peaks.tofile('p.dat')                                                                                                                            
   return maximum1,maximum2
'''
#def FindMinimumCurvaturePoint:
def FindMinimumCurvaturePoint(modelName,centerline,ccaId,icaId):
   point=[0.0]
   maximum=[0,100]
   torsionArray=numpy.zeros(centerline.GetNumberOfPoints())

   for i in range(0,len(torsionArray)):
      point=abs(centerline.GetPointData().GetArray('Curvature').GetComponent(i,1))
      torsionArray[i]=point
      if (i in range (ccaId,icaId+1))and(point<maximum[1]):
         maximum[0]=i
         maximum[1]=point

   return maximum[0]    


def FindTorsionPoint(modelName, centerline, ccaId, icaId):   
   point = [0.0]
   maximum = [0,0]
   torsionArray = numpy.zeros(centerline.GetNumberOfPoints())
   for i in range(0, len(torsionArray)):
      point=abs(centerline.GetPointData().GetArray('Torsion').GetComponent(i, 1))
      torsionArray[i] = point
      if (i in range(ccaId, icaId + 1)) and (point > maximum[1]):
         maximum[0]=i
         maximum[1]=point

   xlabel("pointId")
   ylabel("Torsion")
   title('Maximum Torsion id is ' +str(maximum[0])+' maxCurves:'+str(ccaId)+','+str(icaId))
   plot (torsionArray)

   savefig('/das/pbijari/geometricFactors/figures/torsion/'+modelName+'.png')

   return maximum[0]    


#CutLine cuts a line at specific Id
def CutLine2(inputLine,id1,id2):
   
   points=inputLine.GetPoints()                         
   newPoints=vtk.vtkPoints()
   
   #cutting the centerline at the especific Id
   cutLine=vtk.vtkPolyData()
   line=vtk.vtkCellArray()
   line.InsertNextCell(int(id2)-int(id1)+1)
   print id1,id2 
   for i in range(int(id1),int(id2)+1):
      line.InsertCellPoint(i-int(id1))
      point=[0.0,0.0,0.0]
      points.GetPoint(i,point)
      newPoints.InsertNextPoint(point[0],point[1],point[2])
   cutLine.SetPoints(newPoints)
   cutLine.SetLines(line)   
   return cutLine    


def CutLineKeepArray(inputLine,id1,id2):
   points=inputLine.GetPoints()
   newPoints=vtk.vtkPoints()

   icaCcaRadiusArray=vtk.vtkDoubleArray()
   icaCcaRadiusArray.SetName('MaximumInscribedSphereRadius')
   icaCcaRadiusArray.SetNumberOfComponents(1)
   icaCcaRadiusArray.SetNumberOfTuples(int(id2)-int(id1)+1)
   icaCcaRadiusArray.FillComponent(0,1.0)
   sphereArray=inputLine.GetPointData().GetArray('MaximumInscribedSphereRadius')
   icaCcaRadiusArray.FillComponent(0,1.0)

   #cutting the centerline at the especific Id
   cutLine=vtk.vtkPolyData()
   line=vtk.vtkCellArray()
   line.InsertNextCell(int(id2)-int(id1)+1)
   print id1,id2 
   for i in range(int(id1),int(id2)+1):
      
      line.InsertCellPoint(i-int(id1))
      point=[0.0,0.0,0.0]
      points.GetPoint(i,point)
      newPoints.InsertNextPoint(point[0],point[1],point[2])

      sphere=sphereArray.GetTuple1(i)              
      icaCcaRadiusArray.SetTuple1(i-int(id1),sphere)
#   for i in range(0,int(id2)-int(id1)+1):
#       line.InsertCellPoint(i)

   cutLine.SetPoints(newPoints)
   cutLine.SetLines(line)
   cutLine.GetPointData().AddArray(icaCcaRadiusArray)
   return cutLine    


def CutLineAddArray(inputLine,id1,id2,inArray1,inArray2,inArray3):
   points=inputLine.GetPoints()
   newPoints=vtk.vtkPoints()

   icaCcaRadiusArray=vtk.vtkDoubleArray()
   icaCcaRadiusArray.SetName('cutSectionArea_cutSectionMinR_cutSectionMaxR')
   icaCcaRadiusArray.SetNumberOfComponents(3)
   icaCcaRadiusArray.SetNumberOfTuples(int(id2)-int(id1)+1)
#   icaCcaRadiusArray.FillComponent(0,1.0)
#   sphereArray=inputLine.GetPointData().GetArray('MaximumInscribedSphereRadius')
   #cutting the centerline at the especific Id
   cutLine=vtk.vtkPolyData()
   line=vtk.vtkCellArray()
   line.InsertNextCell(int(id2)-int(id1)+1)
   print id1,id2 
   for i in range(int(id1),int(id2)+1):
      
      line.InsertCellPoint(i-int(id1))
      point=[0.0,0.0,0.0]
      points.GetPoint(i,point)
      newPoints.InsertNextPoint(point[0],point[1],point[2])

      area=inArray1.GetTuple1(i)              
      minSize=inArray2.GetTuple1(i)              
      maxSize=inArray3.GetTuple1(i)              
      icaCcaRadiusArray.SetTuple3(i-int(id1),area,minSize,maxSize)
#   for i in range(0,int(id2)-int(id1)+1):
#       line.InsertCellPoint(i)

   cutLine.SetPoints(newPoints)
   cutLine.SetLines(line)
   cutLine.GetPointData().AddArray(icaCcaRadiusArray)
   return cutLine    


def CutLine(inputLine, Id):
    points = inputLine.GetPoints()                         
    newPoints = vtk.vtkPoints()
   
    # cutting the centerline at the especific Id
    cutLine = vtk.vtkPolyData()
    line = vtk.vtkCellArray()
    line.InsertNextCell(int(Id))
   
    for i in range(0,int(Id)):
        line.InsertCellPoint(i)
        point=[0.0,0.0,0.0]
        points.GetPoint(i,point)
        newPoints.InsertNextPoint(point[0],point[1],point[2])

    cutLine.SetPoints(newPoints)
    cutLine.SetLines(line)   
    return cutLine


#measuring the tortuosity between two points of a line
def CalculateTortuosity(line,pointId1,pointId2):
   id1=pointId1
   id2=pointId2
#   if id1>id2:
#      raw_input("id1 is bigger thaan id2 for calculating tortuosity")              
   point=[0.0,0.0,0.0]
   #the distance between the two end points of the intrested section for toruosity=D in L/D                      
   distance=math.sqrt(vtk.vtkMath.Distance2BetweenPoints(line.GetPoint(id1),line.GetPoint(id2)))
   delta=0
   length=0
   for i in range(id1,id2):

      delta=math.sqrt(vtk.vtkMath.Distance2BetweenPoints(line.GetPoint(i),line.GetPoint(i+1)))
      length=delta+length

#   print 'length two points =',length,'distance=',distance
   tortuosity=(length/distance)-1
   return tortuosity


def CalculateTortuosityIncludingLD(line,pointId1,pointId2):
   id1=pointId1
   id2=pointId2
#   if id1>id2:
#      raw_input("id1 is bigger thaan id2 for calculating tortuosity")              
   point=[0.0,0.0,0.0]
   #the distance between the two end points of the intrested section for toruosity=D in L/D                      
   distance=math.sqrt(vtk.vtkMath.Distance2BetweenPoints(line.GetPoint(id1),line.GetPoint(id2)))
   delta=0
   length=0
   for i in range(id1,id2):

      delta=math.sqrt(vtk.vtkMath.Distance2BetweenPoints(line.GetPoint(i),line.GetPoint(i+1)))
      length=delta+length

   #print 'length two points =',length,'distance=',distance
   tortuosity=(length/distance)-1
   return tortuosity,length,distance


def Length(line, pointId1, pointId2):
   id1=pointId1
   id2=pointId2
   point=[0.0,0.0,0.0]
   #the distance between the two end points of the intrested section for toruosity=D in L/D                      
   distance=math.sqrt(vtk.vtkMath.Distance2BetweenPoints(line.GetPoint(id1),line.GetPoint(id2)))
   delta=0
   length=0
   for i in range(id1,id2):

      delta=math.sqrt(vtk.vtkMath.Distance2BetweenPoints(line.GetPoint(i),line.GetPoint(i+1)))
      length=delta+length

   print 'length two points =',length,'distance=',distance
   tortuosity=(length/distance)-1
   return length



# Calculating the cut section area and the cut section point and the related
# .dat file, the return of this module is the section cut area and the max tube
# radius.  the inputs are the surface FILE NOT VTPXML and the cl FILE and the
# grup Id and the sphereradius.
def ReferenceSystem(surfaceFile,centerlineFile,outdirectory):
   command='vmtksurfacereader -ifile '+surfaceFile+' --pipe vmtkbranchextractor -ifile '+centerlineFile+' -radiusarray MaximumInscribedSphereRadius  --pipe vmtkbifurcationreferencesystems -radiusarray MaximumInscribedSphereRadius -blankingarray Blanking -groupidsarray GroupIds -normalarray Normal -ofile '+outdirectory+'_ref.dat'
   os.system(command)


def CutSection(surfaceFile, centerlineFile, groupId,sphereRadius,outdirectory):
   import numpy
   command='/das/pbijari/mainCodes/pointExtarctor.py'+ ' '+surfaceFile+' '+centerlineFile+' '+ groupId+' '+sphereRadius+' '+outdirectory  
   os.system(command)

   sectionFile = open(outdirectory+'_'+groupId+sphereRadius+'.dat')
   sectionLines = sectionFile.readlines()
   if groupId=='0':
      sectionLine = sectionLines[int(groupId)+1]#the first line is just titles
   if groupId=='3':
      sectionLine = sectionLines[int(groupId)]#the first line is just titles

#we enter this 'if' for the cases with short cca or eca
   if groupId=='2': 
      if len(sectionLines)==2:
   	   sectionLine=sectionLines[1]#for ica5
   if groupId=='2':
      if len(sectionLines)>=3:
   	   sectionLine=sectionLines[2]
   
   sectionArea = float(sectionLine.split()[11])
   sectionRadius = float(sectionLine.split()[12])
   #getting the coordinate of the section points of the cut sections 
   sectionPoint=numpy.array((0.0,0.0,0.0))
   sectionPoint[0]=float(sectionLine.split()[5])
   sectionPoint[1]=float(sectionLine.split()[6])
   sectionPoint[2]=float(sectionLine.split()[7])
   #finding the id of the section point
   cl = ReadPolyData(centerlineFile)
   sectionPointId=FindPoint(cl,sectionPoint)
   return sectionArea, sectionRadius,sectionPointId,sectionPoint

def ExtractCell(inputPoly,cellNumber):
   
   import vtk
   points=inputPoly.GetPoints()
   cell=inputPoly.GetCell(cellNumber)
   cellPointsNo=cell.GetNumberOfPoints()
   cellIds=cell.GetPointIds()
   cellArray=vtk.vtkCellArray()
   
   outputPoly=vtk.vtkPolyData()
   outputPoly.SetPoints(points)

   cellArray.InsertNextCell(cellPointsNo)
   for i in range(0,cellPointsNo):
      cellArray.InsertCellPoint(cellIds.GetId(i))

   outputPoly.SetPolys(cellArray)
   return outputPoly 


