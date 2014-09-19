#!/usr/local/bin/python

import string
import os.path
import math
import time
import sys
sys.path.append('/data/home/alguire/inout_files/python')
import re
import numpy as np
import numpy.linalg as la
import qchem as qc
import random
from datetime import datetime

#fileName = "benzene-freq-631Gs-parallel-cq.out02"
#fileName = "B1A1c2_Toluene-omegaB97XD.out01"
#fileName = "BT1_trans_C2_Toluene.x"
fileName = "emimTf2N.x"
versionString = "0.3"

# Version notes: 0.3
# This code now accomodates for the selected atoms to be
# put into different planes (rather than just the yz
# plane as in previous versions).

# Rotate.py: put this stupid molecule
# in the yz plane

logFile = open("rotate.log","w")
logFile.write("Beginning run of Rotate.py, v"+versionString+". It is ")
timeTuple = datetime.now()
logFile.write(str(timeTuple.year)+"-"+str(timeTuple.month)+"-"+str(timeTuple.day)+" "+str(timeTuple.hour)+":"+str(timeTuple.minute)+":"+str(timeTuple.second)+"\n")
logFile.write("fileName1 = "+fileName+"\n\n")

myFile = open(fileName,"r")
if (fileName[-1]!="x"):
    #Assume a Q-Chem output file
    NAtom, elements, geometry = qc.getGeoData(myFile)
else:
    #Just the Cartesian-formatted geometry information,
    #nothing else in file.
    NAtom, elements, geometry = qc.getGeoDataX(myFile)

def translate(myGeometry,myOffset):
    assert(len(myGeometry)%3==0)
    NAtom = len(myGeometry)/3
    assert(len(myOffset)==3)
    for iii in range(NAtom):
        for mu in range(3):
            myGeometry[iii*3 +mu] += myOffset[mu]
    return myGeometry

def getRotationMatrix(theta,cartToZero,otherCart):
    rotMatrix = np.zeros((3,3))
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    rotMatrix[cartToZero][cartToZero] = cosTheta
    rotMatrix[cartToZero][otherCart] = sinTheta
    rotMatrix[otherCart][cartToZero] = -sinTheta
    rotMatrix[otherCart][otherCart] = cosTheta
    notMyCart = (2*(cartToZero+otherCart))%3
    assert(notMyCart!=cartToZero and notMyCart!=otherCart)
    rotMatrix[notMyCart][notMyCart] = 1.0
    #print rotMatrix
    return rotMatrix

def applyRotationMatrix(myGeometry,rotMatrix):
    NAtom = len(myGeometry)/3
    #print "Before rotation."
    #print qc.geometryToString(elements,geometry)
    for iii in range(NAtom):
        u = myGeometry[iii*3:iii*3+3]
        u = np.dot(u,rotMatrix)
        myGeometry[iii*3:iii*3+3] = u
        #for mu in range(3):
        #    myGeometry[iii*3+mu] = u[mu]
    #print "After rotation."
    #print qc.geometryToString(elements,geometry)
    return myGeometry


def rotate(myGeometry,myAtom,cartToZero,otherCart):
    NAtom = len(myGeometry)/3
    assert(myAtom<NAtom)
    assert(cartToZero!=otherCart)
    if abs(geometry[myAtom*3 +otherCart])>0.0000001:
        theta = np.arctan(geometry[myAtom*3+cartToZero]/geometry[myAtom*3 +otherCart])
    elif geometry[myAtom*3]==0.0:
        theta = 0.0
    else:
        theta = math.pi/2
    rotMatrix = getRotationMatrix(theta,cartToZero,otherCart)
    return applyRotationMatrix(myGeometry,rotMatrix)
    
def getXscore(myGeometry,atomList):
    Xarray = []
    for atomIndex in atomList:
        Xarray.append(myGeometry[atomIndex*3])
    return la.norm(np.array(Xarray))

#If we can get this list of atoms in the yz plane, we'll be all set.
atomList = np.array([0,1,2,3,4,5,6])
#atomList = np.array([35,36,37,38,39])
#atomList = np.array([10,11,12,13,14,15])
#atomList = np.array([8,9,10,11,12,13])
# When we do rotations, which cartesian coordinate
# should be sent to zero?
cartesianCoordToZero = 2 # 0 for x
                         # 1 for y
                         # 2 for z
otherCoords = [0,1,2]
otherCoords.remove(cartesianCoordToZero)
stillRotating = True
NIter = 0
xTotal = getXscore(geometry,atomList)
NConsecutiveConvergences = 0
while(stillRotating):
    #Translate so that the first atom is at the origin
    #Rotate the atoms in the list so that they have x (or y, or z) = 0
    geometry = translate(geometry,(-geometry[atomList[0]*3],-geometry[atomList[0]*3+1],-geometry[atomList[0]*3+2]))
    #for mu in range(1,3):
    for iii in range(len(atomList)):
        atom_i = atomList[iii]
        if iii==0:
            geometry = rotate(geometry,atom_i,0,1)
            geometry = rotate(geometry,atom_i,1,2)
        else:
            mu = otherCoords[random.randint(0,1)]
            geometry = rotate(geometry,atom_i,cartesianCoordToZero,mu)
    xChange = xTotal
    xTotal = getXscore(geometry,atomList)
    xChange = xTotal -xChange
    NIter += 1
    logFile.write("For cycle "+str(NIter)+", xTotal = "+str(xTotal)+", xChange = "+str(xChange)+"\n")
    if (xTotal <0.00001 or abs(xChange)<0.00001 or NIter>1000):
        if(NConsecutiveConvergences>2):
            stillRotating = False
        else:
            NConsecutiveConvergences += 1
    else:
        NConsecutiveConvergences = 0

print qc.geometryToString(elements,geometry)

logFile.write("Script completed.\n")
logFile.close()