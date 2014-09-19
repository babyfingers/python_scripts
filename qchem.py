#!/usr/local/bin/python

#Q-Chem output file anaylsis module

import string
import os.path
import math
import time
import sys
import re
import numpy as np
import numpy.linalg as la

elementMasses = {"H":1, "C":12, "N":14, "O":16, "F":19, "Li":6, "S":32 }
geometryStart_re  = re.compile("Standard Nuclear Orientation")
geometryEnd_re    = re.compile("----------------------------------------------------------------")
atomCoord_re      = re.compile("\s+[1-9][0-9]*\s+[A-Z][a-z]*\s+-?[0-9]+\.[0-9]+\s+-?[0-9]+\.[0-9]+\s+-?[0-9]+\.[0-9]+")
modeRow_re        = re.compile("\s+[A-Z][a-z]?\s+-?[0-9]+\.[0-9]+.*")
num_re            = re.compile("-?[0-9]+\.[0-9]+")
vibBlockStart_re  = re.compile(" Mode:\s+[1-9][0-9]*")
vibBlockEnd_re    = re.compile("\s*TransDip\s+.*")
aimdBlockStart_re = re.compile("TIME STEP #[1-9][0-9]*")
aimdBlockEnd_re   = re.compile("Nuclear Repulsion Energy")
frequency_re      = re.compile("\s*Frequency:\s+.*")
forceConstant_re  = re.compile("\s*Force Cnst:\s+.*")
reducedMass_re    = re.compile("\s*Red\. Mass:\s+.*")

freqString = "frequency"
FCString   = "force_constant"
RMString   = "reduced_mass"
modeString = "vibrational_mode"

A_to_Bohr = 1.889725989 #Number of Bohr radii in an Angstrom
Bohr_to_m = 5.2917721092*math.pow(10,-11)
Har_to_J  = 4.35974434*math.pow(10,-18)
SOL       = 299792458 #Speed of light, in m/s
amu_to_kg = 1.666053892*math.pow(10,-27)
kg_to_amu = 1.0/amu_to_kg
m_to_cm   = 100
cm_to_m   = 0.01
oneNano   = math.pow(10.0,-9)
onePico   = math.pow(10.0,-12)
muX = 0
muY = 1
muZ = 2


quantityList = (freqString,FCString,RMString,modeString)
vibDict = {freqString:frequency_re,FCString:forceConstant_re,RMString:reducedMass_re}
NQuantity = len(quantityList)

def getFilePos(myFile, startPos, myRE, mustFind=False):
    '''Starting from startPos, find the position in a
    file that marks the beginning of a line containing
    a match for myRE. Return the file positions at the
    beginning of this line and at the start of the next
    line in the form of a 2-ple.'''
    doneSearching = False
    myFile.seek(0,2)
    endPos = myFile.tell()
    myFile.seek(startPos)
    currPos = startPos
    while (currPos<endPos and not doneSearching):
        if len(myRE.findall(myFile.readline()))>0:
            doneSearching = True
            lineBeginPos = currPos
        currPos = myFile.tell()
    #If we can't find such a position, return -1
    if (mustFind): assert(doneSearching)
    if doneSearching:
        return (lineBeginPos, currPos)
    else:
        return (-1,-1)


def getGeoFromString(geoString):
    """For a string that's primarily just geometry
    data from a Q-Chem output file (and nothing
    else that might be confused for it), return some
    geometry data. Note that in Q-Chem output files,
    the atom number comes before the element symbol
    and Cartesian coordinates."""
    atomList = atomCoord_re.findall(geoString)
    NAtom = len(atomList)
    elementList = [""]*NAtom
    geometry = np.zeros(3*NAtom)
    for iii in range(NAtom):
        wordList = atomList[iii].split()
        elementList[iii] = wordList[1]
        for mu in range(3):
            geometry[iii*3 +mu] = float(wordList[mu+2])
    return (len(atomList),tuple(elementList),geometry)

def getGeoData(myFile):
    '''Look at an output file and return data about
    the system geometry, including number of atoms
    and geometry vectors.'''
    _, geoStart = getFilePos(myFile, 0, geometryStart_re, True)
    #Find at least one atom coordinate before searching
    # for the end of the geometry section.
    _, geoStarted = getFilePos(myFile, geoStart, atomCoord_re, True)
    endPos, _ = getFilePos(myFile, geoStarted, geometryEnd_re, True)
    geoLength = geoStart -endPos
    myFile.seek(geoStart)
    geoString = myFile.read(geoLength)
    #atomList = atomCoord_re.findall(geoString)
    #NAtom = len(atomList)
    #elementList = [""]*NAtom
    #geometry = np.zeros(NAtom*3)
    #for iii in range(NAtom):
    #    wordList = atomList[iii].split()
    #    elementList[iii] = wordList[1]
    #    for mu in range(3):
    #        geometry[iii*3 +mu] = float(wordList[mu+2])
    #return (len(atomList),tuple(elementList),geometry)
    return getGeoFromString(geoString)
    
def getGeoDataX(myFile):
    '''Look at a file containing only the Cartesian-
    formatted geometry information, then return the
    number of atoms, element list, and geometry
    vector.'''
    geoString = myFile.read()
    geoList = geoString.split()
    assert(len(geoList)%4==0)
    NAtom = len(geoList)/4
    elementList = [""]*NAtom
    geometry = np.zeros(3*NAtom)
    for iii in range(NAtom):
        elementList[iii] = geoList[iii*4]
        for mu in range(3):
            geometry[iii*3 +mu] = float(geoList[iii*4 +mu+1])
    return (NAtom, tuple(elementList), geometry)

def getVibBlockLocations(myFile, NMode):
    '''Return a tuple containing a set of 2-ples
    corresponding to the beginning position and length
    of each vibrational analysis block.'''
    NBlock = (NMode-1)/3 +1
    #print "NBlock = "+str(NBlock)
    blockIndices = [0]*NBlock
    startPos = 0
    for iii in range(NBlock):
        blockPosList = [0]*2
        blockPosList[0], startPos = getFilePos(myFile, startPos, vibBlockStart_re, True)
        startPos, blockPosList[1] = getFilePos(myFile, startPos, vibBlockEnd_re, True)
        blockPosList[1] = blockPosList[1] -blockPosList[0]
        blockIndices[iii] = tuple(blockPosList)
        #print blockIndices[iii]
    return tuple(blockIndices)

def getAIMDblockLocations(myFile,NStep,includeFileInTuple=False):
    '''Return a tuple containing a set of 2-ples
    corresponding to the beginning position and length
    of each AIMD-generated geometry.'''
    blockIndices = [0]*NStep
    startPos = 0
    for IStep in range(NStep):
        if (includeFileInTuple):
            blockPosList = [0]*3
            blockPosList[2] = myFile
        else:
            blockPosList = [0]*2
        blockPosList[0], startPos = getFilePos(myFile, startPos, aimdBlockStart_re, True)
        startPos, blockPosList[1] = getFilePos(myFile, startPos, aimdBlockEnd_re, True)
        blockPosList[1] = blockPosList[1] -blockPosList[0]
        blockIndices[IStep] = tuple(blockPosList)
    return tuple(blockIndices)

def getVibMode(myBlockString,myColumn,NAtom):
    #print myBlockString
    stringList = modeRow_re.findall(myBlockString)
    #print "len(stringList) = "+str(len(stringList))
    #for s in stringList:
    #    print s
    assert(len(stringList)==NAtom)
    myModeVec = np.zeros(NAtom*3)
    for iii in range(NAtom):
        rowWordList = stringList[iii].split()
        for mu in range(3):
            myModeVec[iii*3 +mu] = float(rowWordList[1 +myColumn*3 +mu])
    return myModeVec

def getVibQuantity(myFile,myBlockIndices,myMode,myQuantity,NAtom,NMode):
    quantityIsMode = (myQuantity==modeString)
    assert(myMode>0 and myMode<=NMode)
    modeIndex = myMode-1
    if(not quantityIsMode):
        myRE = vibDict[myQuantity]
    myBlockIndex = modeIndex/3
    myColumn = modeIndex%3
    myFile.seek(myBlockIndices[myBlockIndex][0])
    myBlockString = myFile.read(myBlockIndices[myBlockIndex][1])
    if (not quantityIsMode):
        stringList = myRE.findall(myBlockString)
        assert(len(stringList)==1)
        stringList = num_re.findall(stringList[0])
        assert(len(stringList)>myColumn)
        return float(stringList[myColumn])
    else:
        return getVibMode(myBlockString,myColumn,NAtom)

def getAllScalarVibQuantities(myFile,myBlockIndices,NAtom,NMode):
    myFreq    = []
    myRedMass = []
    myFC      = []
    for iii in range(1,NMode+1):
        myFreq.append(getVibQuantity(myFile,myBlockIndices,iii,freqString,NAtom,NMode))
        myRedMass.append(getVibQuantity(myFile,myBlockIndices,iii,RMString,NAtom,NMode))
        myFC.append(getVibQuantity(myFile,myBlockIndices,iii,FCString,NAtom,NMode))
    return (np.array(myFreq),np.array(myRedMass),np.array(myFC))

def getModeOverlap(mode1, mode2, elements, NAtom):
    NCoord = 3*NAtom
    massWeight = np.zeros((NCoord,NCoord))
    for iii in range(NAtom):
        mass = elementMasses[elements[iii]]
        for mu in range(3):
            massWeight[iii*3 +mu, iii*3 +mu] = mass
    overlap = np.dot(mode1.transpose(),massWeight)
    overlap = np.dot(overlap,mode2)
    assert(overlap.shape==(1,1))
    return overlap[0][0]

def geometryToString(elements, geometry):
    '''Return a string that works as a properly
    formatted Q-Chem geometry section in Cartesian
    coordinates.'''
    NAtom = len(elements)
    assert(3*NAtom==len(geometry))
    geoString = ""
    for iii in range(NAtom):
        geoString += "{:<2} ".format(elements[iii])
        for mu in range(3):
            geoString += "    {:13.10f}".format(geometry[iii*3 +mu])
        geoString += "\n"
    return geoString

def createFile(fileName, fileContentString):
    f = open(fileName,"w")
    f.write(fileContentString)
    f.close()
    return 0

def readFile(fileName):
    fileName = open(fileName,"r")
    fileString = fileName.read()
    fileName.close()
    return fileString

def getCenterOfMass(elements,geometry):
    totalMass = 0.0
    COM = np.zeros(3)
    for iii,element in enumerate(elements):
        mass = elementMasses[element]
        for mu in range(3):
            COM[mu] = COM[mu] +mass*geometry[iii*3 +mu]
        totalMass = totalMass +mass
    return COM/totalMass
    
def geoTranslate(geometry, transVector):
    NAtom = geometry.size/3
    assert(len(transVector)==3)
    translatedGeometry = np.zeros(3*NAtom)
    for iii in range(NAtom):
        for mu in range(3):
            translatedGeometry[iii*3 +mu] = geometry[iii*3 +mu] +transVector[mu]
    return translatedGeometry

def normalize(myBasis):
    '''Given a numpy array of size NVec*Vlen, for NVec
    vectors of length Vlen, return an normalized
    transformation of same. Given a vector of length
    NVec, does the same thing, except for the single
    vector.'''
    NDim = len(myBasis.shape)
    if (NDim==1):
        myNorm = la.norm(myBasis)
        newBasis = np.zeros(myBasis.size)
        if (myNorm>onePico):
            newBasis = myBasis/myNorm
    elif (NDim==2):
        Vlen, NVec = myBasis.shape
        newBasis = np.zeros((Vlen,NVec))
        for iii in range(NVec):
            myNorm = la.norm(myBasis[:,iii])
            if (myNorm>onePico):
                newBasis[:,iii] = myBasis[:,iii]/myNorm
    return newBasis

def projectOut(myBasis,myVector):
    '''Return myBasis after projecting out myVector.'''
    normVector = normalize(myVector)
    proj = np.array([np.dot(myBasis.transpose(),normVector)])
    #print "normVector: ",normVector
    #print "myBasis: ",myBasis
    #print "proj: ",proj
    #print "proj.shape: ",proj.shape
    #print "myBasis.shape: ",myBasis.shape
    subMat = np.dot(np.array([normVector]).transpose(),proj)
    #print "subMat.shape: ",subMat.shape
    return myBasis -subMat

def orthonormalize(myBasis):
    '''Given a numpy array of size NVec*Vlen, for NVec
    vectors of length Vlen, return an orthonormalized
    transformation of same.'''
    assert(len(myBasis.shape)==2)
    Vlen, NVec = myBasis.shape
    newBasis = np.array(myBasis,copy=True)
    vecsToDelete = []
    for iii in range(NVec):
        myNorm = la.norm(newBasis[:,iii])
        #print "myBasis["+str(iii)+"], myNorm = "+str(myNorm)
        myVec = normalize(np.array(newBasis[:,iii],copy=True))
        #print "myVec.shape: ",myVec.shape
        #print "newBasis[:,iii+1:].shape: ",newBasis[:,iii+1:].shape
        newBasis[:,iii] = np.array(myVec,copy=True)
        if (myNorm>oneNano and iii<NVec-1):
            newBasis[:,iii+1:] = projectOut(newBasis[:,iii+1:],myVec)
            #proj = np.array([np.inner(myVec,newBasis[iii+1:,:])])
            #print "proj: ",proj
            #subMat = np.dot(proj.transpose(),np.array([myVec]))
            #print "subMat: ",subMat
            #newBasis[iii+1:,:] = newBasis[iii+1:,:] -subMat
            #print "newBasis: ",newBasis
        elif (myNorm<oneNano):
            vecsToDelete.append(iii)
    #Get rid of any zero vectors
    #print "newBasis: ",newBasis
    if (len(vecsToDelete)>0):
        vecsToDelete.reverse()
        for iii in vecsToDelete:
            newBasis = np.delete(newBasis,iii,0)
    #print "newBasis: ",newBasis
    return normalize(newBasis)

def symmetryTest(M,verbose):
    '''Given square array M, find the difference between Mij
    and Mji and return the largest such value.'''
    nn,mm = M.shape
    assert(nn==mm)
    maxDiff = 0.0
    maxIII = -1
    maxJJJ = -1
    for iii in range(1,nn):
        for jjj in range(iii):
            myDiff = math.fabs(M[iii,jjj]-M[jjj,iii])
            if (myDiff>maxDiff):
                maxDiff = myDiff
                maxIII = iii
                maxJJJ = jjj
    if (verbose):
        if (maxDiff==0.0):
            print "This matrix is perfectly symmetric."
        else:
            print "Maximum deviation from symmetry: "+str(maxDiff)
            print "Occurred at M["+str(maxIII)+","+str(maxJJJ)+"] = "+str(M[maxIII,maxJJJ])
            print "compare to  M["+str(maxJJJ)+","+str(maxIII)+"] = "+str(M[maxJJJ,maxIII])
    return maxDiff
            
def massWeightHessian(Hess,elements,cartToMW):
    '''Convert a Hessian from Cartesian to mass-weighted coordinates,
    or do the reverse, and output the result.'''
    nn,mm = Hess.shape
    assert(nn==mm)
    NAtom = len(elements)    
    assert(NAtom==nn/3)
    outHess = np.array(Hess,copy=True)
    for iii,elementI in enumerate(elements):
        massI = elementMasses[elementI]
        for mui in range(3):
            for jjj, elementJ in enumerate(elements):
                massJ = elementMasses[elementJ]
                massPow = -0.5 if cartToMW else 0.5
                massFac = math.pow(massI*massJ,massPow)
                for muj in range(3):
                    outHess[iii*3+mui,jjj*3+muj] = Hess[iii*3+mui,jjj*3+muj]*massFac
    return outHess

def deleteSmallVectors(myBasis,thresh):
    '''Eliminated vectors from a basis if their
    2-norm is smaller than a threshold value.'''
    vecsToDelete = []
    newBasis = np.array(myBasis,copy=True)
    Vlen, NVec = myBasis.shape
    for iii in range(NVec):
        if la.norm(myBasis[:,iii])<thresh:
            vecsToDelete.append(iii)
    vecsToDelete.reverse()
    for iii in vecsToDelete:
        newBasis = np.delete(newBasis,iii,1)
    return newBasis

def getPrunedBasis(myBasis,NVecToPrune):
    '''Go through basis comprised of array of size NVec*Vlen
    and remove the smallest-norm basis vector, normalizing
    and projecting out the non-small basis vectors as it goes.
    Return the pruned basis set.'''
    Vlen, NVec0 = myBasis.shape
    NDesired = NVec0-NVecToPrune
    #List the vectors from highest norm to lowest-norm
    lB = list(myBasis.transpose())
    lB = sorted(lB,key=la.norm)
    lB.reverse()
    newBasis = np.array(lB).transpose()
    #print "newBasis: ",newBasis
    newBasis = deleteSmallVectors(newBasis,oneNano)
    _, NVec1 = newBasis.shape
    if (NVec1<NVec0):
        print "Pruned "+str(NVec0-NVec1)+" vectors in the initial sweep."
    if (NVec1<=NDesired):
        print "There were enough small vectors as given. Success!"
        return newBasis
    #Project out each vector one at a time from remaining vectors
    NVec2 = NVec1
    for iii in range(NVec1-1):
        if (iii>=NVec2):
            break
        #print "Round "+str(iii)+":"
        newBasis[:,iii+1:] = projectOut(newBasis[:,iii+1:],newBasis[:,iii])
        #print newBasis
        newBasis = deleteSmallVectors(newBasis,oneNano)
        _, NVec2 = newBasis.shape
    #If this has failed, just delete vectors with smallest
    #residual norm.
    _, NVec2 = newBasis.shape
    if (NVec2<NVec1):
        print "Pruned "+str(NVec1-NVec2)+" vectors after projecting out other basis vectors."
    if (NVec2<=NDesired):
        print "Pruned enough vectors out by projecting out other basis vectors."
        return newBasis
    NVecsToDelete = NVec2-NDesired
    lB = list(newBasis)
    lB = sorted(lB,key=la.norm)
    lB.reverse()
    newBasis = np.array(lB)
    #print "newBasis: ",newBasis
    for iii in range(NVecsToDelete):
        print "About to delete column "+str(iii)
        newBasis = np.delete(newBasis,NVec2-1-iii,1)
        #print "newBasis: ",newBasis
    print "Just got rid of the remaining "+str(NVecToPrune)+" smallest vector(s)."
    print "newBasis.shape: ",newBasis.shape
    return newBasis
    
def getMomentOfInertiaTensor(geometry, elements):
    '''Given the center-of-mass geometry vector and
    a list of elements, produce the 3x3 moment of
    inertia tensor, with units of mass in a.m.u.
    and units of length in whatever the geometry
    vector is in.'''
    MOIT = np.zeros((3,3))
    myMasses = [elementMasses[myElement] for myElement in elements]
    for atomI, massI in enumerate(myMasses):
        xcoordI = geometry[atomI*3 +muX]
        ycoordI = geometry[atomI*3 +muY]
        zcoordI = geometry[atomI*3 +muZ]
        MOIT[muX,muX] += massI*(ycoordI*ycoordI +zcoordI*zcoordI)
        MOIT[muY,muY] += massI*(xcoordI*xcoordI +zcoordI*zcoordI)
        MOIT[muZ,muZ] += massI*(xcoordI*xcoordI +ycoordI*ycoordI)
        MOIT[muX,muY] -= massI*xcoordI*ycoordI
        MOIT[muX,muZ] -= massI*xcoordI*zcoordI
        MOIT[muY,muZ] -= massI*ycoordI*zcoordI
        MOIT[muY,muX] += MOIT[muX,muY]
        MOIT[muZ,muX] += MOIT[muX,muZ]
        MOIT[muZ,muY] += MOIT[muY,muZ]        
    return MOIT

class coarseGrainBin:
    """A bin that takes a bunch of quantities for storage
    and averaging, and keeps track of how many it contains."""
    def __init__(self, myCGmax,myTier=0):
        self.coarseGrainMax = myCGmax
        self.NQuant = 0
        self.tier = 0
        self.quant = 0
        self.isFull = False
    def addQuantity(self, myQuant, myTier):
        assert(not self.isFull)
        if (self.NQuant==0):
            self.tier = myTier
            self.quant = myQuant
        else:
            assert (self.tier==myTier)
            self.quant += myQuant
        self.NQuant += 1
        if (self.NQuant==self.coarseGrainMax):
            #If this bin is full, divide contents by
            #number of quantities added and set status
            #to full.
            self.isFull = True

class coarseGrainAverage:
    """Keep track of data while we look through a file
    to build an average quantity by continuousely
    coarse-graining."""
    coarseGrainMax = 1000
    def __init__(self):
        self.avg  = 0
        self.bins = []
    def countOpenBins(self,myTier):
        count = 0
        if (len(self.bins)>myTier):
            for myBin in self.bins[myTier]:
                if (myBin.tier==myTier and not myBin.isFull):
                    count = count +1
        return count
    def getOpenBin(self,myTier):
        if (len(self.bins)>myTier):
            for myBin in self.bins[myTier]:
                if (myBin.tier==myTier and not myBin.isFull):
                    return myBin
        else:
            while(len(self.bins)<=myTier):
                self.bins.append([])
        self.bins[myTier].insert(0,coarseGrainBin(self.coarseGrainMax,myTier))
        return self.bins[myTier][0]
    def getClosedBins(self,myTier):
        if (len(self.bins)>myTier):
            return [myBin for myBin in self.bins[myTier] if myBin.isFull]
        else:
            return []
    def getAllBins(self,myTier):
        return [myBin for myBin in self.bins[myTier]]
    def getMaxTier(self):
        for iTier in range(len(self.bins)-1,-1,-1):
            if(len(self.bins[iTier])>0):
                return iTier
        return -1
    def consolidateBins(self):
        for iTier in range(self.getMaxTier()+1):
            closedBins = self.getClosedBins(iTier)
            for myClosedBin in closedBins:
                self.addQuantity(myClosedBin.quant/self.coarseGrainMax,iTier+1)
                self.bins[iTier].remove(myClosedBin)
    def addQuantity(self,myQuant,myTier=0):
        myBin = self.getOpenBin(myTier)
        myBin.addQuantity(myQuant,myTier)
    def initializeAverage(self):
        NBin = 0
        for iTier in range(self.getMaxTier()+1):
            binList = self.getAllBins(iTier)
            NBin += len(binList)
            if (NBin>0):
                break
        if (NBin==0):
            self.avg = None
            return
        self.consolidateBins()
        myQuant = binList[0].quant
        if (isinstance(myQuant,float)):
            self.avg = 0.0
        elif (isinstance(myQuant,int)):
            self.avg = 0
        elif (isinstance(myQuant,np.ndarray)):
            self.avg = np.zeros(myQuant.shape)
    def getAverage(self):
        self.initializeAverage()
        totalTierWeight = 0.0
        CGMin = 1.0/self.coarseGrainMax
        quantCount = 0
        for iTier in range(self.getMaxTier()+1):
            binsAtThisTier = self.getAllBins(iTier)
            totalTierWeight = totalTierWeight*CGMin
            sumForThisTier = self.avg*CGMin
            for myBin in binsAtThisTier:
                totalTierWeight = totalTierWeight +myBin.NQuant
                sumForThisTier = sumForThisTier +myBin.quant
                quantCount = quantCount +myBin.NQuant*int(math.pow(self.coarseGrainMax,iTier))
            assert(totalTierWeight<self.coarseGrainMax)
            if (totalTierWeight>0.0):
                self.avg = sumForThisTier/totalTierWeight
        #print "Using "+str(quantCount)+" values, the average value is:"
        #print self.avg
        return self.avg
    def report(self):
        for iTier in range(self.getMaxTier()+1):
            print "Looking at tier {}".format(iTier)
            print "\tThis tier has {} bins, containing the following number of quantities:".format(len(self.getAllBins(iTier)))+str([myBin.NQuant for myBin in self.bins[iTier]])
            print "\tThe value of these quantities is: "+str([myBin.quant for myBin in self.bins])
            print "\tThe fill status of these bins is: "+str([myBin.isFull for myBin in self.bins])


