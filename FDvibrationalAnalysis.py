#!/usr/local/bin/python

#*********************************************************
# FDvibrationalAnalysis.py
# From a set of force calculations, calculate the finite
# difference Hessian and do relevant vibrational analysis.
# Note that we assume that the analytic gradients we are
# reading as input are in units of Hartree/a.u., although
# I haven't found any confirmation of that in the
# documentation.
#
# Note: this version is programmed to accept output from
# makehess1.2.py.
#*********************************************************

NExternal = 6
useBook = True#Whether to use book or Gaussian website to construct
               #and remove vibrational modes from the molecular
               #degrees of freedom.
#If we are using the book, the following two aren't necessary.
rotXIndicesSwitch1 = True
rotXIndicesSwitch2 = False
mySystem = 3

if (mySystem==0):
    inputDirectory = "p_fdh2o-test01"
    baseFileName = "fdh2o"
    refFileName = inputDirectory+"/_"+baseFileName+".out"
    freqFileName = "h2o-freq-STO3G.out01"
elif (mySystem==1):
    inputDirectory = "p_fdben_test01"
    baseFileName = "fdben"
    refFileName = inputDirectory+"/_"+baseFileName+".out"
    freqFileName = "benzene-freq-STO3G.out01"
elif (mySystem==2):
    inputDirectory = "p_fdacet-test02"
    baseFileName = "fdacet"
    refFileName = inputDirectory+"/_"+baseFileName+".out"
    freqFileName = "acetaldehyde-freq-STO3G.out01"
elif (mySystem==3):
    inputDirectory = "p_fdforce-venus"
    baseFileName = "fdforcXD"
    refFileName = inputDirectory+"/reference.o"
    freqFileName = "OLD_mock_qchem_output.out"

hh = 0.0001 #Finite difference step size

def optionString(option):
    if isinstance(option,bool):
        if (option):
            return "T"
        else:
            return "F"
    return ""

import numpy as np
from numpy import linalg as la
import math
import re
import os
import sys
sys.path.append('/data/home/alguire/inout_files/python')
import qchem as qc
from datetime import datetime

# Set up the log file for running the script
versionString = "1.0"
logFileNameBase = optionString(useBook)
if (not useBook):
    logFileNameBase += optionString(rotXIndicesSwitch1)+optionString(rotXIndicesSwitch2)
logFile = open(logFileNameBase+".log","w")
logFile.write("Beginning run of FDvibrationalAnalysis.py, v"+versionString+". It is ")
timeT = datetime.now()
logFile.write(str(timeT.year)+"-"+str(timeT.month)+"-"+str(timeT.day)+" ")
logFile.write(str(timeT.hour)+":"+str(timeT.minute)+":"+str(timeT.second)+"\n")
logFile.write("Input directory = "+baseFileName+", base input file name = "+baseFileName+"\n")
logFile.write("Finite difference step size given as "+str(hh)+"\n")
logFile.write("Adjustable parameters of interest:\n")
if (useBook):
    logFile.write("useBook = {}\n".format(useBook))
else:
    logFile.write("useBook = {}, rotXIndicesSwitch1 = {}, rotXIndicesSwitch2 = {}\n".format(useBook,rotXIndicesSwitch1,rotXIndicesSwitch2))

# Gather some information about our system
refGeoOutput = open(refFileName,"r")
NAtom, elements, geometry = qc.getGeoData(refGeoOutput)
refGeoOutput.close()
COM = qc.getCenterOfMass(elements,geometry)
COMgeometry = qc.geoTranslate(geometry,-COM)
logFile.write("Reference file given as "+refFileName+"\n")
logFile.write("This script has determined that NAtom = "+str(NAtom)+"\n")
logFile.write("This script has calculated the center of mass to be "+str(COM)+" A, as written in the reference file.\n")

#Set up initial quantities
muX = 0
muY = 1
muZ = 2
NCoord = 3*NAtom
NMode = NCoord -NExternal
NAstringLen = max(len(str(NAtom)),2)
grad1_re = re.compile("    [123]   [- ][0-9]\.[0-9]{14}$",re.MULTILINE)
grad2_re = re.compile("    [123]   [- ][0-9]\.[0-9]{14}   [- ][0-9]\.[0-9]{14}$",re.MULTILINE)
grad3_re = re.compile("    [123]   [- ][0-9]\.[0-9]{14}   [- ][0-9]\.[0-9]{14}   [- ][0-9]\.[0-9]{14}$",re.MULTILINE)
fourt_re  = re.compile("[- ][0-9]\.[0-9]{14}")
atom_re = re.compile("[a-zA-Z]{1,2}\s+[-]?[0-9]+\.[0-9]{10}\s+[-]?[0-9]+\.[0-9]{10}\s+[-]?[0-9]+\.[0-9]{10}")
ele_re = re.compile("[a-zA-Z]{1,2}")
outFileNameSuffix = "_"+baseFileName+".out"
#os.chdir(inputDirectory)
hhFac = 1.0/(2.0*hh*qc.A_to_Bohr) #Finite difference factor (in a.u.)
FCmat = np.zeros([NCoord,NCoord])
TPI = 2*math.pi

#Construct the Hessian
for atomI in range(NAtom):
    atomS = str(atomI+1).zfill(NAstringLen)
    for muI, muS in enumerate(["x", "y", "z"]):
        tempGradArray = np.zeros(NAtom)
        for PorMS in ["p", "m"]:
            myFileName = inputDirectory+"/"+PorMS +muS +atomS +outFileNameSuffix
            myFile = open(myFileName,"r")
            myFileS = myFile.read()
            myFile.close()
            iGradList1 = grad1_re.findall(myFileS)
            iGradList2 = grad2_re.findall(myFileS)
            iGradList3 = grad3_re.findall(myFileS)
            assert(len(iGradList3)+2*len(iGradList2)/3+len(iGradList1)/3==NAtom)
            atomsLeft = NAtom
            for iii in range((NAtom+2)/3):
                NCol = min(3,atomsLeft)
                if (NCol==3):
                    fGradList = fourt_re.findall(iGradList3[iii*3 +muX])
                    fGradList = fGradList +fourt_re.findall(iGradList3[iii*3 +muY])
                    fGradList = fGradList +fourt_re.findall(iGradList3[iii*3 +muZ])
                elif (NCol==2):
                    fGradList = fourt_re.findall(iGradList2[muX])
                    fGradList = fGradList +fourt_re.findall(iGradList2[muY])
                    fGradList = fGradList +fourt_re.findall(iGradList2[muZ])
                elif (NCol==1):
                    fGradList = fourt_re.findall(iGradList1[muX])
                    fGradList = fGradList +fourt_re.findall(iGradList1[muY])
                    fGradList = fGradList +fourt_re.findall(iGradList1[muZ])
                for jjj in range(NCol):
                    if (PorMS=="p"):
                        FCmat[atomI*3 +muI,iii*9 +jjj*3 +muX] = float(fGradList[jjj +muX*NCol])
                        FCmat[atomI*3 +muI,iii*9 +jjj*3 +muY] = float(fGradList[jjj +muY*NCol])
                        FCmat[atomI*3 +muI,iii*9 +jjj*3 +muZ] = float(fGradList[jjj +muZ*NCol])
                    elif (PorMS=="m"):
                        FCmat[atomI*3 +muI,iii*9 +jjj*3 +muX] = FCmat[atomI*3 +muI,iii*9 +jjj*3 +muX] -float(fGradList[jjj +muX*NCol])
                        FCmat[atomI*3 +muI,iii*9 +jjj*3 +muY] = FCmat[atomI*3 +muI,iii*9 +jjj*3 +muY] -float(fGradList[jjj +muY*NCol])
                        FCmat[atomI*3 +muI,iii*9 +jjj*3 +muZ] = FCmat[atomI*3 +muI,iii*9 +jjj*3 +muZ] -float(fGradList[jjj +muZ*NCol])
                atomsLeft -= NCol

qc.symmetryTest(FCmat,True)

# Apply the finite difference factor and convert everything to SI units
FCmat = FCmat*hhFac #FCmat is now in units of Hartree/Bohr^2
FCmat = FCmat*qc.Har_to_J*math.pow(qc.Bohr_to_m,-2)

# Create mass-weighting matrix MWmat_{ij} = (m_i*m_j)^(-1/2)
MWmat = []
for iii in range(NAtom):
    massFac = math.pow(qc.elementMasses[elements[iii]],-0.5)
    MWmat.append([massFac*math.pow(qc.elementMasses[myEle],-0.5) for myEle in elements])
MWmat = np.array(MWmat)
MWmat = MWmat/qc.amu_to_kg

# Symmetrize force constant matrix and create mass-weighted variant
FCmatCopy = np.array(FCmat,copy=True)
for iii in range(NCoord):
    for jjj in range(NCoord):
        FCmat[iii,jjj] = 0.5*(FCmatCopy[iii,jjj]+FCmatCopy[jjj,iii])
FCMWmat = np.zeros([NCoord,NCoord])
for iii in range(NAtom):
    for jjj in range(NAtom):
        massWeight = MWmat[iii,jjj]
        for mui in range(3):
              for muj in range(3):
                FCMWmat[iii*3 +mui,jjj*3 +muj] = massWeight*FCmat[iii*3 +mui,jjj*3 +muj]

FCMWmat2 = qc.massWeightHessian(FCmat,elements,True)
FCMWmat2 = FCMWmat2/qc.amu_to_kg
diff = FCMWmat2-FCMWmat
#print "FCMWmat : ",FCMWmat
#print "FCMWmat2: ",FCMWmat2
#print "diff: ",diff

#print "FCMWmat:",FCMWmat

#Build external modes (in three steps) in the mass-weighted coordinate system
externalModes = np.zeros((NCoord,6))
# Step A: Build translational modes
for iii, myElement in enumerate(elements):
    massFac = math.sqrt(qc.elementMasses[myElement])
    externalModes[iii*3 +muX,muX] = massFac
    externalModes[iii*3 +muY,muY] = massFac
    externalModes[iii*3 +muZ,muZ] = massFac

# Step B: Build rotational modes
if (useBook):
    for iii, myElement in enumerate(elements):
        massFac = math.sqrt(qc.elementMasses[myElement])
        Xi = COMgeometry[iii*3 +muX]
        Yi = COMgeometry[iii*3 +muY]
        Zi = COMgeometry[iii*3 +muZ]
        externalModes[iii*3 +muZ,3 +muX] =  massFac*Yi
        externalModes[iii*3 +muY,3 +muX] = -massFac*Zi
        externalModes[iii*3 +muX,3 +muY] =  massFac*Zi
        externalModes[iii*3 +muZ,3 +muY] = -massFac*Xi
        externalModes[iii*3 +muY,3 +muZ] =  massFac*Xi
        externalModes[iii*3 +muX,3 +muZ] = -massFac*Yi
else:
    momOfInertiaTensor = qc.getMomentOfInertiaTensor(COMgeometry,elements)
    momOfInertia,rotX = la.eig(momOfInertiaTensor)
    rotP = np.zeros((NAtom,3))
    for iii in range(NAtom):
        for mu1 in range(3):
            for mu2 in range(3):
                if (rotXIndicesSwitch1):
                    rotP[iii,mu1] += rotX[mu2,mu1]*COMgeometry[iii*3 +mu2] #NOT SURE: switch indices on rotX?
                else:
                    rotP[iii,mu1] += rotX[mu1,mu2]*COMgeometry[iii*3 +mu2]
    for iii, myElement in enumerate(elements):
        massFac = math.pow(qc.elementMasses[myElement],-0.5)
        for muj in range(3): #NOT SURE: switch indices on rotX?
            if (rotXIndicesSwitch2):
                externalModes[iii*3 +muj,3 +muX] = massFac*(rotP[iii,muY]*rotX[muj,muZ]-rotP[iii,muZ]*rotX[muj,muY])
                externalModes[iii*3 +muj,3 +muY] = massFac*(rotP[iii,muZ]*rotX[muj,muX]-rotP[iii,muX]*rotX[muj,muZ])
                externalModes[iii*3 +muj,3 +muZ] = massFac*(rotP[iii,muX]*rotX[muj,muY]-rotP[iii,muY]*rotX[muj,muX])
            else:
                externalModes[iii*3 +muj,3 +muX] = massFac*(rotP[iii,muY]*rotX[muZ,muj]-rotP[iii,muZ]*rotX[muY,muj])
                externalModes[iii*3 +muj,3 +muY] = massFac*(rotP[iii,muZ]*rotX[muX,muj]-rotP[iii,muX]*rotX[muZ,muj])
                externalModes[iii*3 +muj,3 +muZ] = massFac*(rotP[iii,muX]*rotX[muY,muj]-rotP[iii,muY]*rotX[muX,muj])

# Step C: orthoormalize external modes
externalModes = qc.orthonormalize(externalModes)
#print "FCMWmat.shape: ",FCMWmat.shape
#print "externalModes.shape: ",externalModes.shape

# Project the external modes out from the mass-weighted Hessian
internalDD = qc.normalize(FCMWmat)
for iii in range(6):
    internalDD = qc.projectOut(internalDD,externalModes[:,iii])
internalDD = qc.getPrunedBasis(internalDD,6)
internalDD = qc.orthonormalize(internalDD)

internalMat = np.dot(internalDD.transpose(),np.dot(FCMWmat,internalDD))

#print "internalMat: ",internalMat
eigV,modeV = la.eig(internalMat)
modeV = modeV.transpose()
eigs = zip(eigV,modeV)
#print "eigs: ",eigs
eigs = sorted(eigs)
eigV,modeV = zip(*eigs)
modeV = np.array(modeV)
modeV = modeV.transpose()
#print "modeV: ",modeV

#Get frequencies in cm-1
freqV = [math.sqrt(myEig)/(TPI*qc.SOL*qc.m_to_cm) if myEig>0 else -math.sqrt(-myEig)/(TPI*qc.SOL*qc.m_to_cm) for myEig in eigV]
freqV = np.array(freqV)
#print "freqV: ",freqV

#Get reduced masses in a.m.u.
cartModes = np.dot(internalDD,modeV)
massFacMat = np.zeros((NCoord,NCoord))
for iii,myElement in enumerate(elements):
    massFac = math.pow(qc.elementMasses[myElement],-0.5)
    for mu in range(3):
        massFacMat[iii*3 +mu,iii*3 +mu] = massFac
cartModes = np.dot(massFacMat,cartModes)
redMass = np.zeros(NMode)
for iii in range(NMode):
    redMass[iii] = 1.0/np.dot(cartModes[:,iii],cartModes[:,iii])
#print "redMass: ",redMass
cartModes = qc.normalize(cartModes)
#Get force constants in mDyn/A
FCk = [math.pow(myFreq*qc.SOL*TPI/qc.cm_to_m,2)*myMass*qc.amu_to_kg*math.pow(10,-2) for myFreq,myMass in zip(freqV,redMass)]
#print "FCk: ",FCk

def getModeError(modeTry,freqFile,myBlockIndices,NAtom,NMode):
    normErrors = []
    maxErrors = []
    minErrors = []
    for iii in range(NMode):
        anlMode = qc.getVibQuantity(freqF,myBlockIndices,iii+1,qc.modeString,NAtom,NMode)
        normError = float("inf")
        for jjj in range(2):
            diff = modeTry[:,iii] +((-1)**jjj)*anlMode
            #print "modeTry: ",modeTry[:,iii]
            #print "anlMode: ",anlMode
            #print "mode {} try {}".format(iii,jjj)
            #print "diff: ",diff
            myNorm = la.norm(diff)
            if (normError>myNorm):
                normError = myNorm
                maxError = np.amax(diff)
                minError = np.amin(diff)
        normErrors.append(normError)
        maxErrors.append(maxError)
        minErrors.append(minError)
    return np.array(normErrors),np.array(maxErrors),np.array(minErrors)

#Check against analytic results
if (freqFileName):
    logFile.write("Comparing finite difference results to analytic frequency analysis calculation.\n")
    logFile.write("Frequency output file: "+freqFileName)
    freqF = open(freqFileName,"r")
    myBlockIndices = qc.getVibBlockLocations(freqF, NMode)
    anlFreq, anlRedMass, anlFC = qc.getAllScalarVibQuantities(freqF,myBlockIndices,NAtom,NMode)
    logFile.write("\nRMS absolute error\n")
    freqError = math.sqrt(((freqV-anlFreq)**2).mean()))
    logFile.write("Frequency:      {: f}\n".format(freqError)
    logFile.write("Reduced mass:   {: f}\n".format(math.sqrt(((redMass-anlRedMass)**2).mean())))
    logFile.write("Force constant: {: f}\n\n".format(math.sqrt(((FCk-anlFC)**2).mean())))
    logFile.write("RMS relative error\n")
    logFile.write("Frequency:      {: f}\n".format(math.sqrt((((freqV-anlFreq)/anlFreq)**2).mean())))
    logFile.write("Reduced mass:   {: f}\n".format(math.sqrt((((redMass-anlRedMass)/anlRedMass)**2).mean())))
    logFile.write("Force constant: {: f}\n\n".format(math.sqrt((((FCk-anlFC)/anlFC)**2).mean())))
    modeError,modeMaxDiff,modeMinDiff = getModeError(cartModes,freqF,myBlockIndices,NAtom,NMode)
    logFile.write("Avg. mode error: {: f}\n".format(modeError.mean()))
    logFile.write("Mode norm error: "+str(modeError))
    logFile.write("\nMode max error:  "+str(modeMaxDiff))
    logFile.write("\nMode min error:  "+str(modeMinDiff))
    freqF.close()
logFile.close()



# Create a mock Q-Chem output file.
def getModeBlock(modeList,elements,myBlockIndex):
    outString = ""
    NAtom = len(elements)
    for iii,iEle in enumerate(elements):
        outString += " {:<2}     ".format(iEle)
        for jjj in range(3):
            item1 = modeList[iii*3   ,myBlockIndex*3 +jjj]
            item2 = modeList[iii*3 +1,myBlockIndex*3 +jjj]
            item3 = modeList[iii*3 +2,myBlockIndex*3 +jjj]
            outString +="   {:6.3f} {:6.3f} {:6.3f}".format(item1,item2,item3)
        outString += "\n"
    return outString

geoString = qc.geometryToString(elements,geometry)
geoStrings = geoString.split("\n")
while ("" in geoStrings):
    geoStrings.remove("")
assert(len(geoStrings)==NAtom)
geoString = ""
for iii,iLine in enumerate(geoStrings):
    geoString += "{:>5}      ".format(iii+1)
    geoString += iLine+"\n"

mockFile = open("mock_Q-chem.out","w")
fileIntro  = "--------------------------------------------------------------\n"
fileIntro += " ----------------------------------------------------------------\n"
fileIntro += "             Standard Nuclear Orientation (Angstroms)\n"
fileIntro += "    I     Atom           X                Y                Z\n"
fileIntro += " ----------------------------------------------------------------\n"
fileIntro += geoString
fileIntro += " ----------------------------------------------------------------\n"

vibIntro  = "\n  JobOver = TRUE\n"
vibIntro += " **********************************************************************\n"
vibIntro += " **                                                                  **\n"
vibIntro += " **                       VIBRATIONAL ANALYSIS                       **\n"
vibIntro += " **                       --------------------                       **\n"
vibIntro += " **                                                                  **\n"
vibIntro += " **        VIBRATIONAL FREQUENCIES (CM**-1) AND NORMAL MODES         **\n"
vibIntro += " **     FORCE CONSTANTS (mDYN/ANGSTROM) AND REDUCED MASSES (AMU)     **\n"
vibIntro += " **                  INFRARED INTENSITIES (KM/MOL)                   **\n"
vibIntro += " **                                                                  **\n"
vibIntro += " **********************************************************************\n\n"

mockFile.write(fileIntro)
mockFile.write(vibIntro)

for iBlock in range(NMode/3):
    modeL = iBlock*3
    modeC = iBlock*3 +1
    modeR = iBlock*3 +2
    blockString  = "\n"
    blockString += " Mode:          {:>8}               {:>8}               {:>8}\n".format(modeL+1,modeC+1,modeR+1)
    blockString += " Frequency:     {:8.2f}               {:8.2f}               {:8.2f}\n".format(freqV[modeL],freqV[modeC],freqV[modeR])
    blockString += " Force Cnst:    {:8.4f}               {:8.4f}               {:8.4f}\n".format(FCk[modeL],FCk[modeC],FCk[modeR])
    blockString += " Red. Mass:     {:8.4f}               {:8.4f}               {:8.4f}\n".format(redMass[modeL],redMass[modeC],redMass[modeR])
    blockString += " IR Active:          YES                    YES                    YES\n"
    blockString += " IR Intens:        0.500                  0.500                  0.500\n"
    blockString += " Raman Active:       YES                    YES                    YES\n"
    blockString += "               X      Y      Z        X      Y      Z        X      Y      Z\n"
    blockString += getModeBlock(cartModes,elements,iBlock)
    blockString += " TransDip   0.000  0.000  0.000    0.000  0.000  0.000    0.000  0.000  0.000\n"
    mockFile.write(blockString)
mockFile.close()
sys.exit()


