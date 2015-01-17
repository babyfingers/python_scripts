#!/usr/local/bin/python

charge  = 0
mult    = 1
states = (1, 2, 3, 4, 5, 6)
#states = (1,2)
NState = len(states)
methodList = ("adia","boys","boysov","beta","betaov","er","ereps","pure")

myBetaREM = 0 #In 0.1 eV^-1
betaFac = 0.1
myBeta = myBetaREM*betaFac
myBetaREM = str(myBetaREM)
dirNameMod = "newDamp-01-12"
#dirNameMod = "wB97X"
#dirNameMod = "figure={}".format(myBeta)
import os
import re
import time

mCHARGE = "xxCHARGExx"
mMULT = "xxMULTxx"
mGEOMETRY = "xxGEOMETRYxx"
mREM_LD_TYPE = "xxREM_LD_TYPExx"
mLD_NSTATE = "xxLD_NSTATExx"
mREM_COMMAND = "xxREM_COMMANDxx"
mLD_COMMENT = "xxLD_COMMENTxx"
mLD_STATES = "xxLD_STATESxx"
mList = (mCHARGE,mMULT,mREM_LD_TYPE,mREM_COMMAND,mLD_COMMENT)
dAdia   = {}
dBoysOV = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"LOC_CIS_OV_SEPARATE TRUE\n",mLD_COMMENT:"Auto-generated BoysOV"}
dBoys   = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"",mLD_COMMENT:"Auto-generated Boys"}
dBetaOV = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"LOC_CIS_OV_SEPARATE TRUE\nLOC_BETA "+myBetaREM+"\n",mLD_COMMENT:"Auto-generated betaOV"}
dBeta   = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"LOC_BETA "+myBetaREM+"\n",mLD_COMMENT:"Auto-generated beta"}
dER     = {mREM_LD_TYPE:"ER_CIS_NUMSTATE",mREM_COMMAND:"",mLD_COMMENT:"Auto-generated ER"}
dEReps  = {mREM_LD_TYPE:"ER_CIS_NUMSTATE",mREM_COMMAND:"LOC_BETA "+myBetaREM+"\nLOC_PEKAR 500\n",mLD_COMMENT:"Auto-generated EReps"}
dPure   = {mREM_LD_TYPE:"PURE_CIS_NUMSTATE",mREM_COMMAND:"LOC_BETA "+myBetaREM+"\n",mLD_COMMENT:"Auto-generated purity"}
dLD     = {"adia":dAdia,"boysov":dBoysOV,"boys":dBoys,"betaov":dBetaOV,"beta":dBeta,"er":dER,"ereps":dEReps,"pure":dPure}
dGeneral = {mCHARGE:str(charge),mMULT:str(mult)}

# User prompts
fileBase = raw_input("Base file name:\n")
promptS = "Choose the LD function:"
for iii,iMethod in enumerate(methodList):
    promptS += " "+iMethod+" ("+str(iii+1)+")"
myMethod = input(promptS+"?\n")
myMethod = methodList[myMethod-1]
myDict = dLD[myMethod]
print "You have chosen: {}".format(myMethod)


#Geometry information
atom_re = re.compile("[A-Za-z]{1,2}\s+-?[0-9]*\.[0-9]*\s+-?[0-9]*\.[0-9]*\s+-?[0-9]*\.[0-9]*")
#os.chdir("PYCM_Cross")
os.chdir("my_geometries")
#os.chdir("my_geo_cross")

fileList = os.listdir(".")
NFile = len(fileList)
geoList = [""]*NFile
for iii, iName in enumerate(fileList):
    f = open(iName,"r")
    geoS = f.read()
    f.close()
    atomList = atom_re.findall(geoS)
    geoList[iii] = "\n".join(atomList)
os.chdir("..")

#Get template information
f = open(fileBase+".TPL","r")
tplS = f.read()
f.close()

#Make changes common to all input files
inS = str(tplS)
for iMark in mList:
    if (iMark in inS):
        if (iMark in myDict):
            iSub = myDict[iMark]
        elif (iMark in dGeneral):
            iSub = dGeneral[iMark]
        else:
            printf("Error- no instance of market string {} in template file.".format(iMark))
            continue
        inS = inS.replace(iMark,iSub)


dirName = "p_{}-{}-NS{}-{}".format(fileBase,myMethod,NState,dirNameMod)
try:
    os.mkdir(dirName)
except:
    pass
os.chdir(dirName)

#Make changes specific to each input file
iStates = " ".join([str(states[iii]) for iii in range(NState)])
baseNames = []
for iii,iGeo in enumerate(geoList):
    iInS = str(inS).replace(mLD_NSTATE,str(NState)).replace(mLD_STATES,iStates).replace(mGEOMETRY,iGeo)
    iBase = fileBase+str(iii).zfill(2)
    baseNames.append(iBase)
    f = open(iBase+".i","w")
    f.write(iInS)
    f.close()

print "Finished creating input files. Now let's run them."
if (fileBase=="adia"):
    qc0 = "ldqchem.csh -t 2 -m 4gb -save lots "
else:
    qc0 = "ldqchem.csh -t 1 -m 2gb "
for iii,iBase in enumerate(baseNames):
    #command = qc0 +" -indir ../p_adia-adia-NS6-w97X/p_{} -in {} -out {} ".format(str(iii).zfill(2),iBase+".i",iBase+".o")
    #command = qc0 +" -indir ../p_adia-adia-NS6-test01/p_{} -in {} -out {} ".format(str(iii).zfill(2),iBase+".i",iBase+".o")
    if (fileBase=="adia"):
        command = qc0 +" -outdir p_{} -in {} -out {} ".format(str(iii).zfill(2),iBase+".i",iBase+".o")
    else:
        command = qc0 +" -indir ../p_adia-adia-NS2-newSave/p_{} -in {} -out {} ".format(str(iii).zfill(2),iBase+".i",iBase+".o")
        #command = qc0 +" -indir ../p_adia-adia-NS6-newTest01/p_{} -in {} -out {} ".format(str(iii).zfill(2),iBase+".i",iBase+".o")
    print command
    os.system(command)
    #print command
    time.sleep(1)
