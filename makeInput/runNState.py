#!/usr/local/bin/python

charge  = 0
mult    = 1
states = (1, 2, 3, 4, 5, 6,7,8,9,10)
NState0 = 2
NState1 = 10
methodList = ("boys","boysov","betaov","er","erov","ereps","pure")

myBeta = "40" #In 10 eV^-1
dirNameMod = "newDamp"
#dirNameMod = "beta4"
import os
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

dBoysOV = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"LOC_CIS_OV_SEPARATE TRUE\n",mLD_COMMENT:"Auto-generated BoysOV"}
dBoys   = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"",mLD_COMMENT:"Auto-generated Boys"}
dBetaOV = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"LOC_CIS_OV_SEPARATE TRUE\nLOC_BETA "+myBeta+"\n",mLD_COMMENT:"Auto-generated betaOV"}
dBeta   = {mREM_LD_TYPE:"BOYS_CIS_NUMSTATE",mREM_COMMAND:"LOC_BETA "+myBeta+"\n",mLD_COMMENT:"Auto-generated beta"}
dER     = {mREM_LD_TYPE:"ER_CIS_NUMSTATE",mREM_COMMAND:"",mLD_COMMENT:"Auto-generated ER"}
dEROV   = {mREM_LD_TYPE:"ER_CIS_NUMSTATE",mREM_COMMAND:"LOC_CIS_OV_SEPARATE TRUE\n",mLD_COMMENT:"Auto-generated EROV"}
dEReps  = {mREM_LD_TYPE:"ER_CIS_NUMSTATE",mREM_COMMAND:"LOC_BETA "+myBeta+"\nLOC_PEKAR 500\n",mLD_COMMENT:"Auto-generated EReps"}
dPure   = {mREM_LD_TYPE:"PURE_CIS_NUMSTATE",mREM_COMMAND:"LOC_BETA "+myBeta+"\n",mLD_COMMENT:"Auto-generated purity"}
dLD     = {"boysov":dBoysOV,"boys":dBoys,"betaov":dBetaOV,"beta":dBeta,"er":dER,"erov":dEROV,"ereps":dEReps,"pure":dPure}
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
f = open(fileBase+".GEO","r")
geoS = f.read()
f.close()

#Get template information
f = open(fileBase+".TPL","r")
tplS = f.read()
f.close()

#Make changes common to all input files
inS = str(tplS).replace(mGEOMETRY,geoS)
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


dirName = "p_{}-{}-NS{}-{}".format(fileBase,myMethod,NState1,dirNameMod)
try:
    os.mkdir(dirName)
except:
    pass
os.chdir(dirName)

#Make changes specific to each input file
iStates = " ".join([str(states[iii]) for iii in range(NState0-1)])
baseNames = []
for iii,iNState in enumerate(range(NState0,NState1+1)):
    iStates += " "+str(states[iii+NState0-1])
    iInS = str(inS).replace(mLD_NSTATE,str(iNState)).replace(mLD_STATES,iStates)
    iBase = fileBase+str(iNState).zfill(2)
    baseNames.append(iBase)
    f = open(iBase+".i","w")
    f.write(iInS)
    f.close()

print "Finished creating input files. Now let's run them."
qc0 = "ldqchem.csh -t 1 -m 2gb -indir ../p_adia-wB97XD_2 "
for iBase in baseNames:
    command = qc0 +" -in {} -out {} ".format(iBase+".i",iBase+".o")
    os.system(command)
    time.sleep(5)
