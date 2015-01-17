#!/usr/local/bin/python

import os
import sys
import time

templateBase = raw_input("What is the name of the template file (without .i extension)?\n")
beta0 = 0
beta1 = 4000#50
NBeta = 11#26
betas = [int(iii*(beta1-beta0)/(NBeta-1)) for iii in range(NBeta)]

f = open(templateBase+".i","r")
templateString = f.read()
f.close()

needNewDir = True
dirInd = 1
while (needNewDir):
    try:
        dirName = "p_beta-"+templateBase+"_"+str(dirInd).zfill(2)
        os.mkdir("p_beta-"+templateBase+"_"+str(dirInd).zfill(2))
        needNewDir = False
        break
    except OSError:
        print "Well, {} aldeady exists. Let's try again.".format(dirName)
        dirInd += 1
os.chdir(dirName)

baseNames = [""]*NBeta

for iii,iBeta in enumerate(betas):
    iString = "#Step {}/{}. Beta range = [{},{}].\n\n".format(iii+1,NBeta,beta0,beta1)
    iString += templateString.replace("xxBETAxx",str(iBeta))
    iBase = templateBase+"-"+str(iii).zfill(2)
    baseNames[iii] = iBase
    f = open(iBase+".i","w")
    f.write(iString)
    f.close()

print "Finished creating input files. Now let's run them."

qc0 = "ldqchem.csh -t 1 -m 2gb -indir ../p_adia-wB97XD "

for iBase in baseNames:
    command = qc0 +" -in {} -out {} ".format(iBase+".i",iBase+".o")
    os.system(command)
    time.sleep(5)
