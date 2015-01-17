#!/usr/local/bin/python

import re
import os

NState = 7
getDiabats = True

myStates = range(NState)
float_re = re.compile("[-]?[0-9]+\.[0-9]+")
EHF_re = re.compile("[-]?[0-9]\.[0-9]+e\+[0-9]+")
myDir = raw_input("Directory:\n")

def getEHF():
    os.system("grep \"EHF = \" *.o > ehf.txt")
    f = open("ehf.txt","r")
    fs = f.read()
    f.close()
    os.system("rm ehf.txt")
    ehf = EHF_re.findall(fs)
    for iii,sss in enumerate(ehf):
        ehf[iii] = float(sss)
    return ehf

os.chdir(myDir)
surfaces = []

for iState in myStates:
    if (getDiabats):
        os.system("grep \"x diabatH({0},{0})\" *.o > tmp.txt".format(iState))
    else:
        os.system("grep \"Total energy for state   {0}\" *.o > tmp.txt".format(iState+1))
    f = open("tmp.txt","r")
    fs = f.read()
    f.close()
    surfaces.append(float_re.findall(fs))
    NPoint = len(surfaces[-1])
os.system("rm tmp.txt")

if (getDiabats):
    EHFlist = getEHF()
else:
    EHFlist = [0.0]*NPoint
print "len(EHFlist) = ",len(EHFlist)

for iPoint in range(NPoint):
    line = ""
    for iState in myStates:
       line += "{:.10f} ".format(float(surfaces[iState][iPoint])+EHFlist[iPoint])
    print line
