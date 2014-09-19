#!/usr/local/bin/python

#Do a 1-D random walk and print the simple random walk
#resulting from keeping a running sum of the value
#of a random variable Z that is -1 or +1 with equal
#probability.

# Random walk: random variable Z is -1 with P = 0.5
#                                        1 with P = 0.5
# Simple random walk: random variable Sn = \sum_{i=1}^{n} Z_i
# so that  <S_n>   = 0
#          <S_n^2> = n


import random
import qchem as qc
import cProfile

NTrial = 10000
NStep = 5000
myFileName = "randomWalk.txt"
getNewData = True

#Print Sn to a file delimeted by carriage returns.
def writeNewData(fileName,NStep,NTrial):
    outputFile = open(fileName,"w")
    iTrial = 0
    outString = ""
    for iTrial in range(NTrial):
        myS = 0
        for iStep in range(NStep):
            myS += 2*random.randint(0,1) -1
        outString += str(myS)+"\n"
        if (iTrial%50==0):
            outputFile.write(outString)
            outString = ""
    outputFile.write(outString)
    outputFile.close()

#Read and analyze the data
def analyzeData(fileName):
    expectS  = qc.coarseGrainAverage()
    expectS2 = qc.coarseGrainAverage()
    dataFile = open(fileName,"r")
    for myLine,line in enumerate(dataFile):
        myS = int(line)
        expectS.addQuantity(myS)
        expectS2.addQuantity(myS*myS)
        if (myLine%10000==0):
            expectS.consolidateBins()
            expectS2.consolidateBins()
            print "myLine = {}".format(myLine)
    ES = expectS.getAverage()
    ES2 = expectS2.getAverage()
    print "ES = {}, ES2 = {}".format(ES,ES2)

if (getNewData):
    writeProfile = cProfile.Profile()
    writeProfile.enable()
    writeNewData(myFileName,NStep,NTrial)
    writeProfile.disable()
    writeProfile.print_stats()
    #cProfile.run(writeNewData(NStep,myFileName))
#cProfile.run(analyzeData(myFileName))
analyzeProfile = cProfile.Profile()
analyzeProfile.enable()
analyzeData(myFileName)
analyzeProfile.disable()
analyzeProfile.print_stats()
