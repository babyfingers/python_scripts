#!/usr/local/bin/python

class FourTensor:
    """A four-tensor with the symmetry properties found in the
    energy-damped purity method. Also contains information about
    coefficients so that a group of such terms can be combined."""
    def __init__(self, coeff, NCos, NSin, myIndices):
        self.prefac = [[coeff,NCos,NSin]]
        self.NCoeff = len(self.prefac)
        self.indices = list(myIndices)
        assert(len(self.indices)==4)
    def __init__(self, masterList):
        self.prefac = [[masterList[0],masterList[1],masterList[2]]]
        self.NCoeff = len(self.prefac)
        self.indices = list(masterList[3:])
        assert(len(self.indices)==4)
    def indexRank(self,in_indices):
        #Return the integer representation of the
        #indices (i.e., [0,1,1,0] -> 110.
        rank = 0
        for iii,iInd in enumerate(in_indices):
            rank += int(pow(10,3-iii))*iInd
        return -rank
    def reorder(self):
        order0 = list(self.indices)
        order1 = [self.indices[1], self.indices[0], self.indices[3], self.indices[2]]
        order2 = [self.indices[2], self.indices[3], self.indices[0], self.indices[1]]
        order3 = [self.indices[3], self.indices[2], self.indices[1], self.indices[0]]
        orders = [order0, order1, order2, order3]
        bestRank = -2000
        for iOrder in orders:
            iRank = self.indexRank(iOrder)
            if (iRank>bestRank):
                bestRank = iRank
                bestOrder = iOrder
        self.indices = list(bestOrder)
    def getString(self):
        myString = ""
        if (self.NCoeff>1):
            myString += "("
        iii=0
        for iii,iPrefac in enumerate(self.prefac):
            iCoeff, iNC, iNS = iPrefac
            if (iCoeff==-1):
                myString += "-"
            elif (iCoeff!=1):
                myString += str(iCoeff)+"*"
            if (iNC>1):
                myString += "ct{}".format(iNC)
            elif (iNC==1):
                myString += "ct"
            if (iNS>1):
                myString += "st{}".format(iNS)
            elif (iNS==1):
                myString += "st"
            if (iii<self.NCoeff-1 and self.prefac[iii+1][0]>0):
                myString += "+"
            iii += 1
        if (self.NCoeff>1):
            myString += ")"
        myString += "*P({}{}{}{})".format(self.indices[0],self.indices[1],self.indices[2],self.indices[3])
        return myString
    def sameIndices(self,otherTensor):
        if (self.indices==otherTensor.indices):
            return True
        return False
    def setToZero(self):
        self.prefac = []
    def addFrom(self,otherTensor):
        assert(self.sameIndices(otherTensor))
        for iCoeff, iNC, iNS in otherTensor.prefac:
            foundMatch = False
            for jPrefac in self.prefac:
                jCoeff, jNC, jNS = jPrefac
                if (jNC==iNC and jNS==iNS):
                    foundMatch = True
                    jPrefac[0] += iCoeff
                    break
            if (not foundMatch):
                self.prefac.append([iCoeff,iNC,iNS])
                self.NCoeff += 1
        self.sweepZeros()
        self.NCoeff = len(self.prefac)
    def sweepZeros(self):
        zeros = []
        for iii,iPrefac in enumerate(self.prefac):
            if (iPrefac[0]==0):
                zeros.append(iii)
        zeros.reverse()
        for iZ in zeros:
            self.prefac.pop(iZ)
    def isZero(self):
        if (len(self.prefac)==0):
            return True
        return False

def expandTerms(inputIndices, outTerms, myPlace):
    myIndex = inputIndices[myPlace]
    if (myIndex==0 or myIndex==1):
        cIndex = inputIndices[myPlace]
        sIndex = (cIndex +1)%2
        cSign = 1
        sSign = 1
        if (cIndex==1): sSign = -1
        newTerms = []
        for iii,iTerm in enumerate(outTerms):
            newTerms.append(list(iTerm))
            outTerms[iii][0] = iTerm[0]*cSign
            newTerms[iii][0] = iTerm[0]*sSign
            outTerms[iii][1] += 1
            newTerms[iii][2] += 1
            #print "outTerms[iii] = ",outTerms[iii]
            outTerms[iii][3+myPlace] = cIndex
            newTerms[iii][3+myPlace] = sIndex
        outTerms.extend(newTerms)
    else:
        for iii,iTerm in enumerate(outTerms):
            outTerms[iii][0] = iTerm[0]
            outTerms[iii][3+myPlace] = myIndex
    if (myPlace<3):
        expandTerms(inputIndices,outTerms,myPlace+1)
    return outTerms

def tensRank(myTensor):
    return myTensor.indexRank(myTensor.indices)

listGood = False
indList = []
while (len(indList)!=4 and not listGood):
    inString = raw_input("Enter tensor indices (4 single-digit indices; 0 will be rotated with 1):\n")
    #inString = "0000"
    indList = []
    for c in inString:
        indList.append(int(c))
    listGood = True
    for ind in indList:
        if (ind!=0 and ind!=1):
            listGood = False

#Prime the pump, then generate all possible combinations.
termList = [[0]*7]
termList[0][0] = 1
termList = expandTerms(indList, termList, 0)
tensList = [FourTensor(iTerm) for iTerm in termList]
for iTens in tensList:
    iTens.reorder()
zeros = []
for iii,iTens in enumerate(tensList):
    for jjj,jTens in enumerate(tensList[iii+1:]):
        if (iTens.sameIndices(jTens)):
            jTens.addFrom(iTens)
            iTens.setToZero()
            zeros.append(iii)
            break
zeros.reverse()
for iZ in zeros:
    tensList.pop(iZ)
tensList = sorted(tensList,key=tensRank)
tensList.reverse()
line = "The rotated expression is:\n"
for iii,iTens in enumerate(tensList):
    line += iTens.getString()
    if (iii<len(tensList)-1):
        line += "\n"
print line
    
    
