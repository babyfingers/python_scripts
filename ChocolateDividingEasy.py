# COmpleted 2014-10-19, 11:00 


# TIMINGS:
# Thinking time: 6 minutes
# Writing time: 36 minutes
# Debug time: 26 minutes
# Total time: 68 minutes

# ANALYSIS:
# I had an idea for a solution quickly enough, although in the execution
# stage I spent some time screwing around with getting numpy on my windows
# machine, which I should have done beforehand. (There was 20 minutes
# beween writing and debug that I didn't count where I was only trying
# to do this, as well as some counted time in the writing stage).
# More practice will probably lead to a faster writing stage.
# Finally, the bugs all had to do with indices and labels. I need to
# be more focused when I try this so that I don't end up cycling
# through a bunch of dumb mistakes.

# Div. 1 success rate: 86.60%
##Mirosz adores sweets. He has just bought a rectangular bar of chocolate. The bar is divided into a grid of square cells. Different cells may have a different quality. You are given the description of the bar in a String[] chocolate. Each character in chocolate is a digit between '0' and '9', inclusive: the quality of one of the cells.
##
##
##
##Mirosz is now going to divide the chocolate into 9 parts: one for him and one for each of his 8 friends. He will do the division by making four cuts: two horizontal and two vertical ones. Each cut must go between two rows or columns of cells. Each of the 9 parts must be non-empty. The quality of a part is the sum of the qualities of all cells it contains.
##
##
##
##Mirosz is well-mannered and he will let his friends choose their pieces first. His friends are even more addicted to chocolate than he is. Therefore, they will certainly choose the pieces with higher quality first, and Mirosz will be left with the worst of the nine pieces.
##
##
##
##You are given the String[] chocolate. Find the optimal places for the four cuts. More precisely, compute and return the largest possible quality of Mirosz's part of the chocolate bar.
##Definition
##    	
##Class:	ChocolateDividingEasy
##Method:	findBest
##Parameters:	String[]
##Returns:	int
##Method signature:	int findBest(String[] chocolate)
##(be sure your method is public)
##    
## 
##Constraints
##-	chocolate will contain between 3 and 50 elements, inclusive. 
##
##-	All elements in chocolate will contain between 3 and 50 characters, inclusive. 
##
##-	All elements in chocolate will contain the same number of characters. 
##
##-	All elements in chocolate will contain only digits ('0'-'9').
import numpy as np

def strToArray(stringList):
    outList = []*len(stringList)
    for s in stringList:
        outList.append([])
        for c in s:
            outList[-1].append(int(c))
    return np.array(outList)

class chocoPart:
    def __init__(self, array):
        self.array = np.array(array,copy=True)
        print "We're in business! chocoPart with"
        print self.array
        self.h1 = 1
        self.h2 = 2
        self.v1 = 1
        self.v2 = 2
        self.hmax, self.vmax = np.shape(self.array)
        #self.hmax -= 1
        #self.vmax -= 1
        print "hmax = {}, vmax = {}".format(self.hmax,self.vmax)
        self.divVals = np.zeros(9)
    def countVals(self):
        hBounds = [0,self.v1,self.v2,self.vmax]
        #print "hBounds",hBounds
        vBounds = [0,self.h1,self.h2,self.hmax]
        #print "vBounds",vBounds
        for col in range(3):
            lB = vBounds[col]
            rB = vBounds[col+1]
            for row in range(3):
                tB = hBounds[row]
                bB = hBounds[row+1]
                #print "subspace: "
                #print self.array[lB:rB,tB:bB]
                self.divVals[col*3 +row] = np.sum(self.array[lB:rB,tB:bB])
                #print "sum: {}".format(self.divVals[col*3 +row])
    def setPart(self,HorV,IorO,val):
        if (HorV==0):
            assert(val>0 and val<=self.hmax)
        else:
            assert(val>0 and val<=self.vmax)
        if (HorV==0 and IorO==0):
            self.h1 = val
        elif (HorV==0 and IorO==1):
            self.h2= val
        elif(HorV==1 and IorO==0):
            self.v1 = val
        elif(IorO==1):
            self.v2 = val
    def getMinVal(self):
        #print "in getMinVal"
        self.countVals()
        #print "divVals: ",self.divVals
        #print "returning {}".format(np.min(self.divVals))
        return np.min(self.divVals)
    def getParts(self):
        return (self.h1, self.h2, self.v1, self.v2)

def doProblem(myString):
    print "********************************************************"
    myArray = strToArray(myString)
    myPart = chocoPart(myArray)
    chocoMax = -1
    maxParts = (0,0,0,0)
    for h1 in range(1,myPart.hmax-1):
        myPart.setPart(0,0,h1)
        for v1 in range(1,myPart.vmax-1):
            myPart.setPart(1,0,v1)
            for h2 in range(h1+1,myPart.hmax):
                myPart.setPart(0,1,h2)
                for v2 in range(v1+1,myPart.vmax):
                    myPart.setPart(1,1,v2)
                    currMin = myPart.getMinVal()
                    if (currMin>chocoMax):
                        chocoMax = currMin
                        maxParts = myPart.getParts()
                        print "New chocoMax = {}, with h1|h2|v1|v2 = {}|{}|{}|{}".format(chocoMax,maxParts[0],maxParts[1],maxParts[2],maxParts[3])
    print "Final answer: {}, with h1|h2|v1|v2 = {}|{}|{}|{}".format(chocoMax,maxParts[0],maxParts[1],maxParts[2],maxParts[3])

s1 = (
"9768",
"6767",
"5313"
)
doProblem(s1)
s2 = (
    "36753562",
    "91270936",
    "06261879",
    "20237592",
    "28973612",
    "93194784"
    )
doProblem(s2)
s3 = (
    "012",
    "345",
    "678"
    )
doProblem(s3)
