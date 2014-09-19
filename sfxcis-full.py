#sfxcis-full.py
#!/usr/local/bin/python
import os.path
import sys
import time
import math
import re

#This script will run inputs for the script sq-full.py for all SF-XCIS excitation combinations.

eols = "\n"

ex_lstl = []
ex_lstr = []
ex1_lstl = [['s1','S2'],['s2','S1'],['s1','S1'],['s2','S2'],['i','S1'],['i','S2'],['s1','A'],['s2','A'],['i','A']]
ex2_lstl = [['i','s1','a','S1'],['i','s2','a','S2'],['i','s1','a','S2'],['i','s2','a','S1'],['I','s2','S1','S2'],['I','s1','S2','S1'],['s1','s2','a','S2'],['s2','s1','a','S1'],['I','s2','A','S2'],['I','s1','A','S1'],['I','s1','A','S2'],['I','s2','A','S1']]
ex3_lstl = [['I','s1','s2','a','S1','S2']]
ex1_lstr = [['s1','S2'],['s2','S1'],['s1','S1'],['s2','S2'],['j','S1'],['j','S2'],['s1','B'],['s2','B'],['j','B']]
ex2_lstr = [['j','s1','b','S1'],['j','s2','b','S2'],['j','s1','b','S2'],['j','s2','b','S1'],['J','s2','S1','S2'],['J','s1','S2','S1'],['s1','s2','b','S2'],['s2','s1','b','S1'],['J','s2','B','S2'],['J','s1','B','S1'],['J','s1','B','S2'],['J','s2','B','S1']]
ex3_lstr = [['J','s1','s2','b','S1','S2']]
ex_lstl.append(ex1_lstl)
ex_lstl.append(ex2_lstl)
ex_lstl.append(ex3_lstl)
ex_lstr.append(ex1_lstr)
ex_lstr.append(ex2_lstr)
ex_lstr.append(ex3_lstr)
nul = len(ex_lstl) #Number of different types of excitations (left side)
nur = len(ex_lstr) #Number of different types of excitations (right side)

# lencheck: Given a list of lists of lists, such
#that the final set of lists are all
#supposed to have the same length, this
#function checks to make sure this is the case.
def lencheck(ulst):
    for l in ulst:
        nli = len(l[0])
        for li in l:
            if (len(li)!=nli):
                print "ERROR! List has inconsistent number of entries.",nli
                print "l = ",l
                break
lencheck(ex_lstl)
lencheck(ex_lstr)

os.system("rm sf-xcis_new.tex")
for uil in range(nul):
    exll = ex_lstl[uil] #The current set of exitations (single, double, triple) on left side.
    nil = len(exll) #The number of different excitations in this set.
    nl = len(exll[0])/2 #What kind of excitations (single -> 1, double -> 2, etc.) on the left side.
    if (nur>uil+2):
        uir_end = uil+2
    else:
        uir_end = nur
    print "uir_end = ",uir_end
    for uir in range(uil,uir_end):
        exlr = ex_lstr[uir] #The current set of excitations (single, double, triple) on right side.
        nir = len(exlr) #The number of different excitations in this set.
        nr = len(exlr[0])/2 #What kind of excitations (single -> 1, double -> 2, etc.) on the right side.
        for il in range(nil):
            indl = exll[il]
            if (uil==uir):
                ir_beg = il #If left and right sides have same number of excitations, don't bother repeating yourself
            else:
                ir_beg = 0 #Otherwise, you aren't repeating yourself.
            for ir in range(ir_beg,nir):
                print "(",il,",",ir,")"
                indr = exlr[ir]
                f = open("tmp.in",'w')
                s = str(nl)+eols+str(nr)+eols+"1"+eols
                for iind in range(2*nl):
                    s += indl[iind] + eols
                for iind in range(2*nr):
                    s += indr[iind] + eols
                s += "1" + eols
                f.write(s)
                f.close()
                os.system("python sq-full.py < tmp.in > tmp.out")


