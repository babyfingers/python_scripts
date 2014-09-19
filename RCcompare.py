#!/usr/local/bin/python

# RCcompare.py
import vectorfun
from vectorfun import atom_re, ele_re, num_re, int_re, dp_re
import string
import os.path
import numpy
import math
import time
import sys
import re

#RCcompare.py
#How well do gradient quantities match up with the reaction coordinate?

qtype_l = ["direrror"]
tabs = "\t"
eols = "\n"
tia_s = "1"
rs = "r"
ws = "w"
outfilenameend = '.o'

magnitude = vectorfun.magnitude
normalize = vectorfun.normalize
doterror = vectorfun.doterror
direrror = vectorfun.direrror
RMSD = vectorfun.RMSD
docmd = vectorfun.docmd
lencheck = vectorfun.lencheck

#dircompvec: compare two sets of vectors
# using direrror, add output
# to list.
def dircompvec(vec_l, nvec):
    v1_l = vec_l[0]
    v2_l = vec_l[1]
    out_l = []
    for iii in range(nvec):
        out_l.append([])
        for jjj in range(nvec):
            out_l[-1].append([])
            v1 = v1_l[iii]
            v2 = v2_l[jjj]
            d = direrror(v1,v2)
            out_l[iii][jjj].append(d)
    return out_l

#dircompR: compare a gradient vector
# to a single vector, e.g. the
# reaction coordinate normal vector.
def dircompR(Rvec, vec_ll, nvec):
    out_l = []
    for grep_i in range(len(vec_ll)):
        vec_l = vec_ll[grep_i]
        out_l.append([])
        for iii in range(nvec):
            out_l[-1].append([])
            for jjj in range(nvec):
                v = vec_l[iii*nvec +jjj]
                d = direrror(Rvec,v)
                out_l[-1][-1].append(d)
    return out_l
            
# Get some user input: data set
#filenamebase1 = "c14eeG" #raw_input("Base input 1 file name (max. 6 characters):\n")
filenamebase1 = "c14eeF"

while (len(filenamebase1)>7):
    print "Sorry, that file name is no good. Please try a different one."
    filenamebase1 = raw_input("Base input 1 file name (max. 6 characters):\n")
dir1 = "p_"+filenamebase1+"_t"+tia_s+"_"
dir1 = dir1 +"full" #raw_input("Directory name addition ("+dir1+"):\n")
#grepm1 = "xxDgradCIS1_wetfxx" #raw_input("Greppable string for target matrix:\n")
grepm1 = "xxEgradxx"
grepm2 = "xxDgradxx"
#transpode1 = 0 #input("Has the target matrix been transposed? (0 no, 1 yes, 2 help):\n")
transpode1 = 1
while (transpode1 != 0 and transpode1 != 1):
    print("Is the column-indexed value the fast variable? Or, maybe more simply,\n" \
          +"is the desired matrix taller than it is wide? If yes, you probably want" \
          +"to choose (0). Otherwise, choose (1).\n")
    transpode1 = input("Has the target matrix been transposed? (0) no, (1) yes:\n")

# Get some user input: universal
numkeep = 2  #input("Number of states (numkeep):\n")
matrixM = 138 #input("Number of columns (fast variable):\n")
matrixN = numkeep*numkeep  #input("Number of rows (slow variable):\n")

#Check out geometry information
f = open(filenamebase1 +".x",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
atom_l = re.findall(atom_re,f_s)
tnatom = len(atom_l)
natom = tnatom/2
ncoord = 3*natom
outfn_l = []
outdn_l = []
ele_l = []

#Check out reaction coordinate information
# Part 1: what points between \zeta=0 and \zeta=1 to use?
f = open(filenamebase1 +".l",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
rc_l = re.findall(num_re,f_s)
nfile = len(rc_l)
# Part 2: what does our reaction coordinate look like in configuration space?
l0 = [] #Initial geometry vector
l1 = [] #Final geometry vector
R_v = [] #Reaction coordinate in configuration space
for jjj in range(natom): #For each atom line...
    ele_l = ele_l +ele_re.findall(atom_l[jjj])
    l0 = l0 +num_re.findall(atom_l[jjj]) #Extract information from the initial geometry
    l1 = l1 +num_re.findall(atom_l[jjj +natom]) #Extract information from the final geometry
    for mu in range(-3,0):
        l0[mu] = float(l0[mu])
        l1[mu] = float(l1[mu])
        R_v.append(l1[mu]-l0[mu])
R_v = numpy.array(R_v)
R_v = normalize(R_v) #Store directional vector for reaction coordinate

docmd("ls " +dir1 +"/*" +outfilenameend +" > tmp1.txt")
f = open("tmp1.txt",rs)
f_s = f.read()
f.close()
file1_l = f_s.split()
file_l = file1_l
mygrads_ll = [0,0]
grepm_l = [grepm1, grepm2]
transpode = transpode1
# gradover_ll[file #][state pair index I][state pair index J][quantity index]
gradover_ll = [] #Store analyzed data here for all files
# gradR_ll = [file #][grep marker index][state index I][state index J]
gradR_ll = []
for file_i in range(nfile):
    #Iterates over the list of files provided; each file represents a different geometry.
    fn1 = file_l[file_i] #File of interest from directory 1
    fn2 = file_l[file_i] #File of interest from directory 2
    fn_l = [fn1, fn2] #List containing both files of interest for this loop
    #Step 1: collect data for this pair of files
    for grep_i in range(2): #Loop over whether we are interested in the file from directory 1 or 2
        grepm = grepm_l[grep_i]
        mygrads_l = []
        for iii in range(matrixN): #Make space in the output array
            mygrads_l.append([])
        fn = fn_l[grep_i]
        docmd("grep " +grepm +" " +fn +" > tmp.txt")
        f = open("tmp.txt",rs)
        f_s = f.read()
        f.close()
        big1_l = dp_re.findall(f_s)
        for i in range(len(big1_l)):
            big1_l[i] = float(big1_l[i])
        if (transpode==0):
            for jjj in range(matrixN):
                for mmm in range(matrixM):
                    mygrads_l[jjj].append(big1_l[mmm*matrixN +jjj])
        elif (transpode==1):
            for jjj in range(matrixN):
                for mmm in range(matrixM):
                    mygrads_l[jjj].append(big1_l[jjj*matrixM +mmm])
        for jjj in range(matrixN):
            mygrads_l[jjj] = numpy.array(mygrads_l[jjj])
        mygrads_ll[grep_i] = mygrads_l
    #Step 2: do analysis for this pair of files
    outdata_l = dircompvec(mygrads_ll,matrixN)
    gradover_ll.append(outdata_l)
    outdata_l = dircompR(R_v, mygrads_ll, numkeep)
    gradR_ll.append(outdata_l)
#Step 3: write output for this pair of files
docmd("rm tmp.txt tmp1.txt")
for iii in range(matrixN):
    for jjj in range(matrixN):
        fname = "dircompout_" +str(iii) +"-" +str(jjj) +".txt"
        ff = open(fname,ws)
        ff.write(tabs)
        for qtype_s in qtype_l:
            ff.write(tabs +qtype_s)
        for file_i in range(nfile):
            ff.write(eols +str(rc_l[file_i]))
            for q_i in range(len(qtype_l)):
                ff.write(tabs +str(gradover_ll[file_i][iii][jjj][q_i]))
        ff.write(eols)
        ff.close()

for grep_i in range(len(grepm_l)):
    grepm = grepm_l[grep_i]
    fname = "dircompR"+grepm+".txt"
    ff = open(fname,ws)
    ff.write(grepm)
    for iii in range(numkeep):
        for jjj in range(numkeep):
            ff.write(tabs +str(iii)+"-"+str(jjj))
    for file_i in range(nfile):
        ff.write(eols +str(rc_l[file_i]))
        for iii in range(numkeep):
            for jjj in range(numkeep):
                ff.write(tabs +str(gradR_ll[file_i][grep_i][iii][jjj]))
    ff.write(eols)
    ff.close()

