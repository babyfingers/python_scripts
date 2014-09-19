#!/usr/local/bin/python

# gradcompare.py
import vectorfun
import string
import os.path
import numpy
import math
import time
import sys
import re

#gradcompare.py: given two directories of linterp.py output files, two separate
# greppable strings, and two separate dimensions for the grepped matrix, (in
# some cases I'll have decided to transpose a matrix before printing), will
# compare matrices from corresponding files based on three criteria:
#     1) Comparing the magnitude of the matrices (root of sum of squares of elements)
#     2) Comparing the direction of the matrices (normalize each, find dot product,
#        square result)
#     3) Comparing RMSD (average the squares of the differences of each element, find
#        the root of the resulting value)
#
# Built from lindc1.4.py.

qtype_l = ["mag1", "mag2", "mag1-mag2", "magerror", "direrror", "doterror", "RMSD"]
tabs = "\t"
eols = "\n"
tia_s = "1"
rs = "r"
ws = "w"
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
outfilenameend = '.o'

magnitude = vectorfun.magnitude
normalize = vectorfun.normalize
doterror = vectorfun.doterror
direrror = vectorfun.direrror
RMSD = vectorfun.RMSD
docmd = vectorfun.docmd
lencheck = vectorfun.lencheck

#compvec: compare two sets of vectors
# using various criteria, add output
# to list.
def compvec(vec_l, out_l, nvec):
    v1_l = vec_l[0]
    v2_l = vec_l[1]
    for jjj in range(nvec):
        out_l.append([])
        v1 = v1_l[jjj]
        v2 = v2_l[jjj]
        m1 = vectorfun.magnitude(v1)
        m2 = vectorfun.magnitude(v2)
        d = direrror(v1,v2)
        dd = doterror(v1,v2)
        r = RMSD(v1,v2)
        out_l[-1].append(m1)
        out_l[-1].append(m2)
        out_l[-1].append(m1-m2)
        out_l[-1].append((m1-m2)/math.fabs(m2))
        out_l[-1].append(d)
        out_l[-1].append(dd)
        out_l[-1].append(r)

# Get some user input: data set 1
#filenamebase1 = "c14eeG" #raw_input("Base input 1 file name (max. 6 characters):\n")
filenamebase1 = "c14eeF"

while (len(filenamebase1)>7):
    print "Sorry, that file name is no good. Please try a different one."
    filenamebase1 = raw_input("Base input 1 file name (max. 6 characters):\n")
dir1 = "p_"+filenamebase1+"_t"+tia_s+"_"
dir1 = dir1 +"full" #raw_input("Directory name addition ("+dir1+"):\n")
#grepm1 = "xxDgradCIS1_wetfxx" #raw_input("Greppable string for target matrix:\n")
grepm1 = "xxEgradCISxx"
#transpode1 = 0 #input("Has the target matrix been transposed? (0 no, 1 yes, 2 help):\n")
transpode1 = 1
while (transpode1 != 0 and transpode1 != 1):
    print("Is the column-indexed value the fast variable? Or, maybe more simply,\n" \
          +"is the desired matrix taller than it is wide? If yes, you probably want" \
          +"to choose (0). Otherwise, choose (1).\n")
    transpode1 = input("Has the target matrix been transposed? (0) no, (1) yes:\n")

# Get some user input: data set 2
filenamebase2 = "c14eeF" #raw_input("Base input 2 file name (max. 6 characters):\n")
while (len(filenamebase2)>7):
        print "Sorry, that file name is no good. Please try a different one."
        filenamebase2 = raw_input("Base input 2 file name (max. 6 characters):\n")
dir2 = "p_"+filenamebase2+"_t"+tia_s+"_"
dir2 = dir2 +"full" #raw_input("Directory name addition ("+dir2+"):\n")
grepm2 = "xxDgradCISxx" #raw_input("Greppable string for target matrix:\n")
transpode2 = 1 #input("Has the target matrix been transposed? (0 no, 1 yes, 2 help):\n")
while (transpode2 != 0 and transpode2 != 1):
    print("Is the column-indexed value the fast variable? Or, maybe more simply,\n" \
          +"is the desired matrix taller than it is wide? If yes, you probably want" \
          +"to choose (0). Otherwise, choose (1).\n")
    transpode2 = input("Has the target matrix been transposed? (0) no, (1) yes:\n")

# Get some user input: universal
numkeep = 2  #input("Number of states (numkeep):\n")
matrixM = 138 #input("Number of columns (fast variable):\n")
matrixN = 4  #input("Number of rows (slow variable):\n")

#Check out geometry information
f = open(filenamebase1 +".x",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
atom_l = re.findall(atom_re,f_s)
tnatoms = len(atom_l)
natoms = tnatoms/2
ncoord = 3*natoms
outfn_l = []
outdn_l = []

#Check out reaction coordinate information
f = open(filenamebase1 +".l",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
rc_l = re.findall(num_re,f_s)

docmd("ls " +dir1 +"/*" +outfilenameend +" > tmp1.txt")
docmd("ls " +dir2 +"/*" +outfilenameend +" > tmp2.txt")
f = open("tmp1.txt",rs)
f_s = f.read()
f.close()
file1_l = f_s.split()
f = open("tmp2.txt",rs)
f_s = f.read()
f.close()
file2_l = f_s.split()
file_ll = [file1_l, file2_l]
grepm_l = [grepm1, grepm2]
mygrads_ll = [0,0]
transpode_l = [transpode1, transpode2]
outdata_ll = [] #Store analyzed data here for all files
for fn_i in range(len(file1_l)):
    #Iterates over the list of files provided; each file represents a different geometry.
    fn1 = file1_l[fn_i]
    fn2 = file2_l[fn_i]
    fn_l = [fn1, fn2]
    #Step 1: collect data for this pair of files
    for file_i in range(2):
        transpode = transpode_l[file_i]
        mygrads_l = []
        for iii in range(matrixN):
            mygrads_l.append([])
        fn = fn_l[file_i]
        grepm = grepm_l[file_i]
        docmd("grep " +grepm +" " +fn +" > tmp.txt")
        f = open("tmp.txt",rs)
        f_s = f.read()
        f.close()
        big1_l = dp_re.findall(f_s)
        for i in range(len(big1_l)):
            big1_l[i] = float(big1_l[i])
        if (transpode==0):
            for jjj in range(matrixN):
                for iii in range(matrixM):
                    mygrads_l[jjj].append(big1_l[iii*matrixN +jjj])
        elif (transpode==1):
            for jjj in range(matrixN):
                for iii in range(matrixM):
                    mygrads_l[jjj].append(big1_l[jjj*matrixM +iii])
        for jjj in range(matrixN):
            mygrads_l[jjj] = numpy.array(mygrads_l[jjj])
        mygrads_ll[file_i] = mygrads_l
    #Step 2: do analysis for this pair of files
    outdata_l = []
    compvec(mygrads_ll,outdata_l,matrixN)
    outdata_ll.append(outdata_l)
    #outdata_ll is indexed by [file #][gradient index][quantity type]
    # and is therefore of size nfiles*numkeep2*6 (as of 2013-09-18)
#Step 3: write output for this pair of files
docmd("rm tmp.txt")
for jjj in range(matrixN):
    fname = "compout_" +str(jjj) +".txt"
    ff = open(fname,ws)
    ff.write(tabs)
    for qtype_s in qtype_l:
        ff.write(tabs +qtype_s)
    for file_i in range(len(outdata_ll)):
        ff.write(eols +str(rc_l[file_i]))
        for q_i in range(len(qtype_l)):
            ff.write(tabs +str(outdata_ll[file_i][jjj][q_i]))
    ff.write(eols)
    ff.close()



