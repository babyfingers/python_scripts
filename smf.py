#!/usr/local/bin/python

# smf.py
import numpy
import string
import os.path
import math
import time
import sys
import re


# smf.py: Shervin Mathematica formatter. Take geometries and vector gradient quantities
# (such as derivative couplings and Hamiltonian gradients) and put them in the format
# that Shervin's .nb file likes so that I can make quiver plots that look nice.

tabs = "\t"
eols = "\n"
tia_s = "1"
rs = "r"
ws = "w"

atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
ele_re = re.compile("[a-zA-Z]{1,2}")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
geo_re = re.compile("\$molecule[^$.]*\$end", re.DOTALL) #Grab the MOLECULE block from input files
outfilenameend = '.o'

#docmd: prints a given string to screen
# and executes it as a command.
def docmd(s):
    #print ">"+s
    os.system(s)

#lencheck: checks that len(mylist)==correct_len,
# otherwise returns 0 and prints a message.
def lencheck(mylist, correct_len, filename_s):
    if (len(mylist) != correct_len):
        print "Error: not seeing the right number of numbers in grepped DC for file name "+filename_s
        print "Skipping this file."
        return 0
    else:
        return 1

# Get some user input
oneormany = input("Do calculation on single file (0) or on output directory from linterp.py (1)?:\n")
if (oneormany==1):
    filenamebase = raw_input("Base input file name (max. 6 characters):\n")
    while (len(filenamebase)>7):
        print "Sorry, that file name is no good. Please try a different one."
        filenamebase = raw_input("Base input file name (max. 6 characters):\n")
    dir = "p_"+filenamebase+"_t"+tia_s+"_"
    dir = dir +raw_input("Directory name addition ("+dir+"):\n")
else:
    filenamebase = raw_input("File name base:\n")
    filename = filenamebase +"." +raw_input("Output file name suffix:\n")
    file_l = [filename]
    rc_l = ["0.0"]
    
numkeep = input("Number of states (numkeep):\n")
numkeep2 = numkeep*numkeep

#Check out geometry information
geo_l = []
if (oneormany==1):
    f = open(filenamebase +".x",rs)
else:
    f = open(filenamebase +".in",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
if (oneormany==0):
    tmp_l = geo_re.findall(f_s)
    atom_l = atom_re.findall(tmp_l[0])
    ele_l = ele_re.findall(tmp_l[0])
    geo_l = list(atom_l)
    natom = len(atom_l)
    tnatom = 2*natom
else:
    ele_l = []
    atom_l = re.findall(atom_re,f_s)
    tnatom = len(atom_l)
    natom = tnatom/2
    for q in range(natom):
        tmp_l = ele_re.findall(atom_l[q])
        ele_l.append(tmp_l[0])
ncoord = 3*natom

ele_d = {"h": 1, "c": 1, "o": 1}
for q in range(natom):
    tmp_ele = ele_l[q].lower()
    ele_l[q] = tmp_ele +ele_d[tmp_ele]
    ele_d[tmp_ele] = ele_d[tmp_ele] +1

outfn_l = []
outdn_l = []
direc_v = []

if (oneormany==1):
    #Obtain derivative coupling vectors
    print "Entering the designated directory, dir = "+dir
    os.chdir(dir)

    docmd("ls *" +outfilenameend+ " > tmp.txt")
    f = open("tmp.txt",rs)
    f_s = f.read()
    f.close()
    outfile_l = f_s.split()
    nout = len(outfile_l)
    print "outfile list = ",outfile_l

    docmd("ls *" +infilenameend+ " > tmp.txt")
    f = open("tmp.txt",rs)
    f_s = f.read()
    f.close()
    infile_l = f_s.split()
    nin = len(infile_l)
    
    for infilename in infile_l:
        f = open(infilename,rs)
        f_s = f.read()
        f.close()
        tmp_l = geo_re.findall(f_s)
        geo_l.append(atom_re.findall(tmp_l[0]))
    print "len(geo_l): "+str(len(geo_l))
    print "# in files: "+str(nin)
    
grads_l = ["xxDerCoup_wetfxx","xxDiaCoup_wetfxx","xxEgradCISxx","xxDgradCISxx"]

#Construct the empty gradient
empty_gradient = []
for q in range(natom):
    empty_gradient.append([])
    for mu in range(3):
        empty_gradient[-1].append(0.0)
    
for num1 in range(numkeep):
    for num2 in range(numkeep):
        for fn in outfile_l:
            for dc in grads_l:
                nact_l = list(empty_gradient) #Storage for the derivative coupling values
                docmd("grep "+dc+" "+fn+" > tmp.txt")
                f = open("tmp.txt",rs)
                f_s = f.read()
                f.close()
                big1_l = dp_re.findall(f_s))
                if (lencheck(big1_l,numkeep2*ncoord, fn)==0):
                    nact_l.append(999)
                else:
                    for i in range(len(big1_l)):
                        big1_l[i] = float(big1_l[i])
                    for q in range(natom):
                        for mu in range(3):
                            if ("grad" in dc):
                                nact_l[q][mu] = big1_l[num1*numkeep*ncoord +num2*ncoord +q*3 +mu] #This is the part that works only for numkeep = 2
                            ekse:
                                nact_l[q][mu] = big1_l[q*numkeep2*3 +mu*numkeep2 +num1*numkeep +num2] #This is the part that works only for numkeep = 2
                docmd("rm tmp.txt")
                fname = fn[0:-2]+"-"+str(num1)+"|"+str(num2)+"_"+dc+".txt"
                ff = open(fname,ws)
                #Write the geometry
                for q in range(natom):
                    ff.write(ele_l[q] +" = {")
                    for mu in range(3):
                        ff.write(str(geo_l[q*3 +mu]))
                        if (mu<2):
                            ff.write(", ")
                    ff.write("};\n")
                ff.write(eols)
                #Write the DC
                for q in range(natom):
                    ff.write("tv" +ele_l[q] +" = {")
                    for mu in range(3):
                        ff.write(str(nact_l[q][mu]))
                        if (mu<2):
                            ff.write(", ")
                    ff.write("};\n")
                ff.write(eols)
                ff.close()


