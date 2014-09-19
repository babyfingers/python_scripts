# lindc.py
#!/usr/local/bin/python
import string
import os.path
import math
import time
import sys
import re


# lindc.py: Calculates the magnitude of derivative couplings from
# a geometry input file (.x) and output files containing the full
# derivative coupling vectors. Meant to integrate with inputs/
# outputs of linterp1.1.py.

eols = "\n"
tia_s = "1"
rs = "r"
ws = "w"
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
outfilenameend = '.o'

#docmd: prints a given string to screen
# and executes it as a command.
def docmd(s):
    print ">"+s
    os.system(s)

# Get some user input
filenamebase = raw_input("Base input file name (max. 6 characters):\n")
while (len(filenamebase)>7):
    print "Sorry, that file name is no good. Please try a different one."
    filenamebase = raw_input("Base input file name (max. 6 characters):\n")
dir = "p_"+filenamebase+"_t"+tia_s+"_"
dir = dir +raw_input("Directory name addition ("+dir+"):\n")

#Check out geometry information
f = open(filenamebase +".x",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
atom_l = re.findall(atom_re,f_s)

tnatoms = len(atom_l)
natoms = tnatoms/2
ncoord = 3*natoms
outfn_l = []
outdn_l = []
direc_v = []

#Create normalized directional vector
normin = 0.0
for j in range(natoms):
    l0 = re.findall(num_re,atom_l[j])
    l1 = re.findall(num_re,atom_l[j+natoms])
    for c in range(3):
        direc_v.append(float(l1[c]) -float(l0[c]))
        normin = normin +direc_v[-1]*direc_v[-1]
normin = 1.0/math.sqrt(normin)
for i in range(ncoord):
    direc_v[i] = direc_v[i]*normin

#Obtain derivative coupling vectors
print "Entering the designated directory, dir = "+dir
os.chdir(dir)

docmd("ls *" +outfilenameend+ " > tmp.txt")
f = open("tmp.txt",rs)
f_s = f.read()
f.close()
file_l = f_s.split()
print "file list = ",file_l

nactm_l = ["xxDerCoup_noetfxx","xxDerCoup_wetfxx","xxDiaCoup_noetfxx","xxDiaCoup_wetfxx"]
nact_l = [] #Where I will store the strings representing the directional DC magnitudes

for fn in file_l:
    for dc in nactm_l:
        docmd("grep "+dc+" "+fn+" > tmp.txt")
        f = open("tmp.txt",rs)
        f_s = f.read()
        f.close()
        big_l = dp_re.findall(f_s)
        if (len(big_l) != 4*ncoord):
            print "Error: not seeing the right number of numbers in grepped DC for file name "+fn
            print "Skipping this file."
            nact_l.append(str(0))
        else:
            for i in range(len(big_l)):
                big_l[i] = float(big_l[i])
            temp = 0.0
            for i in range(ncoord):
                dc_f = big_l[4*i +1]
                temp = temp +dc_f*direc_v[i]
            nact_l.append(str(temp))
docmd("rm tmp.txt")
f = open("DCs.txt",ws)
for dc in nactm_l:
    f.write("\t"+dc+"\t\t")
for fn_i in range(len(file_l)):
    f.write(eols+file_l[fn_i]+":")
    for dc_i in range(len(nactm_l)):
        f.write("\t"+nact_l[fn_i*4 +dc_i])
f.close()



