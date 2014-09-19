#!/usr/local/bin/python

# curve.py
import string
import os.path
import math
import time
import sys
import re

#MYSTERY HACK
# Right now, I want a fast way to compute curvatures around the
# ET minima in the directions of the various gradient vectors.
# Let's use this script to generate that. Print out ".x" files,
# with initial and final geometries 1 a.u. apart, in the direction
# of the gradient vectors. This way, we can repurpose the linter.py
# scripts to run a lot of calcs at once.

tabs = "\t"
eols = "\n"
tia_s = "1"
rs = "r"
ws = "w"
ele_re = re.compile("[a-zA-Z]{1,2}")
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
outfilenameend = '.o'

#print2d
def print2d(v, N, M):
    s = "Running print2d:\n"
    for iii in range(N):
        for jjj in range(M):
            s = s +str(v[iii*M +jjj])+" "
        s = s +eols
    s = s +eols
    print s

#magnitude: returns the 2-norm of a list of
# floats.
def magnitude(v, n):
    norm = 0.0
    for iii in range(n):
        norm = norm +v[iii]*v[iii]
    return math.sqrt(norm)

#normalize: takes list of n doubles, normalizes it
# as if it were an n-vector.
def normalize(v, n):
    normin = 1.0/magnitude(v,n)
    for iii in range(n):
        v[iii] = v[iii]*normin
    

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
oneormany = 1
if (oneormany==1):
    filenamebase = raw_input("Base input file name (max. 6 characters):\n")
    while (len(filenamebase)>7):
        print "Sorry, that file name is no good. Please try a different one."
        filenamebase = raw_input("Base input file name (max. 6 characters):\n")
    dir = "p_"+filenamebase+"_t"+tia_s+"_"
    dir = dir +raw_input("Directory name addition ("+dir+"):\n")
    dot_mode = 2 #input("Compute directional magnitude (1), simple magnitude (2), directional overlap matrices (3), or the works (0)?:\n")
    myfileno = raw_input("File number of interest:\n")
    myfileno = myfileno.zfill(2)
numkeep = input("Number of states (numkeep):\n")
numkeep2 = numkeep*numkeep

#Check out geometry information
if (oneormany==1):
    f = open(filenamebase +".x",rs)
else:
    f = open(filenamebase +".in",rs)
f_s = f.read()
f.close()
f_s = eols +f_s    
atom_l = re.findall(atom_re,f_s) #List of atomic data lines from the ".x" file (for oneormany==1)
if (oneormany==1):
    tnatom = len(atom_l)
    natom = tnatom/2
else:
    natom = len(atom_l)
    tnatom = 2*natom
ncoord = 3*natom

outfn_l = []
outdn_l = []
direc_v = []
ele_l = [] #List of elemental symbols in molecule
l0 = [] #Initial geometry vector
l1 = [] #Final geometry vector
if (oneormany==1):
    #Check out reaction coordinate information
    f = open(filenamebase +".l",rs)
    f_s = f.read()
    f.close()
    f_s = eols +f_s
    rc_l = re.findall(num_re,f_s)

    #Create normalized sdirectional reaction vector, save geometry
    #information.
    for j in range(natom): #For each atom line...
        ele_l = ele_l +ele_re.findall(atom_l[j])
        l0 = l0 +num_re.findall(atom_l[j]) #Extract information from the initial geometry
        l1 = l1 +num_re.findall(atom_l[j+natom]) #Extract information from the final geometry
        for c in range(-3,0):
            l0[c] = float(l0[c])
            l1[c] = float(l1[c])
            direc_v.append(l1[c]-l0[c])
    normalize(direc_v,ncoord) #Store directional vector for reaction coordinate

    #Obtain derivative coupling vectors
    print "Entering the designated directory, dir = "+dir
    os.chdir(dir)

    docmd("ls *" +outfilenameend+ " > tmp.txt")
    f = open("tmp.txt",rs)
    f_s = f.read()
    f.close()
    file_l = f_s.split()
    print "file list = ",file_l
    for file in file_l:
        if myfileno in file:
            file_l = [file]
            break
    print "new file list = ",file_l

nactm_l = ["xxDerCoup_wetfxx"]
nact_v = []
for iii in range(ncoord):
    nact_v.append(0.0)

dotiters = 1
if (dot_mode==0): dotiters = 2
num1range = [0] #For the purposes of this hack, we're only interested in the diagonal
num2range = [1] # element of H_{AB}^{[Q]} for a two-state system.
    
if (dot_mode==1 or dot_mode==2 or dot_mode==0):
    for dotiter in range(dotiters):
        for num1 in num1range:
            for num2 in num2range:
                nact_l = [] #Where I will store the strings representing the directional DC magnitudes
                for fn in file_l:
                    for dc in nactm_l:
                        docmd("grep "+dc+" "+fn+" > tmp.txt")
                        f = open("tmp.txt",rs)
                        f_s = f.read()
                        f.close()
                        big1_l = dp_re.findall(f_s)
                        if (lencheck(big1_l,numkeep2*ncoord, fn)==0):
                            nact_l.append(str(0))
                        else:
                            for i in range(len(big1_l)):
                                big1_l[i] = float(big1_l[i])
                            temp = 0.0
                            for i in range(ncoord):
                                dc_f = big1_l[numkeep2*i +num1*numkeep +num2] #This is the part that works only for numkeep = 2
                                nact_v[i] = float(dc_f)
                                if (dot_mode==1 or (dot_mode==0 and dotiter==0)):
                                    temp = temp +dc_f*direc_v[i]
                                elif (dot_mode==2 or (dot_mode==0 and dotiter==1)):
                                    temp = temp +dc_f*dc_f
                            if (dot_mode==2 or (dot_mode==0 and dotiter==1)): temp = math.sqrt(temp)
                            nact_l.append(str(temp))
                        nact_v_copy = list(nact_v)
                        normalize(nact_v,ncoord)
                #print2d(nact_v,natom,3)
                docmd("rm tmp.txt")
                for minchoice in range(2):
                    if (minchoice == 0):
                        fname = "DCD"+myfileno+".x"
                        geo_v = l1
                    else:
                        fname = "DCA"+myfileno+".x"
                        geo_v = l0
                    ff = open(fname,ws)
                    for coeff in range(2):
                        for q in range(natom):
                            ff.write(ele_l[q])
                            for mu in range(3):
                                ff.write(" " +str(geo_v[q*3 +mu] +coeff*nact_v[q*3 +mu]))
                            ff.write(eols)
                        ff.write(eols+eols)
                    ff.close()

gradm_l = ["xxDgradCISxx"]
grad_v = [] #For the Diabatic coupling gradient vector
for iii in range(ncoord):
    grad_v.append(0.0)
if (dot_mode==1 or dot_mode==2 or dot_mode==0):
    for dotiter in range(dotiters):
        for num1 in num1range:
            for num2 in num2range:
                grad_l = [] #Where I will store the strings representing the directional Hamiltonian gradient magnitudes
                for fn in file_l: #For the purposes of this hack, the file_l only contains one item, selected by the user.
                    for gr in gradm_l: #There is also only one item in this list, should be "xxDgradCISxx"
                        docmd("grep "+gr+" "+fn+" > tmp.txt") #Obtain diabatic Hamiltonian gradients
                        f = open("tmp.txt",rs)
                        f_s = f.read()
                        f.close()
                        big1_l = dp_re.findall(f_s) #Find numbers in this order: left to right, then top to bottom.
                                                    #For old files, such as c14eeF, HABQ has numkeep2 rows, NCoord cols.
                        if (lencheck(big1_l,numkeep2*ncoord,fn)==0):
                            grad_l.append(str(0))
                        else:
                            for i in range(len(big1_l)):
                                big1_l[i] = float(big1_l[i]) #Convert NCoord*numkeep2 gradient values into floats
                            temp = 0.0
                            for i in range(ncoord):
                                grad_v[i] = big1_l[(num1*numkeep +num2)*ncoord +i] #Take the values we are interested in
                                gr_f = grad_v[i]
                                if (dot_mode==1 or (dot_mode==0 and dotiter==0)):
                                    temp = temp +gr_f*direc_v[i]
                                elif (dot_mode==2 or (dot_mode==0 and dotiter==1)): #Make sure dot_mode==2! Should be.
                                    temp = temp +gr_f*gr_f
                            if (dot_mode==2 or (dot_mode==0 and dotiter==1)): temp = math.sqrt(temp)
                            grad_l.append(str(temp))
                        normalize(grad_v,ncoord) #Normalize the HABQ gradient vector
                docmd("rm tmp.txt")
                for minchoice in range(2):
                    if (minchoice == 0):
                        fname = "grD"+myfileno+".x"
                        geo_v = l1
                    else:
                        fname = "grA"+myfileno+".x"
                        geo_v = l0
                    ff = open(fname,ws)
                    for coeff in [0.0, 1.0]: #A cheesy way to write the base geometry, then the modified geometry for the .x
                        for q in range(natom):
                            ff.write(ele_l[q])
                            for mu in range(3):
                                ff.write(" " +str(geo_v[q*3 +mu] +coeff*grad_v[q*3 +mu]))
                            ff.write(eols)
                        ff.write(eols+eols)
                    ff.close()
