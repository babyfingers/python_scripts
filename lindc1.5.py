#!/usr/local/bin/python

# lindc.py
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

# v1.1 added functionality that does the same thing for Hamiltonian
# gradients. Also extends the method to work for all values of numkeep,
# and for all state combinations.

# v1.2 added functionality that allows for the magnitude of the DCs and
# the Hamiltonian gradients to be calculated (i.e., dot these vectors
# with themselves and take the square root of the resulting value).

# v1.3 added functionality that might be able to condense some
# information about the direction of the derivative coupling vector
# into a single-valued function. Just normalize the vectors from
# successive runs and dot them with each other. Did this for the
# gradients too, because why not?

# v1.4 added functionality that does the same as v1.3, but across
# different quantities (i.e., where is the derivative coupling
# pointing w/ respect to the diabatic coupling gradient?).

# v1.5 enable single output processing for full magnitude calculations
# (i.e., obviously can't do directional magnitude calculations without
#  a well-defined initial and final geometry to interpolate between,
#  and can't do directional overlap matrices without multiple files.
#  Well, I guess you could do it between different quantities, but
#  I don't care about that right now).
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
    dot_mode = input("Compute directional magnitude (1), simple magnitude (2), directional overlap matrices (3), or the works (0)?:\n")
else:
    filenamebase = raw_input("File name base:\n")
    filename = filenamebase +"." +raw_input("Output file name suffix:\n")
    dot_mode = 2
    file_l = [filename]
    rc_l = ["0.0"]
    
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
atom_l = re.findall(atom_re,f_s)
if (oneormany==1):
    tnatoms = len(atom_l)
    natoms = tnatoms/2
else:
    natoms = len(atom_l)
    tnatoms = 2*natoms
ncoord = 3*natoms


outfn_l = []
outdn_l = []
direc_v = []

if (oneormany==1):
    #Check out reaction coordinate information
    f = open(filenamebase +".l",rs)
    f_s = f.read()
    f.close()
    f_s = eols +f_s
    rc_l = re.findall(num_re,f_s)

    #Create normalized directional vector
    if (dot_mode==1 or dot_mode==0):
        normin = 0.0
        for j in range(natoms):
            l0 = re.findall(num_re,atom_l[j])
            l1 = re.findall(num_re,atom_l[j+natoms])
            for c in range(3):
                direc_v.append(float(l1[c]) -float(l0[c]))
                normin = normin +direc_v[-1]*direc_v[-1]
        normin = math.sqrt(normin)
        print "distance between endpoints = "+str(normin)+" a.u."
        normin = 1.0/normin
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

dotiters = 1
if (dot_mode==0): dotiters = 2
    
if (dot_mode==1 or dot_mode==2 or dot_mode==0):
    for dotiter in range(dotiters):
        for num1 in range(numkeep):
            for num2 in range(numkeep):
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
                                if (dot_mode==1 or (dot_mode==0 and dotiter==0)):
                                    temp = temp +dc_f*direc_v[i]
                                elif (dot_mode==2 or (dot_mode==0 and dotiter==1)):
                                    temp = temp +dc_f*dc_f
                            if (dot_mode==2 or (dot_mode==0 and dotiter==1)): temp = math.sqrt(temp)
                            nact_l.append(str(temp))
                docmd("rm tmp.txt")
                if (dot_mode==1 or (dot_mode==0 and dotiter==0)):
                    fname = "DCs_"+str(num1)+"-"+str(num2)+".txt"
                elif (dot_mode==2 or (dot_mode==0 and dotiter==1)):
                    fname = "DCM_"+str(num1)+"-"+str(num2)+".txt"
                ff = open(fname,ws)
                ff.write(tabs)
                for dc in nactm_l:
                    ff.write(tabs+dc)
                for fn_i in range(len(file_l)):
                    ff.write(eols+rc_l[fn_i]+" ")
                    for dc_i in range(len(nactm_l)):
                        ff.write(tabs+nact_l[fn_i*numkeep2 +dc_i])
                ff.write(eols)
                ff.close()


gradm_l = ["xxEgradxx","xxEgradCISxx","xxDgradxx","xxDgradCISxx"]
if (dot_mode==1 or dot_mode==2 or dot_mode==0):
    for dotiter in range(dotiters):
        for num1 in range(numkeep):
            for num2 in range(numkeep):
                grad_l = [] #Where I will store the strings representing the directional Hamiltonian gradient magnitudes
                for fn in file_l:
                    for gr in gradm_l:
                        docmd("grep "+gr+" "+fn+" > tmp.txt")
                        f = open("tmp.txt",rs)
                        f_s = f.read()
                        f.close()
                        big1_l = dp_re.findall(f_s)
                        if (lencheck(big1_l,numkeep2*ncoord,fn)==0):
                            grad_l.append(str(0))
                        else:
                            for i in range(len(big1_l)):
                                big1_l[i] = float(big1_l[i])
                            temp = 0.0
                            for i in range(ncoord):
                                gr_f = big1_l[(num1*numkeep +num2)*ncoord +i]
                                if (dot_mode==1 or (dot_mode==0 and dotiter==0)):
                                    temp = temp +gr_f*direc_v[i]
                                elif (dot_mode==2 or (dot_mode==0 and dotiter==1)):
                                    temp = temp +gr_f*gr_f
                            if (dot_mode==2 or (dot_mode==0 and dotiter==1)): temp = math.sqrt(temp)
                            grad_l.append(str(temp))
                docmd("rm tmp.txt")
                if (dot_mode==1 or (dot_mode==0 and dotiter==0)):
                    fname = "grads_"+str(num1)+"-"+str(num2)+".txt"
                elif (dot_mode==2 or (dot_mode==0 and dotiter==1)):
                    fname = "gradM_"+str(num1)+"-"+str(num2)+".txt"
                ff = open(fname,ws)
                ff.write(tabs)
                for grad in gradm_l:
                    ff.write(tabs+grad)
                for fn_i in range(len(file_l)):
                    ff.write(eols+rc_l[fn_i]+" ")
                    for gr_i in range(len(gradm_l)):
                        ff.write(tabs+grad_l[fn_i*numkeep2 +gr_i])
                ff.write(eols)
                ff.close()



# Make a correlation function of the direction part of various quantities -- all possible combinations, if you want.
if (dot_mode==3 or dot_mode==0):
    for qtype in range(4):
        donacts1 = bool(qtype%2) #Quantity 1: do gradients (False) or NACTs (True)
        donacts2 = bool(qtype/2) #Quantity 2: do gradients (False) or NACTs (True)
        if (donacts1):
            tempm1_l = nactm_l
        else:
            tempm1_l = gradm_l
        if (donacts2):
            tempm2_l = nactm_l
        else:
            tempm2_l = gradm_l
        for grp1_s in tempm1_l: #Loop over types of derivative couplings/gradients (grp1_s == grep string)
            for grp2_s in tempm2_l:
                for num1 in range(numkeep2):
                    if (num1%numkeep != num1/numkeep and "Egrad" in grp1_s): continue #Adiabatic Hamiltonian off-diagonals are always 0
                    if (num1%numkeep == num1/numkeep and "DerCoup" in grp1_s): continue #Adiabatic DC off-diagonals are always 0
                    for num2 in range(numkeep2):
                        if (num2%numkeep !=num2/numkeep and "Egrad" in grp2_s): continue
                        if (num2%numkeep == num2/numkeep and "DerCoup" in grp2_s): continue
                        TEMP_do = [] #Where I will store the final matrix of directional overlaps for the DCs/gradients
                        for fn1_i in range(len(file_l)): #Loop over output files (geometry)
                            TEMP_do.append([])
                            fn1 = file_l[fn1_i]
                            docmd("grep "+grp1_s+" "+fn1+" > tmp.txt")
                            f = open("tmp.txt",rs)
                            f1_s = f.read()
                            f.close()
                            big1_l = dp_re.findall(f1_s)
                            for i in range(len(big1_l)):
                                big1_l[i] = float(big1_l[i])
                            for fn2_i in range(len(file_l)):
                                fn2 = file_l[fn2_i]
                                docmd("grep "+grp2_s+" "+fn2+" > tmp.txt")
                                f = open("tmp.txt",rs)
                                f2_s = f.read()
                                f.close()
                                big2_l = dp_re.findall(f2_s)
                                for i in range(len(big2_l)):
                                    big2_l[i] = float(big2_l[i])
                                if (lencheck(big1_l,numkeep2*ncoord,fn1)==0 or lencheck(big2_l,numkeep2*ncoord,fn2)==0):
                                    TEMP_do[fn1_i].append(str(0))
                                else:
                                    for i in range(len(big2_l)):
                                        big2_l[i] = float(big2_l[i])
                                    temp = 0.0
                                    temp1 = 0.0
                                    temp2 = 0.0
                                    tm1_f = 0.0
                                    tm2_f = 0.0 #keep indenting after this line
                                    for i in range(ncoord):
                                        if (donacts1):
                                            tm1_f = big1_l[numkeep2*i +num1]
                                        else:
                                            tm1_f = big1_l[num1*ncoord +i]
                                        if (donacts2):
                                            tm2_f = big2_l[numkeep2*i +num2]
                                        else:
                                            tm2_f = big2_l[num2*ncoord +i]
                                        temp = temp +tm1_f*tm2_f
                                        temp1 = temp1 +tm1_f*tm1_f
                                        temp2 = temp2 +tm2_f*tm2_f
                                    if (temp1 > 0.0 and temp2 > 0.0): temp = temp/math.sqrt(temp1*temp2)
                                    TEMP_do[fn1_i].append(str(temp))
                        docmd("rm tmp.txt")
                        fname = "DO_"+grp1_s[2:-2]+"-"+grp2_s[2:-2]+"_"+str(num1%numkeep)+"-"+str(num1/numkeep)+"_"+str(num2%numkeep)+"-"+str(num2/numkeep)+".txt"
                        if (len(fname) > 0):
                            ff = open(fname,ws)
                            ff.write(tabs)
                            for rc_s in rc_l:
                                ff.write(tabs+" "+rc_s)
                            ff.write(eols)
                            for fn1_i in range(len(file_l)):
                                ff.write(rc_l[fn1_i])
                                for fn2_i in range(len(file_l)):
                                    ff.write(tabs+" "+TEMP_do[fn1_i][fn2_i])
                                ff.write(eols)
                            ff.close()
                    
