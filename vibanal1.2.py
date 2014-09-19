#!/usr/local/bin/python

# vibanal.py
import string
import os.path
import numpy
import math
import time
import sys
import re

#vibanal.py: given an output file for a frequency calculation, extracts normal mode
# and force constant information. Also takes gradient vector information from another
# output file, and combines that with the normal mode information to magical effect.
#
# Built from gradcompare.py

#v1.2: added a fully-quantum treatment to obtain characteristic lengths.

#v1.1: added a bit that does an excited state lifetime approx., as from
# 2012Fatehi2039, and a bit that does a Quantum analysis of the available
# vibrational potential energy.

quantum_treatment = True
#Set to true in order to assign potential energy to the vibrational
#modes as dictated by the quantum harmonic oscillator. Otherwise,
#just use classical equipartition theorem to assign 0.5kBT.

moreinfo = False
#Set to true in order to print out states that have large characteristic
#lengths or overlaps with the gradient vector, and print average energy
#of each mode assuming the quantum harmonic oscillator approximation.

qtype_l = ["mag1", "mag2", "mag1-mag2", "magerror", "direrror", "doterror", "RMSD"]
tabs = "\t"
eols = "\n"
tia_s = "1"
rs = "r"
ws = "w"
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
vib_re = re.compile(" [a-zA-Z]{1,2}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}\s+-?[0-9]\.[0-9]{10}")
force_re = re.compile("Force Cnst:\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+")
mass_re = re.compile("Red\. Mass:\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+")
freq_re = re.compile("Frequency:\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
ele_re = re.compile("[a-zA-Z]{1,2}")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
outfilenameend = '.o'
force_f = 0.064230506330 #Convert force constants from mdyn/A to Hartree/a0^2
mass_f = 1822.88839 #Convert reduced masses from a.m.u. to a.u.
freq_f = 0.000000725163277861 #Converts cm^{-1} to the a.u. frequency unit
tpi = 2*math.pi
#freq_f = 0.000002278167837214 #Convert wavenumbers to Hartree (i.e., quantum harmonic oscillator)
kBT = 0.000943709879 #kBT at room temperature in Hartree
kBTin = 1.0/kBT
Tfactor = 0.5 #Equipartion theorem: each degree of freedom in a Harmonic oscillator has 0.5*kBT energy.
mass_d = {"H":1.00783, "C":12.0, "O":15.99491}
ratiomax = 0.10

#magnitude: calculates the 2-norm of a
# numpy array v.
def magnitude(v):
    return math.sqrt(numpy.inner(v,v))

#normalize: normalizes a numpy array v.
def normalize(v):
    v = v/magnitude(v)
    return v

#doterror: finds the normalized overlap
# between two numpy vectors v1 and v2;
# takes the dot product of v1 and v2,
# then divides by the dot product of
# v2 with itself.
def doterror(v1, v2):
    return numpy.inner(v1,v2)/numpy.inner(v2,v2) -1

#direrror: finds the directional overlap
# between two numpy vectors v1 and v2.
def direrror(v1, v2):
    u1 = v1.copy()
    u2 = v2.copy()
    return 1 -numpy.inner(normalize(u1),normalize(u2))

#RMSD: returns the root mean squared
# deviation between two numpy vectors
# v and u.
def RMSD(v1, v2):
    d = v1-v2
    d = d**2
    return math.sqrt(numpy.average(d))

#Lanalyze: take a gradient, find its
# projection into some normal mode. Use
# force constant information to determine
# a characteristic length. Multiply these
# quantities together to compare to the
# magnitude of the gradient quantity (i.e.
# |V| vs dV/dx -- how does (dV/dX.u)*L
# stack up against |V|?
def Lanalyze(modesv_l,force_l,mass_l,freq_l,grad_v,nact_v,mag,ele_l):
    nmode = len(modesv_l)
    projgrad_l = [] #For each mode, the projection of the gradient onto that mode
    projdirg_l = [] #For each mode, the projection of the DIRECTION of the gradient onto that mode
    DHAB_l = []     #For each mode, the product |dV/dx|*L (\Delta H_{AB})
    ratio_l = []    #For each mode, the ratio |dV/dx|*L/V
    charL_l = []    #For each mode, the characteristic length, L
    mode_l = []     #For each mode, the index (starting with 1); indexed by frequency (low to high)
    lifetime = 0.0  #Total lifetime lower bound, as calculated by the equation from 2012Fatehi2039.
    projmag = 0.0   #Temporary variable that keeps track of the projected gradient magnitude
    gradmag = magnitude(grad_v)
    for mode_i in range(nmode):
        mode_l.append(mode_i+1)
        mode_v = modesv_l[mode_i]
        force = force_l[mode_i]
        projgrad = numpy.inner(mode_v,grad_v)
        projdirg = projgrad/gradmag
        if (quantum_treatment):
            charL = quantumHO(freq_l[mode_i], mass_l[mode_i], mode_i)
        else:
            charL = math.sqrt(2*Tfactor*kBT/force)
        lhs = math.fabs(projgrad*charL)
        projmag = projmag +projgrad*projgrad
        projgrad_l.append(projgrad)
        projdirg_l.append(projdirg)
        charL_l.append(charL)
        DHAB_l.append(lhs)
        ratio_l.append(lhs/math.fabs(mag))
    for q in range(natom):
        mymass_rt = math.sqrt(mass_d[ele_l[q]])
        for mu in range(3):
            mynact = nact_v[q*3 +mu]
            lifetime = lifetime +mynact*mymass_rt
    lifetime = 191.0/lifetime
    projmag = math.sqrt(projmag)
    print "Lifetime upper bound:     " +str(lifetime) +" fs"
    print "Magnitude of projections: " +str(projmag)
    print "Magnitude of gradient:    " +str(gradmag)
    print "Ratio of magnitudes:      " +str(projmag/gradmag)
    for mode_i in range(nmode):
        projbool = (math.fabs(projdirg_l[mode_i]) > 0.15)
        charbool = (math.fabs(charL_l[mode_i]) > 0.75)
        ratibool = (math.fabs(ratio_l[mode_i]) > ratiomax)
        if ((projbool or charbool or ratibool) and moreinfo):
            print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
            if (projbool):
                print "Large overlap for mode " +str(mode_i +1) +"."
            if (charbool):
                print "Large characteristic length for mode " +str(mode_i +1) +"."
            if (ratibool):
                print "Large change in V possible for mode " +str(mode_i +1) +"."
            print "ratio =    " +str(ratio_l[mode_i])
            print "charL =    " +str(charL_l[mode_i])
            print "projgrad = " +str(projgrad_l[mode_i])
            print "projdirg = " +str(projdirg_l[mode_i])
            print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    total_l = []
    DHAB_tot = 0.0
    for mode_i in range(nmode):
        DHAB_tot = DHAB_tot +DHAB_l[mode_i]*DHAB_l[mode_i]
        total_l.append([mode_l[mode_i], projgrad_l[mode_i], charL_l[mode_i], DHAB_l[mode_i], ratio_l[mode_i]])
        #total_l.append([ratio_l[mode_i], charL_l[mode_i], projdirg_l[mode_i], projgrad_l[mode_i], mode_l[mode_i]])
    DHAB_tot = math.sqrt(DHAB_tot)
    ratio_tot = DHAB_tot/math.fabs(mag)
    print "DHAB_tot = "+str(DHAB_tot)+", ratio_tot = "+str(ratio_tot)
    #for i in range(4):
    #    print "before: ",total_l[i]
    total_l = sorted(total_l, key=lambda quantity: -math.fabs(quantity[4])) #Sort by ratio absolute magnitude
    #for i in range(4):
    #    print "after : ",total_l[i]
    mymode = 0
    while (math.fabs(total_l[mymode][4])>ratiomax):
        mymode = mymode +1
    print "There are "+str(mymode)+" out of "+str(nmode)+" modes for which V might change by > "+str(100*ratiomax)+"%."
    mymode = mymode +1
    print "MODE\t\tHABQ\t\t\tDL\t\tDHAB\t\t\tRATIO"
    for mode_i in range(mymode):
        print total_l[mode_i]

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
        m1 = magnitude(v1)
        m2 = magnitude(v2)
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

#grepfls: takes a file name, a string to grep,
# and returns the grepped output.
def grepfls(filename_s, grepmarker_s):
    grepfln = "greptmp.txt"
    docmd("grep " +grepmarker_s +" " +filename_s +" > " +grepfln)
    grepf = open(grepfln,rs)
    grepf_s = grepf.read()
    grepf.close()
    docmd("rm " +grepfln)
    return grepf_s

#gethessdata: build an array of floats
# containing information gleaned from
# Q-Chem's Hessian. Only works for
# quantities for which there is only
# one per mode (e.g., frequency)
def gethessdata(nmode, file_s, conv_f, data_re, data_l):
    temp1_l = data_re.findall(file_s)
    for row_i in range(nmode/3):
        temp2_l = num_re.findall(temp1_l[row_i])
        for col_i in range(3):
            data_l.append(conv_f*float(temp2_l[col_i]))

#quantumHO:
def quantumHO(myfreq_f, mymass_f, mode_i):
    Z = 0.0
    diff = 10.0
    thresh = 0.00000001
    iii = 0
    Pi_l = []
    Ei_l = []
    DXi_l = []
    avgE = 0.0  #average fo the energy (halve for average potential energy)
    avgDX = 0.0 #average of the RMS deviation of X, \sum_i <i|X^2|i>*P_i
    Eterm = tpi*myfreq_f
    DXterm = 1.0/(tpi*mymass_f*myfreq_f)
    while (diff>thresh or iii<5):
        Ei = Eterm*(iii+0.5)
        Ei_l.append(Ei)
        DXi = math.sqrt((iii+0.5)*DXterm)
        DXi_l.append(DXi)
        Pi = math.exp(-Ei*kBTin)
        Pi_l.append(Pi)
        diff = Z
        Z = Z +Pi
        diff = Z -diff
        iii = iii +1
    Zin = 1.0/Z
    maxi = iii #max value +1, which is great, because range(maxi) = [0, ..., maxi-1]
    total_P = 0.0
    for iii in range(maxi):
        Pi = Pi_l[iii]*Zin
        Pi_l[iii] = Pi
        avgE  = avgE  +Pi*Ei_l[iii]
        avgDX = avgDX +Pi*DXi_l[iii]
        total_P = total_P +Pi
    if (moreinfo):
        if ((math.fabs(1-total_P)<0.0000001) and (diff/Z<0.000001)):
            print "---------------MODE "+str(mode_i +1)+"----------------"
            print "avgE =    " +str(avgE)
            print "avgDX =   " +str(avgDX)
        else:
            print "++++ERROR! ERROR! ERROR! Something is wrong!++++"
            print "diff =    " +str(diff)
            print "Z =       " +str(Z)
            print "diff/Z =  " +str(diff/Z)
            print "maxi =    " +str(maxi)
            print "total_P = " +str(total_P)
    return avgDX

#Get some user input:
freqfilename = "c14ee_hess.out02" #raw_input("Full file name for frequency calculation output.")
filenamebase = "c14eeF" #raw_input("Base input 1 file name (max. 6 characters):\n")

while (len(filenamebase)>7):
    print "Sorry, that file name is no good. Please try a different one."
    filenamebase = raw_input("Base input 1 file name (max. 6 characters):\n")
dir = "p_"+filenamebase+"_t"+tia_s+"_"
dir = dir +"full" #raw_input("Directory name addition ("+dir+"):\n")
grepm1 = "xxDgradCISxx" #raw_input("Greppable string for target gradient matrix:\n")
grepm2 = "xxDiaCoup_wetfxx" #raw_input("Greppable string for target derivative coupling matrix:\n")
transpode1 = 1 #input("Has the target gradient matrix been transposed? (0 no, 1 yes, 2 help):\n")
transpode2 = 0 #input("Has the target derivative coupling matrix been transposed? (0 no, 1 yes, 2 help):\n")
while (transpode1 != 0 and transpode1 != 1):
    print("Is the column-indexed value the fast variable? Or, maybe more simply,\n" \
          +"is the desired matrix taller than it is wide? If yes, you probably want" \
          +"to choose (0). Otherwise, choose (1).\n")
    transpode = input("Has the target matrix been transposed? (0) no, (1) yes:\n")

numkeep = 2  #input("Number of states (numkeep):\n")
matrixM = 138 #input("Number of columns (fast variable):\n")
matrixN = 4  #input("Number of rows (slow variable):\n")

#Check out geometry information
f = open(filenamebase +".x",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
atom_l = atom_re.findall(f_s)
tnatom = len(atom_l)
natom = tnatom/2
ncoord = 3*natom
nmode = ncoord -6
outfn_l = []
outdn_l = []

#Check out reaction coordinate information
f = open(filenamebase +".l",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
rc_l = num_re.findall(f_s)

#Step 2: super test time
# Extract data from Hessian output file
f = open(freqfilename,rs)
f_s = f.read()
f.close()

#Extract scalar data
force_l = []
mass_l = []
freq_l = []
gethessdata(nmode, f_s, force_f, force_re, force_l)
gethessdata(nmode, f_s, mass_f, mass_re, mass_l)
gethessdata(nmode, f_s, freq_f, freq_re, freq_l)

# Organize vector data into useful format
#modesv_l[normal mode #] = numpy array of normal mode
vib_l = vib_re.findall(f_s)
modesv_l = []
ele_l = []
for row_i in range(nmode/3):
    tmp_ll = [[], [], []]
    for atom_i in range(natom):
        if (row_i==0):
            ele_l.append((ele_re.findall(vib_l[row_i*natom +atom_i]))[0])
        num_l = num_re.findall(vib_l[row_i*natom +atom_i])
        for col_i in range(3):
            for mu in range(3):
                tmp_ll[col_i].append(float(num_l[col_i*3 +mu]))
    for col_i in range(3):
        modesv_l.append(numpy.array(tmp_ll[col_i]))

#Do mass-weighting to make normal modes orthonormal
for mode_i in range(nmode):
    for atom_i in range(natom):
        ele_f = math.sqrt(mass_d[ele_l[atom_i]])
        for mu in range(3):
            #print "before = " +str( modesv_l[mode_i][atom_i*3 +mu])
            modesv_l[mode_i][atom_i*3 +mu] = ele_f*modesv_l[mode_i][atom_i*3 +mu]
            #print "after  = " +str( modesv_l[mode_i][atom_i*3 +mu])
    modesv_l[mode_i] = normalize(modesv_l[mode_i])
    

#Construct normal mode overlap matrix
moverlap = []
for mode_i in range(nmode):
    modei_v = modesv_l[mode_i]
    moverlap.append([])
    for mode_j in range(nmode):
        modej_v = modesv_l[mode_j]
        moverlap[-1].append(numpy.inner(modei_v,modej_v))

#Test orthogonality
for mode_i in range(nmode-1):
    if (math.fabs(moverlap[mode_i][mode_i+1])>0.0000001):
        print "Large off-diagonal mode overlap; modes " +str(mode_i) +"|" +str(mode_i+1) +": " +str(moverlap[mode_i][mode_i+1])
    #print ['%5.3f' % val for val in moverlap[mode_i]]
#Test normality
for mode_i in range(nmode):
    if (math.fabs(moverlap[mode_i][mode_i]-1)>0.0000001):
        print "Large non-normality of modes; mode " +str(mode_i) +": " +str(moverlap[mode_i][mode_i])

# Extract gradient data

docmd("ls " +dir +"/*" +outfilenameend +" > tmp1.txt")
file_l_full = grepfls("tmp1.txt",".").split()
file_l = [file_l_full[0], file_l_full[-1]]
os.system("rm tmp1.txt")
#mynacts_ll[file #][state pair index] = numpy array of gradient quantity
mynacts_ll = []
#mygrads_ll[file #][state pair index] = numpy array of gradient quantity
mygrads_ll = []
#mag_l[file #] = floating point integer with magnitude of diabatic coupling
mag_l = []
for fn_i in range(len(file_l)):
    print "========================================================"
    print "Now working on file "+file_l[fn_i]+"."
    #Iterates over the list of files provided; each file represents a different geometry.
    fn = file_l[fn_i]
    #Step 1: collect data for this file
    #Step 1a: diabatic coupling magnitude
    # WARNING: only have it set now for 2-state systems #2state
    tmp_l = dp_re.findall(grepfls(fn,"xxHABxx00"))
    mag_l.append(float(tmp_l[1])) #2state
    mygrads_l = []
    mynacts_l = []
    for iii in range(matrixN):
        mygrads_l.append([])
        mynacts_l.append([])
    big1_l = dp_re.findall(grepfls(fn,grepm1))
    big2_l = dp_re.findall(grepfls(fn,grepm2))
    for i in range(len(big1_l)):
        big1_l[i] = float(big1_l[i])
    for i in range(len(big2_l)):
        big2_l[i] = float(big2_l[i])
    for jjj in range(matrixN):
        for iii in range(matrixM):
            if (transpode1==0):
                mygrads_l[jjj].append(big1_l[iii*matrixN +jjj])
            elif (transpode1==1):
                mygrads_l[jjj].append(big1_l[jjj*matrixM +iii])
            if (transpode2==0):
                mynacts_l[jjj].append(big2_l[iii*matrixN +jjj])
            elif (transpode2==1):
                mynacts_l[jjj].append(big2_l[jjj*matrixM +iii])
    for jjj in range(matrixN):
        mygrads_l[jjj] = numpy.array(mygrads_l[jjj])
        mynacts_l[jjj] = numpy.array(mynacts_l[jjj])
    mygrads_ll.append(mygrads_l)
    mynacts_ll.append(mynacts_l)
    Lanalyze(modesv_l,force_l,mass_l,freq_l,mygrads_ll[fn_i][1],mynacts_ll[fn_i][1],mag_l[fn_i],ele_l) #2state
    print "========================================================"

#Step 3: write output for this pair of files
#docmd("rm tmp.txt")
#for jjj in range(matrixN):
#    fname = "compout_" +str(jjj) +".txt"
#    ff = open(fname,ws)
#    ff.write(tabs)
#    for qtype_s in qtype_l:
#        ff.write(tabs +qtype_s)
#    for file_i in range(len(outdata_ll)):
#        ff.write(eols +str(rc_l[file_i]))
#        for q_i in range(len(qtype_l)):
#            ff.write(tabs +str(outdata_ll[file_i][jjj][q_i]))
#    ff.write(eols)
#    ff.close()



