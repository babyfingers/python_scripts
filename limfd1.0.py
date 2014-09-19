# limfd.py
#!/usr/local/bin/python
import string
import os.path
import math
import time
import re

# limfd.py: Given two initial geometries and a template file, this
# script will run a finite-difference derivative coupling calculation
# on each of them (a la mfd.py). This file is a combination of that
# script and linterp.py.

# A note regarding the variable tia_i: it tells the script what value
# to give tia_order, which tells Q-Chem what kind of cross-geometry
# overlap you want. As of 2012-04-24, this is the code:
# tia_order = -3 make cross-geometry S matrix in stvman/mkSTV.C, crash.
#             -2 same as above
#              0 do nothing special, order states by energy
#              1 lite overlap ordering from previous calculation (do not use cross-geo AOs)
#              2 fast overlap (use AO overlap, but not full determinant of MO overlap matrix)
#              3 full overlap ordering from previous calculation 

#Number of points to run along reaction coordinate.
npts = 10
#Finite difference increment
hh = 0.0001

#Other variables
eols = "\n"
tabs = "\t"
spas = " "
uns = "_"
cls = ":"
ws = "w"
ps = "p"
ms = "m"
ss = "s"
Ss = "S"
rs = "r"
ls = "l"
txts  = '.txt'
ppfx = "p_"
mpfx = "m_"
gradpfx = 'grad_'
rds = "read"
fn_m = "fffnnn"
geo_m = "ggeeoo"
scf_m = "sssggg"
tia_m = "ttiiaa"
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
ele_re = re.compile(eols+"[a-zA-Z]{1,2}")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
infilenameend = '.i'
outfilenameend = '.o'
tmp_l = []
x = []
for i in range(npts):
    x.append(i/(npts-1.0))
wid_i = max(2,len(str(npts)))
maxj = 1 #Jobs depend on each other! Have to leave this at 1.
tmpfl = "tmp.txt"
tm = 5.0
stm = 0.3
tia_i = 2 #tia_order variable
tia_s = str(tia_i)
tia_sS = "-"+tia_s
tia_s0 = "0"
narray = []
if (tia_i>1):
    narray = range(5*npts)
else:
    narray = range(3*npts)

#jobwait: takes string and reg. ex. as
# input to stop the program until a
# certain specified job has completed.
def jobwait(waitmsg,qmk_r):
    print waitmsg+" Waiting "+str(stm)+" seconds."
    time.sleep(stm)
    os.system("qstat > "+tmpfl)
    f = open(tmpfl,rs)
    f_s = f.read()
    f.close()
    jobs_running = len(re.findall(qmk_r,f_s))
    os.system("rm "+tmpfl)
    while (jobs_running >= maxj):
        print waitmsg+" Waiting "+str(tm)+" seconds."
        time.sleep(tm)
        os.system("qstat > "+tmpfl)
        f = open(tmpfl,rs)
        f_s = f.read()
        f.close()
        jobs_running = len(re.findall(qmk_r,f_s))
        os.system("rm "+tmpfl)                                                                    

#docmd: prints a given string to screen
# and executes it as a command.
def docmd(s):
    print ">"+s
    os.system(s)

# Get some user input
dir = raw_input("Directory name:\n")
filenamebase = raw_input("Base input file name:\n")
mode1 = input("Just make files (1) or make and run (2)?\n")
f = open(filenamebase+infilenameend,rs)
template = f.read()
f.close()

#Check out geometry information
f = open(filenamebase +".x",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
atom_l = re.findall(atom_re,f_s)
#Make sure the two geometries match up in terms of # of atoms and order
tna = len(atom_l)
na = tna/2
#Check that the total number of atoms is even
if (tna%2 == 1):
    print "tna = "+str(tna)
    print "atom_l"
    print atom_l
    sys.exit("Do not have an even number of atoms (counting both geometries).")
#Check that the first half correspond with the second half
ele_l = []
ele_l2 = []
for i in range(na):
    tmp_l = re.findall(ele_re,atom_l[i])
    ele_l.append((tmp_l[0])[1:])
    tmp_l = re.findall(ele_re,atom_l[i+na])
    ele_l2.append((tmp_l[0])[1:])
    if (ele_l[i] != ele_l2[i]):
        print "atom_l"
        print atom_l
        print "Item "+str(i)+": "+ele_l[i]+" and "+ele_l2[i]
        print atom_l[i]
        print atom_l[i+na]
        sys.exit("Elements of two geometries do not correspond.")

infn_l = []
outfn_l = []
qoutro_0 = []
qoutro_p = []
qoutro_m = []
qoutro_Sp = []
qoutro_Sm = []
cart0_f = []
cart1_f = []
os.system("mkdir "+dir)
os.chdir(dir)
#Determine how much to move each atom for the finite-difference calculation
ssd = 0.0 #sum of squared displacements
for j in range(na):
    l0 = re.findall(num_re,atom_l[j])
    l1 = re.findall(num_re,atom_l[j+na])
    for c in range(3):
        cart0_f.append(float(l0[c]))
        cart1_f.append(float(l1[c]))
        ssd = ssd + (cart0_f[-1]-cart1_f[-1])*(cart0_f[-1]-cart1_f[-1])
ssd = math.sqrt(ssd)
dx = hh/ssd

#Create the input files
for i in range(npts):
    if (i==0):
        scf_s = "gwh"
        qoutro_0.append(" -save diabat -outdir "+ppfx+str(i).zfill(wid_i))
    else:
        scf_s = rds
        qoutro_0.append(" -indir " +ppfx+str(i-1).zfill(wid_i)+" -save diabat -outdir "+ppfx+str(i).zfill(wid_i))
    infn_l.append(ls+uns+str(i).zfill(wid_i)+infilenameend)
    outfn_l.append(ls+uns+str(i).zfill(wid_i)+outfilenameend)
    geo_s0 = ""
    geo_sp = ""
    geo_sm = ""
    for j in range(na):
        geo_s0 = geo_s0 +ele_l[j]
        geo_sp = geo_sp +ele_l[j]
        geo_sm = geo_sm +ele_l[j]
        for c in range(3):
            cartx_f0 = (1-x[i])*cart0_f[3*j +c] +x[i]*cart1_f[3*j +c]
            cartx_fp = (1-x[i]-dx)*cart0_f[3*j +c] +(x[i]+dx)*cart1_f[3*j +c]
            cartx_fm = (1-x[i]+dx)*cart0_f[3*j +c] +(x[i]-dx)*cart1_f[3*j +c]
            geo_s0 = geo_s0 +tabs +str(cartx_f0)
            geo_sp = geo_sp +tabs +str(cartx_fp)
            geo_sm = geo_sm +tabs +str(cartx_fm)
        if (j<na-1):
            geo_s0 = geo_s0 +eols
            geo_sp = geo_sp +eols
            geo_sm = geo_sm +eols
    f = open(infn_l[i],ws)
    f.write(template.replace(geo_m,geo_s0).replace(scf_m,scf_s).replace(tia_m,tia_s0).replace(fn_m,"x = "+str(i)+"/"+str(npts-1)+", central (reference) geometry."))
    f.close()
    f = open(ppfx+infn_l[i],ws)
    f.write(template.replace(geo_m,geo_sp).replace(scf_m,rds).replace(tia_m,tia_s).replace(fn_m,"x = "+str(i)+"/"+str(npts-1)+", +dx geometry."))
    f.close()
    f = open(mpfx+infn_l[i],ws)
    f.write(template.replace(geo_m,geo_sm).replace(scf_m,rds).replace(tia_m,tia_s).replace(fn_m,"x = "+str(i)+"/"+str(npts-1)+", -dx geometry."))
    f.close()
    if (tia_i>1):
        qoutro_Sp.append(" -indir "+ppfx+str(i).zfill(wid_i)+" -save diabat -outdir "+ppfx+Ss+ps+str(i).zfill(wid_i))
        qoutro_Sm.append(" -indir "+ppfx+str(i).zfill(wid_i)+" -save diabat -outdir "+ppfx+Ss+ms+str(i).zfill(wid_i))
        qoutro_p.append(" -indir "+ppfx+Ss+ps+str(i).zfill(wid_i))
        qoutro_m.append(" -indir "+ppfx+Ss+ms+str(i).zfill(wid_i))
        f = open(Ss+ppfx+infn_l[i],ws)
        f.write(template.replace(geo_m,geo_s0+eols+geo_sp).replace(scf_m,rds).replace(tia_m,tia_sS).replace(fn_m,"Cross-geometry AO overlap calculation for x = "+str(i)+"/"+str(npts-1)+", +dx geometry."))
        f.close()
        f = open(Ss+mpfx+infn_l[i],ws)
        f.write(template.replace(geo_m,geo_s0+eols+geo_sm).replace(scf_m,rds).replace(tia_m,tia_sS).replace(fn_m,"Cross-geometry AO overlap calculation for x = "+str(i)+"/"+str(npts-1)+", -dx geometry."))
        f.close()
    else:
        qoutro_p.append(" -indir "+ppfx+str(i).zfill(wid_i))
        qoutro_m.append(" -indir "+ppfx+str(i).zfill(wid_i))

#Run files
qintro = "mmqchem.csh -in "
qmezzo = " -out "
if (mode1==2):
    qmk_r = re.compile(ls+uns+"[0-9]+"+outfilenameend + ".* [ERQ] batch")
    docmd(qintro +infn_l[0] +qmezzo +outfn_l[0] +qoutro_0[0])
    if (tia_i>1):
        jobwait("",qmk_r)
        docmd(qintro +Ss+ppfx+infn_l[0] +qmezzo +Ss+ppfx+outfn_l[0] +qoutro_Sp[0])
        docmd(qintro +Ss+mpfx+infn_l[0] +qmezzo +Ss+mpfx+outfn_l[0] +qoutro_Sm[0])    
    jobwait("",qmk_r)
    docmd(qintro +ppfx+infn_l[0] +qmezzo +ppfx+outfn_l[0] +qoutro_p[0])
    docmd(qintro +mpfx+infn_l[0] +qmezzo +mpfx+outfn_l[0] +qoutro_m[0])
    for i in range(1,npts):
        docmd(qintro +infn_l[i] +qmezzo +outfn_l[i] +qoutro_0[i])
        if (tia_i>1):
            jobwait("",qmk_r)
            docmd(qintro +Ss+ppfx+infn_l[i] +qmezzo +Ss+ppfx+outfn_l[i] +qoutro_Sp[i])
            docmd(qintro +Ss+mpfx+infn_l[i] +qmezzo +Ss+mpfx+outfn_l[i] +qoutro_Sm[i])
        jobwait("",qmk_r)
        docmd(qintro +ppfx+infn_l[i] +qmezzo +ppfx+outfn_l[i] +qoutro_p[i])
        docmd(qintro +mpfx+infn_l[i] +qmezzo +mpfx+outfn_l[i] +qoutro_m[i])
    jobwait("",qmk_r)
    #The jobs are done running. Do the analysis.
    os.system("grep xxMOxOvr_stdxx *"+outfilenameend+" > MOxOvr_std.txt")
    os.system("grep xxMOxOvr_Newxx *"+outfilenameend+" > MOxOvr_New.txt")
    quant = ["MOxOvr_std","MOxOvr_New"]
    quantfn = ["MOxOvr_std.txt","MOxOvr_New.txt"]
    numq = len(quant)
    # Get names of states and save to array
    os.system("grep xxSta *"+outfilenameend+" > states")
    f = open("states",rs)
    states_raw = f.readline()
    f.close()
    os.system("rm states")
    int_re = re.compile("\s-?[0-9]+")
    states = re.findall(int_re,states_raw)
    numkeep = len(states)
    numkeep2 = numkeep*numkeep
#    print "numkeep = " +str(numkeep)
    for i in range(numkeep):
        s = states[i]
        states[i] = s[1:]
    print "Will generate and save gradients to .txt files."
    quantgfn = []
    for q in quant:
        quantgfn.append(gradpfx+q+txts)
    #Loop over quantities (in this case, all are measurements of state overlap)
    for q in range(numq):
        print "Now working on the quantity " + str(quant[q])
        fq = open(quantfn[q],rs)
        fq_s = fq.read()
        fq.close()
        fq_s = eols +fq_s
#        print "fq_s = " +fq_s
        gradq = []
        # Extract relevant values as floats and compute gradient
        #Loop over all the different hs (this calculation can do a finite-difference
        # approximation with many different values of h--the finite difference-- as
        # parameters. Each might give a slightly different result until h is
        # sufficiently small as to converge on the actual derivative.)
        thinc = 1/(1.889725989*2*dx) #dx is in A, but we want our gradients in Bohr
        print "dx = " +str(dx)
        for i in range(npts):
            #Get a list of all "plus" (x+h) and "minus" (x-h) values for overlap from the
            #output file we are currently interested in.
            p_re = re.compile(eols+ppfx+outfn_l[i]+".*")
            m_re = re.compile(eols+mpfx+outfn_l[i]+".*")
            p_raw = re.findall(p_re,fq_s)
            m_raw = re.findall(m_re,fq_s)
#            print "len(p_raw) = " +str(len(p_raw))
            #For each state of interest, print out the 'p' and 'm' overlap values;
            # currently, these values are still in string form and in groups of numkeep
            # (notice how we only have to loop over I, not just over I and J).
            for I in range(len(p_raw)):
                print "p_raw[" +states[I] +"] = " +p_raw[I]
                print "m_raw[" +states[I] +"] = " +m_raw[I]
            #For each state of interest, extract individual numbers for overlap values.
            for I in range(numkeep):
                print "Now doing state I(R) = " +states[I]
                p_ar = re.findall(dp_re,p_raw[I])
                m_ar = re.findall(dp_re,m_raw[I])
                l_p = len(p_ar)
                l_m = len(m_ar)
                #Catch any discrepancy between the number of values obtained and what we expect.
                if (l_p != l_m):
                    print "Problem with " + quant[q]
                    print "l_p = " + str(l_p) + ", l_m = " + str(l_m)
                    sys.exit("Error: output array lengths are not the same.")
                elif (l_p != numkeep):
                    print "Problem with " + quant[q]
                    print "l_p = l_m = " + str(l_p) + ", numkeep = " + str(numkeep)
                    sys.exit("Error: output array lengths are not equal to the number of states considered in diabatization.")
                #Loop over J (within loop over I) to calculate numkeep x numkeep matrix of our approximate
                # derivative couplings.
                for J in range(numkeep):
                    #Set strings to floating-point numberes
                    pfloat_J = float(p_ar[J])
                    mfloat_J = float(m_ar[J])
                    #Perform finite-difference calculation and save to gradq
                    gradq.append(thinc*(pfloat_J-mfloat_J))
                    print "<" +states[I] +"|d" +states[J] +"> = " +str(gradq[-1])
        # Write it all to a .txt file (note that we are still in the quantities loop, so
        #each quantity will get its own text file).
#        print "Total number of gradients = ",len(gradq)
#        print "q"+str(q)
        fq = open(quantgfn[q],ws)
        k = 0
        fq.write("states:")
        for I in states:
            for J in states:
                fq.write(spas +"<" +I +"|d" +J +">")
        fq.write(eols)
        for i in range(npts):
            fq.write(str(i) +cls)
            for I in states:
                for J in states:
#                    print "i"+str(i)+"I"+str(I)+"J"+str(J)
                    fq.write(spas +str(gradq[k]))
                    k += 1
            fq.write(eols)
        fq.close()
    
