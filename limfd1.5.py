# limfd.py
#!/usr/local/bin/python
import string
import os.path
import math
import time
import re

# Version 1.5: Minor change: added a control variable that allows you
#              to run jobs in such a way that the 'base' jobs (i.e.,
#              neither +h nor -h) do not depend on each other, allowing
#              them all to be run simultaneously if desired. This boolean
#              variable is called sequential. Also, had the script label
#              relevant input files with the value of hh, for later
#              reference (previously, if you looked at a set of input files,
#              it was impossible to tell what hh was).

# Version 1.4: IMPORTANT CHANGE! Corrected the code so that it now
#              calculates accurate finite-difference directional
#              derivative approximations. Explanation: we have dx = hh/ssd,
#              where hh is the finite difference step, ssd is the sum
#              of the squares of the difference x1-x0, where x0 is the
#              initial geometry in the interpolation, and x1 is the final
#              geometry. When computing the finite-difference, we want
#              NACT = ((f(x)+dx*(x1-x0))-(f(x)-dx*(x1-x0)))/(2*hh).
#              Previously, we had a denominator of 2*dx, which is incorrect.

# Version 1.3: Added functionality to allow for external generation of
#              the points along the linear interpolation route. In other
#              words, now have the ability to read different values of x
#              from a file rather than generating them within the code.
#              Also, added code to incorporate the CIS_GUESS_DISK REM
#              variable functionality (CIS amplitudes are read from job
#              to job.)
# Version 1.2: Altered the queuing system to improve job parallelization.
#              Also added extraction of different types of overlap and
#              calculation of respective gradients.
# Version 1.1: Minor change to the way files are named.

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

#Number of points to run along reaction coordinate (if auto-generating
# linear interpolation points).
npts = 10
#Finite difference increment
hh = 0.00001

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
cis_m = "ccggdd"
tia_m = "ttiiaa"
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
ele_re = re.compile(eols+"[a-zA-Z]{1,2}")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
infilenameend = '.i'
outfilenameend = '.o'
tmp_l = []
maxj = 80 #Maximum number of jobs to run at once.
tmpfl = "tmp.txt"
tm = 10.0
stm = 0.3
sequential = False #Run jobs in sequence (True) or in series (False)?
tia_i = 1 #tia_order variable
tia_s = str(tia_i)
tia_sS = "-"+tia_s
tia_s0 = "0"

#limq: takes a bunch of arrays as arguments,
# each with information about the jobs that
# need to be run and their dependencies. Allows
# you to run a bunch of jobs in parallel with
# an arbitrary set of dependencies among them.
#   cmd:   (string) submission command
#   qmk_r: (regular expression) to check for once
#          job is in submission queue.
#   djid:  (integer) non-unique job dependency
#          ID; i.e., the job ID of the job which
#          which needs to be completed before
#          this job can be submitted.
#   n:     (integer) total number of jobs
#   maxjobs: (integer) maximum number of jobs we
#            want simultaneously running.
#   tm:    (float) number of seconds to wait when
#          taking a break.
# Note: the ID of a job is equivalent to its index
# in any one of these arrays.
def limq(cmd,qmk_r,djid,n,tm):
    nj = 0 #The number of jobs we have already submitted
    queue_l = [] #List of indices of jobs (JIDs) waiting
                 #to be submitted.
    queue_l.extend(newq(djid,-1)) #Add jobs that can be run right away to queue
                                  #(no dependencies).
    while (nj<n):
        taketm(tm)
        for qj in queue_l:
            if (qcheck(qmk_r[qj])==0 and qcheck(qmk_r_g)<maxj):
                docmd(cmd[qj])
                nj = nj +1
                queue_l.remove(qj)
                queue_l.extend(newq(djid,qj))
                break
    while (qcheck(qmk_r_g)>0):
        taketm(tm)

#taketm: instructs script to sleep for tm
# seconds, and announces it to standard
# output.
def taketm(tm):
    print "Waiting "+str(tm)+" seconds."
    time.sleep(tm)

#newq: given a list of integers and a
# potentially matching integer, return
# a tuple containing the indices of the list
# for which the contents of the list match
# the integer. Ex.: int_list = [1,2,0,1],
# int_match = 1. newq returns [0,3] because
# int_list[0]=int_list[3]=int_match.
def newq(int_list,int_match):
    match = []
    for i in range(len(int_list)):
        if (int_list[i]==int_match):
            match.append(i)
    return match

#qcheck: taking a regular expression as its
# argument, pings qstat and checks for it.
# Returns the number of results for a search
# of the qstat output for the regular
# expression.
def qcheck(qmk_r):
    os.system("qstat > "+tmpfl)
    f = open(tmpfl,rs)
    f_s = f.read()
    f.close()
    return len(re.findall(qmk_r,f_s))

#jobwait: takes string and reg. ex. as
# input to stop the program until a
# certain specified job has completed.
def jobwait(waitmsg,qmk_r):
    print waitmsg+" Waiting "+str(stm)+" seconds."
    time.sleep(stm)
    jobs_running = qcheck(qmk_r)
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
filenamebase = raw_input("Base input file name (max. 6 characters):\n")
while (len(filenamebase)>7):
    print "Sorry, that file name is no good. Please try a different one."
    filenamebase = raw_input("Base input file name (max. 6 characters):\n")
qmk_r_g = re.compile(filenamebase+".* [ERQ] batch") #Generic to test total # of jobs
mode1 = input("Just make files (1), make and run (2), or just do analysis (3)?\n")
mode_x = input("Auto-generate linear interpolation points (1), or read from file (2)?\n")
f = open(filenamebase+infilenameend,rs)
template = f.read()
f.close()

#Generate/read linear interpolation points
x = []
if (mode_x==1):
    for i in range(npts):
        x.append(i/(npts-1.0))
elif (mode_x==2):
    f = open(filenamebase+".l",rs)
    f_s = f.read()
    f.close()
    tmp_l = re.findall(num_re,f_s)
    for s in tmp_l:
        x.append(float(s))
    npts = len(x)
print x
wid_i = max(2,len(str(npts)))
if (tia_i>1):
    njobs = 5*npts
else:
    njobs = 3*npts

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
    if (i==0 or not sequential):
        scf_s = "gwh"
        cis_s = "false"
        qoutro_0.append(" -save diabat -outdir "+ppfx+str(i).zfill(wid_i))
    else:
        scf_s = rds
        cis_s = "true"
        qoutro_0.append(" -indir " +ppfx+str(i-1).zfill(wid_i)+" -save diabat -outdir "+ppfx+str(i).zfill(wid_i))
    infn_l.append(filenamebase+uns+str(i).zfill(wid_i)+infilenameend)
    outfn_l.append(filenamebase+uns+str(i).zfill(wid_i)+outfilenameend)
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
    f.write(template.replace(geo_m,geo_s0).replace(scf_m,scf_s).replace(cis_m,cis_s).replace(tia_m,tia_s0).replace(fn_m,"x = "+str(i)+"/"+str(npts-1)+", central (reference) geometry."))
    f.close()
    f = open(ppfx+infn_l[i],ws)
    f.write(template.replace(geo_m,geo_sp).replace(scf_m,rds).replace(cis_m,"true").replace(tia_m,tia_s).replace(fn_m,"x = "+str(i)+"/"+str(npts-1)+", +dx geometry. h = "+str(hh)))
    f.close()
    f = open(mpfx+infn_l[i],ws)
    f.write(template.replace(geo_m,geo_sm).replace(scf_m,rds).replace(cis_m,"true").replace(tia_m,tia_s).replace(fn_m,"x = "+str(i)+"/"+str(npts-1)+", -dx geometry. h = "+str(hh)))
    f.close()
    if (tia_i>1):
        qoutro_Sp.append(" -indir "+ppfx+str(i).zfill(wid_i)+" -save diabat -outdir "+ppfx+Ss+ps+str(i).zfill(wid_i))
        qoutro_Sm.append(" -indir "+ppfx+str(i).zfill(wid_i)+" -save diabat -outdir "+ppfx+Ss+ms+str(i).zfill(wid_i))
        qoutro_p.append(" -indir "+ppfx+Ss+ps+str(i).zfill(wid_i))
        qoutro_m.append(" -indir "+ppfx+Ss+ms+str(i).zfill(wid_i))
        f = open(Ss+ppfx+infn_l[i],ws)
        f.write(template.replace(geo_m,geo_s0+eols+geo_sp).replace(scf_m,rds).replace(cis_m,"true").replace(tia_m,tia_sS).replace(fn_m,"Cross-geometry AO overlap calculation for x = "+str(i)+"/"+str(npts-1)+", +dx geometry. h = "+str(hh)))
        f.close()
        f = open(Ss+mpfx+infn_l[i],ws)
        f.write(template.replace(geo_m,geo_s0+eols+geo_sm).replace(scf_m,rds).replace(cis_m,"true").replace(tia_m,tia_sS).replace(fn_m,"Cross-geometry AO overlap calculation for x = "+str(i)+"/"+str(npts-1)+", -dx geometry. h = "+str(hh)))
        f.close()
    else:
        qoutro_p.append(" -indir "+ppfx+str(i).zfill(wid_i))
        qoutro_m.append(" -indir "+ppfx+str(i).zfill(wid_i))

#Run files
if (mode1==2):
    qintro = "mmqchem.csh -in "
    qmezzo = " -out "
    cmd_l   = []
    qmk_r_l = []
    djid_l  = []
    #Do 'base' jobs first; not plus and minus jobs for finite difference.
    qmk_r_l.append(re.compile("This is dummy text that will never match anything that comes up in qstat."))
    cmd_l.append(qintro +infn_l[0] +qmezzo +outfn_l[0] +qoutro_0[0]) #First job to submit
    djid_l.append(-1)
    for i in range(1,npts):
        cmd_l.append(qintro +infn_l[i] +qmezzo +outfn_l[i] +qoutro_0[i])
        if (sequential):
            djid_l.append(i-1)
            qmk_r_l.append(re.compile("\s"+outfn_l[i-1]+".* [ERQ] batch"))
        else:
            djid_l.append(-1)
            qmk_r_l.append(re.compile("This is dummy text that will never match anything that comes up in qstat."))
    #Now do the offshoot jobs
    for i in range(0,npts):
        if (tia_i>1):
            #First, add arrays for the overlap ("S") calculations
            cmd_l.append(qintro +Ss+ppfx+infn_l[i] +qmezzo +Ss+ppfx+outfn_l[i] +qoutro_Sp[i])
            qmk_r_l.append(re.compile("\s"+outfn_l[i]+".* [ERQ] batch"))
            djid_l.append(i)
            cmd_l.append(qintro +Ss+mpfx+infn_l[i] +qmezzo +Ss+mpfx+outfn_l[i] +qoutro_Sm[i])
            qmk_r_l.append(re.compile("\s"+outfn_l[i]+".* [ERQ] batch"))
            djid_l.append(i)
            #Then, add certain arrays for the p+m calculations
            qmk_r_l.append(re.compile("\s"+Ss+ppfx+outfn_l[i]+".* [ERQ] batch"))
            djid_l.append(npts -1 +4*i)
            qmk_r_l.append(re.compile("\s"+Ss+mpfx+outfn_l[i]+".* [ERQ] batch"))
            djid_l.append(npts +4*i)
        else:
            #If the p+m calculations don't depend on the overlap calculations, things are simpler.
            qmk_r_l.append(re.compile("\s"+outfn_l[i]+".* [ERQ] batch"))
            djid_l.append(i)
            qmk_r_l.append(re.compile("\s"+outfn_l[i]+".* [ERQ] batch"))
            djid_l.append(i)
        cmd_l.append(qintro +ppfx+infn_l[i] +qmezzo +ppfx+outfn_l[i] +qoutro_p[i])
        cmd_l.append(qintro +mpfx+infn_l[i] +qmezzo +mpfx+outfn_l[i] +qoutro_m[i])
    #We've constructed the necessary arrays. Now run all the jobs.
    limq(cmd_l,qmk_r_l,djid_l,njobs,tm)

if (mode1==2 or mode1==3):
    #Do the analysis.
    quant = ["MOxOvr_std","MOxOvr_New","aMOxOvr_stdaa","aMOxOvr_Newaa","aMOxOvr_stdab","aMOxOvr_Newab"]
    numq = len(quant)
    # Get names of states and save to array
    os.system("grep xxSta *"+outfilenameend+" > states")
    f = open("states",rs)
    states_raw = f.readline()
    f.close()
    os.system("rm states")
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
    print "hh = " +str(hh)
    print "dx = " +str(dx)
    for q in range(numq):
        os.system("grep xx"+quant[q]+"xx *"+outfilenameend+" > "+quant[q]+txts)
        print "Now working on the quantity " +quant[q]
        fq = open(quant[q]+txts,rs)
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
        thinc = 1/(1.889725989*2*hh) #hh is in A, but we want our gradients in Bohr
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
#            for I in range(len(p_raw)):
#                print "p_raw[" +states[I] +"] = " +p_raw[I]
#                print "m_raw[" +states[I] +"] = " +m_raw[I]
            #For each state of interest, extract individual numbers for overlap values.
            for I in range(numkeep):
#                print "Now doing state I(R) = " +states[I]
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
#                    print "<" +states[I] +"|d" +states[J] +"> = " +str(gradq[-1])
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
    
