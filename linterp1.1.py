# linterp1.1.py
#!/usr/local/bin/python
import string
import os.path
import math
import time
import sys
import re


# linterp1.1.py: Given two initial geometries and a template file, this
# script will run a number of files with geometries in between those
# two initial points, either at directed points by a .l file, or
# by automatically interpolating a number of evenly-spaced points.
# Derived from the script limfd1.7.py.

#Number of points to run along reaction coordinate (if auto-generating
# linear interpolation points).
npts = 26
savetype = "dc"
sequential = False #Run jobs in sequence (True) or in series (False)?
tia_i = 1 #tia_order variable

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
gradpfx = 'grad_'
rds = "read"
fn_m = "fffnnn"  #comment at top of input file
cha_m = "ccchhh" #charge
mul_m = "mmmuuu" #multiplicity
geo_m = "ggeeoo" #full geometry
scf_m = "sssggg" #argument for scf_guess REM variable
cis_m = "ccggdd" #argument for cis_guess_disk REM variable
tia_m = "ttiiaa" #argument for tia_order REM variable
sin_re = re.compile("^cis_singlets[\s=]+true\s?.*", re.I | re.M)
tri_re = re.compile("^cis_triplets[\s=]+true\s?.*", re.I | re.M)
unr_re = re.compile("^unrestricted[\s=]+true\s?.*", re.I | re.M)
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
ele_re = re.compile(eols+"[a-zA-Z]{1,2}")
mul_re = re.compile(mul_m+"\s*[1-9][0-9]*")
cha_re = re.compile(cha_m+"\s*-?[0-9]+")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
infilenameend = '.i'
outfilenameend = '.o'
tmp_l = []
maxj = 80 #Maximum number of jobs to run at once.
tmpfl = "tmp.txt"
tm = 4.0
stm = 0.5
tia_s = str(tia_i)
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
filenamebase = raw_input("Base input file name (max. 6 characters):\n")
while (len(filenamebase)>7):
    print "Sorry, that file name is no good. Please try a different one."
    filenamebase = raw_input("Base input file name (max. 6 characters):\n")
dir = "p_"+filenamebase+"_t"+tia_s+"_"
dir = dir +raw_input("Directory name addition ("+dir+"):\n")
qmk_r_g = re.compile(filenamebase+".* [ERQ] batch") #Generic to test total # of jobs
mode1 = input("Just make files (1), just run (2), or make and run (3)?\n")
mode_x = input("Auto-generate linear interpolation points (1), or read from file (2)?\n")
if (mode1==2 or mode1==3):
    mode_p = input("Start calculations fresh (1), or read SCF and CIS data from save files that already exist (2)?\n")
    if (mode_p == 2):
        sequential = False
else:
    mode_p = 1
f = open(filenamebase+infilenameend,rs)
template = f.read()
f.close()

#Generate/read linear interpolation points
x = []
xs = []
if (mode_x==1):
    for i in range(npts):
        x.append(i/(npts-1.0))
        xs.append(str(x[i]))
elif (mode_x==2):
    f = open(filenamebase+".l",rs)
    f_s = f.read()
    f.close()
    tmp_l = re.findall(num_re,f_s)
    for s in tmp_l:
        x.append(float(s))
        xs.append(s)
    npts = len(x)
wid_i = max(2,len(str(npts)))
njobs = npts

#Check out geometry information
f = open(filenamebase +".x",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
atom_l = re.findall(atom_re,f_s)
cha_l = re.findall(cha_re,f_s)
mul_l = re.findall(mul_re,f_s)
cha_l = re.findall(int_re,cha_l[0])
mul_l = re.findall(int_re,mul_l[0])
cha_i = int(cha_l[0])
mul_i = int(mul_l[0])
cha_s = str(cha_i)
mul_s = str(mul_i)

#Check out the kind of excited states
sin_i = len(re.findall(sin_re,template))
tri_i = len(re.findall(tri_re,template))
unr_i = len(re.findall(unr_re,template))

if ((sin_i>0 and tri_i>0) or unr_i>0):
    cgdt_i = 1
elif (sin_i==0 and tri_i>0):
    cgdt_i = 0
elif (sin_i>0 and tri_i==0):
    cgdt_i = 2
#elif (sin_i>1 or tri_i>1):
#    print "Error-- too many \"cis_singlets\" or \"cis_triplets\" commands."
#    sys.exit()
elif (sin_i==0 or tri_i==0):
    print "Error-- too few \"cis_singlets\" or \"cis_triplets\" commands."
    sys.exit()
cgdt_s = str(cgdt_i)

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
outdn_l = []
qoutro_0 = []
cart0_f = []
cart1_f = []
os.system("mkdir "+dir)
os.chdir(dir)

for j in range(na):
    l0 = re.findall(num_re,atom_l[j])
    l1 = re.findall(num_re,atom_l[j+na])
    for c in range(3):
        cart0_f.append(float(l0[c]))
        cart1_f.append(float(l1[c]))

#Create the input files and lists for file names and input commands
scf_st = "read\nscf_algorithm gdm"
cis_st = "true\ncis_guess_disk_type "+cgdt_s
for i in range(npts):
    outdn_l.append(ppfx+str(i).zfill(wid_i))
    if (mode_p==1 and (i==0 or not sequential)):
        scf_s = "gwh"
        cis_s = "false"
        qoutro_0.append(" -save "+savetype+" -outdir "+outdn_l[i])
    else:
        scf_s = scf_st
        cis_s = cis_st
        if (mode_p==1):
            qoutro_0.append(" -indir " +outdn_l[i-1]+" -save "+savetype+" -outdir "+outdn_l[i])
        else:
            qoutro_0.append(" -indir " +outdn_l[i]+" -save "+savetype+" -outdir "+outdn_l[i])
    infn_l.append(filenamebase+uns+str(i).zfill(wid_i)+infilenameend)
    outfn_l.append(filenamebase+uns+str(i).zfill(wid_i)+outfilenameend)
    geo_s0 = ""
    for j in range(na):
        geo_s0 = geo_s0 +ele_l[j]
        for c in range(3):
            cartx_f0 = (1-x[i])*cart0_f[3*j +c] +x[i]*cart1_f[3*j +c]
            geo_s0 = geo_s0 +tabs +str(cartx_f0)
        if (j<na-1):
            geo_s0 = geo_s0 +eols
    if (mode1==1 or mode1==3):
        f = open(infn_l[i],ws)
        f.write(template.replace(geo_m,geo_s0).replace(scf_m,scf_s).replace(cis_m,cis_s).replace(tia_m,tia_s0).replace(mul_m,mul_s).replace(cha_m,cha_s).replace(fn_m,"x = "+str(i)+"/"+str(npts-1)+", l = "+xs[i]+"."))
        f.close()

#Run files
if (mode1==2 or mode1==3):
    qintro = "dcqchem.csh -in "
    qmezzo = " -out "
    cmd_l   = []
    qmk_r_l = []
    djid_l  = []
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
    #We've constructed the necessary arrays. Now run all the jobs.
    limq(cmd_l,qmk_r_l,djid_l,njobs,tm)

