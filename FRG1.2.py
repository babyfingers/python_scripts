# Snqchem.py
#!/usr/local/bin/python
import os.path
import time
import re

# v1.2: modified to account for new use of tia_order.

# v1.1: add markers, replacement variables for eps, temp. This way, the
# template only need be changed for different systems, and more control
# is added to the script. In addition, some variable names are changed
# to be more informative.

# Snqchem.csh: updated for use of cross-geometry overlap matrices and
# calculations, in particular to facilitate the use of the newest
# version of the ER-eps code, which uses both cross-geometry overlaps
# and Jacobi sweeps for function optimization.

#This program takes a single input file and increments the values of a single
# parameter (from 'pta' to 'ptb' in increments of 'inc'; the parameter is marked
# by the initial string in its line given by 'poi') and generates a ton of input
# files and optionally runs them and extracts certain info.

# v1.1: adds new functionality. Give opportunity to do a bunch of calculations at once
#       instead of just starting where the last one left off. Furthermore, can now
#       accomodate calls to myqchem.csh or solqchem.csh (must be changed manually!)

# Variables: you probably will only need to change things on the following few lines.
maxjobs = 16 #If we just do all the jobs at once, how many should we run at a time?
pta  = pppta  #Start the position value where?
ptb  = ppptb  #End the position value where?
inc  = iinncc # Position value increment
veps = seeppss #Dielectric coefficient
vTmp = sTTmmpp #Temperature
ttia = str(ttovlp) #The type of overlap calculation to perform
meps = 'eeppss'
seps = str(veps)
mTmp = 'TTmmpp'
sTmp = str(vTmp)
mzz1  = 'zz1'  #Some marker for the position variable
mzz2  = 'zz2'  #Some marker for the position variable +1
msgf = 'sssggg' #some marker for the scf_guess variable
mtia = 'tttooo' #some marker for the tia_order variable
tmpfl = 'tempout.txt'
os.system("pwd > "+tmpfl)
f = open(tmpfl,'r')
fileloc = f.read()
fileloc = fileloc[0:-1]
f.close()
os.system("rm "+tmpfl)
infilenameend = '.i'
outfilenameend = '.o'
sz  = '0' #Just a string of 0
so  = '1'
su  = '_'
ss  = ' '
sS  = 'S'
nfiles = int(round(1+((ptb-pta)/inc)))
narray = range(nfiles)
vals = []
szz1 = []
ssg = []
stia = []
rd = 'read'
tm = 2.0 #Seconds between issuing qstat
def chars_before_period(str):
    for i in range(len(str)):
        if str[i]=='.':
            return i
    return len(str)

def chars_after_period(str):
    for i in range(len(str)):
        if str[i]=='.':
            return len(str)-i-1
    return 0

def jobwait(waitmsg,qmk_r,Mscfg):
    jobs_running = maxjobs +1
    while ((jobs_running>0 and Mscfg!=3) or (jobs_running>maxjobs and Mscfg==3)):
        print waitmsg+" Waiting "+str(tm)+" seconds."
        time.sleep(tm)
        os.system("qstat > "+tmpfl)
        f = open(tmpfl,'r')
        f_s = f.read()
        f.close()
        jobs_running = len(re.findall(qmk_r,f_s))
        os.system("rm "+tmpfl)

os.chdir(fileloc)
dir = raw_input("Directory name:\n")
filenamebase = raw_input("Base input file name:\n")
Mmkrun = input("Just make files (1) or make and run (2)?\n")
Mstab = input("Run normal calculation (1) or stability analysis (2)?\n")
if (Mstab==1):
    Mscfg = input("Input SCF guess from previous calculation (1), from individual set of SCF guess folders (2), or start fresh each time (3)?\n")
else:
    Mscfg = 1
if (Mstab==1 and Mscfg==1):
    # Although in principal we can do this for Mstab==2, it would require a little reworking of the code.
    Mstart = input("Start fresh with scf_guess = gwh (1) or start where previous calculation left off (2)?\n")
else:
    Mstart = 1
if (Mscfg==1 and Mstart==1):
    ssg.append('gwh')
    stia.append(sz)
elif (Mscfg != 3):
    ssg.append(rd)
    stia.append(ttia)

if (Mscfg!=3):
    for i in range(nfiles-1):
        ssg.append(rd)
        stia.append(ttia)
else:
    for i in range(nfiles):
        ssg.append('gwh')
        stia.append(ttia)

#Since we are labelling files based on the zzz variable,
# let's make sure they each have the same number of zeros
# before and after the decimal so that it is neat and
# all the files order themselves well.
max_cbp_s = 1
max_cap_s = 2
for i in narray:
    vals.append(pta + inc*i)
    szz1.append(str(vals[i]))
    szz1.append(str(vals[i]))
    cbp = chars_before_period(szz1[-1])
    cap = chars_after_period(szz1[-1])
    if cbp>max_cbp_s:
        max_cbp_s = cbp
    if cap>max_cap_s:
        max_cap_s = cap
szz1.pop(-1)
for i in range(2*nfiles-1):
    qmk_r = re.compile("[a-z]")
    if (len(re.findall(qmk_r,szz1[i])) == 0):
        while chars_before_period(szz1[i])<max_cbp_s:
            szz1[i] = sz + szz1[i]
        while chars_after_period(szz1[i])<max_cap_s:
            szz1[i] = szz1[i] + sz

#Now open up the template file and replace the relevant markers
# with their respective variables to create each input file.
template = open(filenamebase+infilenameend,'r')
Stemplate = open(sS+filenamebase+infilenameend,'r')
template_s = template.read()
Stemplate_s = Stemplate.read()
os.system("mkdir "+dir)
os.chdir(dir)
infile = []
outfile = []
for i in narray:
    infile.append(filenamebase+szz1[2*i]+infilenameend)
    outfile.append(szz1[2*i]+filenamebase+outfilenameend)
    print "Writing "+infile[2*i]
    f_i = open(infile[2*i],'w')
    f_i.write(template_s.replace(mzz1,szz1[2*i]).replace(msgf,ssg[i]).replace(mtia,stia[i]).replace(meps,seps).replace(mTmp,sTmp))
    f_i.close()
    if (i < nfiles-1):
        infile.append(sS+filenamebase+szz1[2*i+1]+infilenameend)
        outfile.append(sS+szz1[2*i+1]+filenamebase+outfilenameend)
        print "Writing "+infile[2*i+1]
        f_i = open(infile[2*i+1],'w')
        f_i.write(Stemplate_s.replace(mzz1,szz1[2*i+1]).replace(mzz2,szz1[2*i+2]))
        f_i.close()                
template.close()
Stemplate.close()

if Mmkrun==2:
    #Then we want to run these files, too. Since each job will depend on the
    # last for its scf_guess input, we want to make sure that the last one
    # has finished before we start the next.
    qmk_r = re.compile(filenamebase+outfilenameend+".* [ERQ] batch") #What to look for in qstat to indicate a running job
    #If we're not picking up where an old calculation left off, we have to do job 0 differently
    if (Mstart==1 and Mscfg!=3):
        if (Mstab==1 and Mscfg==1):
            sv = " -outdir scf_read -save diabat"
        elif (Mstab==2):
            sv = " -outdir scf_read"+szz1[0]+" -save diabat"
        elif (Mscfg==2):
            sv = " -indir scf_read"+szz1[0]
        command = "snqchem.csh -in "+infile[0]+" -out "+outfile[0]+sv
        print "Executing command: "+command
        os.system(command)
        read_job_start = 1
    else:
        read_job_start = 0
    #Proceed with the remaining jobs
    for i in range(read_job_start,2*nfiles-1):
        jobs_running = maxjobs+1 #When this drops below 1, the job has completed
        if (Mstab==1 and Mscfg==1):
            sv = " -indir scf_read -outdir scf_read -save diabat"
        elif (Mstab==2):
            sv = " -indir scf_read"+szz1[i-1]+" -outdir scf_read"+szz1[i]+" -save diabat"
        elif (Mscfg==2):
            sv = " -indir scf_read"+szz1[i]
        elif (Mscfg==3):
            sv = ''
        waitmsg = ''
        jobwait('',qmk_r,Mscfg)
        command = "snqchem.csh -in "+infile[i]+" -out "+outfile[i]+sv
        print "Executing command: "+command
        os.system(command)
    jobs_running = 1
    waitmsg = "All jobs submitted."
    jobwait(waitmsg,qmk_r,Mscfg)
    os.system("grep xxAd *"+filenamebase+outfilenameend+" > adiabatic.txt")
    os.system("grep xxDi *"+filenamebase+outfilenameend+" > diabatic.txt")
    os.system("grep xxDCo *"+filenamebase+outfilenameend+"> dcoupling.txt")
    os.system("grep xxIni *"+filenamebase+outfilenameend+" > idist.txt")
    os.system("grep xxInt *"+filenamebase+outfilenameend+" > mdist.txt")
    os.system("grep xxFi *"+filenamebase+outfilenameend+" > fdist.txt")
    os.system("grep xxRo *"+filenamebase+outfilenameend+" > rotmatrix.txt")
    os.system("grep xxXMATi *"+filenamebase+outfilenameend+" > xmati.txt")
    os.system("grep xxYMATi *"+filenamebase+outfilenameend+" > ymati.txt")
    os.system("grep xxZMATi *"+filenamebase+outfilenameend+" > zmati.txt")
    os.system("grep xxXMATm *"+filenamebase+outfilenameend+" > xmatm.txt")
    os.system("grep xxYMATm *"+filenamebase+outfilenameend+" > ymatm.txt")
    os.system("grep xxZMATm *"+filenamebase+outfilenameend+" > zmatm.txt")
    os.system("grep xxXMATf *"+filenamebase+outfilenameend+" > xmatf.txt")
    os.system("grep xxYMATf *"+filenamebase+outfilenameend+" > ymatf.txt")
    os.system("grep xxZMATf *"+filenamebase+outfilenameend+" > zmatf.txt")
    os.system("grep xxXDIPi *"+filenamebase+outfilenameend+" > xdipi.txt")
    os.system("grep xxYDIPi *"+filenamebase+outfilenameend+" > ydipi.txt")
    os.system("grep xxZDIPi *"+filenamebase+outfilenameend+" > zdipi.txt")
    os.system("grep xxXDIPm *"+filenamebase+outfilenameend+" > xdipm.txt")
    os.system("grep xxYDIPm *"+filenamebase+outfilenameend+" > ydipm.txt")
    os.system("grep xxZDIPm *"+filenamebase+outfilenameend+" > zdipm.txt")
    os.system("grep xxXDIPf *"+filenamebase+outfilenameend+" > xdipf.txt")
    os.system("grep xxYDIPf *"+filenamebase+outfilenameend+" > ydipf.txt")
    os.system("grep xxZDIPf *"+filenamebase+outfilenameend+" > zdipf.txt")
    os.system("grep xxDC2s  *"+filenamebase+outfilenameend+" > dc2sum.txt")
    os.system("grep xxS2i *"+filenamebase+outfilenameend+" > s2i.txt")
    os.system("grep xxS2m *"+filenamebase+outfilenameend+" > s2m.txt")
    os.system("grep xxS2f *"+filenamebase+outfilenameend+" > s2f.txt")
