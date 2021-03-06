# mfrog.py
#!/usr/local/bin/python
import os.path
import time
import re

#This program takes a single input file and increments the values of a single
# parameter (from 'pta' to 'ptb' in increments of 'inc'; the parameter is marked
# by the string 'zzz') and generates a ton of input files and optionally runs them
# and extracts certain info.

# v1.1: adds new functionality. Give opportunity to do a bunch of calculations at once
#       instead of just starting where the last one left off. Furthermore, can now
#       accomodate calls to myqchem.csh or solqchem.csh (must be changed manually!)

# Variables: you probably will only need to change things on the following few lines.
maxjobs = 1  #If we just do all the jobs at once, how many should we run at a time?
pta  = 7.00  #Start the value where?
ptb  = 8.0  #End the value where?
inc  = 0.01
m1  = 'zzz'  #Some marker for the position variable
m2 = 'sssggg' #some marker for the scf_guess variable
m3 = 'tttooo' #some marker for the tia_order variable
tmpfl = 'tempout.txt'
os.system("pwd > "+tmpfl)
f = open(tmpfl,'r')
fileloc = f.read()
fileloc = fileloc[0:-1]
f.close()
os.system("rm "+tmpfl)
infilenameend = '.in'
outfilenameend = '.out'
sz  = '0' #Just a string of 0
so  = '1'
su  = '_'
ss  = ' '
nfiles = int(round(1+((ptb-pta)/inc)))
narray = range(nfiles)
vals = []
s1 = []
s2 = []
s3 = []
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

os.chdir(fileloc)
dir = raw_input("Directory name:\n")
filenamebase = raw_input("Base input file name:\n")
mode1 = input("Just make files (1) or make and run (2)?\n")
mode2 = input("Run normal calculation (1) or stability analysis (2)?\n")
if (mode2==1):
    mode3 = input("Input SCF guess from previous calculation (1), from individual set of SCF guess folders (2), or start fresh each time (3)?\n")
else:
    mode3 = 1
if (mode2==1 and mode3==1):
    # Although in principal we can do this for mode2==2, it would require a little reworking of the code.
    mode4 = input("Start fresh with scf_guess = gwh (1) or start where previous calculation left off (2)?\n")
else:
    mode4 = 1
if (mode3==1 and mode4==1):
    s2.append('gwh')
    s3.append(sz)
elif (mode3 != 3):
    s2.append(rd)
    s3.append(so)

if (mode3!=3):
    for i in range(nfiles-1):
        s2.append(rd)
        s3.append(so)
else:
    for i in range(nfiles):
        s2.append('gwh')
        s3.append(so)

#Since we are labelling files based on the zzz variable,
# let's make sure they each have the same number of zeros
# before and after the decimal so that it is neat and
# all the files order themselves well.
max_cbp_s = 2
max_cap_s = 2
for i in narray:
    vals.append(pta + inc*i)
    s1.append(str(vals[i]))
    cbp = chars_before_period(s1[i])
    cap = chars_after_period(s1[i])
    if cbp>max_cbp_s:
        max_cbp_s = cbp
    if cap>max_cap_s:
        max_cap_s = cap
for i in narray:
    while chars_before_period(s1[i])<max_cbp_s:
        s1[i] = sz + s1[i]
    while chars_after_period(s1[i])<max_cap_s:
        s1[i] = s1[i] + sz

#Now open up the template file and replace the relevant markers
# with their respective variables to create each input file.
template = open(filenamebase+infilenameend,'r')
template_s = template.read()
os.system("mkdir "+dir)
os.chdir(dir)
infile = []
for i in narray:
    infile.append(filenamebase+s1[i]+infilenameend)
    print "Writing "+infile[i]
    f_i = open(infile[i],'w')
    f_i.write(template_s.replace(m1,s1[i]).replace(m2,s2[i]).replace(m3,s3[i]))
    f_i.close()
template.close()

if mode1==2:
    #Then we want to run these files, too. Since each job will depend on the
    # last for its scf_guess input, we want to make sure that the last one
    # has finished before we start the next.
    qmk_r = re.compile(filenamebase+outfilenameend+".* [ERQ] batch") #What to look for in qstat to indicate a running job
    #If we're not picking up where an old calculation left off, we have to do job 0 differently
    if (mode4==1 and mode3!=3):
        outfile = s1[0]+filenamebase+outfilenameend
        if (mode2==1 and mode3==1):
            sv = " -outdir scf_read -save diabat"
        elif (mode2==2):
            sv = " -outdir scf_read"+s1[0]+" -save diabat"
        elif (mode3==2):
            sv = " -indir scf_read"+s1[0]
        command = "solqchem.csh -in "+infile[0]+" -out "+outfile+sv
        print "Executing command: "+command
        os.system(command)
        read_job_start = 1
    else:
        read_job_start = 0
    #Proceed with the remaining jobs
    for i in range(read_job_start,nfiles):
        jobs_running = maxjobs+1 #When this drops below 1, the job has completed
        if (mode2==1 and mode3==1):
            sv = " -indir scf_read -outdir scf_read -save diabat"
        elif (mode2==2):
            sv = " -indir scf_read"+s1[i-1]+" -outdir scf_read"+s1[i]+" -save diabat"
        elif (mode3==2):
            sv = " -indir scf_read"+s1[i]
        elif (mode3==3):
            sv = ''
        outfile = s1[i]+filenamebase+outfilenameend
        while ((jobs_running>0 and mode3!=3) or (jobs_running>maxjobs and mode3==3)):
            print "Waiting "+str(tm)+" seconds."
            time.sleep(tm)
            os.system("qstat > "+tmpfl)
            f = open(tmpfl,'r')
            f_s = f.read()
            f.close()
            jobs_running = len(re.findall(qmk_r,f_s))
            os.system("rm "+tmpfl)
        command = "solqchem.csh -in "+infile[i]+" -out "+outfile+sv
        print "Executing command: "+command
        os.system(command)
    jobs_running = 1
    while (jobs_running>0):
        print "All jobs submitted. Waiting "+str(tm)+" seconds."
        time.sleep(tm)
        os.system("qstat > "+tmpfl)
        f = open(tmpfl,'r')
        f_s = f.read()
        f.close()
        jobs_running = len(re.findall(qmk_r,f_s))
        os.system("rm "+tmpfl)
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
