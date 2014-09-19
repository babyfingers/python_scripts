# makeinput.py
#!/usr/local/bin/python
import os.path
import sys
import time
import math
import re

#This program is designed to take a specially formatted template input
# file and use it as the basis for constructing a full atomic Hessian
# via the finite difference method. Variables you should be aware of:
#
# natom: # of atoms you want to move in this system
# h: the size of the nuclear coordinate incrementation for the FDM calculation (Angstroms)
# maxj: max number of jobs this program will have running at once
# tm: time that the program rests before checking qstat to see how many jobs it has running. (seconds)
# startjob: default = 1. The job that this program will run first. Can be used to restart the program
#           from where it left off if the script fails unexpectedly.

#User-adjusted variables (see descriptions above)
natom    = 2
h        = 0.00001
maxj     = 5
tm       = 2
startjob = 0

#Other variables
ndf    = natom*3 #Number of degrees of freedom
nfiles = 2*ndf*ndf +1  #Number of input files to generate/jobs to be run
ifl    = 0 #Integer to keep track of the files as we write/run them
narray = range(natom)
as  = 'a'
rs  = 'r'
ws  = 'w'
xs  = 'x'
ys  = 'y'
zs  = 'z'
eqs = '='
zes = '0'
uns = '_'
sps = ' '
pps = 'pp'
pms = 'pm'
mps = 'mp'
mms = 'mm'
eols = '\n'
tmpfl = 'tempout.txt'
# Options for scf_guess:
#  the reference case should be gwh, the rest should read off of that.
gwhs = 'gwh'
rds  = 'read'


# These are the markers: these variables store the strings that are currently in the
# template that mark the areas we will want to replace when we actually create the
# input files we will run. There is one for each degree of freedom and another for
# the scf_guess rem variable.
msg = "sssggg"
mx = []
my = []
mz = []
for i in narray:
    mx.append(str(i))
infilenameend = '.in'
outfilenameend = '.out'
th = 2*h

#Useful functions

#chars_before_period: takes string as input.
# Returns an integer equal to the number of
# characters before a period appears. If no
# period is encountered, it returns the length
# of the string.
def chars_before_period(str):
    for i in range(len(str)):
        if str[i]=='.':
            return i
    return len(str)

#chars_after_period: takes string as input.
# Returns an integer equal to the number of
# characters counted in a string after a period
# is encountered. If no period is encountered,
# it returns zero.
def chars_after_period(str):
    for i in range(len(str)):
        if str[i]=='.':
            return len(str)-i-1
    return 0

#Begin the program proper

# Get some user input
dir = raw_input("Directory name:\n")
filenamebase = 'temporary_string'
lfnb = len(filenamebase)
while (lfnb>8):
    filenamebase = raw_input("Base input file name (max eight characters):\n")
    lfnb = len(filenamebase)
    if (lfnb>8):
        print "Your input, \'" + filenamebase + "\', has too many characters! Try again."
mode1 = input("Just make files (1), just run files (2), or make and run (3)?\n")

# Make a list of strings comprising the necessary input and output file names
max_cbpx = 2
for i in narray:
    cbpx = chars_before_period(mx[i])
    if (cbpx>max_cbpx):
	max_cbpx = cbpx
for i in narray:
    while (chars_before_period(mx[i])<max_cbpx):
	mx[i] = zes + mx[i]
    mx[i] = mx[i]
for i in narray:
    my.append(mx[i])
    mz.append(mx[i])
    mx[i] = xs + mx[i]
    my[i] = ys + my[i]
    mz[i] = zs + mz[i]
m = mx +my +mz

infile = []
outfile = []
#The initial (reference) run, without modification
infile.append(filenamebase+uns+infilenameend)
outfile.append(uns+filenamebase+outfilenameend)
ifl += 1
#The off-diagonal elements
for dfi in range(ndf):
    for dfj in range(dfi+1,ndf):
        dfi_s  = m[dfi]
        dfj_s  = m[dfj]
	inflbs = filenamebase+uns+dfi_s+dfj_s
	outflbs = dfi_s+dfj_s+uns+filenamebase
	infile.append(inflbs+pps+infilenameend)
	infile.append(inflbs+mps+infilenameend)
	infile.append(inflbs+pms+infilenameend)
	infile.append(inflbs+mms+infilenameend)
	outfile.append(pps+outflbs+outfilenameend)
	outfile.append(mps+outflbs+outfilenameend)
	outfile.append(pms+outflbs+outfilenameend)
	outfile.append(mms+outflbs+outfilenameend)
        ifl += 4
for dfi in range(ndf):
    dfi_s = m[dfi]
    inflbs = filenamebase+uns+dfi_s
    outflbs = dfi_s+uns+filenamebase
    infile.append(inflbs+pps+infilenameend)
    infile.append(inflbs+mms+infilenameend)
    outfile.append(pps+outflbs+outfilenameend)
    outfile.append(mms+outflbs+outfilenameend)
    ifl += 2
# Open up the formatted template file and get to work
#Part 1: writing the input files
if (mode1==1 or mode1==3):
        ifl = 0 #integer to keep track of the number of files written
	template = open(filenamebase+infilenameend,rs)
	template_s = template.read()
	os.system("mkdir "+dir)
	os.chdir(dir)

	# Now we have generated our full list of markers for every variable. Each has the form "[c][n]=",
	# where [c] is a cartesian axis ('x','y', or 'z'), and [n] is the index of the variable. We have used
	# the function chars_before_period to make sure that all strings [n] have the same length.
	dec_re  = re.compile("[0-9]*\.[0-9]+") #Regular expression for a number containing a decimal place
	#First: do reference case
	f = open(infile[0],ws)
	f.write(template_s.replace(msg,gwhs))
	f.close()
        ifl += 1
	#Second: do the off-diagonal case
	for dfi in range(ndf):
	    for dfj in range(dfi+1,ndf):
		intro  = ''
		bridge = ''
		outro  = ''
		dfi_s  = m[dfi]
		dfj_s  = m[dfj]
		#Search for the lines containing the flag variables
		# We will use regular expressions. Each line where the variables are defined must have the
		# format startline-variablename-equals sign-blank spaces-variable value-endline.
		# For example, it might look like x01=      4.2
		li_re = re.compile(dfi_s + eqs + "\s*[0-9]*\.[0-9]+")
		lj_re = re.compile(dfj_s + eqs + "\s*[0-9]*\.[0-9]+")
		tmplist_i = re.findall(li_re,template_s)
		tmplist_j = re.findall(lj_re,template_s)
		li = tmplist_i[0]
		lj = tmplist_j[0]
		# Take these lines and extract the values of the variables.
		tmplist_i = re.findall(dec_re,li)
		tmplist_j = re.findall(dec_re,lj)
		dfi_s = tmplist_i[0]
		dfj_s = tmplist_j[0]
		dfi_f = float(dfi_s)
		dfj_f = float(dfj_s)
		# Increment/decrement as necessary, revert to strings and complete the line
		lip = li.replace(dfi_s,str(dfi_f + h))
		lim = li.replace(dfi_s,str(dfi_f - h))
		ljp = lj.replace(dfj_s,str(dfj_f + h))
		ljm = lj.replace(dfj_s,str(dfj_f - h))
		# Open the four relevant input files
		fpp = open(infile[ifl],ws)
		fmp = open(infile[ifl+1],ws)
		fpm = open(infile[ifl+2],ws)
		fmm = open(infile[ifl+3],ws)
		# Write to them with relevant replacements
		fpp.write(template_s.replace(li,lip).replace(lj,ljp).replace(msg,rds))
		fmp.write(template_s.replace(li,lim).replace(lj,ljp).replace(msg,rds))
		fpm.write(template_s.replace(li,lip).replace(lj,ljm).replace(msg,rds))
		fmm.write(template_s.replace(li,lim).replace(lj,ljm).replace(msg,rds))
		# Close the four input files
		fpp.close()
		fmp.close()
		fpm.close()
		fmm.close()
                ifl += 4
	#Third: do the diagonal case (those in which we increment a single degree of freedom
	# twice, instead of incrementing two degrees of freedom once.)
	for dfi in range(ndf):
	    dfi_s = m[dfi]
	    #Search for the lines containing the variable marker of interest
	    li_re = re.compile(dfi_s + eqs + "\s*[0-9]*\.[0-9]+")
	    tmplist_i = re.findall(li_re,template_s)
	    li = tmplist_i[0]
	    # Take this line and extract the values of the variables.
	    tmplist_i = re.findall(dec_re,li)
	    dfi_s = tmplist_i[0]
	    dfi_f = float(dfj_s)
	    # Increment/decrement as necessary, revert to strings and complete the line
	    lip = li.replace(dfi_s,str(dfi_f + th))
	    lim = li.replace(dfi_s,str(dfi_f - th))
	    # Open the two relevant input files
	    fpp = open(infile[ifl],ws)
	    fmm = open(infile[ifl+1],ws)
	    # Write to them with relevant replacements
	    fpp.write(template_s.replace(li,lip).replace(msg,rds))
	    fmm.write(template_s.replace(li,lim).replace(msg,rds))
	    # Close the two input files
	    fpp.close()
	    fmm.close()
            ifl += 2
if (mode1==1):
    print str(ifl) + " files made! My work here is done."
elif (mode1==3):
    print str(ifl) + " files made! Now I will run them."
#Part 2: running the files
# Since we are running a lot of jobs, we want to make sure that we don't monopolize the cluster,
# so we will only run as many as maxj jobs at a time. Consequently, we need to have the script
# constantly checking how many jobs it has running and it will only submit more when this number
# drops below maxj. We also want to keep track of how many jobs we have already submitted in
# case this script unexpectedly fails; then we can start it back up again from where we left off.
# We will want to save this number to some log file.
if (mode1==2):
    #If we are running from previously made files, move into the
    # desired directory.
    os.chdir(dir)
if (mode1==2 or mode1==3):
    #Declare some variables
    ifl = 0 #integer to keep track of files run
    logn = filenamebase + ".slog" #Name of script log file
    firsttime = 0.0 #Used to store the time of the first event
    lasttime = 0.0; #Used to store the time of the previous event
    thistime = 0.0; #Used to store the time of the most recent event
    qstat = "qstat > " + tmpfl
    qsbintro = "myqchem.csh -in "
    qsbmezzo = " -out "
    qsboutro0 = " -outdir scf_read -save standard"
    qsboutro1 = " -indir scf_read"
    qmk_r = re.compile(uns + filenamebase + outfilenameend + ".* R batch")
    #Run all of the jobs. Run the reference job first to get a good input for
    # scf_guess= read. After that, run at most maxj jobs at a time.
    for ifl in range(startjob,nfiles):
        if (ifl>0):
            os.system(qstat)
            f = open(tmpfl,rs)
            fs = f.read()
            f.close()
            nrunning = len(re.findall(qmk_r,fs))
            while (nrunning > maxj):
                print "Waiting " + str(tm) + " seconds."
                time.sleep(tm)
                os.system(qstat)
                f = open(tmpfl,rs)
                fs = f.read()
                f.close()
                nrunning = len(re.findall(qmk_r,fs))
            command = qsbintro + infile[ifl] + qsbmezzo + outfile[ifl] + qsboutro1
            os.system(command)
            thistime = time.time()
            flog = open(logn,as)
            flog.write(str(thistime-lasttime)+" seconds have passed since the last command."+eols)
            flog.write("Running command: "+command+eols)
            flog.write("Job number: " +str(ifl) + eols)
            flog.close()
            lasttime = thistime
            print "Running command:\n>" + command
        else:
            command = qsbintro + infile[ifl] + qsbmezzo + outfile[ifl] + qsboutro0
            print "Running command:\n>" + command
            os.system(command)
            firsttime = time.time()
            lasttime = firsttime
            flog = open(logn,ws)
            flog.write("Running command: "+command+eols)
            flog.write("Job number: " + str(ifl) + eols)
            flog.close()
            os.system(qstat)
            f = open(tmpfl,rs)
            fs = f.read()
            f.close()
            nrunning = len(re.findall(qmk_r,fs))
            while (nrunning>0):
                print "Waiting " + str(tm) + " seconds."
                time.sleep(tm)
                os.system(qstat)
                f = open(tmpfl,rs)
                fs = f.read()
                f.close()
                nrunning = len(re.findall(qmk_r,fs))
    thistime = time.time()
    os.system("rm " + tmpfl)
    flog = open(logn,as)
    flog.write(eols+"All jobs submitted. Total time = " + str(thistime-firsttime) + " seconds.")
    flog.close()
    os.system("grep xxAd *"+filenamebase+outfilenameend+" > adiabatic.txt")
    os.system("grep xxDi *"+filenamebase+outfilenameend+" > diabatic.txt")
    os.system("grep xxDCo *"+filenamebase+outfilenameend+" > dcoupling.txt")
    os.system("grep xxIn *"+filenamebase+outfilenameend+" > idist.txt")
    os.system("grep xxFi *"+filenamebase+outfilenameend+" > fdist.txt")
    os.system("grep xxRo *"+filenamebase+outfilenameend+" > rotmatrix.txt")
    os.system("grep xxXMATi *"+filenamebase+outfilenameend+" > xmati.txt")
    os.system("grep xxYMATi *"+filenamebase+outfilenameend+" > ymati.txt")
    os.system("grep xxZMATi *"+filenamebase+outfilenameend+" > zmati.txt")
    os.system("grep xxXMATf *"+filenamebase+outfilenameend+" > xmatf.txt")
    os.system("grep xxYMATf *"+filenamebase+outfilenameend+" > ymatf.txt")
    os.system("grep xxZMATf *"+filenamebase+outfilenameend+" > zmatf.txt")
    os.system("grep xxXDIPi *"+filenamebase+outfilenameend+" > xdipi.txt")
    os.system("grep xxYDIPi *"+filenamebase+outfilenameend+" > ydipi.txt")
    os.system("grep xxZDIPi *"+filenamebase+outfilenameend+" > zdipi.txt")
    os.system("grep xxXDIPf *"+filenamebase+outfilenameend+" > xdipf.txt")
    os.system("grep xxYDIPf *"+filenamebase+outfilenameend+" > ydipf.txt")
    os.system("grep xxZDIPf *"+filenamebase+outfilenameend+" > zdipf.txt")
    os.system("grep xxDC2s  *"+filenamebase+outfilenameend+" > dc2sum.txt")

if (mode1==2 or mode1==3):
    print "Files run successfully. My work here is done."
