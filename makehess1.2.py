# Makeinput.py
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
# startjob: default = 0. The job that this program will run first. Can be used to restart the program
#           from where it left off if the script fails unexpectedly.

#User-adjusted variables (see descriptions above)
natom    = 46
h        = 0.0001
maxj     = 128
tm       = 120
startjob = 0

#Ask user for the kind of calculation to do: gradient or Hessian?
goh = input("Do a gradient (1) or Hessian (2) calculation?\n")
ndf = natom*3 #Number of degrees of freedom
if (goh==1):
    nfiles = 2*ndf +1  #Number of input files to generate/jobs to be run for gradient calculation
elif (goh==2):
    nfiles = 2*ndf*ndf + 1 #Number of input files to generate/jobs to be run for Hessian calculation

#Other variables
ifl    = 0 #Integer to keep track of the files as we write/run them
narray = range(natom)
as  = 'a'
rs  = 'r'
ws  = 'w'
xs  = 'x'
ys  = 'y'
zs  = 'z'
Ss  = 'S'
sps = ' '
eqs = '='
zes = '0'
uns = '_'
sps = ' '
tbs = '\t'
cls = ':'
ps  = 'p'
ms  = 'm'
pps = 'pp'
pms = 'pm'
mps = 'mp'
mms = 'mm'
eols  = '\n'
txts  = '.txt'
grads = 'grad_'
hesss = 'hess_'
dim = [xs,ys,zs]
ndim = len(dim)
tmpfl = 'tempout.txt'
# Regular expressions
dec_re  = re.compile("-?[0-9]*\.[0-9]+") #Regular expression for a double containing a decimal place
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs

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
    mx.append(str(i+1))

#Other variables
infilenameend = '.in'
outfilenameend = '.out'
th = 2*h
thin = 1/th
thin2 = thin*thin

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
if (goh==1):
    max_char = 8
elif (goh==2):
    max_char = 3
while (lfnb>max_char):
    filenamebase = raw_input("Base input file name (max " + str(max_char) + " characters):\n")
    lfnb = len(filenamebase)
    if (lfnb>max_char):
        print "Your input, \'" + filenamebase + "\', has too many characters! Try again."
mode1 = input("Just make files (1), just run files (2), make and run (3), or just do analysis (4)?\n")

# Make a string containing information necessary for constructing file names
#and identifying relevant variables.
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
# Create file for the initial (reference) run, without modification
infile.append(filenamebase+uns+infilenameend)
outfile.append(uns+filenamebase+outfilenameend)
ifl += 1
if (goh==1):
	#Create file names necessary to calculate gradient
	for dfi in range(ndf):
             dfi_s = m[dfi]
             inflbs = filenamebase+uns+dfi_s
             outflbs = dfi_s+uns+filenamebase
             infile.append(inflbs+ps+infilenameend)
             infile.append(inflbs+ms+infilenameend)
             outfile.append(ps+outflbs+outfilenameend)
             outfile.append(ms+outflbs+outfilenameend)
             ifl += 2
elif (goh==2):
	#Create list of file names necessary to calculate Hessian
	# The off-diagonal elements
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
	# The diagonal elements
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
	#First: do reference case
	f = open(infile[0],ws)
	f.write(template_s.replace(msg,gwhs))
	f.close()
        ifl += 1
        if (goh==1):
        	#Write all of the files necessary to calculate the gradient
		for dfi in range(ndf):
                    dfi_s = m[dfi]
                    #Search for the lines containing the variable marker of interest
                    li_re = re.compile(dfi_s + eqs + "\s*-?[0-9]*\.[0-9]+")
                    tmplist_i = re.findall(li_re,template_s)
                    li = tmplist_i[0]
                    # Take this line and extract the values of the variables.
                    tmplist_i = re.findall(dec_re,li)
                    dfi_s = tmplist_i[0]
                    dfi_f = float(dfi_s)
                    # Increment/decrement as necessary, revert to strings and complete the line
                    lip = li.replace(dfi_s,str(dfi_f + h))
                    lim = li.replace(dfi_s,str(dfi_f - h))
                    # Open the two relevant input files
                    fp = open(infile[ifl],ws)
                    fm = open(infile[ifl+1],ws)
                    # Write to them with relevant replacements
                    fp.write(template_s.replace(li,lip).replace(msg,rds))
                    fm.write(template_s.replace(li,lim).replace(msg,rds))
                    # Close the two input files
                    fp.close()
                    fm.close()
                    ifl += 2
        elif (goh==2):
		#Write all of the files necessary to calculate the Hessian
		#Step 1: do the off-diagonal case
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
			li_re = re.compile(dfi_s + eqs + "\s*-?[0-9]*\.[0-9]+")
			lj_re = re.compile(dfj_s + eqs + "\s*-?[0-9]*\.[0-9]+")
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
		#Step 2: do the diagonal case (those in which we increment a single degree of freedom
		# twice, instead of incrementing two degrees of freedom once.)
		for dfi in range(ndf):
		    dfi_s = m[dfi]
		    #Search for the lines containing the variable marker of interest
		    li_re = re.compile(dfi_s + eqs + "\s*-?[0-9]*\.[0-9]+")
		    tmplist_i = re.findall(li_re,template_s)
		    li = tmplist_i[0]
		    # Take this line and extract the values of the variables.
		    tmplist_i = re.findall(dec_re,li)
		    dfi_s = tmplist_i[0]
		    dfi_f = float(dfi_s)
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
if (mode1==2 or mode1==4):
    #If we are running from previously made files, move into the
    # desired directory.
    os.chdir(dir)
if (mode1==2 or mode1==3):
    #Declare/re-initialize some variables
    ifl = 0 #integer to keep track of files run
    firsttime = 0.0 #Used to store the time of the first event
    lasttime = 0.0; #Used to store the time of the previous event
    thistime = 0.0; #Used to store the time of the most recent event
    qstat = "qstat > " + tmpfl
    qsbintro = "myqchem.csh -in "
    qsbmezzo = " -out "
    qsboutro0 = " -outdir scf_read -save standard"
    qsboutro1 = " -indir scf_read"
    logn = filenamebase + ".slog" #Name of script log file
    qmk_r = re.compile(uns + filenamebase + outfilenameend + ".* [RQ] batch")
    #Run all of the jobs. Run the reference job first to get a good input for
    # scf_guess= read. After that, run at most maxj jobs at a time.
    for ifl in range(startjob,nfiles):
        if (ifl>0):
            time.sleep(2)
            os.system(qstat)
            f = open(tmpfl,rs)
            fs = f.read()
            f.close()
            nrunning = len(re.findall(qmk_r,fs))
            while (nrunning >= maxj):
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
        elif (ifl==0):
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
    flog.write(eols+"All jobs finished. Total time = " + str(thistime-firsttime) + " seconds.")
    flog.close()
    os.system("grep xxAd *"+filenamebase+outfilenameend+" > adiabatic.txt")
    os.system("grep xxDi *"+filenamebase+outfilenameend+" > diabatic.txt")
    os.system("grep xxDCo *"+filenamebase+outfilenameend+" > dcoupling.txt")
    os.system("grep xxIn *"+filenamebase+outfilenameend+" > idist.txt")
    os.system("grep xxFi *"+filenamebase+outfilenameend+" > fdist.txt")
    os.system("grep xxXMATi *"+filenamebase+outfilenameend+" > xmati.txt")
    os.system("grep xxYMATi *"+filenamebase+outfilenameend+" > ymati.txt")
    os.system("grep xxZMATi *"+filenamebase+outfilenameend+" > zmati.txt")
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
    os.system("grep xxUxx *"+filenamebase+outfilenameend+" > U.txt")
    flog = open(logn,as)
    flog.write(eols+eols+"Quantities grepped and piped to .txt files."+eols)
    flog.close()
if (mode1==2 or mode1==3):
    print "Files run successfully."

if (mode1==2 or mode1==3 or mode1==4):
    print "Now analyzing the output."
    os.system("ls *txt > "+tmpfl)
    f = open(tmpfl,rs)
    quant = []
    quantfn = []
    for line in f:
        quant.append(line[0:-5])
        quantfn.append(line[0:-1])
    f.close()
    # Remove a couple troublesome items
    quant.remove('U')
    quantfn.remove('U.txt')
    quant.remove('xmati')
    quantfn.remove('xmati.txt')
    quant.remove('ymati')
    quantfn.remove('ymati.txt')
    quant.remove('zmati')
    quantfn.remove('zmati.txt')
    quant.remove('xmatf')
    quantfn.remove('xmatf.txt')
    quant.remove('ymatf')
    quantfn.remove('ymatf.txt')
    quant.remove('zmatf')
    quantfn.remove('zmatf.txt')
    quant.remove('dc2sum')
    quantfn.remove('dc2sum.txt')
    quant.remove('s2i')
    quantfn.remove('s2i.txt')
    quant.remove('s2f')
    quantfn.remove('s2f.txt')
    numq = len(quant)
    os.system("rm "+tmpfl)
    
    # Get names of states and save to array
    os.system("grep xxSta *"+filenamebase+outfilenameend+" > states")
    f = open("states",rs)
    states_raw = f.readline()
    f.close()
    os.system("rm states")
    int_re = re.compile("\s-?[0-9]+")
    states = re.findall(int_re,states_raw)
    numkeep = len(states)
    print "numkeep = " +str(numkeep)
    for i in range(numkeep):
        s = states[i]
        states[i] = s[1:]
    if (goh==1):
        print "Will generate and save gradients to .txt files."
        quantgfn = []
        for q in quant:
            quantgfn.append(grads+q+txts)
        for q in range(numq):
            print "Now working on the quantity " + str(quant[q])
            fq = open(quantfn[q],rs)
            fq_s = fq.read()
            fq.close()
            fq_s = eols+fq_s
            gradq = []
            # Extract relevant values as floats and compute gradient
            for dfi in range(ndf):
                dfi_s = m[dfi]
		p_re = re.compile(eols+ps+dfi_s+uns+".*\n")
		m_re = re.compile(eols+ms+dfi_s+uns+".*\n")
		p_raw = re.findall(p_re,fq_s)
		m_raw = re.findall(m_re,fq_s)
		p_ar = re.findall(dp_re,p_raw[0])
		m_ar = re.findall(dp_re,m_raw[0])
		l_p = len(p_ar)
		l_m = len(m_ar)
		if (l_p != l_m):
		    print "Problem with " + quant[q]
		    print "l_p = " + str(l_p) + ", l_m = " + str(l_m)
		    sys.exit("Error: output array lengths are not the same.")
		elif (l_p != numkeep):
		    print "Problem with " + quant[q]
		    print "l_p = l_m = " + str(l_p) + ", numkeep = " + numkeep
		    sys.exit("Error: output array lengths are not equal to the number of states considered in diabatization.")
		for j in range(l_p):
		    pfloat_j = float(p_ar[j])
		    mfloat_j = float(m_ar[j])
		    gradq.append(thin*(pfloat_j-mfloat_j))
            # Write it all to a .txt file
	    fq = open(quantgfn[q],ws)
            k = 0
            fq.write("states:")
            for j in states:
                fq.write(sps + j)
            fq.write(eols)
            for dfi in range(ndf):
                dfi_s = m[dfi]
                fq.write(dfi_s + cls)
                for j in states:
                    fq.write(sps + str(gradq[k]))
                    k += 1
                fq.write(eols)
            fq.close()
    elif (goh==2):
        print "Will generate and save Hessian matrices to .txt files."
        quanthfn = []
        for q in quant:
            for j in states:
                quanthfn.append(hesss+q+uns+j+txts)
        # Initialize Hessian matrix
        hessq = []
        for dfi in range(ndf):
            hessq.append([])
            for dfj in range(ndf):
                hessq[dfi].append([])
                for k in range(numkeep):
                    hessq[dfi][dfj].append(0.0)
        # Now, extract necessary quantities from output files
        for q in range(numq):
            print "Now working on the quantity " + str(quant[q])
            fq = open(quantfn[q],rs)
            fq_s = fq.read()
            fq.close()
            fq_s = eols + fq_s
            gradq = []
            # Do off-diagonal case
            for dfi in range(ndf-1):
                dfi_s = m[dfi]
		for dfj in range(dfi+1,ndf):
                    #print "dfi, dfj = " + str(dfi) + ", " + str(dfj)
		    dfj_s = m[dfj]
		    pp_re = re.compile(eols+pps+dfi_s+dfj_s+uns+".*\n")
		    pm_re = re.compile(eols+pms+dfi_s+dfj_s+uns+".*\n")
		    mp_re = re.compile(eols+mps+dfi_s+dfj_s+uns+".*\n")
		    mm_re = re.compile(eols+mms+dfi_s+dfj_s+uns+".*\n")
		    pp_raw = re.findall(pp_re,fq_s)
		    pm_raw = re.findall(pm_re,fq_s)
		    mp_raw = re.findall(mp_re,fq_s)
		    mm_raw = re.findall(mm_re,fq_s)
		    pp_ar = re.findall(dp_re,pp_raw[0])
		    pm_ar = re.findall(dp_re,pm_raw[0])
		    mp_ar = re.findall(dp_re,mp_raw[0])
		    mm_ar = re.findall(dp_re,mm_raw[0])
		    l_pp = len(pp_ar)
		    l_pm = len(pm_ar)
		    l_mp = len(mp_ar)
		    l_mm = len(mm_ar)
		    if (l_pp != l_pm)|(l_pp != l_mp)|(l_pp != l_mm)|(l_pm != l_mp)|(l_pm != l_mm)|(l_mp != l_mm):
			print "Problem with " + quant[q]
			print "l_pp = " + str(l_pp) + ", l_pm = " + str(l_pm) + ", l_mp = " + str(l_mp) + ", l_mm = " + str(l_mm)
			sys.exit("Error: output array lengths are not the same.")
		    elif (l_pp != numkeep):
			print "Problem with " + quant[q]
			print "l_pp = l_pm = etc. = " + str(l_pp) + ", numkeep = " + numkeep
			sys.exit("Error: output array lengths are not equal to the number of states considered in diabatization.")
		    for k in range(numkeep):
			ppfloat_k = float(pp_ar[k])
			pmfloat_k = float(pm_ar[k])
			mpfloat_k = float(mp_ar[k])
			mmfloat_k = float(mm_ar[k])
                        hessq[dfi][dfj][k] = thin2*(ppfloat_k-pmfloat_k-mpfloat_k+mmfloat_k)
                        hessq[dfj][dfi][k] = hessq[dfi][dfj][k]
            # Do diagonal case
            #  First, get term from reference configuration (this will be used in the
            #   calculation of all the diagonal Hessian terms.)
            ref_re = re.compile(eols + uns+filenamebase+".*\n")
            ref_raw = re.findall(ref_re,fq_s)
            ref_ar = re.findall(dp_re,ref_raw[0])
            l_ref = len(ref_ar)
            if (l_ref != numkeep):
                print "Problem with " + quant[q]
                print "l_ref = " + str(l_ref) + ", numkeep = " + str(numkeep)
                sys.exit("Error: output array lengths are not equal to the number of states considered in diabatization.")
            ref_float = []
            for k in range(numkeep):
                ref_float.append(float(ref_ar[k]))
            #  And now the rest
            for dfi in range(ndf):
                dfi_s = m[dfi]
                pp_re = re.compile(eols + pps+dfi_s+uns+filenamebase+".*\n")
                mm_re = re.compile(eols + mms+dfi_s+uns+filenamebase+".*\n")
                pp_raw = re.findall(pp_re,fq_s)
                mm_raw = re.findall(mm_re,fq_s)
                pp_ar = re.findall(dp_re,pp_raw[0])
                mm_ar = re.findall(dp_re,mm_raw[0])
                l_pp = len(pp_ar)
                l_mm = len(mm_ar)
                if (l_pp != l_mm):
                    print "Problem with " + quant[q]
                    print "l_pp = " + str(l_pp) + ", l_mm = " + str(l_mm)
                    sys.exit("Error: output array lengths are not the same.")
                elif (l_pp != numkeep):
                    print "Problem with " + quant[q]
                    print "l_pp = l_mm = " + str(l_pp) + ", numkeep = " + numkeep
                    sys.exit("Error: output array lengths are not equal to the number of states considered in diabatization.")
                for k in range(numkeep):
                    ppfloat_k = float(pp_ar[k])
                    mmfloat_k = float(mm_ar[k])
                    hessq[dfi][dfi][k] = thin2*(ppfloat_k-2*ref_float[k]+mmfloat_k)
            # Write it up to a .txt file
            for k in range(numkeep):
                fqk = open(quanthfn[numkeep*q +k],ws)
                fqk.write(Ss + states[k] + cls)
                for dfi in range(ndf):
                    dfi_s = m[dfi]
                    fqk.write(tbs + dfi_s)
                fqk.write(eols)
                for dfi in range(ndf):
                    dfi_s = m[dfi]
                    fqk.write(dfi_s + cls)
                    for dfj in range(ndf):
                        fqk.write(tbs + str(hessq[dfi][dfj][k]))
                    fqk.write(eols)
                fqk.close()
