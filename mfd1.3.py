# mfd.py
#!/usr/local/bin/python
import os.path
import time
import math
import sys
import re

# v1.3: modified for the new use of the tia_order variable. As of
# 2011-10-31, this is the code:
# tia_order = -2 cross-geometry overlap calculation
#              0 no overlap (energy ordering)
#              1 'light' overlap (only CIS coefficient overlap)
#              2 full overlap (includes MO information
# Also, finally added in the jobwait subroutine and a few other
# small changes.

# Quick + dirty program to find finite-difference values for state overlap
# between different states at different nuclear geometries. Modified from
# makehess.py. It's actually quite restricted in its scope-- meant to be
# used with a certain modified version of Q-Chem, which calculates a
# numkeep x numkeep matrix of state overlaps between the current and
# previous calculations. (Hence the special Q-Chem script, mmqchem.csh.

#User-adjusted variables
hh       = 0.0001 #Incremental h value
ha       = 0.0001 #Starting h value
hb       = 0.0009 #Ending h value
maxj     = 1 #For now, only works for 1
jobstart = 0
tm       = 1.0
dtia     = 2 #What kind of overlap?

ttia     = str(dtia)
goh = 1 # Needs to be this (cruft from makehess.py)
nh = int((hb-ha)/hh +0.5)+1
print "(hb-ha)/hh  +0.5 = " +str((hb-ha)/hh +0.5)
print "int((hb-ha)/hh +0.5) = " +str(int((hb-ha)/hh +0.5))
print "hh = " +str(hh)
nfiles = 2*nh +1  #Number of input files to generate/jobs to be run
                             # Includes initial point and all possible values
                             # of h times 2 (for +h and -h)
hrange = []
hsrange = []
narray = range(nh)
for h in narray:
    hrange.append(ha + h*hh)
    hsrange.append(str(hrange[h]))
    print "hc = " + hsrange[h]
    
#Other variables
ifl    = 0 #Integer to keep track of the files as we write/run them
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
ons = '1'
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
ng2 = '-2'
eols  = '\n'
txts  = '.txt'
grads = 'grad_'
hesss = 'hess_'
zzzs = 'zzz='
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
m1 = "sssggg"
m2 = "tttooo"

#Other variables
infilenameend = '.i'
outfilenameend = '.o'

#Useful functions

#jobwait: takes string and reg. ex. as
# input to stop the program until a
# certain specified job has completed.
def jobwait(waitmsg,qmk_r):
    jobs_running = maxj+1
    while (jobs_running>=maxj):
        print waitmsg+" Waiting "+str(tm)+" seconds."
        time.sleep(tm)
        os.system("qstat > "+tmpfl)
        f = open(tmpfl,'r')
        f_s = f.read()
        f.close()
        jobs_running = len(re.findall(qmk_r,f_s))
        os.system("rm "+tmpfl)                                                                    

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
max_char = 5
while (lfnb>max_char):
    filenamebase = raw_input("Base input file name (max " + str(max_char) + " characters):\n")
    lfnb = len(filenamebase)
    if (lfnb>max_char):
        print "Your input, \'" + filenamebase + "\', has too many characters! Try again."
mode1 = input("Just make files (1), just run files (2), make and run (3), or just do analysis (4)?\n")

# Make a string containing information necessary for constructing file names
#and identifying relevant variables.
max_cbpx = 1
max_capx = 2
#Check to make sure no variables are in scientific notation before adjusting
# the number of zeroes in the file names.
nsn = True
for s in hsrange:
    if ('e' in s):
        nsn = False
        break
if (nsn):
    for i in narray:
        cbpx = chars_before_period(hsrange[i])
        if (cbpx>max_cbpx):
            max_cbpx = cbpx
    for i in narray:
        while (chars_before_period(hsrange[i])<max_cbpx):
            hsrange[i] = zes + hsrange[i]

    for i in narray:
        capx = chars_after_period(hsrange[i])
        if (capx>max_capx):
            max_capx = capx
    for i in narray:
        while (chars_after_period(hsrange[i])<max_capx):
            hsrange[i] = hsrange[i] + zes
    for i in narray:
        print "hsrange[" +str(i)+ "] = " +hsrange[i]

infile = []
outfile = []
if (dtia==2):
    Sinfile = []
    Soutfile = []
# Create file for the initial (reference) run, without modification
infile.append(filenamebase+uns+infilenameend)
outfile.append(uns+filenamebase+outfilenameend)
if (dtia==2):
    Sinfile.append(Ss+filenamebase+uns+infilenameend) #Just a placeholder--
                                                      # won't actually use
    Soutfile.append(Ss+uns+filenamebase+outfilenameend) #Ditto
ifl += 1
if (goh==1):
    #Create file names necessary to calculate gradient
    for hi in narray:
        hc = hrange[hi]
        hc_s = hsrange[hi]
        inflbs = filenamebase+uns+hc_s
        outflbs = hc_s+uns+filenamebase
        infile.append(inflbs+ps+infilenameend)
        infile.append(inflbs+ms+infilenameend)
        outfile.append(ps+outflbs+outfilenameend)
        outfile.append(ms+outflbs+outfilenameend)
        if (dtia==2):
            Sinfile.append(Ss+inflbs+ps+infilenameend)
            Sinfile.append(Ss+inflbs+ms+infilenameend)
            Soutfile.append(Ss+ps+outflbs+outfilenameend)
            Soutfile.append(Ss+ms+outflbs+outfilenameend)
        ifl += 2
print "Done creating file names, ifl = " + str(ifl) + ", nfiles = " +str(nfiles)

# Open up the formatted template file and get to work
#Part 1: writing the input files
if (mode1==1 or mode1==3):
        ifl = 0 #integer to keep track of the number of files written
	template = open(filenamebase+infilenameend,rs)
	template_s = template.read()
        if (dtia==2):
            Stemplate = open(Ss+filenamebase+infilenameend,rs)
            Stemplate_s = Stemplate.read()
	os.system("mkdir "+dir)
	os.chdir(dir)

	#First: do reference case
	f = open(infile[0],ws)
	f.write(template_s.replace(m1,gwhs).replace(m2,zes))
	f.close()
        ifl += 1
        if (goh==1):
        	#Write all of the files necessary to calculate the gradient
                for hi in narray:
                    hc = hrange[hi]
                    hc_s = hsrange[hi]
                    #Search for the lines containing the variable marker of interest
                    li_re = re.compile(zzzs + "\s*-?[0-9]*\.[0-9]+")
                    tmplist_i = re.findall(li_re,template_s)
                    li = tmplist_i[0]
                    # Take this line and extract the values of the variables.
                    tmplist_i = re.findall(dec_re,li)
                    zzz_s = tmplist_i[0]
                    zzz_f = float(zzz_s)
                    # Increment/decrement as necessary, revert to strings and complete the line
                    lip = li.replace(zzz_s,str(zzz_f +hc))
                    lim = li.replace(zzz_s,str(zzz_f -hc))
                    # Open the two relevant input files
                    fp = open(infile[ifl],ws)
                    fm = open(infile[ifl+1],ws)
                    if (dtia==2):
                        Sfp = open(Sinfile[ifl],ws)
                        Sfm = open(Sinfile[ifl+1],ws)
                    # Write to them with relevant replacements
                    fp.write(template_s.replace(li,lip).replace(m1,rds).replace(m2,ttia))
                    fm.write(template_s.replace(li,lim).replace(m1,rds).replace(m2,ttia))
                    if (dtia==2):
                        Sfp.write(Stemplate_s.replace(li,lip).replace(m1,rds).replace(m2,ng2))
                        Sfm.write(Stemplate_s.replace(li,lim).replace(m1,rds).replace(m2,ng2))
                    # Close the two input files
                    fp.close()
                    fm.close()
                    if (dtia==2):
                        Sfp.close()
                        Sfm.close()
                        ifl += 2
print "Done writing files, ifl = " + str(ifl) + ", nfiles = " +str(nfiles)
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
    qsbintro = "mmqchem.csh -in "
    qsbmezzo = " -out "
    qsboutro0 = " -outdir scf_read -save diabat"
    qsboutro1 = " -indir scf_read"
    Sqsboutro = " -indir scf_read -outdir tmp_dir -save diabat"
    logn = filenamebase + ".slog" #Name of script log file
    qmk_r = re.compile(uns + filenamebase + outfilenameend + ".* [ERQ] batch")
    Sqmk_r = re.compile(Ss +".*"+ uns + filenamebase + outfilenameend + ".* [ERQ] batch")
    #Run all of the jobs. Run the reference job first to get a good input for
    # scf_guess= read. After that, run at most maxj jobs at a time.
    for ifl in range(jobstart,nfiles):
        if (ifl>0):
            #Part one: run the overlap file
            if (dtia==2):
                jobwait('',qmk_r)
                command = qsbintro + Sinfile[ifl] + qsbmezzo + Soutfile[ifl] + Sqsboutro
                print ">" + command
                os.system(command)
                thistime = time.time()
                flog = open(logn,as)
                flog.write(str(thistime-lasttime)+" seconds have passed since the last command."+eols)
                flog.write("Running command: "+command+eols)
                flog.write("Job number: " +str(ifl) + eols)
                flog.close()
                lasttime = thistime
            #Part two: the actual file
            if (dtia==2):
                jobwait('',Sqmk_r)
            else:
                jobwait('',qmk_r)
            #command = "rm " + Soutfile
            #os.system(command)
            #print ">" + command
            command = "mv tmp_dir/1213.0 scf_read/1213.0"
            print ">" + command
            os.system(command)
            command = "rm -r tmp_dir"
            print ">" + command
            os.system(command)
            command = qsbintro +infile[ifl] +qsbmezzo +outfile[ifl] +qsboutro1
            print ">" + command
            os.system(command)
            thistime = time.time()
            flog = open(logn,as)
            flog.write(str(thistime-lasttime)+" seconds have passed since the last command."+eols)
            flog.write("Running command: "+command+eols)
            flog.write("Job number: " +str(ifl) + eols)
            flog.close()
            lasttime = thistime
        elif (ifl==0):
            command = qsbintro +infile[ifl] +qsbmezzo +outfile[ifl] +qsboutro0
            print ">" + command
            os.system(command)
            firsttime = time.time()
            lasttime = firsttime
            flog = open(logn,ws)
            flog.write("Running command: "+command+eols)
            flog.write("Job number: "+str(ifl)+eols)
            flog.close()
            jobwait('',qmk_r)
    jobwait('All jobs submitted. ',qmk_r)
    thistime = time.time()
    os.system("rm " + tmpfl)
    flog = open(logn,as)
    flog.write(eols+"All jobs finished. Total time = " + str(thistime-firsttime) + " seconds.")
    flog.close()
if (mode1==2 or mode1==3):
    print "Files run successfully."

if (mode1==2 or mode1==3 or mode1==4):
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
    os.system("grep xxMOxOvr_stdxx *"+filenamebase+outfilenameend+" > MOxOvr_std.txt")
    os.system("grep xxMOxOvr_Newxx *"+filenamebase+outfilenameend+" > MOxOvr_New.txt")
    os.system("grep xxaMOxOvr_stdaaxx *"+filenamebase+outfilenameend+" > aMOxOvr_stdaa.txt")
    os.system("grep xxaMOxOvr_stdabxx *"+filenamebase+outfilenameend+" > aMOxOvr_stdab.txt")
    os.system("grep xxaMOxOvr_Newaaxx *"+filenamebase+outfilenameend+" > aMOxOvr_Newaa.txt")
    os.system("grep xxaMOxOvr_Newabxx *"+filenamebase+outfilenameend+" > aMOxOvr_Newab.txt")
    print "Now analyzing the output."
    quant = ["MOxOvr_std","aMOxOvr_stdaa","aMOxOvr_stdab","MOxOvr_New","aMOxOvr_Newaa","aMOxOvr_Newab"]
    quantfn = ["MOxOvr_std.txt","aMOxOvr_stdaa.txt","aMOxOvr_stdab.txt","MOxOvr_New.txt","aMOxOvr_Newaa.txt","aMOxOvr_Newab.txt"]
    numq = len(quant)
    # Get names of states and save to array
    os.system("grep xxSta *"+filenamebase+outfilenameend+" > states")
    f = open("states",rs)
    states_raw = f.readline()
    f.close()
    os.system("rm states")
    int_re = re.compile("\s-?[0-9]+")
    states = re.findall(int_re,states_raw)
    numkeep = len(states)
    numkeep2 = numkeep*numkeep
    print "numkeep = " +str(numkeep)
    for i in range(numkeep):
        s = states[i]
        states[i] = s[1:]
    if (goh==1):
        print "Will generate and save gradients to .txt files."
        quantgfn = []
        for q in quant:
            quantgfn.append(grads+q+txts)
        #Loop over quantities (in this case, all are measurements of state overlap)
        for q in range(numq):
            print "Now working on the quantity " + str(quant[q])
            fq = open(quantfn[q],rs)
            fq_s = fq.read()
            fq.close()
            fq_s = eols +fq_s
            print "fq_s = " +fq_s
            gradq = []
            # Extract relevant values as floats and compute gradient
            #Loop over all the different hs (this calculation can do a finite-difference
            # approximation with many different values of h--the finite difference-- as
            # parameters. Each might give a slightly different result until h is
            # sufficiently small as to converge on the actual derivative.)
            for hi in narray:
                hc = hrange[hi]
                hc_s = hsrange[hi]
                thinc = 1/(1.889725989*2*hc) #hc is in A, but we want our gradients in Bohr
                print "hc = " +hc_s
                print "thinc = " +str(thinc)
                #Get a list of all "plus" (x+h) and "minus" (x-h) values for overlap from the
                #output file we are currently interested in.
		p_re = re.compile(eols+ps+hc_s+uns+".*")
		m_re = re.compile(eols+ms+hc_s+uns+".*")
		p_raw = re.findall(p_re,fq_s)
		m_raw = re.findall(m_re,fq_s)
                print "len(p_raw) = " +str(len(p_raw))
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
            print "Total number of gradients = ",len(gradq)
	    fq = open(quantgfn[q],ws)
            k = 0
            fq.write("states:")
            for I in states:
                for J in states:
                    fq.write(sps +"<" +I +"|d" +J +">")
            fq.write(eols)
            for hi in narray:
                hc = hrange[hi]
                hc_s = hsrange[hi]
                fq.write(hc_s +cls)
                for I in states:
                    for J in states:
                        fq.write(sps +str(gradq[k]))
                        k += 1
                fq.write(eols)
            fq.close()
