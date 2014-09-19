# makeinput.py
#!/usr/local/bin/python
import os.path
import time
#This program takes a single input file and increments the values of a single
# parameter (from 'pta' to 'ptb' in increments of 'inc'; the parameter is marked
# by the initial string in its line given by 'poi') and generates a ton of input
# files and optionally runs them and extracts certain info.

# Variables: you probably will only need to change things on the following few lines.
poi1 = 'F'  #Some marker for the command line of interest
pta  = 1.00  #Start the value where?
ptb  = 10.0  #End the value where?
inc  = 0.05
poi2 = 'scf_guess'
fileloc = '/data/home/alguire/inout_files/LiF/'
infilenameend = '.in'
outfilenameend = '.out'
sz       = '0' #Just a string of 0
undrscr  = '_'
spc      = ' '
nfiles = int(round(1+((ptb-pta)/inc)))
narray = range(nfiles)
vals = []
s1 = []
s2 = ['gwh']
rd = 'read'
tmpfl = 'tempout.txt'
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

def cutup(sloi,spoi):
    l_loi = len(sloi)
    l_poi = len(spoi)
    loi_a = '';
    loi_b = '';
    mk = -1
    i = 0
    for j in range(l_loi):
        if (mk == -1):
            if (sloi[i:i+l_poi]==spoi):
                mk = 0
                i += l_poi
        elif (mk == 0):
            if (sloi[i]!=spc):
                loi_a = sloi[:i]
                mk = 1
        elif (mk == 1):
            if (sloi[i]==spc or sloi[i]=='\n'):
                loi_b = sloi[i:]
                break
        i += 1
    return [loi_a,loi_b]

os.chdir(fileloc)
dir = raw_input("Directory name:\n")
filenamebase = raw_input("Base input file name:\n")
mode1 = input("Just make files (1) or make and run (2)?\n")
for i in range(nfiles-1):
    s2.append(rd)
max_cbp_s = 2
max_cap_s = 2
l_poi1 = len(poi1)
l_poi2 = len(poi2)
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
template = open(filenamebase+infilenameend,'r')
intro1 = ''
intro2 = ''
loi1   = ''
loi2   = ''
loi1_a = ''
loi1_b = ''
loi2_a = ''
loi2_b = ''
outro  = ''
reached_poi1 = False
reached_poi2 = False
for line in template:
    if not reached_poi1:
        if poi1 in line:
            loi1 = line
            reached_poi1 = True
        else:
            intro1 += line
    elif not reached_poi2:
        if poi2 in line:
            loi2 = line
            reached_poi2 = True
        else:
            intro2 += line
    else:
        outro += line
template.close()
list1  = cutup(loi1,poi1)
list2  = cutup(loi2,poi2)
loi1_a = list1[0]
loi1_b = list1[1]
loi2_a = list2[0]
loi2_b = list2[1]
os.system("mkdir "+dir)
os.chdir(dir)
infile = []
for i in narray:
    infile.append(filenamebase+s1[i]+infilenameend)
    print "Writing "+infile[i]
    f = open(infile[i],'w')
    f.write(intro1)
    f.write(loi1_a+s1[i]+loi1_b)
    f.write(intro2)
    f.write(loi2_a+s2[i]+loi2_b)
    f.write(outro)
    f.close()
if mode1==2:
    qmk = "C batch" #What to look for in qstat to indicate a finished job
    for i in range(nfiles):
        job_running = -1 #When this changes, the job has completed
        sv = " -indir scf_read -outdir scf_read -save standard"
        outfile = s1[i]+filenamebase+outfilenameend
        command = "myqchem.csh -in "+infile[i]+" -out "+outfile+sv
        print "Executing command: "+command
        os.system(command)
        while (job_running == -1):
            print "Waiting " + str(tm) + " seconds."
            time.sleep(tm)
            os.system("qstat | grep " + outfile + " > "+tmpfl)
            f = open(tmpfl,'r')
            for line in f:
                job_running = line.find(qmk)
            f.close()
            os.system("rm "+tmpfl)

    os.system("grep xxAd *"+filenamebase+outfilenameend+" > adiabatic.txt")
    os.system("grep xxDi *"+filenamebase+outfilenameend+" > diabatic.txt")
    os.system("grep xxDCo *"+filenamebase+outfilenameend+"> dcoupling.txt")
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

