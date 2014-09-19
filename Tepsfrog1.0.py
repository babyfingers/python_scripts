# Tepsfrog.py
#!/usr/local/bin/python
import os.path
import time
import re

#Tepsfrog.py: a modified version of Snqchem.csh that will be used to run
# multiple copies of Snqchem.csh simultaneously at different settings for
# epsilon. Where Snqchem.csh would take a template of a Q-Chem input file
# and make copies for different calculations, this should do the same for
# a template file of Snqchem.csh (titled FRG.py for now).

# Variables: you probably will only need to change things on the following few lines.
maxjobs = 2 #If we just do all the jobs at once, how many should we run at a time?
#Start and end values
#z (position) values
zta  = 4.0 #Start the z values where?
ztb  = 5.0 #End the z values where?
zinc = 1.0 #Increment by what value?
#eps (dielectric coefficient) values
eta  = 10
etb  = 11
einc = 1
#T (temperature) values
Tta  = 10
Ttb  = 11
Tinc = 1

szta = str(zta)
sztb = str(ztb)
szinc = str(zinc)
mta = 'pppta' #String to replace in template
mtb = 'ppptb' #For range of z values
minc = 'iinncc' #For z increment
meps = 'seeppss' #Some marker for the epsilon variable
mTmp = 'sTTmmpp' #Some marker for the temperature variable
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
eols = "\n"
tm = 3.0 #Seconds between issuing qstat

def mkarr(pta,ptb,inc):
    #Given range and increment, make array of values
    n = int(round(1+((ptb-pta)/inc)))
    arr = []
    for i in range(n):
        s_i = str(pta +i*inc)
        arr.append(s_i)
    return arr

def jobwait_longer(qmk_r,maxjobs):
    #Waits until job slot is clear after two checks in a row
    jobs_running = maxjobs
    lowjobct = 0
    while (lowjobct<2):
        time.sleep(tm)
        os.system("qstat | grep lg > "+tmpfl)
        f = open(tmpfl,'r')
        f_s = f.read()
        f.close()
        jobs_running = len(re.findall(qmk_r,f_s))
        if (jobs_running<maxjobs):
            lowjobct += 1
        else:
            lowjobct = 0
        os.system("rm "+tmpfl)

#Arrays of variables to replace
earr = mkarr(eta,etb,einc)
Tarr = mkarr(Tta,Ttb,Tinc)
for Ts in Tarr:
    for es in earr:
        if (len(Ts)+len(es)>6):
            print "Sorry, the sum of the lengths of temparature and epsilon strings cannot exceed 6."
            sys.exit()

os.chdir(fileloc)
dir = raw_input("Directory name:\n")
infilenamebase_template = raw_input("Base input file name:\n")
Mmkrn = input("Just make files (1) or make and run (2)?\n")
#Stage 1: make directories sorted by T. Into each place a copy of the input
# file template. Then, for each epsilon (within each directory) make a new
# version of the script (for each T, eps).
infilenameend = ".i"
outfilenameend = ".o"
scriptname_template = "FRG"
scriptnameend = ".py"
template_scrp = open(scriptname_template+scriptnameend,'r')
template_scrp_s = template_scrp.read()
os.system("mkdir "+dir)
os.chdir(dir)
dirprfx = "T="
for Ts in Tarr:
    Tdir = dirprfx + Ts
    os.system("mkdir "+Tdir)
    os.chdir(Tdir)
    for es in earr:
        inname_i = "T"+Ts+"e"+es+infilenameend
        Sinname_i = sS+"T"+Ts+"e"+es+infilenameend
        scriptname_i = scriptname_template+"T"+Ts+"e"+es+scriptnameend
        os.system("cp ../../"+infilenamebase_template+infilenameend+ss+inname_i)
        os.system("cp ../../"+sS+infilenamebase_template+infilenameend+ss+Sinname_i)
        f_i = open(scriptname_i,'w')
        f_i.write(template_scrp_s.replace(mta,szta).replace(mtb,sztb).replace(minc,szinc).replace(meps,es).replace(mTmp,Ts))
        f_i.close()
    os.chdir("../")
template_scrp.close()

if Mmkrn==2:
    #...then we want to run these files, too. Since each job will depend on the
    # last for its scf_guess input, we want to make sure that the last one
    # has finished before we start the next.
    cmdname = "cmd.tmp"
    cmdend = "2\n1\n1\n1\n"
    qmk_r = re.compile("T[1-9][0-9]{0,3}e[1-9][0-9]{0,2}"+outfilenameend+".* [ERQ] batch") #What to look for in qstat to indicate a running job
    for Ts in Tarr:
        print "Working on T = "+Ts
        Tdir = dirprfx + Ts
        os.chdir(Tdir)
        for es in earr:
            print "Working on eps = " +es
            dirname = "eps="+es
            inname_i = "T"+Ts+"e"+es
            scriptname_i = scriptname_template+"T"+Ts+"e"+es+scriptnameend
            f_i = open(cmdname,'w')
            f_i.write(dirname+eols+inname_i+eols+cmdend)
            f_i.seek(0,0)
            jobs_running = maxjobs+1 #When this drops below maxjobs, submit another job
            jobwait_longer(qmk_r,maxjobs)
            command = "nohup python " +scriptname_i + " < " +cmdname + " &"
            f_i.close()
            print "Executing command: "+command
            os.system(command)
            os.system("rm "+cmdname)
        os.chdir("../")
    jobwait_longer(qmk_r,1)
    print "This script has reached its end. All jobs have been completed."
