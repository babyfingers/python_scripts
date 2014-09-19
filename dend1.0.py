# dend.py
#!/usr/local/bin/python
import os.path
import time
import re

#This program is a modified version of Snqchem.csh. Designed to address the problem
# proposed by Brad Veldkamp. If we have a DBA system, one method of figuring out which
# of a number of donor groups is superior is to figure out where the electron and hole
# density localizes after the donation. The distance between these will give us an idea
# of the overlap; shorter distances -> greater overlap -> faster electron transfer.
# I figure this can be obtained by calculating the electric dipole moments of the
# neutral, anionic and cationic states of the molecule (the extra electron will localize
# on the acceptor, the extra hole will localize on the donor) and then find the
# differences between the charged states and the neutral state dipoles (method 1). What
# remains will be the dipole moments of the lone electron, lone hole-- which can be easily
# translated into the centroid of each. Another way would be to look at the centroid of the
# HOMO of the anionic calculation vs. the centroid of the LUMO of the cationic state
# (method 2).

# Variables: you probably will only need to change things on the following few lines.
mch = 'ccc' #Marker string for charge
msm = 'sss' #Marker string for spin multiplicity
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
st  = '2'
su  = '_'
ss  = ' '
tm = 2.0 #seconds to wait before calling qstat

def jobwait(waitmsg,qmk_r):
    jobs_running = 1
    while jobs_running>0:
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
dir_full = fileloc+"/"+dir
filenamebase = raw_input("Base input file name:\n")
Mmkrn = input("Just make files (1), make and run (2), or just do analysis (3)?\n")
infl = []
oufl = []
ch = [sz,so,'-1']
sm = [so,st,st]
infl.append(filenamebase+"_nu"+infilenameend)
infl.append(filenamebase+"_ca"+infilenameend)
infl.append(filenamebase+"_an"+infilenameend)
oufl.append("nu_"+filenamebase+outfilenameend)
oufl.append("ca_"+filenamebase+outfilenameend)
oufl.append("an_"+filenamebase+outfilenameend)
qmk_r = re.compile(filenamebase+outfilenameend+".* [ERQ] batch") #What to look for in qstat to indicate a running job

if (Mmkrn==1 or Mmkrn==2):
    #Open up the template file and replace the relevant markers
    # with their respective variables to create each input file.
    template = open(filenamebase+infilenameend,'r')
    template_s = template.read()
    os.system("mkdir "+dir)
    os.chdir(dir)
    for i in range(3):
        f = open(infl[i],'w')
        f.write(template_s.replace(mch,ch[i]).replace(msm,sm[i]))
        f.close()
    template.close()

if Mmkrn==2:
    #Run the jobs.
    os.chdir(dir_full)
    for i in range(3):
        command = "myqchem.csh -in "+infl[i]+" -out "+oufl[i]
        os.system(command)
        time.sleep(2.0)
    jobwait("All jobs submitted.",qmk_r)

if (Mmkrn==2 or Mmkrn==3):
    os.chdir(dir_full)
    #Analysis: extract dipole moments, do relevant arithmetic.
    dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs
    os.system("grep xxDIPgxx *"+filenamebase+outfilenameend+" > dipg.txt")
    os.system("grep xxDIP_HOMOxx *"+filenamebase+outfilenameend+" > diphomo.txt")
    os.system("grep xxDIP_LUMOxx *"+filenamebase+outfilenameend+" > diplumo.txt")
    #Open up text file with all dipoles, extract them as list of 9 strings
    f = open("dipg.txt",'r');
    f_s = f.read()
    list1 = re.findall(dp_re,f_s)
    f.close()
    #Convert dipole values into sets of 3 floats
    anlist = []
    calist = []
    nulist = []
    for i in range(3):
        anlist.append(float(list1[i]))
        calist.append(float(list1[i+3]))
        nulist.append(float(list1[i+6]))
    #Convert dipoles into centroids of electron, hole
    h_v = []
    e_v = []
    distance = 0.0
    for i in range(3):
        h_v.append(calist[i]-nulist[i])
        e_v.append(nulist[i]-anlist[i])
        distance += (h_v[i]-e_v[i])**2
    distance = distance**0.5
    print "RESULTS: METHOD 1"
    print "Location of hole:",h_v[0],h_v[1],h_v[2]
    print "Location of electron:",e_v[0],e_v[1],e_v[2]
    print "Distance between two:",distance
    #Try a different method: use HOMO and LUMO dipoles
    f = open("diphomo.txt",'r')
    f_s = f.read()
    list1 = re.findall(dp_re,f_s)
    f.close()
    f = open("diplumo.txt",'r')
    f_s = f.read()
    list2 = re.findall(dp_re,f_s)
    f.close()
    h_v = []
    e_v = []
    distance = 0.0
    for i in range(3):
        e_v.append(-float(list1[i]))
        h_v.append(-float(list2[i+3]))
        distance += (h_v[i]-e_v[i])**2
    distance = distance**0.5
    print "RESULTS: METHOD 2"
    print "Location of hole:",h_v[0],h_v[1],h_v[2]
    print "Location of electron:",e_v[0],e_v[1],e_v[2]
    print "Distance between two:",distance
