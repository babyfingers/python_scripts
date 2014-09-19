# autorun.py
#!/usr/local/bin/python
import string
import os.path
import time
import re


#Given a list of geometries separated by spaces, will fill a template file with each of the
# geometries and run them in turn.

#Other variables
eols = "\n"
tabs = "\t"
spas = " "
uns = "_"
ws = "w"
ss = "s"
rs = "r"
ls = "l"
fn_m = "fffnnn"
geo_m = "ggeeoo"
scf_m = "sssggg"
block_re = re.compile("xxx.*?xxx",re.S)
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]*\.[0-9]*e?-?[0-9]*\s+-?[0-9]*\.[0-9]*e?-?[0-9]*\s+-?[0-9]*\.[0-9]*e?-?[0-9]*")
ele_re = re.compile(eols+"[a-zA-Z]{1,2}")
num_re = re.compile("-?[0-9]*\.[0-9]*e?-?[0-9]*")
infilenameend = '.i'
outfilenameend = '.o'
tmp_l = []
x = []
n = 10
for i in range(n):
    x.append(i/(n-1.0))
wid_i = max(2,len(str(n)))
maxj = 1
tmpfl = "tmp.txt"
tm = 5.0
stm = 0.3

#jobwait: takes string and reg. ex. as
# input to stop the program until a
# certain specified job has completed.
def jobwait(waitmsg,qmk_r):
    print waitmsg+" Waiting "+str(stm)+" seconds."
    time.sleep(stm)
    os.system("qstat > "+tmpfl)
    f = open(tmpfl,rs)
    f_s = f.read()
    f.close()
    jobs_running = len(re.findall(qmk_r,f_s))
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
filenamebase = raw_input("Base input file name:\n")
dir = raw_input("Directory name:\n")
mode1 = input("Just make files (1) or make and run (2)?\n")
if (mode1==2):
    mode2 = input("Don't save plots (1) or save plots (2)?\n")
else:
    mode2 = 1
f = open(filenamebase+infilenameend,rs)
template = f.read()
f.close()

#Obtain and consolidate geometry information
f = open(filenamebase +".x",rs)
f_s = f.read()
f.close()
f_s = eols +f_s
block_l = re.findall(block_re,f_s)
ng = len(block_l)
for i in range(ng):
    block_l[i] = block_l[i].replace("xxx\n","")
    block_l[i] = block_l[i].replace("\nxxx","")

#Make sure the geometries match up in terms of # of atoms and order
#Check that the geometries correspond
ele_l = re.findall(ele_re,block_l[0])
ele_l2 = []
na = len(ele_l)
for i in range(1,ng):
    ele_l2 = re.findall(ele_re,block_l[i])
    for e in range(na):
        if (ele_l[e] != ele_l2[e]):
            print "Item ",e,": "+ele_l[e]+" and "+ele_l2[e]
            sys.exit("Elements of two geometries do not correspond.")

#Make input files.
infn_l = []
outfn_l = []
os.system("mkdir "+dir)
os.chdir(dir)
for i in range(ng):
    if (i==0): scf_s = "gwh"
    else: scf_s = "read"
    infn_l.append(ls+uns+str(i).zfill(wid_i)+infilenameend)
    outfn_l.append(ls+uns+str(i).zfill(wid_i)+outfilenameend)
    f = open(infn_l[i],ws)
    f.write(template.replace(geo_m,block_l[i]).replace(scf_m,scf_s).replace(fn_m,"x = "+str(i)+"/"+str(n-1)))
    f.close()

#Run files
qintro = "myqchem.csh -in "
qmezzo = " -out "
qoutro = []
if (mode2==1):
    qoutro.append(" -outdir tmp_dir -save standard")
    for i in range(1,ng):
        qoutro.append(" -indir scf_dir -save standard -outdir tmp_dir")
elif (mode2==2):
    qoutro.append(" -outdir "+ls+uns+str(0).zfill(wid_i)+" -save plots")
    for i in range(1,ng):
        qoutro.append(" -indir "+ls+uns+str(i-1).zfill(wid_i)+" -outdir "+ls+uns+str(i).zfill(wid_i)+" -save plots")
if (mode1==2):
    qmk_r = re.compile(ls+uns+"[0-9]+"+outfilenameend + ".* [ERQ] batch")
    docmd(qintro +infn_l[0] +qmezzo +outfn_l[0] +qoutro[0])
    jobwait("",qmk_r)
    for i in range(1,n):
        docmd("rm -r scf_dir")
        docmd("mv tmp_dir scf_dir")
        docmd(qintro +infn_l[i] +qmezzo +outfn_l[i] +qoutro[i])
        jobwait("",qmk_r)
