# linterp.py
#!/usr/local/bin/python
import string
import os.path
import time
import re


#Given two geometries and some other information, this program will create and sumbit
# a number of input files that represent geometries linearly interpolated between
# the two initial geometries.

#Number of points to run along reaction coordinate.
n = 11

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
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]*\.[0-9]*\s+-?[0-9]*\.[0-9]*\s+-?[0-9]*\.[0-9]*")
ele_re = re.compile(eols+"[a-zA-Z]{1,2}")
num_re = re.compile("-?[0-9]*\.[0-9]*")
infilenameend = '.i'
outfilenameend = '.o'
tmp_l = []
x = []
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
f = open(filenamebase+infilenameend,rs)
template = f.read()
f.close()

#Check out geometry information
f = open(filenamebase +".x",rs)
f_s = f.read()
f.close()
f_s = eols + f_s
atom_l = re.findall(atom_re,f_s)
#Make sure the two geometries match up in terms of # of atoms and order
tna = len(atom_l)
na = tna/2
#Check that the total number of atoms is even
if (tna%2 == 1):
    print "tna = "+str(tna)
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
        print "Item ",i,": "+ele_l[i]+" and "+ele_l2[i]
        sys.exit("Elements of two geometries do not correspond.")

#Make input files.
infn_l = []
outfn_l = []
os.system("mkdir "+dir)
os.chdir(dir)
for i in range(n):
    if (i==0): scf_s = "gwh"
    else: scf_s = "read"
    infn_l.append(ls+uns+str(i).zfill(wid_i)+infilenameend)
    outfn_l.append(ls+uns+str(i).zfill(wid_i)+outfilenameend)
    geo_s = ""
    for j in range(na):
        cart0_s = re.findall(num_re,atom_l[j])
        cart1_s = re.findall(num_re,atom_l[j+na])
        geo_s = geo_s +ele_l[j]
        for c in range(3):
            cartx_f = (1-x[i])*float(cart0_s[c]) +x[i]*float(cart1_s[c])
            geo_s = geo_s +tabs +str(cartx_f)
        if (j<na-1): geo_s = geo_s +eols
    f = open(infn_l[i],ws)
    f.write(template.replace(geo_m,geo_s).replace(scf_m,scf_s).replace(fn_m,"x = "+str(i)+"/"+str(n-1)))
    f.close()

#Run files
qintro = "myqchem.csh -in "
qmezzo = " -out "
qoutro = []
qoutro.append(" -outdir tmp_dir -save standard")
qoutro.append(" -indir scf_dir -save standard -outdir tmp_dir")
if (mode1==2):
    qmk_r = re.compile(ls+uns+"[0-9]+"+outfilenameend + ".* [ERQ] batch")
    docmd(qintro +infn_l[0] +qmezzo +outfn_l[0] +qoutro[0])
    jobwait("",qmk_r)
    for i in range(1,n):
        docmd("rm -r scf_dir")
        docmd("mv tmp_dir scf_dir")
        docmd(qintro +infn_l[i] +qmezzo +outfn_l[i] +qoutro[1])
        jobwait("",qmk_r)
