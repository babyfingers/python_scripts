#!/usr/local/bin/python2.7
# makintern.py
import string
import os.path
import math
import time
import sys
import re
import numpy


# makintern.py: takes two geometries of the same molecule and applies
# translation and rotation operations iteratively until they are as
# close to each other as possible, as defined by the 2-norm.

rot_conv = 0.000000001 #Tolerance for converging one rotational degree of freedom within goldenss
Rdist_conv = 0.00000001 #Tolerance for converging iterations of rotational minimization
dist_conv = 0.00000001 #Tolerance for converging iterations of translational + rotational minimization
max_rot = 200 #Maximum number of rotational iterations per run
max_run = 200 #Maximum number of total iterations
tPI = 2.0*math.pi
npoints = 20 #Number of points used in the initial 'brute force' minimization of rotation angles
dt = tPI/npoints #Increment for above, in radians
golR = 0.6180339887498950
golC = 1.0 -golR
tabs = "\t"
eols = "\n"
tia_s = "1"
rs = "r"
ws = "w"
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")

#docmd: prints a given string to screen
# and executes it as a command.
def docmd(s):
    #print ">"+s
    os.system(s)

#distance: calculates the 2-norm of the
# difference between two similarly-sized arrays.
def distance(v1, v2):
    v12 = v1-v2
    v12 = numpy.append(v12[:,0],v12[:,1:3])
    return math.sqrt(numpy.inner(v12,v12))

#makeRx: returns a rotation matrix for rotating
# about the x (0), y (1) or z (2) axes by
# theta radians.
def makeRx(theta,rotDirection):
    Rx = numpy.zeros((3,3))
    if (rotDirection==0):
        Rx[0][0] = 1.0
        Rx[1][1] = math.cos(theta)
        Rx[2][2] = math.cos(theta)
        Rx[1][2] = math.sin(theta)
        Rx[2][1] = -math.sin(theta)
    elif (rotDirection==1):
        Rx[1][1] = 1.0
        Rx[0][0] = math.cos(theta)
        Rx[2][2] = math.cos(theta)
        Rx[2][0] = math.sin(theta)
        Rx[0][2] = -math.sin(theta)
    elif (rotDirection==2):
        Rx[2][2] = 1.0
        Rx[0][0] = math.cos(theta)
        Rx[1][1] = math.cos(theta)
        Rx[0][1] = math.sin(theta)
        Rx[1][0] = -math.sin(theta)
    return Rx

def minbrac(g0_f, g1_f, kkk):
    dist_l = []
    #Brute force: look at distances over all possible ranges of theta
    for th_i in range(npoints):
        th_f = th_i*dt
        Rx = makeRx(th_f,kkk)
        rot_f = numpy.dot(cartVector1,Rx) #Rotate the 2nd endpoint
        dist_l.append(distance(g0_f,rot_f)) #Keep track of distances
    return dist_l.index(min(dist_l))*dt

def shft3(l, iii, jjj, kkk, x):
    l[iii] = l[jjj]
    l[jjj] = l[kkk]
    l[kkk] = x

def shft2(l, iii, jjj, x):
    l[iii] = l[jjj]
    l[jjj] = x

def shftarr3(l, iii, jjj, kkk, x):
    l[iii] = l[jjj].copy()
    l[jjj] = l[kkk].copy()
    l[kkk] = x.copy()

def shftarr2(l, iii, jjj, x):
    l[iii] = l[jjj].copy()
    l[jjj] = x.copy()
        
def golupdate(t_l, R_l, g_l, d_l, g0_f, g1_f, iii, kkk):
    R_l[iii] = makeRx(t_l[iii],kkk)     #Rotation matrix in directi
    g_l[iii] = numpy.dot(g1_f,R_l[iii]) #New, modified geometry (g1 rotated by R)
    d_l[iii] = distance(g0_f, g_l[iii])

def goldenss(t_l, g0_f, g1_f, kkk, rot_conv):
    if (math.fabs(t_l[3]-t_l[1])>math.fabs(t_l[1]-t_l[0])):
        t_l[2] = t_l[1] +golC*(t_l[3]-t_l[1])
    else:
        t_l[2] = t_l[1]
        t_l[1] -= golC*(t_l[1]-t_l[0])
    R_l = [0.0, 0.0, 0.0, 0.0] #Will store rotation matrices (zeroes to hold place for now)
    g_l = [0.0, 0.0, 0.0, 0.0] #Will store geometry matrices (zeroes to hold place for now)
    d_l = [0.0, 0.0, 0.0, 0.0] #Will store distances
    for iii in range(4):
        golupdate(t_l, R_l, g_l, d_l, g0_f, g1_f, iii, kkk)
    while(math.fabs(t_l[3]-t_l[0]) > rot_conv*(math.fabs(t_l[1])+math.fabs(t_l[2]))):
        if (d_l[2]<d_l[1]):
            shft3(t_l,0,1,2,golC*t_l[2] +golR*t_l[3])
            shftarr3(R_l,0,1,2,numpy.zeros((1,1)))
            shftarr3(g_l,0,1,2,numpy.zeros((1,1)))
            shft3(d_l,0,1,2,0.0)
            golupdate(t_l, R_l, g_l, d_l, g0_f, g1_f, 2, kkk)
        else:
            shft3(t_l,3,2,1,golC*t_l[1] +golR*t_l[0])
            shftarr3(R_l,3,2,1,numpy.zeros((1,1)))
            shftarr3(g_l,3,2,1,numpy.zeros((1,1)))
            shft3(d_l,3,2,1,0.0)
            golupdate(t_l, R_l, g_l, d_l, g0_f, g1_f, 1, kkk)
    return g_l[1].copy()

def printGeometry(elementList,cartVector):
    printString = ""
    NAtom = len(elementList)
    for iii in range(NAtom):
        printString += elementList[iii]
        for mu in range(3):
            printString += " "+str(cartVector[iii,mu])
        printString += "\n"
    print printString

# Get some user input
inputFileName = raw_input("Input file name:\n")

#Check out geometry information
f = open(inputFileName,rs)
f_s = f.read()
f.close()
f_s = eols +f_s
atom_l = re.findall(atom_re,f_s)
tNAtom = len(atom_l)
NAtom = tNAtom/2

#Store the geometry vectors
elementList = []
geo0_v = [[],[],[]]
geo1_v = [[],[],[]]

for iAtom in range(NAtom):
    str0 = atom_l[iAtom]
    str1 = atom_l[iAtom +NAtom]
    words0 = str0.split()
    words1 = str1.split()
    element0 = words0[0]
    element1 = words1[0]
    assert(element0==element1)
    elementList.append(element0)
    for mu in range(3):
        geo0_v[mu].append(float(words0[mu+1]))
        geo1_v[mu].append(float(words1[mu+1]))
elementList = tuple(elementList)
cartVector0 = (numpy.array(geo0_v)).transpose()
cartVector1 = (numpy.array(geo1_v)).transpose()

dist_init = distance(cartVector0,cartVector1)
print "Initial distance = "+str(dist_init)+" Angstroms.\n"

i = 0
converged = False
dist_prev = 0.0
dist_now = distance(cartVector0,cartVector1)
conv_count = 0 #Number of times the convergence requirement has been reached
Rdist_prev = 0.0
Rdist_now = dist_now
Rconv_count = 0 #Number of times the convergence requirement has been reached
Tdist_prev = 0.0
Tdist_now = dist_now
while (not converged):
    #First, minimize translational distance
    #  (translate by the negative of the average values
    #   along each Cartesian coordinate.)
    geodiff_f = cartVector0-cartVector1
    avg_f = numpy.zeros((3,1))
    for ccc in range(3):
        avg_f[ccc] = numpy.average(geodiff_f[:,ccc])
    temp_f = numpy.zeros((3,3))
    temp_f[0][0] = avg_f[0]
    temp_f[1][1] = avg_f[1]
    temp_f[2][2] = avg_f[2]
    trans_f = numpy.ones((NAtom,3))
    trans_f = numpy.dot(trans_f,temp_f)
    cartVector1 += trans_f
    Tdist_prev = Tdist_now
    Tdist_now = distance(cartVector0,cartVector1)
    Tdist_diff = Tdist_now -Tdist_prev
    #print "Tdist: "+str(Tdist_prev)+" -> "+str(Tdist_now)+". Difference = "+str(Tdist_diff)+"."
    #Second, minimize rotational distance
    Rconverged = False
    j = 0
    while (not Rconverged):
        #print "Rotational optimization: i = "+str(i)+", j = "+str(j)
        #Optimize rotational angle over x, y and z axes
        l_t = [0.0, 0.0, 0.0, 0.0]
        for kkk in range(3):
            l_t[1] = minbrac(cartVector0, cartVector1, kkk)
            l_t[0] = l_t[1] -dt
            l_t[3] = l_t[1] +dt
            #Perform minimization within bounds by using a golden section search
            cartVector1 = goldenss(l_t, cartVector0, cartVector1, kkk, rot_conv)
        Rdist_prev = Rdist_now
        Rdist_now = distance(cartVector0,cartVector1)
        Rdist_diff = Rdist_now -Rdist_prev
        #print "Rdist: "+str(Rdist_prev)+" -> "+str(Rdist_now)+". Difference = "+str(Rdist_diff)+"."
        j += 1
        if (math.fabs(Rdist_diff)<Rdist_conv):
            Rconv_count += 1
        else:
            Rconv_count = 0
        if (j>max_rot or Rconv_count>1): Rconverged = True
        #if (Rdist_diff<Rdist_conv): Rconverged = True
    #Evaluate convergence criterion
    dist_prev = dist_now
    dist_now = distance(cartVector0,cartVector1)
    dist_diff = dist_now -dist_prev
    print "dist: "+str(dist_prev)+" -> "+str(dist_now)+". Difference = "+str(dist_diff)+"."
    i += 1
    if (math.fabs(dist_diff)<dist_conv):
        conv_count += 1
    else:
        conv_count = 0
    if (i>max_run or conv_count>1): converged = True
    #Print distance
    print "Distance on iteration "+str(i)+": "+str(dist_now)+" Angstroms.\n"
#Output
print "First geometry- unchanged:"
printGeometry(elementList,cartVector0)
print "Second geometry- after rotation and translation:"
printGeometry(elementList,cartVector1)
print "Initial distance apart = "+str(dist_init)+" Angstroms, final distance apart = "+str(dist_now)+" Angstroms."
Atoa0 = 1.889725989
print "Initial distance apart = "+str(dist_init*Atoa0)+" a.u.,      final distance apart = "+str(dist_now*Atoa0)+" a.u."
