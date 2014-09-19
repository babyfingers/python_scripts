#vectorfun.py

# Module containing functions I've defined to process
# gradient quantities from Q-Chem outputs.
# Created 2013-10-17

import string
import os.path
import numpy
import math
import time
import sys
import re

#Some helpful regular expressions
eols = "\n"
atom_re = re.compile(eols+"[a-zA-Z]{1,2}\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*\s+-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
ele_re = re.compile("[a-zA-Z]{1,2}")
num_re = re.compile("-?[0-9]+\.?[0-9]*[eEdD]?[-+]?[0-9]*")
int_re = re.compile("\s-?[0-9]+")
dp_re   = re.compile("-?[0-9]\.[0-9]{14}e[+-][0-9]{2,3}") # Reg ex for the scientific notation for these outputs 

#magnitude: calculates the 2-norm of a
# numpy array v.
def magnitude(v):
    return math.sqrt(numpy.inner(v,v))

#normalize: normalizes a numpy array v.
def normalize(v):
    if (magnitude(v) > 0):
        v = v/magnitude(v)
    return v

#doterror: finds the normalized overlap
# between two numpy vectors v1 and v2;
# takes the dot product of v1 and v2,
# then divides by the dot product of
# v2 with itself.
def doterror(v1, v2):
    return numpy.inner(v1,v2)/numpy.inner(v2,v2) -1

#direrror: finds the directional overlap
# between two numpy vectors v1 and v2.
def direrror(v1, v2):
    u1 = v1.copy()
    u2 = v2.copy()
    return 1 -numpy.inner(normalize(u1),normalize(u2))

#RMSD: returns the root mean squared
# deviation between two numpy vectors
# v and u.
def RMSD(v1, v2):
    d = v1-v2
    d = d**2
    return math.sqrt(numpy.average(d))

#docmd: prints a given string to screen
# and executes it as a command.
def docmd(s):
    #print ">"+s
    os.system(s)

#lencheck: checks that len(mylist)==correct_len,
# otherwise returns 0 and prints a message.
def lencheck(mylist, correct_len, filename_s):
    if (len(mylist) != correct_len):
        print "Error: not seeing the right number of numbers in grepped DC for file name "+filename_s
        print "Skipping this file."
        return 0
    else:
        return 1
    
#grepfls: takes a file name, a string to grep,
# and returns the grepped output.
def grepfls(filename_s, grepmarker_s):
    grepfln = "greptmp.txt"
    docmd("grep " +grepmarker_s +" " +filename_s +" > " +grepfln)
    grepf = open(grepfln,rs)
    grepf_s = grepf.read()
    grepf.close()
    docmd("rm " +grepfln)
    return grepf_s


