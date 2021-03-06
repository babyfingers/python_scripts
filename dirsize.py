#!/usr/local/bin/python2.7

# Use "du" to determine the size on disk of a given directory.

import re
import os
import sys
import time

def readfile(filename):
    f = open(filename,"r")
    fs = f.read()
    f.close()
    return fs

tmpfile = "dirsize.tmp"
duSize_re = re.compile("([0-9]+)\W+total")

def getSizeOnDisk(path):
    #Return size (on disk) of the item at path in bytes
    os.system("du -sc {} > {}".format(path,tmpfile))
    fs = readfile(tmpfile)
    os.remove(tmpfile)
    l = duSize_re.findall(fs)
    return int(l[0])*1024

#target = raw_input("What directory would you like the size of?\n")
#print "size = {} bytes".format(getSizeOnDisk(target))

def getPaths(pathfile,user):
    pathList = []
    f = open(pathfile,"r")
    for line in f:
        #Comment mark, ignore everything after on this line
        hashIndex = line.find("#")
        if (hashIndex!=-1):
            line = line[:hashIndex]
        line = line.strip()
        if (len(line)==0):
            continue
        #Check if the path is in this user's directory
        if (line.find("/data/home/{}/".format(user))!=0):
            print "*************************************"
            print "Problem with path = {}".format(line)
            print "This does not appear to be an absolute path pointing"
            print "to a location in your home directory. Skipping."
            print "*************************************"
            continue
        #Check if the path is legit
        if (not os.path.exists(line)):
            print "*************************************"
            print "Problem with path = {}".format(line)
            print "I can't find it. Skipping."
            print "*************************************"
            continue
        #Add to list
        pathList.append(line)
    f.close()
    return pathList

pathfile = "/data/home/alguire/pathfile.txt"
paths = getPaths(pathfile,"alguire")
total = 0
for path in paths:
    pathSize = getSizeOnDisk(path)
    total += pathSize
    print "Path {} has size {}.".format(path,pathSize)
print "Total size of paths = {}".format(total)

