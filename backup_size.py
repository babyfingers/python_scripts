#!/usr/local/bin/python2.7

# backup_size.py
# Given a path list file of the format used by offsite_backup.py,
# ("backup_please"), tell the user how large their requested
# backup will be (in bytes).

vers = 0.1

import re
import os
import sys
from datetime import datetime as dt

# Define some functions
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
            print "    ++ERROR++"
            print "Problem with path = {}".format(line)
            print "This does not appear to be an absolute path pointing"
            print "to a location in your home directory. Skipping."
            print "    +++++++++"
            continue
        #Check if the path is legit
        if (not os.path.exists(line)):
            print "    ++ERROR++"
            print "Problem with path = {}".format(line)
            print "I can't find it. Skipping."
            print "    +++++++++"
            continue
        #Add to list
        pathList.append(line)
    f.close()
    return pathList

def openFile(filePath,willWrite):
    #Check if there is a log file at the chosen path.
    #If there is a log file, append at the end. If
    #there is not, create one. Otherwise, throw an
    #error.
    file_exists = os.path.exists(filePath)
    file_is_file = os.path.isfile(filePath)
    if (file_exists and file_is_file):
        if (willWrite):
            outFile = open(filePath,"a")
            outFile.seek(0,2)
        else:
            outFile = open(filePath,"r")
    elif(file_exists and not file_is_file):
        print "Error! Something exists at {} and it is not a file.".format(filePath)
        return "NOT_A_FILE"
    elif(not file_exists):
        if (willWrite):
            outFile = open(filePath,"w")
        else:
            return "NOT_A_FILE"
    return outFile

def getBackupSize(bupFilePath,username):
    print "Contents of {}:".format(bupFilePath)
    print "--------------------------------------------------------"
    print readfile(bupFilePath)
    print "--------------------------------------------------------"
    paths = getPaths(bupFilePath,username)
    mem_total = 0
    for path in paths:
        pathSize = getSizeOnDisk(path)
        mem_total += pathSize
        print "Path {} has size {}".format(path,pathSize)
    print "Total size of paths = {} bytes = {} GB".format(mem_total,mem_total/1000000000.0)

orig_path = os.getcwd()
username = os.getlogin()
reqFilePath = "/data/home/{}/backup_please".format(username)
bupFile = openFile(reqFilePath,False)
if (bupFile == "NOT_A_FILE"):
    print "Cannot find backup request file at \n\t{}".format(reqFilePath)
    print "Cannot proceed with backup if request file cannot be read."
    result_str = "aborted"
else:
    result_str = ""
    bupFile.close()
  
if (result_str != "aborted"):
    getBackupSize(reqFilePath,username)
    bupFile.close()

