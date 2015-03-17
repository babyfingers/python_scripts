#!/usr/local/bin/python2.7

# offsite_backup.py
# Go into each user's home directory and look for a file called
# backup_please. It should contain a list of directories. After
# checking that each directory is legit and the total size does
# not exceed some designated amount.

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
script_path = os.getcwd()

def getSizeOnDisk(path):
    #Return size (on disk) of the item at path in bytes
    os.system("du -sc {} > {}".format(path,tmpfile))
    fs = readfile(tmpfile)
    os.remove(tmpfile)
    l = duSize_re.findall(fs)
    return int(l[0])*1024

#target = raw_input("What directory would you like the size of?\n")
#print "size = {} bytes".format(getSizeOnDisk(target))

def writeList(write_string, write_to):
    has_error = False
    in_type = type(write_to)
    if (in_type is file):
        write_to.write(write_string)
    elif (write_to=="PRINT"):
        print write_string
    elif (in_type is tuple or in_type is list):
        for wt in write_to:
            if (type(wt) is file):
                wt.write(write_string)
            elif(wt=="PRINT"):
                print write_string
            else:
                has_error = True
                break
    else:
        has_error = True
    if (has_error):
        print "The 'write_to' input variable is unacceptable."
        print "It must either be a file, the string \"PRINT\", or a"
        print "list of files or the string \"PRINT\"."
        assert(False)

def getPaths(pathfile,user,write_list):
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
            writeList("    ++ERROR++\n",write_list)
            writeList("Problem with path = {}\n".format(line),write_list)
            writeList("This does not appear to be an absolute path pointing\n",write_list)
            writeList("to a location in your home directory. Skipping.\n",write_list)
            writeList("    +++++++++\n",write_list)
            continue
        #Check if the path is legit
        if (not os.path.exists(line)):
            writeList("    ++ERROR++\n",write_list)
            writeList("Problem with path = {}\n".format(line),write_list)
            writeList("I can't find it. Skipping.\n",write_list)
            writeList("    +++++++++\n",write_list)
            continue
        #Add to list
        pathList.append(line)
    f.close()
    return pathList

def backupPaths(paths,username):
    os.chdir("/data/home")
    for path in paths:
        path = path.replace("/data/home/","")
        #Given a path that you want to back up [...]
        #command = "rsync -avR --log-file='/data/home/amber/backup/rsync.log' {} backup@romulus:'/Volumes/My\ Passport\ for\ Mac/amber'".format(myPath)
        #command = "rsync -avR --log-file='/data/home/amber/backup/rsync.log' {} backup@romulus:''".format(path)
        #command = "rsync -avR {} /data/home/alguire/backup_test/{}".format(path,username)
        command = "rsync -avR --log-file='{}' {} backup@romulus:'diabatical'".format(os.path.join(script_path,"rsync.log"),path)
        print command
        os.system(command)
    os.chdir(username)

def openFile(filePath,write_out,willWrite):
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
        writeList("Error! Something exists at {} and it is not a file.\n".format(filePath),write_out)
        return "NOT_A_FILE"
    elif(not file_exists):
        if (willWrite):
            outFile = open(filePath,"w")
        else:
            return "NOT_A_FILE"
    return outFile

# Set up log files
os.chdir("/data/home/")

tm0 = dt.now()
mainLogPath = os.path.join(script_path,"main_backup_log_{}-{}-{}.txt".format(tm0.year,tm0.month,tm0.day))

mainLogFile = openFile(mainLogPath,"PRINT",True)
if (mainLogFile=="NOT_A_FILE"):
    print "Error! We can't write the log file to {}"
    print "because there is some non-file in the way."
    print "This must be fixed before backup can proceed."
    sys.exit()

writeList("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",mainLogFile)
writeList("~~~~~~~ STARTING BACKUP LOOP, v{} {:02d}-{:02d} {:02d}:{:02d}:{:02d} ~~~~~~~~\n".format(vers,tm0.year,tm0.month,tm0.day,tm0.hour,tm0.minute,tm0.second),mainLogFile)

def do_backup(bupFilePath,write_list,username,mem_max):
    writeList("Maximum backup size for user {} is {} bytes.\n".format(username,mem_max),write_list)
    writeList("Contents of {}:\n".format(bupFilePath),write_list)
    writeList("--------------------------------------------------------\n",write_list)
    writeList(readfile(bupFilePath),write_list)
    writeList("--------------------------------------------------------\n",write_list)
    paths = getPaths(bupFilePath,username,write_list)
    mem_total = 0
    for path in paths:
        pathSize = getSizeOnDisk(path)
        mem_total += pathSize
        writeList("Path {} has size {}.\n".format(path,pathSize),write_list)
    writeList("Total size of paths = {} bytes\n".format(mem_total),write_list)
    if (mem_total>mem_max):
        writeList("Sorry, your requested backups are too large (max size = {} bytes).\n".format(mem_max),write_list)
        result_str = "aborted"
    else:
        backupPaths(paths,username)
        result_str = "complete"
    return result_str

# (start loop over users)
user_list = [u for u in os.listdir("/data/home/") if os.path.isdir(u)]
#user_list = ["alguire"]
user_limits = {}
for u in user_list:
    user_limits[u] = 10000000000 #10 GB

for username in user_list:
    mem_max = user_limits[username]
    user_home = "/data/home/" +username
    bupFilePath = user_home+"/backup_please"
    backup_success = False
    
    write_list = mainLogFile
    writeList("\n\n",write_list)
    writeList("Now working on user {}\n".format(username),write_list)

    os.chdir(user_home)

    #Set up local log file
    locLogPath = user_home +"/offsite_backup.log"
    locLogFile = openFile(locLogPath,write_list,True)
    t0 = dt.now()
    if (locLogFile == "NOT_A_FILE"):
        writeList("********************************************************\n",write_list)
        writeList("******* offsite_backup v{}, {}-{:02d}-{:02d} {:02d}:{:02d}:{:02d} *******\n".format(vers,t0.year,t0.month,t0.day,t0.hour,t0.minute,t0.second),write_list)
        writeList("Something other than a text file is where the local log file should be.\n",write_list)
        writeList("\t"+locLogPath,write_list)
        writeList("Please move this file if you would like a local log of the backup.\n",write_list)
    else:
        #write_list = [mainLogFile,locLogFile]
        #writeList("\n\n",locLogFile)
        writeList("********************************************************\n",write_list)
        writeList("******* offsite_backup v{}, {}-{:02d}-{:02d} {:02d}:{:02d}:{:02d} *******\n".format(vers,t0.year,t0.month,t0.day,t0.hour,t0.minute,t0.second),write_list)
    os.system("chmod o+r "+locLogPath)

    #Look for the file that will tell us what to back up
    bupFile = openFile(bupFilePath,write_list,False)
    if (bupFile == "NOT_A_FILE"):
        writeList("Cannot find backup request file at \n\t{}\n".format(bupFilePath),write_list)
        writeList("Cannot proceed with backup if request file cannot be read.\n",write_list)
        result_str = "aborted"
    else:
        result_str = ""
        bupFile.close()
        
    if (result_str != "aborted"):
        result_str = do_backup(bupFilePath,write_list,username,mem_max)
    t1 = dt.now()
    writeList("******* backup {} on {}-{:02d}-{:02d} {:02d}:{:02d}:{:02d} *******\n".format(result_str,t1.year,t1.month,t1.day,t1.hour,t1.minute,t1.second),write_list)
    t10 = t1-t0
    writeList("Elapsed time: {} seconds.\n".format(t10.total_seconds()),write_list)
    writeList("********************************************************\n",write_list)

    if (locLogFile != "NOT_A_FILE"):
        locLogFile.close()
        write_list = mainLogFile
tm1 = dt.now()
tm10 = tm1-tm0
writeList("~~~~~~~~~ ENDING BACKUP LOOP, {}-{:02d}-{:02d} {:02d}:{:02d}:{:02d} ~~~~~~~~\n".format(tm1.year,tm1.month,tm1.day,tm1.hour,tm1.minute,tm1.second),mainLogFile)
writeList("Elapsed time: {} seconds.\n".format(tm10.total_seconds()),write_list)
writeList("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n",mainLogFile)

mainLogFile.close()

