#!/usr/local/bin/python

#Useful functions for submitting and monitoring
#queue jobs.

import time
import sys
import re
import os
# A lot of these functions revolve around the use of the jobData
# array, which is produced from getJobData by looking at the
# output of "qstat -n1". The first dimension is the number of
# jobs of interest. The second describes the types of information
# about those jobs. Right now there are four types:
#jobData[i][0] : Jobname
#jobData[i][1] : Requested memory
#jobData[i][2] : Status
#jobData[i][3] : Node/processor list (E.g. "n5/0+n5/1")
# Maybe at some point it will make sense to organize this into
# a dictionary.

GJIfileName = "getJobInfo.txt"

def getJobData(userName,jobNameBase=""):
    '''Get information about a subset of jobs
    in the queue belonging to a certain user
    and with a specific job name.'''
    assert(len(jobNameBase)<16)
    os.system("qstat -n1 | grep \""+userName+".*"+jobNameBase+"\" > "+GJIfileName)
    qstatFile = open(GJIfileName,"r")
    jobData = []
    for line in qstatFile:
        columns = line.split()
        jobData.append((columns[3],columns[7],columns[9],columns[11]))
    qstatFile.close()
    os.system("rm "+GJIfileName)
    return tuple(jobData)

def incompleteJobs(jobData):
    '''Look through a jobData tuple and count
    the number of jobs with status Q or R.'''
    NIncomplete = 0
    for dataTup in jobData:
        if dataTup[2]=="Q" or dataTup[2]=="R":
            NIncomplete += 1
    return NIncomplete

def timeToSubmitJobs(maxJobs,userName,jobNameBase=""):
    '''Return True if the current number of submitted
    jobs is below the maximum desired.'''
    jobData = getJobData(userName,jobNameBase)
    if len(jobData)<maxJobs:
        return True
    else:
        return False

def submitCommand(submitStrings,inputName,outputName=""):
    '''Return a string to use as the command
    to run a job.'''
    NInputString = len(submitStrings)
    assert(NInputString>0 and NInputString<4)
    if (NInputString==1):
        command = submitStrings[0] +" "+inputName +" " +outputName
    elif (NInputString==2):
        command = submitStrings[0] +" "+inputName +submitStrings[1] +outputName
    elif (NInputString==3):
        command = submitStrings[0] +" "+inputName +submitStrings[1] +outputName +submitStrings[2]
    os.system(command)
    return command

def manageJobs(submitStrings,userName,jobList,maxJobs,waitTime,jobNameBase=""):
    '''As long as there are fewer jobs in the
    queue than maxJobs, keep submitted jobs
    from jobList.'''
    assert(len(jobList[0])==len(jobList[1]))
    totalJobs = len(jobList[0])
    currentJob = 0
    while (currentJob<totalJobs):
        if timeToSubmitJobs(maxJobs,userName,jobNameBase):
            inputName = jobList[0][currentJob]
            if (len(jobList)>1):
                outputName = jobList[1][currentJob]
                print "Running command:" +submitCommand(submitStrings,inputName,outputName)
            else:
                print "Running command:" +submitCommand(submitStrings,inputName)
            currentJob += 1
        print "Waiting "+str(waitTime)+" seconds."
        time.sleep(waitTime)
    print "\nAll jobs submitted.\n"

def templateInputFormatting(templateString,replacementList):
    '''Takes a template string containing markers
    for replacement, as well as a list of tuples
    that contain those markers as well as their
    replacements. This function should return
    a string describing an input file that is
    ready to run.'''
    NReplacement = len(replacementList)
    finalString = templateString
    for iii in range(NReplacement):
        finalString = finalString.replace(replacementList[iii][0],replacementList[iii][1])
    return finalString

