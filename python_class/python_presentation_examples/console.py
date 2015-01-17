#!/usr/local/bin/python
import os

os.system('qstat | grep alguire > temp.txt')
f = open('temp.txt','r')
fs = f.read()
f.close()

print fs

if "R batch" in fs:
    print "Job has not yet completed."
else:
    print "Job has completed."

