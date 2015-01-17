#!/usr/local/bin/python
import os
import re

os.system('qstat | grep alguire > temp.txt')
f = open('temp.txt','r')
fs = f.read()
f.close()

print fs

rgx = re.compile("[RQ] batch")
num_matches = len(rgx.findall(fs))

if num_matches > 0:
    print "Job has not yet completed."
else:
    print "Job has completed."

