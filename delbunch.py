#!usr/local/bin/python2.7

# delbunch.py
import os.path

# delbunch.py: kill a bunch of jobs all at once.

jobid_0 = 568681
jobid_1 = 568776

for jjj in range(jobid_0, jobid_1+1):
    os.system("qdel "+str(jjj))

