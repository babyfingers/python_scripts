#!/usr/local/bin/python

# delbunch.py
import os.path

# delbunch.py: kill a bunch of jobs all at once.

jobid_0 = 593600
jobid_1 = 593624

for jjj in range(jobid_0, jobid_1+1):
    os.system("qdel "+str(jjj))


