#!/usr/local/bin/python
# delall.py
#import string
import os.path

os.system("pwd")
os.system("ls ../ > tmp1.txt")
f = open("tmp1.txt","r")
f_s = f.read()
f.close()
file_l = f_s.split()
os.system("rm tmp1.txt")
for s in file_l:
    os.system("svn delete "+s)
    #print "Hey, it's "+s
