# mser.py
#!/usr/local/bin/python
import os.path
import time
import re

#Just a little script to submit every input file in a directory
tmpfl = 'tempout.txt'
os.system("ls *.in > "+tmpfl) #List all files in the directory that end with ".in" and write them to a file called "tempout.txt"
f = open(tmpfl,'r')
file_s = f.read() #Copy contents of file to a string variable file_s
f.close()
os.system("rm "+tmpfl)
print file_s
file_re = re.compile(".*\.in\n")
file_list = re.findall(file_re,file_s) #
s1_A = "eeppss"
s1_B = "100"
s2_A = "tteemmpp"
s2_B = "200"
for file in file_list:
    basefile = file[0:-4].replace(s1_A,s1_B).replace(s2_A,s2_B) #Replace certain strings (s1_A, s2_A) with other strings (s1_B, s2_B)
    os.system("myqchem.csh -in "+basefile+".in -out " +basefile +".out") #Submit the file to Q-Chem with this new name.
    time.sleep(1.0) #Wait one second
    
