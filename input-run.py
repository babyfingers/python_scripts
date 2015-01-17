#!/usr/local/bin/python

# inout-copy.py
import os
import time

#functional_list = ["hf","b3","om","ox"]
functional_list = ["hf"]
atom_list = ["C","O","H","h"]
direction_list = ["x","y","z"]

command_file_name = "command.txt"

input_end = ".i"

command_end = "\nscript4\n3\n1\n"

for func in functional_list:
    for atom in atom_list:
        for dir in direction_list:
            f = open(command_file_name,"w")
            f_s = "F"+func+atom+dir+command_end
            f.write(f_s)
            f.close()
            print f_s
            command = "./limfd-tdhf.py < "+command_file_name+" &"
            os.system(command)
            print command
            time.sleep(5)
os.system("rm "+command_file_name)
