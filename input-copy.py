#!/usr/local/bin/python

# inout-copy.py
import os

#functional_list = ["hf","b3","om","ox"]
functional_list = ["hf"]
atom_list = ["C","O","H","h"]
direction_list = ["x","y","z"]

input_end = ".i"

for func in functional_list:
    for atom in atom_list:
        for dir in direction_list:
            command = "cp "+func+input_end+" F"+func+atom+dir+input_end
            #print command
            os.system(command)
