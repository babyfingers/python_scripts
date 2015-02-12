#!/usr/local/bin/python2.7

import os
import sys
import argparse

# This script is a wrapper for the Quantum Espresso
# planewave (pw.x) submission script, ~/bin/pw.csh.
# It should read information from the input file
# and use that information to set important environment
# variables before running the submission script.

#Define arguments (and the flags for passing on the arguments)
parser = argparse.ArgumentParser()
flagDict = {}
parser.add_argument("-i",default=None,help="Input file name")
flagDict["i"] = "-in"
parser.add_argument("-o",default=None,help="Output file name")
flagDict["o"] = "-out"
parser.add_argument("-w",default=None,nargs="?",help="Wall time")
flagDict["w"] = "-w"
parser.add_argument("-p",default=None,nargs="?",help="Number of processors")
flagDict["p"] = "-p"
parser.add_argument("-m",default=None,nargs="?",help="Amount of memory")
flagDict["m"] = "-m"
parser.add_argument("-save",default=None,nargs="?",help="Save flag")
flagDict["save"] = "-save"
parser.add_argument("-outdir",default=None,nargs="?",help="Output directory")
flagDict["outdir"] = "-outdir"
parser.add_argument("-indir",default=None,nargs="?",help="Input directory")
flagDict["indir"] = "-indir"

#Process arguments
args = vars(parser.parse_args())
if (args["i"]==None):
    print "An template file must be specified with the option -i \"[template file]\"."
    sys.exit()
if (args["i"][-4:]!=".tmp"):
    print "The template file must have the extension '.tmp'."
    sys.exit()
if (args["o"]==None):
    print "An output file must be specified with the option -o \"[output file]\"."
    sys.exit()
if (not os.path.isfile(args["i"])):
    print "The specified input file does not exist."
    sys.exit()
if (os.path.isfile(args["o"]) or os.path.isdir(args["o"])):
    print "The specified output file exists already, please use a different file name."
    sys.exit()

thisJustIn = False
#Change the input file to reflect the output directory specified by
#"outdir". If not specified, use default based on output file name.
if (args["indir"]!=None and not os.path.isdir(args["indir"])):
    print "Cannot locate specified input directory."
    sys.exit()
if (args["indir"]!=None and args["outdir"]==None and args["save"]==None):
    #In Quantum Espresso, we can't specify an input directory distinctly
    #from an output directory. So if we don't need to save anything, copy
    #the contents of the input directory to a separate scratch directory
    #that we won't copy back from scratch.
    thisJustIn = True
    args["outdir"] = args["indir"] +"_out"
    print "Output directory name automatically set to {}".format(args["outdir"])
if (args["outdir"]==None):
    args["outdir"] = "p_"+args["o"].replace(".out","")
    print "Output directory name automatically set to {}".format(args["outdir"])
if (os.path.isdir(args["outdir"]) or os.path.isfile(args["outdir"])):
    print "The specified output directory exists already, please use a different name."
    sys.exit()
if (thisJustIn):
    print "cp -r {} {}".format(args["indir"],args["outdir"])
    os.system("cp -r {} {}".format(args["indir"],args["outdir"]))

#Create the proper input file from the template file
myTemp = args["i"]
myIn = args["i"][:-4]+".in"
f = open(myTemp,"r")
fs0 = f.read()
f.close()
fs1 = fs0.replace("xxOUTDIRxx",args["outdir"])
if (fs1==fs0):
    print "Please add flag {} to input file.".format("outdir='xxOUTDIRxx'")
    sys.exit()
f = open(myIn,"w")
f.write(fs1)
f.close()
args["i"] = myIn

#Put arguments together into command line for shell script
cmd = "/data/home/alguire/bin/pw.csh"
for myVar,myVal in args.iteritems():
    if (myVal!=None):
        myFlag = flagDict[myVar]
        cmd += " {} {}".format(myFlag,myVal)

print cmd
os.system(cmd)

