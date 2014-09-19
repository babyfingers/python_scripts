# bnzmk.py
#!/usr/local/bin/python
import math
import os.path
import time
import re

# bnzmk.py: Given a C-C length and a C-H length, construct cartesian
# coordinates for a benzene molecule.

tPI = math.pi/3.0
Rcc = input("RADIUS of C atom from center of ring?\n")
Rch = input("C-H bond length?\n")
Rch = Rch + Rcc
zc = 0.0
zh = 0.0
for i in range(6):
    t = i*tPI
    ct = math.cos(t)
    st = math.sin(t)
    xc = Rcc*ct
    xh = Rch*ct
    yc = Rcc*st
    yh = Rch*st
    print "C %10.7f %10.7f %10.7f" % (xc, yc, zc)
    print "H %10.7f %10.7f %10.7f" % (xh, yh, zh)
