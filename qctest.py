#!/usr/local/bin/python

import numpy as np
from numpy import random as rand
import qchem as qc

v = rand.ranf((3,2))
#v[0,0] = 1
#v[0,1] = 0
print "init. v: ",v
v = qc.orthonormalize(v)
print "final v: ",v
print "vTv: "
print np.dot(v.transpose(),v)


B = rand.ranf((4,3))
B[:,0] = np.zeros(4)
B[0,0] = 1
B[:,1] = np.zeros(4)
B[1,1] = 1
B[:,2] = B[:,2]/5.0
print "B: ",B
B = qc.getPrunedBasis(B,1)
print "pruned B: ",B
