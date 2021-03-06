#sq.py
#!/usr/local/bin/python
import os.path
import sys
import time
import math
import re

# Prior to this: v1.1

# 2012-06-21: Fixed up the way the delta function
# information is transferred to the final, formatted
# output. In addition, standardized the way in which
# indices are replaced based on delta function
# information, ensuring that indices that are coupled
# via two or more other indices are handled
# appropriately. (I.e., if we have D(a)(r)D(b)(r),
# all rs and bs are replaced with as. If we have
# D(a)(r)D(i)(r), the term is deleted.

# 2012-06-20: Adding an option to consider when
# the one-electron operator is diagonal (i.e.,
# probably only ever the Fock matrix.) Also,
# removed redundancies in summation symbols in
# the output.

# Everything before this: v1.0

# 2011-12-07: master terms complete for the
# following overlap pairs: <1|1>, <1|2>,  <1|3>,
# <2|2>, <2|3>, <3|3>. Made minor fixes to style
# and made it so a summation symbol is not
# included if there are no indices over which to
# sum. Now prints to a file ("sf-xcis.tex"),
# and running the script "sfxcis.py" should
# run this script with all possible overlaps
# for SF-XCIS.

# 2011-12-06: LaTeX formatting is complete. To
# indicate a beta orbital, use a capital letter
# in the index. This will be converted to a
# lowercase index with a bar over it in the
# output. Any numerals will be converted to
# subtext, but they will not be included in any
# summation lists because it will be assumed
# that they are specific orbitals (i.e., s_{1}).

# 2011-12-05: as of the current version of this
# code, only expressions for operators between
# single-single, single-double, double-double,
# and double-triple excitations can be handled.
# To expand this list, a new master expression
# must be written, and the occ_lst and vrt_lst
# lists must be modified.

# sq.py: Program to automate some of this second
# quantization stuff. Hopefully we wlil be able
# to input a series of characters and get out a
# situation-specific set of delta functions (and
# maybe ultimately a proper set of coefficients
# and operators.)

#Index groups
occ_lsta = ['i','j','k','l','m','n']
vrt_lsta = ['a','b','c','d','e','f']
occ_lstb = ['I','J','K','L','M','N']
vrt_lstb = ['A','B','C','D','E','F']
wld_lst = ['r','s']
#Regular expressions
#Recognize a full term composed of delta functions
trm_re = re.compile("[+\-]D\([a-zA-Z0-9()]*")
#Recognize a single delta function
dlt_re = re.compile("D\([a-zA-Z0-9]+\)\([a-zA-Z0-9]+\)")
#Extract an index from that delta function
ind_re = re.compile("\([a-zA-Z0-9]+\)")
#Capital letter
cp_re = re.compile("[A-Z]")
#Numeral
nt_re = re.compile("[0-9]+")
#Master strings for delta function expressions (snm, tnm)
# s:  second quantization operators expanded into delta function terms. The logic happens here.
# t:  more for formatting; when the script prints out a solution, it will print out a modified
#     version of this string.
# nm: these numbers refer to the number of excitations on the left (n) and the right (m) of the operator
s11 = "+D(r)(s)D<a><b>D<i><j> -D<a><b>D<i>(s)D<j>(r) +D<i><j>D<a>(r)D<b>(s)"
t11 = "t_{<i>}^{I<a>}X_{<r><s>}t_{<j>}^{J<b>}"
s12 = "+D<a><b>D<c>(s)D<i><k>D<j>(r) -D<a><b>D<c>(s)D<i><j>D<k>(r) -D<a><c>D<b>(s)D<i><k>D<j>(r) +D<a><c>D<b>(s)D<i><j>D<k>(r)"
t12 = "t_{<i>}^{I<a>}X_{<r><s>}t_{<j><k>}^{J<b><c>}"
s13 = "0"
t13 = "t_{<i>}^{I<a>}X_{<r><s>}t_{<j><k><l>}^{J<b><c><d>}"
s22 = "+D<i><k>D<j><l>D<b>(r)D<a><c>D<d>(s) -D<i><k>D<j><l>D<b>(r)D<a><d>D<c>(s) -D<i><k>D<j><l>D<a>(r)D<b><c>D<d>(s) +D<i><k>D<j><l>D<a>(r)D<b>(s)D<c><d> -D<i><l>D<j><k>D<b>(r)D<a><c>D<d>(s) +D<i><l>D<j><k>D<b>(r)D<a><d>D<c>(s) +D<i><l>D<j><k>D<a>(r)D<b><c>D<d>(s) -D<i><l>D<j><k>D<a>(r)D<b>(s)D<c><d> +D<a><c>D<b><d>D<i><k>D<j><l>D(r)(s) -D<a><c>D<b><d>D<i><l>D<j><k>D(r)(s) -D<a><c>D<b><d>D<i><k>D<j>(s)D<l>(r) +D<a><c>D<b><d>D<i>(s)D<j><k>D<l>(r) +D<a><c>D<b><d>D<i><l>D<j>(s)D<k>(r) -D<a><c>D<b><d>D<i>(s)D<j><l>D<k>(r) -D<a><d>D<b><c>D<i><k>D<j><l>D(r)(s) +D<a><d>D<b><c>D<i><l>D<j><k>D(r)(s) +D<a><d>D<b><c>D<i><k>D<j>(s)D<l>(r) -D<a><d>D<b><c>D<i>(s)D<j><k>D<l>(r) -D<a><d>D<b><c>D<i><l>D<j>(s)D<k>(r) +D<a><d>D<b><c>D<i>(s)D<j><l>D<k>(r)"
t22 = "t_{<i><j>}^{I<a><b>}X_{<r><s>}t_{<k><l>}^{J<c><d>}"
s23 = "+D<a><d>D<b><e>D<c>(s)D<i><l>D<j><m>D<k>(r) -D<a><d>D<b><e>D<c>(s)D<i><m>D<j><l>D<k>(r) -D<a><e>D<b><d>D<c>(s)D<i><l>D<j><m>D<k>(r) +D<a><e>D<b><d>D<c>(s)D<i><m>D<j><l>D<k>(r) +D<a><d>D<b><e>D<c>(s)D<i><m>D<j><j>D<l>(r) -D<a><d>D<b><e>D<c>(s)D<i><k>D<j><m>D<l>(r) -D<a><e>D<b><d>D<c>(s)D<i><m>D<j><k>D<l>(r) +D<a><e>D<b><d>D<c>(s)D<i><k>D<j><m>D<l>(r) +D<a><d>D<b><e>D<c>(s)D<i><k>D<j><l>D<m>(r) -D<a><d>D<b><e>D<c>(s)D<i><l>D<j><k>D<m>(r) -D<a><e>D<b><d>D<c>(s)D<i><k>D<j><l>D<m>(r) +D<a><e>D<b><d>D<c>(s)D<i><l>D<j><k>D<m>(r) +D<a><e>D<b><c>D<d>(s)D<i><l>D<j><m>D<k>(r) -D<a><e>D<b><c>D<d>(s)D<i><m>D<j><l>D<k>(r) -D<a><c>D<b><e>D<d>(s)D<i><l>D<j><m>D<k>(r) +D<a><c>D<b><e>D<d>(s)D<i><m>D<j><l>D<k>(r) +D<a><e>D<b><c>D<d>(s)D<i><m>D<j><j>D<l>(r) -D<a><e>D<b><c>D<d>(s)D<i><k>D<j><m>D<l>(r) -D<a><c>D<b><e>D<d>(s)D<i><m>D<j><k>D<l>(r) +D<a><c>D<b><e>D<d>(s)D<i><k>D<j><m>D<l>(r) +D<a><e>D<b><c>D<d>(s)D<i><k>D<j><l>D<m>(r) -D<a><e>D<b><c>D<d>(s)D<i><l>D<j><k>D<m>(r) -D<a><c>D<b><e>D<d>(s)D<i><k>D<j><l>D<m>(r) +D<a><c>D<b><e>D<d>(s)D<i><l>D<j><k>D<m>(r) +D<a><c>D<b><d>D<e>(s)D<i><l>D<j><m>D<k>(r) -D<a><c>D<b><d>D<e>(s)D<i><m>D<j><l>D<k>(r) -D<a><d>D<b><c>D<e>(s)D<i><l>D<j><m>D<k>(r) +D<a><d>D<b><c>D<e>(s)D<i><m>D<j><l>D<k>(r) +D<a><c>D<b><e>D<d>(s)D<i><m>D<j><j>D<l>(r) -D<a><c>D<b><e>D<d>(s)D<i><k>D<j><m>D<l>(r) -D<a><d>D<b><c>D<e>(s)D<i><m>D<j><k>D<l>(r) +D<a><d>D<b><c>D<e>(s)D<i><k>D<j><m>D<l>(r) +D<a><c>D<b><e>D<d>(s)D<i><k>D<j><l>D<m>(r) -D<a><c>D<b><e>D<d>(s)D<i><l>D<j><k>D<m>(r) -D<a><d>D<b><c>D<e>(s)D<i><k>D<j><l>D<m>(r) +D<a><d>D<b><c>D<e>(s)D<i><l>D<j><k>D<m>(r)"
t23 = "t_{<i><j>}^{I<a><b>}X_{<r><s>}t_{<k><l><m>}^{J<c><d><e>}"
s33 = "+D<a><d>D<b><e>D<c><f>D<i><n>D<j><l>D<k><m>D(r)(s) -D<a><d>D<b><e>D<c><f>D<i><n>D<j><m>D<k><l>D(r)(s) -D<a><d>D<b><e>D<c><f>D<i><m>D<j><l>D<k><n>D(r)(s) +D<a><d>D<b><e>D<c><f>D<i><m>D<j><n>D<k><l>D(r)(s) +D<a><d>D<b><e>D<c><f>D<i><l>D<j><m>D<k><n>D(r)(s) -D<a><d>D<b><e>D<c><f>D<i><l>D<j><n>D<k><m>D(r)(s) -D<a><d>D<b><e>D<c><f>D<i>(s)D<j><l>D<k><m>D<n>(r) +D<a><d>D<b><e>D<c><f>D<i>(s)D<j><m>D<k><l>D<n>(r) +D<a><d>D<b><e>D<c><f>D<i><m>D<j><l>D<k>(s)D<n>(r) -D<a><d>D<b><e>D<c><f>D<i><m>D<j>(s)D<k><l>D<n>(r) -D<a><d>D<b><e>D<c><f>D<i><l>D<j><m>D<k>(s)D<n>(r) +D<a><d>D<b><e>D<c><f>D<i><l>D<j>(s)D<k><m>D<n>(r) +D<a><d>D<b><e>D<c><f>D<i>(s)D<j><l>D<k><n>D<m>(r) -D<a><d>D<b><e>D<c><f>D<i>(s)D<j><n>D<k><l>D<m>(r) -D<a><d>D<b><e>D<c><f>D<i><n>D<j><l>D<k>(s)D<m>(r) +D<a><d>D<b><e>D<c><f>D<i><n>D<j>(s)D<k><l>D<m>(r) +D<a><d>D<b><e>D<c><f>D<i><l>D<j><n>D<k>(s)D<m>(r) -D<a><d>D<b><e>D<c><f>D<i><l>D<j>(s)D<k><n>D<m>(r) -D<a><d>D<b><e>D<c><f>D<i>(s)D<j><m>D<k><n>D<l>(r) +D<a><d>D<b><e>D<c><f>D<i>(s)D<j><n>D<k><m>D<l>(r) +D<a><d>D<b><e>D<c><f>D<i><n>D<j><m>D<k>(s)D<l>(r) -D<a><d>D<b><e>D<c><f>D<i><n>D<j>(s)D<k><m>D<l>(r) -D<a><d>D<b><e>D<c><f>D<i><m>D<j><n>D<k>(s)D<l>(r) +D<a><d>D<b><e>D<c><f>D<i><m>D<j>(s)D<k><n>D<l>(r) -D<a><d>D<b><f>D<c><e>D<i><n>D<j><l>D<k><m>D(r)(s) +D<a><d>D<b><f>D<c><e>D<i><n>D<j><m>D<k><l>D(r)(s) +D<a><d>D<b><f>D<c><e>D<i><m>D<j><l>D<k><n>D(r)(s) -D<a><d>D<b><f>D<c><e>D<i><m>D<j><n>D<k><l>D(r)(s) -D<a><d>D<b><f>D<c><e>D<i><l>D<j><m>D<k><n>D(r)(s) +D<a><d>D<b><f>D<c><e>D<i><l>D<j><n>D<k><m>D(r)(s) +D<a><d>D<b><f>D<c><e>D<i>(s)D<j><l>D<k><m>D<n>(r) -D<a><d>D<b><f>D<c><e>D<i>(s)D<j><m>D<k><l>D<n>(r) -D<a><d>D<b><f>D<c><e>D<i><m>D<j><l>D<k>(s)D<n>(r) +D<a><d>D<b><f>D<c><e>D<i><m>D<j>(s)D<k><l>D<n>(r) +D<a><d>D<b><f>D<c><e>D<i><l>D<j><m>D<k>(s)D<n>(r) -D<a><d>D<b><f>D<c><e>D<i><l>D<j>(s)D<k><m>D<n>(r) -D<a><d>D<b><f>D<c><e>D<i>(s)D<j><l>D<k><n>D<m>(r) +D<a><d>D<b><f>D<c><e>D<i>(s)D<j><n>D<k><l>D<m>(r) +D<a><d>D<b><f>D<c><e>D<i><n>D<j><l>D<k>(s)D<m>(r) -D<a><d>D<b><f>D<c><e>D<i><n>D<j>(s)D<k><l>D<m>(r) -D<a><d>D<b><f>D<c><e>D<i><l>D<j><n>D<k>(s)D<m>(r) +D<a><d>D<b><f>D<c><e>D<i><l>D<j>(s)D<k><n>D<m>(r) +D<a><d>D<b><f>D<c><e>D<i>(s)D<j><m>D<k><n>D<l>(r) -D<a><d>D<b><f>D<c><e>D<i>(s)D<j><n>D<k><m>D<l>(r) -D<a><d>D<b><f>D<c><e>D<i><n>D<j><m>D<k>(s)D<l>(r) +D<a><d>D<b><f>D<c><e>D<i><n>D<j>(s)D<k><m>D<l>(r) +D<a><d>D<b><f>D<c><e>D<i><m>D<j><n>D<k>(s)D<l>(r) -D<a><d>D<b><f>D<c><e>D<i><m>D<j>(s)D<k><n>D<l>(r) -D<a><e>D<b><d>D<c><f>D<i><n>D<j><l>D<k><m>D(r)(s) +D<a><e>D<b><d>D<c><f>D<i><n>D<j><m>D<k><l>D(r)(s) +D<a><e>D<b><d>D<c><f>D<i><m>D<j><l>D<k><n>D(r)(s) -D<a><e>D<b><d>D<c><f>D<i><m>D<j><n>D<k><l>D(r)(s) -D<a><e>D<b><d>D<c><f>D<i><l>D<j><m>D<k><n>D(r)(s) +D<a><e>D<b><d>D<c><f>D<i><l>D<j><n>D<k><m>D(r)(s) +D<a><e>D<b><d>D<c><f>D<i>(s)D<j><l>D<k><m>D<n>(r) -D<a><e>D<b><d>D<c><f>D<i>(s)D<j><m>D<k><l>D<n>(r) -D<a><e>D<b><d>D<c><f>D<i><m>D<j><l>D<k>(s)D<n>(r) +D<a><e>D<b><d>D<c><f>D<i><m>D<j>(s)D<k><l>D<n>(r) +D<a><e>D<b><d>D<c><f>D<i><l>D<j><m>D<k>(s)D<n>(r) -D<a><e>D<b><d>D<c><f>D<i><l>D<j>(s)D<k><m>D<n>(r) -D<a><e>D<b><d>D<c><f>D<i>(s)D<j><l>D<k><n>D<m>(r) +D<a><e>D<b><d>D<c><f>D<i>(s)D<j><n>D<k><l>D<m>(r) +D<a><e>D<b><d>D<c><f>D<i><n>D<j><l>D<k>(s)D<m>(r) -D<a><e>D<b><d>D<c><f>D<i><n>D<j>(s)D<k><l>D<m>(r) -D<a><e>D<b><d>D<c><f>D<i><l>D<j><n>D<k>(s)D<m>(r) +D<a><e>D<b><d>D<c><f>D<i><l>D<j>(s)D<k><n>D<m>(r) +D<a><e>D<b><d>D<c><f>D<i>(s)D<j><m>D<k><n>D<l>(r) -D<a><e>D<b><d>D<c><f>D<i>(s)D<j><n>D<k><m>D<l>(r) -D<a><e>D<b><d>D<c><f>D<i><n>D<j><m>D<k>(s)D<l>(r) +D<a><e>D<b><d>D<c><f>D<i><n>D<j>(s)D<k><m>D<l>(r) +D<a><e>D<b><d>D<c><f>D<i><m>D<j><n>D<k>(s)D<l>(r) -D<a><e>D<b><d>D<c><f>D<i><m>D<j>(s)D<k><n>D<l>(r) +D<a><e>D<b><f>D<c><d>D<i><n>D<j><l>D<k><m>D(r)(s) -D<a><e>D<b><f>D<c><d>D<i><n>D<j><m>D<k><l>D(r)(s) -D<a><e>D<b><f>D<c><d>D<i><m>D<j><l>D<k><n>D(r)(s) +D<a><e>D<b><f>D<c><d>D<i><m>D<j><n>D<k><l>D(r)(s) +D<a><e>D<b><f>D<c><d>D<i><l>D<j><m>D<k><n>D(r)(s) -D<a><e>D<b><f>D<c><d>D<i><l>D<j><n>D<k><m>D(r)(s) -D<a><e>D<b><f>D<c><d>D<i>(s)D<j><l>D<k><m>D<n>(r) +D<a><e>D<b><f>D<c><d>D<i>(s)D<j><m>D<k><l>D<n>(r) +D<a><e>D<b><f>D<c><d>D<i><m>D<j><l>D<k>(s)D<n>(r) -D<a><e>D<b><f>D<c><d>D<i><m>D<j>(s)D<k><l>D<n>(r) -D<a><e>D<b><f>D<c><d>D<i><l>D<j><m>D<k>(s)D<n>(r) +D<a><e>D<b><f>D<c><d>D<i><l>D<j>(s)D<k><m>D<n>(r) +D<a><e>D<b><f>D<c><d>D<i>(s)D<j><l>D<k><n>D<m>(r) -D<a><e>D<b><f>D<c><d>D<i>(s)D<j><n>D<k><l>D<m>(r) -D<a><e>D<b><f>D<c><d>D<i><n>D<j><l>D<k>(s)D<m>(r) +D<a><e>D<b><f>D<c><d>D<i><n>D<j>(s)D<k><l>D<m>(r) +D<a><e>D<b><f>D<c><d>D<i><l>D<j><n>D<k>(s)D<m>(r) -D<a><e>D<b><f>D<c><d>D<i><l>D<j>(s)D<k><n>D<m>(r) -D<a><e>D<b><f>D<c><d>D<i>(s)D<j><m>D<k><n>D<l>(r) +D<a><e>D<b><f>D<c><d>D<i>(s)D<j><n>D<k><m>D<l>(r) +D<a><e>D<b><f>D<c><d>D<i><n>D<j><m>D<k>(s)D<l>(r) -D<a><e>D<b><f>D<c><d>D<i><n>D<j>(s)D<k><m>D<l>(r) -D<a><e>D<b><f>D<c><d>D<i><m>D<j><n>D<k>(s)D<l>(r) +D<a><e>D<b><f>D<c><d>D<i><m>D<j>(s)D<k><n>D<l>(r) +D<a><f>D<b><d>D<c><e>D<i><n>D<j><l>D<k><m>D(r)(s) -D<a><f>D<b><d>D<c><e>D<i><n>D<j><m>D<k><l>D(r)(s) -D<a><f>D<b><d>D<c><e>D<i><m>D<j><l>D<k><n>D(r)(s) +D<a><f>D<b><d>D<c><e>D<i><m>D<j><n>D<k><l>D(r)(s) +D<a><f>D<b><d>D<c><e>D<i><l>D<j><m>D<k><n>D(r)(s) -D<a><f>D<b><d>D<c><e>D<i><l>D<j><n>D<k><m>D(r)(s) -D<a><f>D<b><d>D<c><e>D<i>(s)D<j><l>D<k><m>D<n>(r) +D<a><f>D<b><d>D<c><e>D<i>(s)D<j><m>D<k><l>D<n>(r) +D<a><f>D<b><d>D<c><e>D<i><m>D<j><l>D<k>(s)D<n>(r) -D<a><f>D<b><d>D<c><e>D<i><m>D<j>(s)D<k><l>D<n>(r) -D<a><f>D<b><d>D<c><e>D<i><l>D<j><m>D<k>(s)D<n>(r) +D<a><f>D<b><d>D<c><e>D<i><l>D<j>(s)D<k><m>D<n>(r) +D<a><f>D<b><d>D<c><e>D<i>(s)D<j><l>D<k><n>D<m>(r) -D<a><f>D<b><d>D<c><e>D<i>(s)D<j><n>D<k><l>D<m>(r) -D<a><f>D<b><d>D<c><e>D<i><n>D<j><l>D<k>(s)D<m>(r) +D<a><f>D<b><d>D<c><e>D<i><n>D<j>(s)D<k><l>D<m>(r) +D<a><f>D<b><d>D<c><e>D<i><l>D<j><n>D<k>(s)D<m>(r) -D<a><f>D<b><d>D<c><e>D<i><l>D<j>(s)D<k><n>D<m>(r) -D<a><f>D<b><d>D<c><e>D<i>(s)D<j><m>D<k><n>D<l>(r) +D<a><f>D<b><d>D<c><e>D<i>(s)D<j><n>D<k><m>D<l>(r) +D<a><f>D<b><d>D<c><e>D<i><n>D<j><m>D<k>(s)D<l>(r) -D<a><f>D<b><d>D<c><e>D<i><n>D<j>(s)D<k><m>D<l>(r) -D<a><f>D<b><d>D<c><e>D<i><m>D<j><n>D<k>(s)D<l>(r) +D<a><f>D<b><d>D<c><e>D<i><m>D<j>(s)D<k><n>D<l>(r) -D<a><f>D<b><e>D<c><d>D<i><n>D<j><l>D<k><m>D(r)(s) +D<a><f>D<b><e>D<c><d>D<i><n>D<j><m>D<k><l>D(r)(s) +D<a><f>D<b><e>D<c><d>D<i><m>D<j><l>D<k><n>D(r)(s) -D<a><f>D<b><e>D<c><d>D<i><m>D<j><n>D<k><l>D(r)(s) -D<a><f>D<b><e>D<c><d>D<i><l>D<j><m>D<k><n>D(r)(s) +D<a><f>D<b><e>D<c><d>D<i><l>D<j><n>D<k><m>D(r)(s) +D<a><f>D<b><e>D<c><d>D<i>(s)D<j><l>D<k><m>D<n>(r) -D<a><f>D<b><e>D<c><d>D<i>(s)D<j><m>D<k><l>D<n>(r) -D<a><f>D<b><e>D<c><d>D<i><m>D<j><l>D<k>(s)D<n>(r) +D<a><f>D<b><e>D<c><d>D<i><m>D<j>(s)D<k><l>D<n>(r) +D<a><f>D<b><e>D<c><d>D<i><l>D<j><m>D<k>(s)D<n>(r) -D<a><f>D<b><e>D<c><d>D<i><l>D<j>(s)D<k><m>D<n>(r) -D<a><f>D<b><e>D<c><d>D<i>(s)D<j><l>D<k><n>D<m>(r) +D<a><f>D<b><e>D<c><d>D<i>(s)D<j><n>D<k><l>D<m>(r) +D<a><f>D<b><e>D<c><d>D<i><n>D<j><l>D<k>(s)D<m>(r) -D<a><f>D<b><e>D<c><d>D<i><n>D<j>(s)D<k><l>D<m>(r) -D<a><f>D<b><e>D<c><d>D<i><l>D<j><n>D<k>(s)D<m>(r) +D<a><f>D<b><e>D<c><d>D<i><l>D<j>(s)D<k><n>D<m>(r) +D<a><f>D<b><e>D<c><d>D<i>(s)D<j><m>D<k><n>D<l>(r) -D<a><f>D<b><e>D<c><d>D<i>(s)D<j><n>D<k><m>D<l>(r) -D<a><f>D<b><e>D<c><d>D<i><n>D<j><m>D<k>(s)D<l>(r) +D<a><f>D<b><e>D<c><d>D<i><n>D<j>(s)D<k><m>D<l>(r) +D<a><f>D<b><e>D<c><d>D<i><m>D<j><n>D<k>(s)D<l>(r) -D<a><f>D<b><e>D<c><d>D<i><m>D<j>(s)D<k><n>D<l>(r) +D<a>(r)D<b><d>D<c><e>D<f>(s)D<i><l>D<j><m>D<k><n> -D<a>(r)D<b><e>D<c><d>D<f>(s)D<i><l>D<j><m>D<k><n> -D<a>(r)D<b><d>D<c><f>D<e>(s)D<i><l>D<j><m>D<k><n> +D<a>(r)D<b><f>D<c><d>D<e>(s)D<i><l>D<j><m>D<k><n> +D<a>(r)D<b><e>D<c><f>D<d>(s)D<i><l>D<j><m>D<k><n> -D<a>(r)D<b><f>D<c><e>D<d>(s)D<i><l>D<j><m>D<k><n> -D<a><d>D<b>(r)D<c><e>D<f>(s)D<i><l>D<j><m>D<k><n> +D<a><e>D<b>(r)D<c><d>D<f>(s)D<i><l>D<j><m>D<k><n> +D<a><d>D<b>(r)D<c><f>D<e>(s)D<i><l>D<j><m>D<k><n> -D<a><f>D<b>(r)D<c><d>D<e>(s)D<i><l>D<j><m>D<k><n> -D<a><e>D<b>(r)D<c><f>D<d>(s)D<i><l>D<j><m>D<k><n> +D<a><f>D<b>(r)D<c><e>D<d>(s)D<i><l>D<j><m>D<k><n> +D<a><d>D<b><e>D<c>(r)D<f>(s)D<i><l>D<j><m>D<k><n> -D<a><e>D<b><d>D<c>(r)D<f>(s)D<i><l>D<j><m>D<k><n> -D<a><d>D<b><f>D<c>(r)D<e>(s)D<i><l>D<j><m>D<k><n> +D<a><f>D<b><d>D<c>(r)D<e>(s)D<i><l>D<j><m>D<k><n> +D<a><e>D<b><f>D<c>(r)D<d>(s)D<i><l>D<j><m>D<k><n> -D<a><f>D<b><e>D<c>(r)D<d>(s)D<i><l>D<j><m>D<k><n> -D<a>(r)D<b><d>D<c><e>D<f>(s)D<i><m>D<j><l>D<k><n> +D<a>(r)D<b><e>D<c><d>D<f>(s)D<i><m>D<j><l>D<k><n> +D<a>(r)D<b><d>D<c><f>D<e>(s)D<i><m>D<j><l>D<k><n> -D<a>(r)D<b><f>D<c><d>D<e>(s)D<i><m>D<j><l>D<k><n> -D<a>(r)D<b><e>D<c><f>D<d>(s)D<i><m>D<j><l>D<k><n> +D<a>(r)D<b><f>D<c><e>D<d>(s)D<i><m>D<j><l>D<k><n> +D<a><d>D<b>(r)D<c><e>D<f>(s)D<i><m>D<j><l>D<k><n> -D<a><e>D<b>(r)D<c><d>D<f>(s)D<i><m>D<j><l>D<k><n> -D<a><d>D<b>(r)D<c><f>D<e>(s)D<i><m>D<j><l>D<k><n> +D<a><f>D<b>(r)D<c><d>D<e>(s)D<i><m>D<j><l>D<k><n> +D<a><e>D<b>(r)D<c><f>D<d>(s)D<i><m>D<j><l>D<k><n> -D<a><f>D<b>(r)D<c><e>D<d>(s)D<i><m>D<j><l>D<k><n> -D<a><d>D<b><e>D<c>(r)D<f>(s)D<i><m>D<j><l>D<k><n> +D<a><e>D<b><d>D<c>(r)D<f>(s)D<i><m>D<j><l>D<k><n> +D<a><d>D<b><f>D<c>(r)D<e>(s)D<i><m>D<j><l>D<k><n> -D<a><f>D<b><d>D<c>(r)D<e>(s)D<i><m>D<j><l>D<k><n> -D<a><e>D<b><f>D<c>(r)D<d>(s)D<i><m>D<j><l>D<k><n> +D<a><f>D<b><e>D<c>(r)D<d>(s)D<i><m>D<j><l>D<k><n> -D<a>(r)D<b><d>D<c><e>D<f>(s)D<i><l>D<j><n>D<k><m> +D<a>(r)D<b><e>D<c><d>D<f>(s)D<i><l>D<j><n>D<k><m> +D<a>(r)D<b><d>D<c><f>D<e>(s)D<i><l>D<j><n>D<k><m> -D<a>(r)D<b><f>D<c><d>D<e>(s)D<i><l>D<j><n>D<k><m> -D<a>(r)D<b><e>D<c><f>D<d>(s)D<i><l>D<j><n>D<k><m> +D<a>(r)D<b><f>D<c><e>D<d>(s)D<i><l>D<j><n>D<k><m> +D<a><d>D<b>(r)D<c><e>D<f>(s)D<i><l>D<j><n>D<k><m> -D<a><e>D<b>(r)D<c><d>D<f>(s)D<i><l>D<j><n>D<k><m> -D<a><d>D<b>(r)D<c><f>D<e>(s)D<i><l>D<j><n>D<k><m> +D<a><f>D<b>(r)D<c><d>D<e>(s)D<i><l>D<j><n>D<k><m> +D<a><e>D<b>(r)D<c><f>D<d>(s)D<i><l>D<j><n>D<k><m> -D<a><f>D<b>(r)D<c><e>D<d>(s)D<i><l>D<j><n>D<k><m> -D<a><d>D<b><e>D<c>(r)D<f>(s)D<i><l>D<j><n>D<k><m> +D<a><e>D<b><d>D<c>(r)D<f>(s)D<i><l>D<j><n>D<k><m> +D<a><d>D<b><f>D<c>(r)D<e>(s)D<i><l>D<j><n>D<k><m> -D<a><f>D<b><d>D<c>(r)D<e>(s)D<i><l>D<j><n>D<k><m> -D<a><e>D<b><f>D<c>(r)D<d>(s)D<i><l>D<j><n>D<k><m> +D<a><f>D<b><e>D<c>(r)D<d>(s)D<i><l>D<j><n>D<k><m> +D<a>(r)D<b><d>D<c><e>D<f>(s)D<i><n>D<j><l>D<k><m> -D<a>(r)D<b><e>D<c><d>D<f>(s)D<i><n>D<j><l>D<k><m> -D<a>(r)D<b><d>D<c><f>D<e>(s)D<i><n>D<j><l>D<k><m> +D<a>(r)D<b><f>D<c><d>D<e>(s)D<i><n>D<j><l>D<k><m> +D<a>(r)D<b><e>D<c><f>D<d>(s)D<i><n>D<j><l>D<k><m> -D<a>(r)D<b><f>D<c><e>D<d>(s)D<i><n>D<j><l>D<k><m> -D<a><d>D<b>(r)D<c><e>D<f>(s)D<i><n>D<j><l>D<k><m> +D<a><e>D<b>(r)D<c><d>D<f>(s)D<i><n>D<j><l>D<k><m> +D<a><d>D<b>(r)D<c><f>D<e>(s)D<i><n>D<j><l>D<k><m> -D<a><f>D<b>(r)D<c><d>D<e>(s)D<i><n>D<j><l>D<k><m> -D<a><e>D<b>(r)D<c><f>D<d>(s)D<i><n>D<j><l>D<k><m> +D<a><f>D<b>(r)D<c><e>D<d>(s)D<i><n>D<j><l>D<k><m> +D<a><d>D<b><e>D<c>(r)D<f>(s)D<i><n>D<j><l>D<k><m> -D<a><e>D<b><d>D<c>(r)D<f>(s)D<i><n>D<j><l>D<k><m> -D<a><d>D<b><f>D<c>(r)D<e>(s)D<i><n>D<j><l>D<k><m> +D<a><f>D<b><d>D<c>(r)D<e>(s)D<i><n>D<j><l>D<k><m> +D<a><e>D<b><f>D<c>(r)D<d>(s)D<i><n>D<j><l>D<k><m> -D<a><f>D<b><e>D<c>(r)D<d>(s)D<i><n>D<j><l>D<k><m> +D<a>(r)D<b><d>D<c><e>D<f>(s)D<i><m>D<j><n>D<k><l> -D<a>(r)D<b><e>D<c><d>D<f>(s)D<i><m>D<j><n>D<k><l> -D<a>(r)D<b><d>D<c><f>D<e>(s)D<i><m>D<j><n>D<k><l> +D<a>(r)D<b><f>D<c><d>D<e>(s)D<i><m>D<j><n>D<k><l> +D<a>(r)D<b><e>D<c><f>D<d>(s)D<i><m>D<j><n>D<k><l> -D<a>(r)D<b><f>D<c><e>D<d>(s)D<i><m>D<j><n>D<k><l> -D<a><d>D<b>(r)D<c><e>D<f>(s)D<i><m>D<j><n>D<k><l> +D<a><e>D<b>(r)D<c><d>D<f>(s)D<i><m>D<j><n>D<k><l> +D<a><d>D<b>(r)D<c><f>D<e>(s)D<i><m>D<j><n>D<k><l> -D<a><f>D<b>(r)D<c><d>D<e>(s)D<i><m>D<j><n>D<k><l> -D<a><e>D<b>(r)D<c><f>D<d>(s)D<i><m>D<j><n>D<k><l> +D<a><f>D<b>(r)D<c><e>D<d>(s)D<i><m>D<j><n>D<k><l> +D<a><d>D<b><e>D<c>(r)D<f>(s)D<i><m>D<j><n>D<k><l> -D<a><e>D<b><d>D<c>(r)D<f>(s)D<i><m>D<j><n>D<k><l> -D<a><d>D<b><f>D<c>(r)D<e>(s)D<i><m>D<j><n>D<k><l> +D<a><f>D<b><d>D<c>(r)D<e>(s)D<i><m>D<j><n>D<k><l> +D<a><e>D<b><f>D<c>(r)D<d>(s)D<i><m>D<j><n>D<k><l> -D<a><f>D<b><e>D<c>(r)D<d>(s)D<i><m>D<j><n>D<k><l> -D<a>(r)D<b><d>D<c><e>D<f>(s)D<i><n>D<j><m>D<k><l> +D<a>(r)D<b><e>D<c><d>D<f>(s)D<i><n>D<j><m>D<k><l> +D<a>(r)D<b><d>D<c><f>D<e>(s)D<i><n>D<j><m>D<k><l> -D<a>(r)D<b><f>D<c><d>D<e>(s)D<i><n>D<j><m>D<k><l> -D<a>(r)D<b><e>D<c><f>D<d>(s)D<i><n>D<j><m>D<k><l> +D<a>(r)D<b><f>D<c><e>D<d>(s)D<i><n>D<j><m>D<k><l> +D<a><d>D<b>(r)D<c><e>D<f>(s)D<i><n>D<j><m>D<k><l> -D<a><e>D<b>(r)D<c><d>D<f>(s)D<i><n>D<j><m>D<k><l> -D<a><d>D<b>(r)D<c><f>D<e>(s)D<i><n>D<j><m>D<k><l> +D<a><f>D<b>(r)D<c><d>D<e>(s)D<i><n>D<j><m>D<k><l> +D<a><e>D<b>(r)D<c><f>D<d>(s)D<i><n>D<j><m>D<k><l> -D<a><f>D<b>(r)D<c><e>D<d>(s)D<i><n>D<j><m>D<k><l> -D<a><d>D<b><e>D<c>(r)D<f>(s)D<i><n>D<j><m>D<k><l> +D<a><e>D<b><d>D<c>(r)D<f>(s)D<i><n>D<j><m>D<k><l> +D<a><d>D<b><f>D<c>(r)D<e>(s)D<i><n>D<j><m>D<k><l> -D<a><f>D<b><d>D<c>(r)D<e>(s)D<i><n>D<j><m>D<k><l> -D<a><e>D<b><f>D<c>(r)D<d>(s)D<i><n>D<j><m>D<k><l> +D<a><f>D<b><e>D<c>(r)D<d>(s)D<i><n>D<j><m>D<k><l>"
t33 = "t_{<i><j><k>}^{I<a><b><c>}X_{<r><s>}t_{<l><m><n>}^{J<d><e><f>}"
#Useful functions

# sort_sumindex: takes a list of strings, organizes alphabeticaly
#and removes redundand entries. Intended for use when printing
#summation symbols to the output, but can be used for other things.
def sort_sumindex(l):
    for s in l:
        while (l.count(s)>1):
            l.remove(s)
    l.sort()
    return l

# dlt: a delta function that sorts input based on
#the index notation we've established.
def dlt(x,y):
    #Case 1: x is the same as y
    if (x==y):
        return ''
    #Case 2: x, y are in the same group, or one or both is wild
    if ((x in occ_lsta) and (y in occ_lsta)) or ((x in vrt_lsta) and (y in vrt_lsta)) or ((x in occ_lstb) and (y in occ_lstb)) or ((x in vrt_lstb) and (y in vrt_lstb)) or (x in wld_lst) or (y in wld_lst):
        return 'D('+x+')('+y+')'
    #Case 3: x, y can never be the same
    return '0'

# findindex: a little space-saver. Given a well-formatted
#delta function string, will extract the indices and return
#them as a list.
def findindex(s):
    lst = ind_re.findall(s)
    l = lst[0]
    l = l[1:-1]
    lst[0] = l
    l = lst[1]
    l = l[1:-1]
    lst[1] = l
    return lst

# tpose: takes a 2-list and switches their elements. Returns
#the new 2-list.
def tpose(l):
    l1 = l[0]
    l[0] = l[1]
    l[1] = l1
    return l

# indecide: takes a 2-list of indices of a delta function
#and chooses the one that should be eliminateds in the final
#expression. Only designed to handle indices that can be
#found in the same group, which should be fine based on the
#way this program has been designed. Returns a 2-list again,
#but with the index to be removed in the second position.
def indecide(l):
    #r will store the rank of the indices of l. Return index
    #with higher rank.
    print "indecide: l = ",l
    r = []    
    #Check for wild indices (if present, move to last position)
    if (l[0] in wld_lst) and (l[1] in wld_lst):
        r.append(wld_lst.index(l[0]))
        r.append(wld_lst.index(l[1]))
    elif (l[0] in wld_lst):
        l = tpose(l)
        return l
    elif (l[1] in wld_lst):
        return l
    #Find appropriate group
    for x in l:
        if x in occ_lsta:
            r.append(occ_lsta.index(x))
        elif x in vrt_lsta:
            r.append(vrt_lsta.index(x))
        elif x in occ_lstb:
            r.append(occ_lstb.index(x))
        elif x in vrt_lstb:
            r.append(vrt_lstb.index(x))
    print "indecide: r = ",r
    if (r[0]>r[1]):
        l = tpose(l)
    print "indecide: final l = ",l
    return l

# indformat: takes a proposed index and prepares it in various
#ways for LaTeX. First, if a capital letter is detected, it is
#demoted to lowercase and the capitalized part of the index is
#enclosed in a \bar{}, so that when entering indices a capital
#may stand for the index of a beta orbital. Second, if any
#numbers are present, they are enclosed within subscript
#notation (_{}).
def indformat(s):
    cp_lst = cp_re.findall(s)
    nt_lst = nt_re.findall(s)
    if (len(cp_lst)>0):
        cp = cp_lst[0]
        s = s.replace(cp,"\\bar{"+cp.lower()+"}")
    if (len(nt_lst)>0):
        nt = nt_lst[0]
        s = s.replace(nt,"_{"+nt+"}")
    return s
#Step 0: get user input on type of expression
#Get user input on the type of overlap (singlet-singlet, singlet-doublet, etc.)
ln = input("Number of excitations on the left (1, 2, or 3):\n")
while (ln not in [1,2,3]):
    ln = input("Try again (1, 2, or 3):\n")
rn = input("Number of excitations on the right (1, 2, or 3):\n")
while (rn not in [1,2,3] or rn<ln):
    rn = input("Try again (1, 2, or 3, and no less than first number):\n")
#Use this information to select appropriate master expression
if (ln==1 and rn==1):
    si = s11
    ti = t11
elif (ln==1 and rn==2):
    si = s12
    ti = t12
elif (ln==1 and rn==3):
    si = s13
    ti = t13
elif (ln==2 and rn==2):
    si = s22
    ti = t22
elif (ln==2 and rn==3):
    si = s23
    ti = t23
elif (ln==3 and rn==3):
    si = s33
    ti = t33
#Get user input on the specific indices to be used
fn = ln + rn
locc = []
rocc = []
lvrt = []
rvrt = []
for i in range(ln):
    q = "Occupied orbital "+str(i)+" on the left side:\n"
    s = raw_input(q)
    locc.append(s)
for i in range(ln):
    q = "Virtual orbital "+str(i)+" on the left side:\n"
    s = raw_input(q)
    lvrt.append(s)
for i in range(rn):
    q = "Occupied orbital "+str(i)+" on the right side:\n"
    s = raw_input(q)
    rocc.append(s)
for i in range(rn):
    q = "Virtual orbital "+str(i)+" on the right side:\n"
    s = raw_input(q)
    rvrt.append(s)
focc = list(locc+rocc)
fvrt = list(lvrt+rvrt)
forb = list(focc+fvrt)
#Special case: operator is diagonal (i.e. Fock matrix)
op_diag = input("Is one-electron operator X diagonal (0) or not (1)?\n")
if (op_diag==0):
    si = si.replace("(s)","(r)")
    ti = ti.replace("<s>","<r>")
    wld_lst.remove("s")

#Step 1: Do situation-specific replacements
print "template si:",si
print "template ti:",ti
for i in range(fn):
    print "Situational replacement: "+occ_lsta[i]+" -> "+focc[i]
    print "Situational replacement: "+vrt_lsta[i]+" -> "+fvrt[i]
    occi = "<"+occ_lsta[i]+">"
    vrti = "<"+vrt_lsta[i]+">"
    occo = "("+focc[i]+")"
    vrto = "("+fvrt[i]+")"
    si = si.replace(occi,occo).replace(vrti,vrto)
    ti = ti.replace(occi,occo).replace(vrti,vrto)
    print "pre-situational ti:",ti
ti = ti.replace("(","<").replace(")",">")
print "situational si:",si
print "situational ti:",ti
#sit_lst: master list of all indices we will start with for this situation
# (i.e., neither the template indices nor the final, delta-function augmented
#  indices).
sit_lst = sort_sumindex(forb)+wld_lst
print "sit_lst:",sit_lst
#Step 2: Clean up delta function expressions 
so = ''
trm_list = trm_re.findall(si)
for t in trm_list:
    dlt_list = dlt_re.findall(t)
    st = t[0]
    for di in dlt_list:
        print "di = "+di
        ind = findindex(di)
        do = dlt(ind[0],ind[1]) #Determines whether to delete function (D = 1) or term (D = 0) or leave it alone.
        print "do = "+do
        if (do=='0'):
            st = ""
            break
        st += do
    # Look for wild card matches
    twld = list(wld_lst) #Copy wld_list over to a temporary list that we might alter
    dlt_list = dlt_re.findall(st)
    ind_lst = []
    for di in dlt_list:
        ind_lst.append(findindex(di))
    for l in ind_lst:
        ind = indecide(l)
        if (ind[1] in wld_lst):
            print "consult twld: wild index found: ",l
            #If the second index matches a wild card, do one of two things:
            # (1) if the second index matches something in the temporary wild
            #     card list, replace it with the first index.
            # (2) otherwise, the matching index has already been replaced
            #     within the temporary wild card list. If this is the case,
            #     create a delta function between the two (as long as the
            #     delta function would be nontrivial, i.e. not D(a)(a)
            s0 = ind[0]
            s1 = twld[wld_lst.index(ind[1])]
            if (ind[1] in twld):
                print "Replacing temporary wild card list from overlap."
                print "(before) twld = ",twld
                for i in range(len(twld)):
                    twld[i] = twld[i].replace(ind[1],ind[0])
                print "(after) twld = ",twld
            else:
                # Create delta function if nontrivial
                if (s0 != s1):
                    print "Adding delta function from wild card overlap."
                    print "(before) st = ",st
                    s2 = "D("+s0+")("+s1+")"
                    st = st+s2
                    print "(after)  st = ",st
    if (st!=''):
        so += st + ' '

#Step 3: Convert delta function terms into proper expression
# with excitation coefficients and operators and print in
# LaTeX format.
#Preamble
if (op_diag==0):
    latout = "\\begin{flalign*}\n&\\textrm{(diag) }\langle\Phi_{"
else:
    latout = "\\begin{flalign*}\n&\langle\Phi_{"
for s in locc:
    latout += indformat(s)
latout += "}^{"
for s in lvrt:
    latout += indformat(s)
latout += "}\\vert\hat{X}\\vert\Phi_{"
for s in rocc:
    latout += indformat(s)
latout += "}^{"
for s in rvrt:
    latout += indformat(s)
latout += "}\\rangle = "
#Term-by-term
trm_list = trm_re.findall(so)
#No terms? Whole thing is zero.
if (len(trm_list)==0):
    latout += "0"
else:
    #Look at delta function expressions term by term to construct string to be printed
    nterms = 1
    for t in trm_list:
        print "NEW TERM!"
        print "t = "+t
        #Assign sign
        latout += t[0]
        #Extract relevant information from delta functions
        #Get list of delta functions
        dlt_list = dlt_re.findall(t)
        ind_lst = []
        #Transform list of delta functions into list of 2-tuples
        for d in dlt_list:
            print "d = "+d
            ind = findindex(d)
            print "ind = ",ind
            ind_lst.append(ind)
        print "ind_lst = ",ind_lst
        #Assign list of indices over which to sum
        #Get full list of possible indices
        sumind = []
        for o in focc:
            nt_lst = nt_re.findall(o)
            if (len(nt_lst)==0):
                sumind.append(o)
        for v in fvrt:
            nt_lst = nt_re.findall(v)
            if (len(nt_lst)==0):
                sumind.append(v)
        tocc = list(focc) #Eventually replace duplicates for full expression
        tvrt = list(fvrt) # (ditto)
        twld = list(wld_lst) # (ditto)
        #Then remove unnecessary indices
        #Make list of indices to remove from summation in final
        #expression (also reorder so that the one to be removed
        #is always second in its entry in ind_lst).
        indremove = []
        for i in range(len(ind_lst)):
            ind = indecide(ind_lst[i])
            print "indremove: examining delta function pair ind = ",ind
            ind_lst[i] = ind
            if (ind[1] not in indremove):
                indremove.append(ind[1])
                print "Adding ",ind[1],"to indremove."
        #Copy master output term and make appropriate replacements
        tt = str(ti)
        #Replace remaining indices (for placement in expression)
        # Take a copy of the list of occupied and virtual indicies. If we have delta
        # function information from this term that says one of them should be replaced
        # with another, do so.
        rep_lst = [] #List of lists: each list contains the indices that have been matched
                     # (in this term) from delta functions to other indices. Each list
                     # corresponds positionally to matches to the master list of indices,
                     # sit_list.
                     # (i.e., if the first index of sit_list is 'A', then the first list in
                     #  rep_lst contains all the indices that have been mached to 'A' via
                     #  delta functions for this term.)
        for sit in sit_lst:
            rep_lst.append([])
        for l in ind_lst:
            print "Delta function pair: l =",l
            rep_lst[sit_lst.index(l[1])].append(l[0])
        #Now that we've created a list of all indices coupled via delta function,
        # make sure they are all compatible.
        #  Step (1): Go through full list of lists and append the list of matches
        #            of the higher-ranked one with those of the lower-ranked one.
        print "(before) rep_lst:",rep_lst
        print "sit_lst =",sit_lst,", len(sit_lst) =",len(sit_lst)
        for si in range(len(sit_lst)-1,0,-1): # Start from the last (lowest-ranked), skip the first (highest-ranked) index
            print "si =",si
            print "l =",l
            l = rep_lst[si] #So for every replacement list...
            for rsi in range(len(l)): # ...compare every index within that replacement list...
                rin = l[rsi]
                print "rsi =",rsi,", rin =",rin
                for csi in range(si): #...to every index in the full list PRECEEDING the one we're on now.
                    print "csi =",csi,", sit_lst[csi] =",sit_lst[csi]
                    if (rin==sit_lst[csi]): #If they match...
                        print "MATCH FOUND!"
                        rep_lst[csi] += l # Add all matches of the current index to those of the higher-ranking index
                        rep_lst[csi].append(sit_lst[si]) # Also, add the lower-ranking index itself to list of the higher-ranking index
                        rep_lst[si] = []  # and remove all matches from the current list.
                        break
        print "(after)  rep_lst:",rep_lst
        #  Step (2): organize each list of index matches
        for li in range(len(rep_lst)):
            if (len(rep_lst[li])>1):
                print "rep_lst["+str(li)+"] before sort:",rep_lst[li]
                rep_lst[li].sort()
                print "rep_lst["+str(li)+"] after sort :",rep_lst[li]
        #  Step (3): Check within each index (all things that match r, then all things that match s, etc.)
        #            for full negations of the term (alpha to beta, virtual to occupied, etc.)
        for li in range(len(rep_lst)): #Look at all replacement lists
            l = rep_lst[li]
            if len(l)>0: #Actually, only look at the ones that have entries.
                if (dlt(sit_lst[li],l[0])=='0'):
                    # If the first element in the replacement list can't be equal to the
                    # corresponding index, the whole term is zero.
                    print "NEGATION! Term eliminated. One of these is not like the other; sit_lst[li] = ",sit_lst[li],"l = ",l
                    tt = ''
                    break
                
                for ri in range(len(l)-1):
                    # If two elements in the same replacement list can't be equal to each
                    # other, the whole term is zero.
                    if (dlt(l[ri],l[ri+1])=='0'):
                        print "NEGATION! Term eliminated. One of these is not like the other; l = ",l
                        tt = ''
                        break
                for ri in range(len(l)):
                    #No longer need to index over any legit matches. As long as the replacing index is not
                    # itself being replaced, add any indices in the replacement list to the list of indices
                    # NOT to sum over in the final expression.
                    if (l[ri] != sit_lst[li] and (l[ri] not in indremove)):
                        indremove.append(l[ri])
                        print "(match) Adding ",l[ri],"to indremove."
        #Remove unnecessary indices from subscript of summation symbol
        for ind in indremove:
            if ind in sumind:
                sumind.remove(ind)
        #Take those temporary lists of indices that have been modified
        # according to delta function information (see above). Then
        # format the string we plan on printing (tt, in this case-- the single
        # term we've been working on in this loop)
        print "Full delta term:"
        for d in dlt_list:
            print d
        print "(before replacement) tt = ",tt
        #Replace indices in the initial expression with others, according to the delta
        #function information.
        for si in range(len(sit_lst)): #Loop over every index we start with.
            l = rep_lst[si]
            indo = indformat(sit_lst[si]) #Replace every index in this replacement list with its corresponding high-ranking index.
            if (len(l)>0):
                print "indo = "+indo
                for rsi in range(len(l)):
                    indi = "<"+l[rsi]+">"
                    print "indi = "+indi
                    tt = tt.replace(indi,indo)
            else:
                #If there is nothing in the replacement list, then there's a chance
                #that this index has the wrong formatting. Replace it with itself to
                #give it the correct formatting.
                indi = "<"+sit_lst[si]+">"
                tt = tt.replace(indi,indo)
        print "(after replacement)  tt = ",tt
        print "indremove = ",indremove
        #Write summation symbol (if necessary) and then finish up
        if (len(sumind)>0):
            sumind = sort_sumindex(sumind)
            latout += "\sum_{"
            for x in sumind:
                latout += indformat(x)
            latout += "}"
        tt = tt.replace("X_{rr}","X_{GS}") #If nothing matches r, then it just sums over all occupied orbitals.
        if (tt != ""):
            latout += tt + " "
            nterms += 1
        if (nterms%5 == 0):
            latout += "&\\\\\n&"
if (latout[-1] != "&"):
    latout += "&"
latout += "\n\end{flalign*} \n"
print latout
f = open("sf-xcis.tex",'a')
if ("0" not in latout):
    f.write(latout)
f.close()

