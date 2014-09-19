#sq-full.py
#!/usr/local/bin/python
import os.path
import sys
import time
import math
import re

#sq-full: this script is intended to take all of the guesswork
#         out of second quantization calculations by doing all
#         of the logical operations for you and returning an
#         easily understandible output. Unlike sq.py, this
#         script does not rely on hand-computed template delta
#         functions, but rather generates them anew for each
#         calculation. In addition, this script accommodates
#         more types of specialized outputs, and is hopefully
#         better-coded and more reliable in general.

# v1.0: Appears to be mostly bug-free at this time. Concerns:
#       always assumes that you sum over all 'wild' operator
#       indices (i.e., if you use a one-electron operator
#       X_{rs}c_r^+c_s, it assumes you always sum over r and
#       s. This may need to be changed if you only want the
#       result for a single term in the operator matrix.

#CUSTOM SETTINGS
lstate_s = "t_{<i>}^{I<a>}"
rstate_s = "t_{<j>}^{J<b>}"
op_s = "X_{<r><s>}"
oop_s = "X"

terms_per_line = 5
filenamebase = "sf-xcis_new"

#STEP ONE: define a few functions and variables that will
#be used to set the parameters of the problem. Instead of using a
#bunch of if/elif blocks, I'm going to try to define a dictionary
#which references some functions.
tia_re = re.compile("t_\{(?:\<[A-z][0-9]?\>)+\}\^\{[IJ](?:\<[A-z][0-9]?\>)+\}")
gia_re = re.compile("g_\{(?:\<[A-z][0-9]?\>)+\}\^\{(?:\<[A-z][0-9]?\>)+\}")
op_re  = re.compile("[Xni]_\{(?:\<[A-z][0-9]?\>)+\}")
ssq_re = re.compile("\(([A-z][0-9]?[+\-])\)")
tri_re = re.compile("\<([A-z][0-9]?)\>")
triwb_re = re.compile("\<[A-z][0-9]?\>")
cp_re = re.compile("[A-Z]") #Capital letter
nt_re = re.compile("[0-9]+") #Numeral
ind_re = re.compile("[A-z][0-9]?") #Just an index (any character, maybe followed by a single digit
occa_lst = ['i','j','k','l','m','n','s1','s2'] #Occupied, alpha indices
vrta_lst = ['a','b','c','d','e','f']           #Virtual,  alpha indices
occb_lst = ['I','J','K','L','M','N']           #Occupied, beta  indices
vrtb_lst = ['A','B','C','D','E','F','S1','S2'] #Virtual,  beta  indices
occ_lst = list(occa_lst + occb_lst)  #Occupied        indices
vrt_lst = list(vrta_lst + vrtb_lst)  #Virtual         indices
alp_lst = list(occa_lst + vrta_lst)  #Alpha           indices
bet_lst = list(occb_lst + vrtb_lst)  #Beta            indices
wld_lst = ['r','s','t','u']          #Wild card       indices
uni_lst = ['s1','s2','S1','S2']      #Unique indices-- do not sum over these
ind_lst = list(occ_lst+vrt_lst+wld_lst) #All          indices
ind_lst.sort()
key_lst = [] #List of all indices that appear in our initial expression
sin_lst = [] #The list of indices over which to sum in the final output
spr_lst = [] #Same as above, but for prime term (just subtract wild indices)
sus_prm = "" #For the primary term, string to actually put underneath the sum symbol in the final output
sus_lst = [] #For each term, string to actually put underneath the sum symbol in the final output
locc_lst = [] #Lists of primary excitation indices
rocc_lst = []
lvrt_lst = []
rvrt_lst = []
slocc_lst = [] #Equivalent lists for secondary excitation indices
srocc_lst = []
slvrt_lst = []
srvrt_lst = []
#get_ind: General function to ask the user about which
#         indices to use when constructing the function.
#         ex_n:  number of excitations to worry about
#         ex_s:  string that specifies which excitation
#                we're talking about
def get_ind(ex_n,ex_s):
    l = []
    for e in range(ex_n):
        print "Please type the index for occupied orbital "+str(e+1)+" on the "+ex_s+"."
        l.append(raw_input(""))
    for e in range(ex_n):
        print "Please type the index for virtual orbital "+str(e+1)+" on the "+ex_s+"."
        l.append(raw_input(""))
    for ind in l:
        if (ind not in key_lst):
            key_lst.append(ind)
    return l

#construct_state: General function to take index lists
#                 and information about the state to
#                 construct the state term for LaTeX
#                 template creation. Returns a string.
#                 lor:    left (0) or right (1)
#                 pre_ex: number of previous excitations.
#                 o_lst:  list of occupied indices
#                 v_lst:  list of virtual indices
def construct_state(lor,pre_ex,o_lst,v_lst):
    if (pre_ex==0):
        s = "t_{"
    elif (pre_ex==1):
        s = "g_{"
    for ind in o_lst:
        s += "<"+ind+">"
    s += "}^{"
    if (pre_ex==0):
        if (lor==0):
            s += "I"
        elif (lor==1):
            s += "J"
    for ind in v_lst:
        s += "<"+ind+">"
    s += "}"
    return s

#lt_to_sq: Take information from the previously-constructed
#          LaTeX template string, use that to build an appropriate
#          second quantization expression.
def lt_to_sq(lt_s):
    s = ""
    op_lst  = re.findall(op_re ,lt_s)
    tia_lst = re.findall(tia_re,lt_s)
    gia_lst = re.findall(gia_re,lt_s)
    #Left side, first excitations
    temp_lst = re.findall(tri_re,tia_lst[0])
    hnt = len(temp_lst)/2
    for i in range(hnt):
        s += "("+temp_lst[i]    +"+)"
        s += "("+temp_lst[i+hnt]+"-)"
    #Left side, second excitations
    if (len(gia_lst)>0):
        temp_lst = re.findall(tri_re,gia_lst[0])
        hnt = len(temp_lst)/2
        for i in range(hnt):
            s += "("+temp_lst[i]    +"+)"
            s += "("+temp_lst[i+hnt]+"-)"
    #Operator
    temp_lst = re.findall(tri_re,op_lst[0])
    hnt = len(temp_lst)/2
    for i in range(hnt):
        s += "("+temp_lst[i]    +"+)"
    for i in range(hnt):
        s += "("+temp_lst[i+hnt]+"-)"
    #Right side, second excitations
    if (len(gia_lst)>0):
        temp_lst = re.findall(tri_re,gia_lst[1])
        hnt = len(temp_lst)/2
        for i in range(hnt):
            s += "("+temp_lst[i+hnt]+"+)"
            s += "("+temp_lst[i]    +"-)"
    #Right side, second excitations
    temp_lst = re.findall(tri_re,tia_lst[1])
    hnt = len(temp_lst)/2
    for i in range(hnt):
        s += "("+temp_lst[i+hnt]+"+)"
        s += "("+temp_lst[i]    +"-)"
    return s


#STEP TWO: Get user input, use functions and variables defined in
#step one to fully define the problem to be solved.

#key_sort: Organize key_lst so that occupied
#          orbitals come first, followed by
#          virtual, followed by wild card, and
#          the remainder of the sorting is done
#          alphabetically.
def key_sort(l):
    l.sort()
    temp_lst = []
    for ind in l:
        if (ind in occ_lst):
            temp_lst.append(ind)
    for ind in l:
        if (ind in vrt_lst):
            temp_lst.append(ind)
    for ind in l:
        if (ind in wld_lst):
            temp_lst.append(ind)
    return temp_lst

print "Welcome to sq-full, the second-quantization calculator."
print "What kind of states would you like to combine?"
lstate_mode = input("Left state : singles (1), doubles (2), triples (3), singles on singles (4), or custom (0)?\n")
print "You have chosen: ",lstate_mode
rstate_mode = input("Right state: singles (1), doubles (2), triples (3), singles on singles (4), or custom (0)?\n")
print "You have chosen: ",rstate_mode
print "What kind of operator would you like to use?"
op_mode = input("One-electron operator X (1), two-electron operator Pi (2), Fock matrix F (3), or custom (0)?\n")
print "You have chosen: ",op_mode
if (lstate_mode==0 or rstate_mode==0 or op_mode==0):
    print "Please note that all custom inputs are assumed to be already in the script."
    print "If you haven't done this, please exit the script now and edit sq-full.py."

#Set the operator string
if (op_mode==1):
    #Arbitrary one-electron operator
    op_s = "X_{<r><s>}"
    oop_s = "X"
    key_lst.append("r")
    key_lst.append("s")
elif (op_mode==2):
    #Arbitrary two-electron operator
    op_s = "\Pi_{<r><s><t><u>}"
    oop_s = "\Pi"
    key_lst.append("r")
    key_lst.append("s")
    key_lst.append("t")
    key_lst.append("u")
elif (op_mode==3):
    #Fock Matrix (diagonal)
    op_s = "\epsilon_{<r><r>}"
    oop_s = "F"
    key_lst.append("r")

#Set the state coefficient strings
if (lstate_mode>0 and lstate_mode<4):
    temp_l = get_ind(lstate_mode,"left side")
    locc_lst = list(temp_l[:lstate_mode])
    lvrt_lst = list(temp_l[lstate_mode:])
    lstate_s = construct_state(0,0,locc_lst,lvrt_lst)
elif (lstate_mode>3):
    temp_l = get_ind(lstate_mode-3,"left side, initial excitation")
    locc_lst = list(temp_l[:lstate_mode-3])
    lvrt_lst = list(temp_l[lstate_mode-3:])
    lstate_s = construct_state(0,0,locc_lst,lvrt_lst)
    temp_l = get_ind(1,"left side, second excitation")
    slocc_lst += list(temp_l[0])
    slvrt_lst += list(temp_l[1])
    lstate_s += construct_state(0,1,slocc_lst,slvrt_lst)
if (rstate_mode>0 and rstate_mode<4):
    temp_l = get_ind(rstate_mode,"right side")
    rocc_lst = list(temp_l[:rstate_mode])
    rvrt_lst = list(temp_l[rstate_mode:])
    rstate_s = construct_state(1,0,rocc_lst,rvrt_lst)
elif (rstate_mode>3):
    temp_l = get_ind(rstate_mode-3,"right side, initial excitation")
    rocc_lst = list(temp_l[:rstate_mode-3])
    rvrt_lst = list(temp_l[rstate_mode-3:])
    rstate_s = construct_state(1,0,rocc_lst,rvrt_lst)
    temp_l = get_ind(1,"right side, second excitation")
    srocc_lst += list(temp_l[0])
    srvrt_lst += list(temp_l[1])
    rstate_s += construct_state(1,1,srocc_lst,srvrt_lst)

#Write full LaTeX template expression
lt_tmp = lstate_s+op_s+rstate_s

print "lt_tmp:",lt_tmp
input_string = lt_to_sq(lt_tmp)
print "input_string:",input_string
#Figure out which indices to sum over
key_lst = key_sort(key_lst)
print "Currently, this program plans to print the results which sum over the following indices:\n",key_lst
ind_mode = input("Would you like to sum over all of these indices (1) or would you like to sum over a set of indices which you will specify now (2)?\n")
print "You have chosen: ",ind_mode
if (ind_mode==1):
    sin_lst = list(key_lst)
elif (ind_mode==2):
    ind_input = raw_input("Please enter the indices that you would like to sum over, then press return when you are finished:\n")
    print "You have entered: ",ind_input
    sin_lst = re.findall(ind_re,ind_input)
    sin_lst = key_sort(sin_lst)

#Clear out wild card indices from the list of indices to sum
#over (if we end up summing over one of these, it will
#just be over all occupied states, so we'll take care of
#that later (X_{rr} -> X_{GS}). Also, clear out unique
#indices (s1,s2,S1,S2).
#NOTE: This actually gets complicated with two-electron
#      operators. Don't bother.
spr_lst = list(sin_lst)
for indi in range(len(spr_lst)-1,-1,-1):
    if (spr_lst[indi] in wld_lst):
        spr_lst.pop(indi)

#Clear out unique indices (s1,s2,S1,S2).
for indi in uni_lst:
    if (indi in sin_lst):
        sin_lst.remove(indi)

#lt_tmp: The template for each term you ultimately want in the
#        final output. Each index should be surrounded by
#        triangular brackets. In the final output, each of
#        these indices with new ones according to each delta
#        function term ultimately extracted from the
#        input_string second quantization expression.
#        Key: t    = CIS amplitude (single, double, etc.
#             g    = further excitation
#             X    = generic one-electron operator
#             \Pi  = generic two-electron operator
#             \epsilon = Fock operator (diagonal in space of MOs)
#lt_tmp = "t_{<i>}^{I<a>}X_{<r><s>}t_{<j>}^{J<b>}"



#Input_string: this is the second quantization expression that
#you want evaluated. Put terms in the desired order, with each
#operator indexed by a single character following group
#indexing conventions to indicate occupied/virtual/wild orbitals
#(see occ_lst, vrt_lst, wld_lst for details). Also, note that
#capital indices indicate beta orbitals, and lowercase indices
#indicate alpha orbitals. Finally, the formatting for each
#operator should be as follows, surrounded by parentheses:
# (a+)
#where the index is given by a single character, and the type of
#operator is given by a single "+" (creation) or "-"
#(annihilation) character. A full expression might look like
# input_string = "(i+)(a-)(r+)(s-)(b+)(j-)"
#If lt_tmp is defined, the operators will be obtained from that
#expression.
#input_string = "(i+)(a-)(r+)(s-)(b+)(j-)"
#NOTE: the preceeding is obsolete. The variable input_string
#      will now be constructed based on the LaTeX string
#      template variable, lt_tmp.

#STEP THREE: resolve the second-quantization operators, using
#appropriate algebraic operations, into delta function terms.

#This code is intended to automate laborious and tricky second 
# quantization calculations. A few notes:
# (1) Any time these functions talk about a list-formatted term,
#     they are referring to something of the form
#             l = [[1,[[['a','b'],['i','j']],["r+","s-"]]]
#     COUNT THE NUMBER OF BRACKETS! They are there for a reason.
#     All of the terms in the main program are stored in those
#     brackets, and when the list is submitted to a function, it
#     is given its own set of extraneous brackets even though
#     it's on its own, so that if the function needs to return
#     additional terms it is able to. Beyond that, we have three
#     items: (i) (l[0])[0], the sign, an integer, either
#     1 or -1. (ii) (l[0])[1], the delta function list, a list
#     of tuples that represent the delta function coefficients
#     of the term. (iii) (l[0])[2], the creation/annihilation
#     operators, represented as a list of two-character strings,
#     such as "a+", where the first character represents the index
#     and the second character denotes creation/annihilation.
# (2) Any time these functions talk about a list-formatted delta
#     function, it takes the form dl = [[1,[[a,b],[i,j],[r,s]]]]
#     First, there is an outer shell bracket (as with the term
#     list). Then, there is the function itself, which has an
#     integer coefficient and then a list of tuples. The tuples
#     each represent a delta function between those indices.


# is_cre: Takes second quantization operator string as input.
#         Returns 1 (true) if the operator is a creation
#         operator, returns 0 (false) if it is a annihilation
#         operator, and returns -1 otherwise.
def is_cre(s):
    if (s[-1]=="+"):
        return 1
    elif (s[-1]=="-"):
        return 0
    else:
        return -1

# orb_type: Takes second quantization operator string as
#           input. Returns 0 if index is 'wild' (occupied
#           or virtual), 1 if index is occupied, and 2 if
#           index is virtual. Otherwise, return 3 (I guess).
def orb_type(s):
    if (s[:-1] in wld_lst):
        return 0
    elif (s[:-1] in occ_lst):
        return 1
    elif (s[:-1] in vrt_lst):
        return 2
    else:
        print "orb_type: could not identify orbital s:",s
        sys.exit() #Crash
        return 3 #I guess

# term: Takes a standard list-formatted second 
#       quantization term and returns the same (acts as
#       a constructor-- otherwise, trying to copy these
#       lists will just create a pointer to the same
#       list, not an actual copy).
def term(l):
    l2 = []
    for t in l:
        l2.append([0,[],[]])
        (l2[-1])[0] = t[0]
        for d in t[1]:
            ((l2[-1])[1]).append(list(d))
        for s in t[2]:
            ((l2[-1])[2]).append(str(s))
    return l2

# ind_key: Takes an index (string that appears in the
#          list of indices) as an argument, returns
#          its rank in the master list of indices for
#          this particular run (includes only indices
#          used for this calculation, sorts
#          occupied, then virtual, then wild card, each
#          set of which is in alphabetical order.
def ind_key(s):
    print "ind_key: s: ",s
    print "ind_key = ",key_lst.index(s)
    return key_lst.index(s)

# no_repeat: Takes a list, removes all redundant entries.
def no_repeat(l):
    nr = len(l)
    while (nr>0):
        nr = 0
        for item in l:
            n_item = l.count(item)
            if (n_item>1):
                nr += 1                
                for i in range(1,n_item):
                    l.remove(item)
                break
        if (nr == len(l)):
            nr = 0
    return l

# stol: Takes a properly-formatted creation/annihilation
#       operator string as input, returns same in the list
#       format that will be used throughout the script.
def stol(s):
    sq_list = re.findall(ssq_re,s)
    l = []
    l.append([1,[],[]])
    for sq in sq_list:
        ((l[0])[2]).append(sq)
    print "l=",l
    return l

# choose_move: Takes a standard list-formatted second 
#              quantization term and returns the index
#              that should next be switched with its
#              next-nearest neighbor to the right.
#              Failing to find an occupied creation
#              or virtual annihilation to move to the
#              right, or an occupied annihilation or
#              virtual creation to move to the left,
#              the algorithm will return -1.
def choose_move(l):
    #Find something you want to move to the right
    for i in range(len((l[0])[2])-2,-1,-1):
        o = str(((l[0])[2])[i])
        ot = orb_type(o)
        ic = is_cre(o)
        if (ot==1 and ic==1) or (ot==2 and ic==0):
            print "choose_move: ",i
            return i
    #Failing that, find something you want to
    #move to the left.
    for i in range(1,len((l[0])[2])):
        o = str(((l[0])[2])[i])
        ot = orb_type(o)
        ic = is_cre(o)
        if (ot==1 and ic==0) or (ot==2 and ic==1):
            print "choose_move: ",i-1
            return i-1
    #Failing that, crash
    print "choose_move: failed to find a good move"
    sys.exit()
    return -1

# d_clean: Takes a standard list-formatted delta
#          function term and returns an empty
#          list if there are any impossible
#          direct couplings (occ-virt, alp-bet,
#          same excitations acting on each other
#          (single excitation from i followed by
#          another single excitation from i, for
#          example).
def d_clean(dl):
    for d in dl:
        if (d[0] in occ_lst) or (d[1] in occ_lst):
            if (d[0] in vrt_lst) or (d[1] in vrt_lst):
                print "d_clean: cleaning! occupied to virtual! d = ",d
                return []
        if (d[0] in alp_lst) or (d[1] in alp_lst):
            if (d[0] in bet_lst) or (d[1] in bet_lst):
                print "d_clean: cleaning! alpha to beta! d = ",d
                return []
        if (d[0] in uni_lst) or (d[1] in uni_lst):
            if (d[0]!=d[1] and (d[0] not in wld_lst and d[1] not in wld_lst)):
                print "d_clean: cleaning! non-matching unique! d = ",d
                return []
        if (d[0] in locc_lst) or (d[1] in locc_lst):
            if (d[0] in slocc_lst) or (d[1] in slocc_lst):
                print "d_clean: cleaning! double vacancy! d = ",d," slocc_lst = ",slocc_lst
                return []
        if (d[0] in rocc_lst) or (d[1] in rocc_lst):
            if (d[0] in srocc_lst) or (d[1] in srocc_lst):
                print "d_clean: cleaning! double vacancy! d = ",d," srocc_lst = ",srocc_lst
                return []
        if (d[0] in lvrt_lst) or (d[1] in lvrt_lst):
            if (d[0] in slvrt_lst) or (d[1] in slvrt_lst):
                print "d_clean: cleaning! double vacancy! d = ",d," slvrt_lst = ",slvrt_lst
                return []
        if (d[0] in rvrt_lst) or (d[1] in rvrt_lst):
            if (d[0] in srvrt_lst) or (d[1] in srvrt_lst):
                print "d_clean: cleaning! double vacancy! d = ",d," srvrt_lst = ",srvrt_lst
                return []
    return dl

# cleanup: Takes a standard list-formatted second 
#          quantization term and applies the
#          following clean-up operations (in order)
#          (1) deletes terms with occupied creations
#          or virtual annihilations at the rightmost
#          positions, or the corresponding type at the
#          leftmost position (2) deletes terms with
#          impossible delta function coefficients, (3)
#          evaluates terms which only have wild card
#          creation/annihilation pairs, creating a
#          wild-card delta function. Returns 'cleaned-up'
#          term, or empty brackets if term has been
#          deleted.
def cleanup(l):
    # (0) Secret bonus check: is the number of creation
    # operators equal to the number of annihilation operators?
    # If not, delete the term. Also, do the cre/ann operator
    # strings not end with either a "+" or a "-"? Delete.
    nc = 0
    nd = 0
    for o in (l[0])[2]:
        ic = is_cre(o)
        if (ic==1):
            nc += 1
        elif (ic==0):
            nd += 1
        else:
            print "cleanup: can't tell whether it is creation or annihilation operator"
            return []
    if (nc != nd):
        print "cleanup: non-equal numbers of creation and annihilation operators."
        return []
    # (1) Delete terms according to creation/annihilation
    # position (check first and last).
    if (nc>0):
        o = str(((l[0])[2])[-1])
        ot = orb_type(o)
        ic = is_cre(o)
        if (ot==1 and ic==1) or (ot==2 and ic==0):
            print "cleanup: occ. creation or virt. annihilation at rightmost position."
            return []
        o = str(((l[0])[2])[0])
        ot = orb_type(o)
        ic = is_cre(o)
        if (ot==1 and ic==0) or (ot==2 and ic==1):
            print "cleanup: occ. annihilation or virt. creation at leftmost position."
            return []
    # (2) Delete terms with impossible (virt-occ, alpha-
    # beta) direct delta-function coupling. Don't worry
    # about indirect coupling (occ-wild-virt) for now.
    dl = (l[0])[1]
    if (dl!=[]):
        dl = d_clean(dl)
        if (dl==[]):
            return []
    # (3) If only wild-card creation/annihilations
    # remain, create a wild-card delta function to
    # the coefficient.
    if (ca_count(l)==2):
        o0 = str(((l[0])[2])[0])
        o1 = str(((l[0])[2])[1])
        if (not orb_type(o0)) and (not orb_type(o1)):
            if (is_cre(o0) and not is_cre(o1)):
                ((l[0])[1]).append([o0[:-1],o1[:-1]])
                (l[0])[2] = []
    elif (ca_count(l)==4):
        o0 = str(((l[0])[2])[0])
        o1 = str(((l[0])[2])[1])
        o2 = str(((l[0])[2])[2])
        o3 = str(((l[0])[2])[3])
        if (not orb_type(o0) and not orb_type(o1) and not orb_type(o2) and not orb_type(o3)):
            if (is_cre(o0) and is_cre(o1) and not is_cre(o2) and not is_cre(o3)):
                l += term(l) #Make a copy of l (now we have two terms)
                ((l[0])[1]).append([o0[:-1],o3[:-1]]) #Add delta functions
                ((l[0])[1]).append([o1[:-1],o2[:-1]])
                ((l[1])[1]).append([o0[:-1],o2[:-1]])
                ((l[1])[1]).append([o1[:-1],o3[:-1]])
                (l[1])[0] = -(l[1])[0] #The second term changes sign
                (l[0])[2] = [] #Remove all creation/annihilation operators
                (l[1])[2] = []
    return l

# sq_switch: Takes a standard list-formatted second 
#            quantization term, and
#            attempts to switch the two creation/
#            annihilation operators at sensible
#            positions, and all this entails
#            (creating additional terms, altering
#            sign, and ultimately making the switch).
def sq_switch(l):
    n = choose_move(l)
    if (n==-1):
        return l
    o0 = str(((l[0])[2])[n])
    o1 = str(((l[0])[2])[n+1])
    #Case 1: both are creation/annihilation operators.
    #        Flip the sign of the term and exchange
    #        operator positions.
    if (is_cre(o0)==is_cre(o1)):
        (l[0])[0] = -(l[0])[0]
        ((l[0])[2])[n] = o1
        ((l[0])[2])[n+1] = o0
    #Case 2: both are not the same type of operator.
    #        Do the same thing as before, but also
    #        create a second term in which the two
    #        operators are removed and the coefficient
    #        contains an extra delta function.
    else:
        l += term(l) #Make a copy of l (now we have two terms)
        ((l[0])[1]).append([o0[:-1],o1[:-1]])
        ((l[0])[2]).pop(n+1)
        ((l[0])[2]).pop(n)
        (l[1])[0] = -(l[0])[0]
        ((l[1])[2])[n] = o1
        ((l[1])[2])[n+1] = o0
    return l

# ca_count: Takes a standard list-formatted second 
#           quantization term, and counts the number
#           of creation/annihilation operators still
#           remaining. Returns this number.
def ca_count(l):
    return len((l[0])[2])

# hi_coup: Takes a standard delta function list,
#          checks for higher-order couplings, e.g.
#          D(a)(r)D(r)(b) (a is indirectly coupled
#          to b). Returns modified delta function
#          list (e.g., earlier example would become
#          D(a)(b)D(a)(r)). Also removes redundant
#          delta functions and eliminates terms if
#          appropriate (e.g., for D(a)(r)D(j)(r)).
def hi_coup(dl):
    print "hi_coup: dl (1): ",dl
    print "hi_coup: ind_lst: ",ind_lst
    # (1) Convert delta function list into a non-
    #     pairwise list so we can sort everything
    #     out.
    mat_llst = []
    for d in (dl[0])[1]:
        matched = False
        for li in range(len(mat_llst)):
            if ((d[0] in mat_llst[li]) and (d[1] in mat_llst[li])):
                matched = True
                print "hi_coup: full match found. d:",d,", mat_llst[li]:",mat_llst[li]
                break
            elif (d[0] in mat_llst[li]):
                print "hi_coup: partial match found (0). d:",d,", mat_llst[li]:",mat_llst[li]
                mat_llst[li].append(d[1])
                matched = True
                break
            elif (d[1] in mat_llst[li]):
                print "hi_coup: partial match found (1). d:",d,", mat_llst[li]:",mat_llst[li]
                mat_llst[li].append(d[0])
                matched = True
                break
        if (not matched):
            print "hi_coup: no match found. d:",d,", mat_llst:",mat_llst
            mat_llst.append(d)
    print "mat_llst (1):",mat_llst
    # (1.5) Loop over lists to look for further matches
    #       until no more are found.
    nmatches = 1
    while (nmatches>0):
        nmatches = 0
        print "hi_coup: mat_llst:",mat_llst
        for li in range(len(mat_llst)):
            for ind in mat_llst[li]:
                print "hi_coup: ind:",ind
                for mli in range(li+1,len(mat_llst)-1):
                    if (ind in mat_llst[mli]):
                        nmatches += 1
                        print "hi_coup: further match found. mat_llst[li]:",mat_llst[li],", mat_llst[mli]:",mat_llst[mli]
                        mat_llst[li] += list(mat_llst[mli])
                        mat_llst.pop(mli)
                        break
                if (nmatches>0): break
            if (nmatches>0): break
    # (2) Clean up: remove redundancies and get lists
    #     in proper order.
    for li in range(len(mat_llst)):
        mat_llst[li] = no_repeat(mat_llst[li])
        mat_llst[li].sort(key=ind_key)
    print "mat_llst (2):",mat_llst
    # (3) Use this list to create
    #     standardized set of delta functions (if
    #     we have [a,b,r,s], create delta
    #     function list [[a,b],[a,s],[a,r]]).
    (dl[0])[1] = []
    for i in range(len(mat_llst)):
        iind = (mat_llst[i])[0]
        print "hi_coup: mat_llst[i]:",mat_llst[i]
        print "hi_coup: iind:",iind
        for j in range(1,len(mat_llst[i])):
            jind = (mat_llst[i])[j]
            print "hi_coup: jind:",jind
            ((dl[0])[1]).append([iind,jind])
            print "hi_coup: appending ["+iind+","+jind+"]"
    print "hi_coup: dl (2): ",dl
    # (4) Apply standard cleaning procedure to
    #     the new list of delta functions. The
    #     list format should remove any
    #     problematic indirect couplings.
    (dl[0])[1] = d_clean((dl[0])[1])
    print "hi_coup: dl (3): ",dl
    return dl

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
    for i in range(len(cp_lst)):
        s = s.replace(cp_lst[i],"\\bar{"+cp_lst[i].lower()+"}")
    for i in range(len(nt_lst)):
        s = s.replace(nt_lst[i],"_{"+nt_lst[i]+"}")
    return s

# compile_sumind: take a list of indices, compile them into
#                 the actual string that will appear under
#                 the summation symbol in the final output.
#                 Also, takes an integer argument describing
#                 the type of sum; 0 for non-prime (i.e.,
#                 right side of the equals sign) and 1 for
#                 the prime term (left side of equals sign).
def compile_sumind(l,isp):
    s = ""
    s_final = "\substack{"
    for indi in range(len(l)):
        ind = l[indi]
        s += ind
        if (not isp):
            if (ind in locc_lst and len(slocc_lst)>0):
                s += "\\neq "
                for sind in slocc_lst:
                    s += sind
            elif (ind in rocc_lst and len(srocc_lst)>0):
                s += "\\neq "
                for sind in srocc_lst:
                    s += sind
            elif (ind in lvrt_lst and len(slvrt_lst)>0):
                s += "\\neq "
                for sind in slvrt_lst:
                    s += sind
            elif (ind in rvrt_lst and len(srvrt_lst)>0):
                s += "\\neq "
                for sind in srvrt_lst:
                    s += sind
        s += ","
        if (len(s)>8 and indi<len(l)-1):
            s_final += s +"\\\\"
            s = ""
    s_final += s
    if (s_final[-1]==","):
        s_final = s_final[:-1]
    return s_final+"}"
    

# Take input string, format it as a proper list (generalized formatting
# that will be used for the rest of the program).
mysqt_l = stol(input_string)
print "mysqt_l: ",mysqt_l
print "len(mysqt_l): ",len(mysqt_l)
delta_l = [] #Final list of delta function terms
# Perform second quantization logic operations until the creation and
# annihilation operators have all been converted into delta functions.
while(len(mysqt_l)>0):
    #Apply a logic operation to each term, append the results to the
    #end of the list and then delete the term.
    for mysq_i in range(len(mysqt_l)-1,-1,-1):
        mysq_t = mysqt_l[mysq_i]
        print "mysqt_l before switch: ",mysqt_l
        mysqt_l += sq_switch([mysq_t])
        mysqt_l.pop(mysq_i)
        print "mysqt_l after switch:  ",mysqt_l
    #Clean up the resulting terms
    for mysq_i in range(len(mysqt_l)-1,-1,-1):
        mysq_t = mysqt_l[mysq_i]
        print "mysqt_l before cleanup: ",mysqt_l
        mysqt_l += cleanup([mysq_t])
        mysqt_l.pop(mysq_i)
        print "mysqt_l after cleanup:  ",mysqt_l
    #Finally, remove terms which no longer have creation/annihilation
    #operators, and move the delta functions to the delta function list.
    for mysq_i in range(len(mysqt_l)-1,-1,-1):
        mysq_t = mysqt_l[mysq_i]
        if (ca_count([mysq_t])==0):
            delta_l.append(mysq_t[0:2])
            mysqt_l.pop(mysq_i)

# At this point, all of the creation/annihilation operators have been
# resolved into delta functions.
print "(before) delta_l: ",delta_l

# Check for higher-order delta function couplings for each term, and
# standardize delta function format.
for ti in range(len(delta_l)-1,-1,-1):
    delta_l += hi_coup([delta_l[ti]])
    delta_l.pop(ti)
delta_l.reverse()
print "(after)  delta_l: ",delta_l

# Combine standard output template with delta function and sign
# information to create list of standard outputs.
# out_lst format: list of two lists. The first sub-list has sign
#                 information for the terms. The second list has
#                 the remainder of the properly-formatted LaTeX
#                 terms (take LaTeX term template, substitute
#                 using delta functions).
out_lst = [[],[]]
print "lt_tmp:",lt_tmp
sus_prm = compile_sumind(spr_lst,1) #Compile the string that will
                                    #appear underneath the summation
                                    #symbol for the primary term.
for t in delta_l:
    tmp_sin_lst = list(sin_lst)
    (out_lst[0]).append(t[0])
    tlt_tmp = str(lt_tmp)
    for d in t[1]:
        print "tlt_tmp (before): d: ",d,", tlt_tmp:",tlt_tmp
        tlt_tmp = tlt_tmp.replace("<"+d[1]+">",indformat(d[0]))
        print "tlt_tmp (after) : d: ",d,", tlt_tmp:",tlt_tmp
        if (d[1] in tmp_sin_lst):
            tmp_sin_lst.remove(d[1])
    sus_lst.append(compile_sumind(tmp_sin_lst,0)) #Compile the string that will
                                                #appear underneath the summation
                                                #symbol for this term.
    tmp_lst = re.findall(triwb_re,tlt_tmp)
    tmp_lst = no_repeat(tmp_lst)
    for s in tmp_lst:
        sl = re.findall(tri_re,s)
        sp = sl[0]
        print "tlt_tmp (before): s: ",s,", tlt_tmp:",tlt_tmp
        tlt_tmp = tlt_tmp.replace(s,indformat(sp))
        print "tlt_tmp: replacing ",s," with ",indformat(sp)
        print "tlt_tmp (after) : s: ",s,", tlt_tmp:",tlt_tmp
    (out_lst[1]).append(str(tlt_tmp))

for ti in range(len(out_lst[1])):
    (out_lst[1])[ti] = re.sub(r"\\eps","\eps",(out_lst[1])[ti])
    (out_lst[1])[ti] = re.sub(r"\\Pi","\Pi",(out_lst[1])[ti])
print "out_lst (1):",out_lst
#Check for matching terms. If they are found, put them together
#by changing their numerical coefficients.
for ti in range(len(out_lst[1])-1,0,-1):
    for cti in range(ti-1,-1,-1):
        if ((out_lst[1])[ti]==(out_lst[1])[cti]):
            (out_lst[0])[cti] = (out_lst[0])[cti] + (out_lst[0])[ti]
            (out_lst[0])[ti] = 0
            print "Found like terms! Putting them together. (out_lst[1])[ti] =",(out_lst[1])[ti],"(out_lst[1])[cti] =",(out_lst[1])[cti]
print "out_lst (2):",out_lst
#Eliminate terms with a coefficient of zero.
for ti in range(len(out_lst[1])-1,0,-1):
    if ((out_lst[0])[ti]==0):
        print "Pop goes the zero term."
        (out_lst[0]).pop(ti)
        (out_lst[1]).pop(ti)
        sus_lst.pop(ti)
print "out_lst (3):",out_lst


#STEP FOUR: Take these delta function terms and compile them into
#a LaTeX-formatted string for convenient reading.

latout = "\\begin{flalign*}\n&\sum_{"+indformat(sus_prm)+"}"+lstate_s+"\langle\Phi_{"
for ind in locc_lst:
    latout += indformat(ind)
if (len(slocc_lst)>0):
    latout += ";"
    for ind in slocc_lst:
        latout += indformat(ind)
latout += "}^{"
for ind in lvrt_lst:
    latout += indformat(ind)
if (len(slvrt_lst)>0):
    latout += ";"
    for ind in slvrt_lst:
        latout += indformat(ind)
latout += "}\\vert "+oop_s+"\\vert\Phi_{"
for ind in rocc_lst:
    latout += indformat(ind)
if (len(srocc_lst)>0):
    latout += ";"
    for ind in srocc_lst:
        latout += indformat(ind)
latout += "}^{"
for ind in rvrt_lst:
    latout += indformat(ind)
if (len(srvrt_lst)>0):
    latout += ";"
    for ind in srvrt_lst:
        latout += indformat(ind)
latout += "}\\rangle "+rstate_s+" ="
latout = latout.replace(">","").replace("<","")
if (len(out_lst[0])==0):
    print ("len(out_lst[0])==0. out_lst:",out_lst)
    latout += "0"
else:
    nterms = 1
    for ti in range(len(out_lst[0])):
        nterms += 1
        if (nterms>2 and (out_lst[0])[ti] > 0):
            latout += " +"
        elif ((out_lst[0])[ti] < 0):
            latout += " -"
        if (abs((out_lst[0])[ti]) != 1):
            latout += str(abs((out_lst[0])[ti]))
        if (len(sus_lst[ti])>0):
            latout += "\sum_{"+indformat(sus_lst[ti])+"}"
        latout += str((out_lst[1])[ti])
        if (nterms%terms_per_line == 0):
            latout += "&\\\\\n&"
if (latout[-1] != "&"):
    latout += "&"
latout += "\n\end{flalign*} \n"
#for ind in wld_lst:
#    latout = latout.replace("{"+ind+ind+"}","{GS}")
for ind in key_lst:
    latout = latout.replace("epsilon_{"+ind+ind+"}","epsilon_{"+ind+"}")
print latout

f = open(filenamebase+'.tex','a')
if ("0" not in latout):
    f.write(latout)
f.close()
