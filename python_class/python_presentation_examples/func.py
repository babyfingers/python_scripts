#!/usr/local/bin/python

def is_bigger_word(w1,w2):
    print "Word 1 = ",w1
    print "Word 2 = ",w2
    if (len(w1)>len(w2)):
        return True
    else:
        return False

el1 = "einsteinium"
el2 = "rutherfordium"
if (is_bigger_word(el1,el2)):
    print "Einstein is the smartest"
else:
    print "Einstein is pretty smart"



