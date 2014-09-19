#!/usr/bin/env python
# Simple Daikon-style invariant checker
# Andreas Zeller, May 2012
# Complete the provided code around lines 28 and 44
# Do not modify the __repr__ functions.
# Modify only the classes Range and Invariants,
# if you need additional functions, make sure
# they are inside the classes.

import sys
import math
import random

def square_root(x, eps = 0.00001):
    assert x >= 0
    y = math.sqrt(x)
    assert abs(square(y) - x) <= eps
    return y
    
def square(x):
    return x * x

# The Range class tracks the types and value ranges for a single variable.
class Range:
    def __init__(self):
        self.min  = None  # Minimum value seen
        self.max  = None  # Maximum value seen
    # Invoke this for every value
    def track(self, value):
        if (self.min == None):
            self.min = value
        elif (self.min > value):
            self.min = value
        if (self.max == None):
            self.max = value
        elif (self.max < value):
            self.max = value
    def __repr__(self):
        return repr(self.min) + ".." + repr(self.max)

class Var_type:
    def __init__(self):
        self.example  = None  # Just a single example        
    # Invoke this for every value
    def track(self, value):
        self.example = value
    def __repr__(self):
        return repr(self.example)

class Set_pattern:
    def __init__(self):
        self.set  = []  # All possible values stored here
    # Invoke this for every value
    def track(self, value):
        if (value not in self.set):
            self.set.append(value)
        self.set = sorted(self.set, key=abs)
    def __repr__(self):
        return repr(self.set)

class Relations:
    def __init__(self):
        self.diffrange = Range()
        self.diffmin = None
        self.diffmax = None
    # Invoke this for every pair of values
    def track(self, v1, v2):
        if (v1 != None and v2 != None):
            diff = v1 -v2
            self.diffrange.track(diff)
            self.diffmin = self.diffrange.min
            self.diffmax = self.diffrange.max
    def __repr__(self):
        return repr(self.diffrange) +".."+ repr(self.diffrange)
                                                
    
# The Invariants class tracks all Ranges for all variables seen.
class Invariants:
    def __init__(self):
        # Mapping (Function Name) -> (Event type) -> (Variable Name)
        # e.g. self.vars["sqrt"]["call"]["x"] = Range()
        # holds the range for the argument x when calling sqrt(x)
        self.vars = {}
        
    def track(self, frame, event, arg):
        event_l = ["call", "return"]
        vars_dict = frame.f_locals.copy()
        if (event == "return"):
            vars_dict["ret"] = arg
        codename = frame.f_code.co_name
        if (event in event_l):
            if (codename not in self.vars):
                self.vars[codename] = {}
            if (event not in self.vars[codename]):
                self.vars[codename][event] = {}
            for var1 in vars_dict:
                if (var1 not in self.vars[codename][event]):
                    self.vars[codename][event][var1] = {}
                    self.vars[codename][event][var1][var1] = {}
                    self.vars[codename][event][var1][var1]["type"] = Var_type()
                    self.vars[codename][event][var1][var1]["range"] = Range()
                    self.vars[codename][event][var1][var1]["set_pattern"] = Set_pattern()
                    for var2 in vars_dict:
                        if (var1 != var2 and var2 not in self.vars[codename][event][var1]):
                            self.vars[codename][event][var1][var2] = {}
                            self.vars[codename][event][var1][var2]["relations"] = Relations()
            for var1,val1 in vars_dict.iteritems():
                self.vars[codename][event][var1][var1]["type"].track(val1)
                self.vars[codename][event][var1][var1]["range"].track(val1)
                self.vars[codename][event][var1][var1]["set_pattern"].track(val1)
                for var2,val2 in vars_dict.iteritems():
                    if (var1 != var2):
                        self.vars[codename][event][var1][var2]["relations"].track(val1,val2)
    def __repr__(self):
        # Return the tracked invariants
        s = ""
        prop_l = ["type", "set_pattern", "range"]
        for function, events in self.vars.iteritems():
            for event, var1s in events.iteritems():
                s += event + " " + function + ":\n"
                for var1 in var1s:
                    #for pname, property in var1s[var1][var1].iteritems():
                    for pname in prop_l:
                        s += "    assert "
                        property = var1s[var1][var1][pname]
                        if (pname == "range"):
                            if property.min == property.max:
                                s += var1 + " == " + repr(property.min)
                            else:
                                s += repr(property.min) + " <= " + var1 + " <= " + repr(property.max)
                        elif (pname == "set_pattern"):
                            s += var1 +" in set(["
                            for i in range(len(property.set)):
                                s += str(property.set[i])
                                if (i < len(property.set)-1):
                                    s += ", "
                                else:
                                    s += "])"
                        elif (pname == "type"):
                            s += "isinstance("+var1+", type(" +repr(property.example) +"))"
                        s += "\n"
                    for var2 in var1s[var1]:
                        if (var1 != var2 and "relations" in var1s[var1][var2]):
                            diffmin = var1s[var1][var2]["relations"].diffmin
                            diffmax = var1s[var1][var2]["relations"].diffmax
                            if (diffmin == diffmax == 0.0):
                                s += "    assert " +var1 +" == " +var2 +"\n"
                            elif (diffmin >= 0.0 and diffmax >= 0.0):
                                s += "    assert " +var1 +" >= " +var2 +"\n"
                            elif (diffmin <= 0.0 and diffmax <= 0.0):
                                s += "    assert " +var1 +" <= " +var2 +"\n"
        return s

invariants = Invariants()
    
def traceit(frame, event, arg):
    invariants.track(frame, event, arg)
    return traceit

sys.settrace(traceit)
# Tester. Increase the range for more precise results when running locally
eps = 0.000001
#test_l = [3, 0, -10]
for i in range(1, 10):
    r = int(random.random() * 1000) # An integer value between 0 and 999.99
#for i in range(0,3):
#    r = test_l[i]
    z = square_root(r, eps)
    z = square(z)
sys.settrace(None)
print invariants


# Example sample of a correct output:
#
#return double:
#    assert isinstance(x, type(-293438.402))
#    assert x in set([9.348, -293438.402, 34.6363])
#    assert -293438.402 <= x <= 34.6363
#    assert x <= ret
