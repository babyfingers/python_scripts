#!/usr/local/bin/python

f = open('text.in','r')
for line in f:
    print line,
f.seek(0)
f_string = f.read()
f.close()

f_string = f_string.replace('element','helium')
print "****************"
print f_string

f = open('text.out','w')
f.write(f_string)
f.close()
