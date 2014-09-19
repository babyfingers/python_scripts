# test.py
#!usr/local/bin/python

# test.py: testing out a few things.

def sayhi():
    """Just say hi to the nice folks."""
    print "Hi, folks!"

sayhi()
print "sayhi.__name__ = "+sayhi.__name__
if (sayhi.__module__):
    print "sayhi.__module__ = "+sayhi.__module__
if (sayhi.__doc__):
    print "sayhi.__doc__ = "+sayhi.__doc__
if (sayhi.func_code):
    print "sayhi.func_code = "+str(sayhi.func_code)
exec sayhi.func_code
