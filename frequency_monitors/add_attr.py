#!/usr/bin/env python
import h5py

# Modified from https://gist.github.com/dideler/2395703
def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0] == '-fname':  # Found a "-name value" pair.
            opts['fname'] = argv[1]  # Add key and value to the dictionary.
            argv = argv[1:]
        if argv[0] == '-attr':
            opts[argv[1]] =  float(argv[2])
            argv = argv[2:]
           # Reduce the argument list by copying it starting from index 1.
        else:
            argv = argv[1:]
    return opts

def add_attributes(fname, attributes):
	# Attributes is a dictionary
	try:
		f = h5py.File(fname, 'r+')
	except IOError:
		print('Error: Tried to open non-existent file.')
		return

	for attr in attributes:
		f.attrs.create(attr, attributes[attr])

	f.close()

if __name__ == '__main__':
    from sys import argv
    myargs = getopts(argv)

    if 'fname' in myargs:
        fname = myargs['fname']
    else:
        print('Error: No input file provided!')

    myargs.pop('fname', None)
    add_attributes(fname, myargs)
