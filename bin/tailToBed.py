#!/usr/bin/env python
import sys

if __name__ == '__main__':
    fh = open(sys.argv[1],'r')
    fh.readline()#args
    h = fh.readline()[1:].strip()
    header = {}
    for pos, item in enumerate(h.split('\t')):
        exec("%s=%d" % (item, pos))
    
    for line in fh.readlines():
        data = line.strip().split()
        print "{0}\t{1}\t{2}\tHT{3}.{4}".format(data[uRef], \
                        data[uBreak], data[dBreak], \
                        data[id], data[annot])
    fh.close()
