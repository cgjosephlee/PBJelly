#!/usr/bin/env python
import sys

if __name__ == '__main__':
    fh = open(sys.argv[1],'r')
    fh.readline()#args
    h = fh.readline()[1:]
    header = {}
    for pos, item in enumerate(h.split(' ')):
        header[item] = pos
        
    for line in fh.readlines():
        data = line.strip().split()
        print "{0}\t{1}\t{2}\tHT{3}".format(data[header["chrom"]], \
                        data[header["uBreak"]], data[header["dBreak"]], \
                        data[header["id"]])
    fh.close()
