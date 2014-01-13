#!/usr/bin/env python
import sys
import argparse

USAGE = "Turn .spots results into a .bed"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=USAGE, \
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", metavar="SPOTS", \
                        help="Results to convert")
    parser.add_argument("-b", "--brief", action="store_true",\
                        help="Only output columns 1-4 of .bed")
    
    args = parser.parse_args()
    
    fh = open(args.input,'r')
       
    for line in fh.readlines():
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            h = line.strip()[1:]
            header = {}
            for pos, item in enumerate(h.split('\t')):
                header[item] = pos
            continue
        chrom, os, s, ins, ine, e, oe, info = line.strip().split('\t')
        data = {}
        for i in info.split(';'):
            key,val = i.split('=')
            data[key] = val
        d = filter(lambda x: x.startswith('label'), info.split(';'))[0].split('=')[1]
        
        if args.brief:
            print "\t".join([chrom, os, oe, d])
        else:
            print "\t".join([chrom, os, oe, d, '0', '+', ins, ine])
    fh.close()
            
            

