#!/usr/bin/env python
import sys, math, json
from collections import defaultdict
from optparse import OptionParser   

USAGE = """Usage: %prog <input.bed> <output.wig> [--options]
Converts a bed file into a wig coverage plot.
Use '-' as input.bed to read standard in"""

if __name__ == '__main__':
    parser = OptionParser(usage=USAGE)

    parser.add_option("-t", "--thick", \
        help="Only Use Thick Starts Coordinates for coverage", 
        action="store_true")

    opts, args = parser.parse_args()

    try:
        input = args[0]
        output = args[1]
    except IndexError:
        parser.error("Invalid number of arguments!")
    
    fh = open(input,'r') if input != '-' else sys.stdin
    bounds = {}#[minstart, minend]
    coords = defaultdict(list)#[(start,end),...]
    
    for line in fh.readlines():
        data = line.split('\t')
        key = data[0]

        if opts.thick == True:
            start = int(data[6]); stop = int(data[7])
        else:
            start = int(data[1]); stop = int(data[2])
        
        if start < 0 or stop < 0:
            continue
        
        coords[key].append((start,stop))
        
        if bounds.has_key(key):
            bounds[key][0] = min(start, bounds[key][0])
            bounds[key][1] = max(stop,  bounds[key][1])
        else:
            bounds[key] = [start, stop]
    fh.close()
    
    coverage = {}
    for key in bounds:
        #abs is wrong
        coverage[key] = [0]*(abs(bounds[key][1]) - abs(bounds[key][0]))
        lBound = bounds[key][0]
        for start,end in coords[key]:
            start -= lBound; end -= lBound
            coverage[key][start:end] = \
                [coverage[key][i]+1 for i in xrange(start, end)]
    
    del(coords) #don't need
    print json.dumps(coverage)
    fout = open(output, 'w') if output != '-' else sys.stdout
    #collapsing
    for key in coverage:
        start = bounds[key][0]
        end = start + 1
        curval = coverage[key][0]

        for base in coverage[key][1:]:
            if base != curval:
                fout.write("%s\t%d\t%d\t%.4f\n" % (key, start, end, curval))
                start = end
                end = start + 1
                curval = base
            else:
                end += 1
        
    fout.close()
