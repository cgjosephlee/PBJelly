#!/usr/bin/env python
import sys, argparse
from collections import defaultdict
from pbsuite.utils.FileHandlers import FastqFile


USAGE = "Replace single pass reads from a ZMW with CCS read when created"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument("filtered_subreads", type=str, \
                        help="Fastq of single pass reads")
    parser.add_argument("ccs_reads", type=str, \
                        help="Fastq of ccs reads")
    
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output fastq file (STDOUT)")
    args = parser.parse_args()

    sub = FastqFile(args.filtered_subreads)
    ccs = FastqFile(args.ccs_reads)
    
    cKeys = ccs.keys()
    if args.output != None:
        output = open(args.output,'w')
    else:
        output = sys.stdout
        
    #name: numBases
    ccsReads = defaultdict(int)
    subReads = {}
    ccsPases = defaultdict(int)
    for read in sub:
        ccsKey = "/".join(read.split('/')[:2])
        try:
            seq = ccs[ccsKey]
            ccsPases[ccsKey] += 1
            subReads[read] = len(sub[read].seq)
            if ccsReads[ccsKey] > 0:
                continue
            ccsReads[ccsKey] = len(sub[read].seq)
        except KeyError:
            seq = sub[read]
        output.write(seq.toString())
    
    output.close()

    def medianLen(lst):
        lst.sort()
        p = len(lst)/2
        return lst[p]
    a = len(ccsReads.keys())
    sys.stderr.write("+CCS reads    : %d\n" % (a))
    b = sum(ccsReads.values())
    sys.stderr.write("+CCS bases    : %d\n" % (b))
    c = len(subReads.keys())
    sys.stderr.write("-SUB reads    : %d\n" % (c))
    d = sum(subReads.values())
    sys.stderr.write("-SUB bases    : %d\n" % (d))
    sys.stderr.write("Avg SUB Len   : %d\n" % (int(sum(subReads.values()))/float(len(subReads))))
    sys.stderr.write("Avg CCS Len   : %d\n" % (int(sum(ccsReads.values()))/float(len(ccsReads))))
    sys.stderr.write("Avg SUB / CCS : %.3f\n" % (sum(ccsPases.values())/float(len(ccsPases))))
    
    #sys.stderr.write("Med CCS RdLen : %.2f\n" % (medianLen(ccsReads.values())))
    #sys.stderr.write("Med Sub RdLen : %.2f\n" % (medianLen(subReads.values())))
