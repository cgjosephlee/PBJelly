#!/usr/bin/env python
import sys
import pysam
from pbsuite.honey.Force import *

fh = open(sys.argv[1])
fh.readline(); fh.readline() #Header

bam = pysam.Samfile(sys.argv[2])

#400bp around bp must be spanned
def bpCheck(bam, chrom, point, BUFFER=100):
    nSpan = 0
    nCount = 0
    for read in bam.fetch(reference=chrom, start = max(0, point-BUFFER), end=point+BUFFER):
        if read.pos <= point-BUFFER and point+BUFFER <= read.aend:
            nSpan += 1
        nCount += 1
        #elif read.pos < point and point < read.aend:
        #nCount += 1
    return nSpan, nCount
    
print "id\tnReads\tnZMWs\tavgCov\tpctVar\tSpan1\tCoverage1\tPct1\tSpan2\tCoverage2\tPct2"
for line in fh.readlines():
    data = line.strip().split("\t")
    chr1 = data[2]
    bp1 = int(data[3])
    chr2 = data[5]
    bp2 = int(data[6])

    nReads = int(data[10])
    nZMWs = int(data[11])
    
    nSpan1, nCount1 = bpCheck(bam, chr1, bp1)
    nSpan2, nCount2 = bpCheck(bam, chr2, bp2)
    
    avgCov = (nCount1 + nCount2) / 2
    
    
    print "\t".join([str(x) for x in [data[0], \
                                      nReads, \
                                      nZMWs, \
                                      avgCov, \
                                      "%.3f" % (float(nReads)/avgCov), \
                                      nSpan1, \
                                      nCount1, \
                                      "%.3f" % (float(nSpan1)/nCount1), \
                                      nSpan2, \
                                      nCount2, \
                                      "%.3f" % (float(nSpan2)/nCount2)]])
    
