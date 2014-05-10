#!/usr/bin/env python
import sys, argparse, re
from collections import defaultdict, namedtuple, Counter
import pysam
import pbsuite.honey.TGraf as tails
import pbsuite.honey.HSpots as spots
from pbsuite.utils.setupLogging import *

USAGE = """
Checks if there are any reads in the given .bam that support predicted SVs

Takes a .bed with structure
    chrom  start  end  name  svtype size

- svtype must be one of DEL, INS, MIS
- size is what the SV's size is estimated to be.  
- DEL and MIS size should equal the sv's span (end - start). INS is the 
    number of inserted bases

RegionBuffer is the +- space in which you consider reads for support 
    around predicted sv
Half of region buffer gives the area where reads must land to support 
    predicted sv.
"""

class Variant():
    def __init__(self, chrom, start, end, svtype, size):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.svtype = svtype
        self.size = size
        #sequence eventually
    
    def __str__(self):
        return "%s:%d-%d(%s)%d" % (self.chrom, self.start, self.end, self.svtype, self.size)
    
def parseArgs(args):
    parser = argparse.ArgumentParser(prog="Honey.py force", description=USAGE, \
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", metavar="BAM", type=str, \
                        help="Assembled Contigs Bam")
    parser.add_argument("bed", metavar="BED", type=str, \
                        help="Bed of locations to force SV Calls")
    parser.add_argument("-s", "--sizebuffer", type=float, default=0.35, \
                        help=("Buffer of estimated sv size to "
                              "create match (%(default)s)"))
    parser.add_argument("-r", "--regionbuffer", type=int, default=400, \
                        help="Buffer sv region prediction (%(default)s)")
    parser.add_argument("-o", "--overlapbuffer", type=float, default=0.50, \
                        help="Percent overlap required from calls to tails (%(default)s)")
    parser.add_argument("-m", "--minIns", type=int, default=5, \
                        help="Minimum insertion error size to consider (%(default)s)")
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    args = parser.parse_args(args)
    setupLogging(args.debug)
    
    return args

def tailsSearch(bam, bed, args):
    """
    Populate the answer dictionary by looking for tails 
    through the bam

    Returns a list of pbsuite.honey.TGraf.Bnode that support
    """
    
    chrom, start, end, name, svtype, estsize = bed
    
    fetchS = max(0, start - args.regionbuffer)
    fetchE = min(end + args.regionbuffer, bam.lengths[bam.references.index(chrom)])
    points = tails.makeBreakReads(bam.fetch(chrom, fetchS, fetchE), getrname = bam.getrname)
    reads = []
    for key in points:
        #eventually will need tloc work
        if key.split('_')[0] != chrom:
            continue
        #eventually will need a reference allele check for tails
        for read in points[key]:
            anno = read.annotate()
            if anno in ['TLOC', 'INV']:
                anno = 'MIS'
            
            #TLOCs...
            if chrom != bam.getrname(read.read.tid):
                continue
            
            if anno != svtype:
                #Not perfect..
                continue
            
            #within reciprocal ovl
            maxStart = max(start, read.uBreak)
            minEnd   = min(end, read.dBreak)
            if minEnd <= maxStart: #No overlap
                continue
            maxSpan  = max(end-start, read.dBreak - read.uBreak)
            recipOvl = end-start / float(maxSpan)
            if recipOvl < args.overlapbuffer:#not enough overlap
                continue
            
            reads.append(read)
    ret = Counter()
    for i in reads:
        ret[i.annotate()] += 1
    k = ret.keys()
    k.sort()
    return ["%s:t:%d" % (x, ret[x]) for x in k]

    #return reads

def spotsSearch(bam, bed, args):
    """
    take a pysam.Samfile and fetch reads in chrom/start/end region
    see if any reads support the call

    But this doesn't take into account that I have specific groupIds to use...
    """
    chrom, start, end, name, svtype, size = bed
    
    fetchS = max(0, start - args.regionbuffer)
    fetchE = min(end + args.regionbuffer, bam.lengths[bam.references.index(chrom)])
    ret = Counter()
    for read in bam.fetch(chrom, fetchS, fetchE):
        if read.pos > start or read.aend < end:
            #Not spanning our region
            logging.debug("%s doesn't span region" % (read.qname))
            continue
        
        #I'm going to need md if I get good a MIS
        cigar = spots.expandCigar(read.cigar)
        regionStart = start - (args.regionbuffer / 2)
        regionEnd = end + (args.regionbuffer /2)
        readPosition = read.pos
    
        numDel = 0
        numRef = 0
        before = False #need to see if we span
        after = False

        numIns = 0
        pins = False
        pinsz = 0
        curMd = 0
        def pinsLoad(start, size):
            if not pins:
                return False, 0, 0
            rs = size if size >= args.minIns else 0
            return False, 0, rs
        
        for code in cigar:
            #make sure we're spanning
            if readPosition < regionStart:
                before = True
                if code != 1: 
                    readPosition += 1
                continue
            if readPosition >= regionEnd:
                after = True
                if code != 1: 
                    readPosition += 1
                continue
            #we're in the span, check the code
            if code == 0:
                readPosition += 1
                numRef += 1
                pins, pinsz, t = pinsLoad(readPosition, pinsz)
                numIns += t
            elif code == 1: #ins
                pins = True
                pinsz += 1
            elif code == 2: #del
                numDel += 1
                readPosition += 1
                pins, pinsz, t = pinsLoad(readPosition, pinsz)
                numIns += t

        #what's the +- difference to validate the size
        leeway = size * args.sizebuffer
        logging.debug(("regionStart, regionEnd, estSize, leeway, "
                       "estSize+leeway, estSize-leeway, numIns, numDel"))
        logging.debug("%d %d %d %d %d %d %d %d" % (regionStart, regionEnd, \
                        size, leeway, size+leeway, size-leeway, numIns, numDel))
        
        if svtype == 'DEL':
            if size + leeway >= numDel >= size - leeway:
                ret["DEL"] += 1
                #ret.append("DEL " + read.qname)
            else:
                ret["REF"] += 1
                #ret.append("REF " + read.qname)
        elif svtype == 'INS':
            if size + leeway >= numIns >= size - leeway:
                ret["INS"] += 1
                #ret.append("INS " + read.qname)
            else:
                ret["REF"] += 1
                #ret.append("REF " + read.qname)
        else:
            #ret["WAIT"] += 1
            pass#there is no small inverson fixing
            #ret.append("WAIT " + read.qname)
    k = ret.keys()
    k.sort()
    return ["%s:s:%d" % (x, ret[x]) for x in k]

def run(args):
    args = parseArgs(args)
    bam = pysam.Samfile(args.bam)

    #putative caller
    fh = open(args.bed)
    
    tails.BUFFER = 100 #it's okay if these don't combine into one. -- maybe even prefered
    
    vtypes = ['INS', 'DEL', 'MIS']#, 'UNK'] #unk eventually for -find any-
    i = 0
    for line in fh.readlines():
        myentry = line.strip().split('\t')
        if myentry[4] not in vtypes:
            if myentry[4] == 'UNK':
                logging.warning("Bed Entry %s is UNK and can't be forced... skipping" % (strmyentry))
                sys.stdout.write(line.strip() + "\t.\n")
                continue
            else:
                logging.error("Bed Entry %s svtype column isn't one of %s" % (str(myentry), str(vtypes)))
                exit(1)
        if myentry[0] not in bam.references:
            logging.error("Invalid Chromosome %s" % myentry[0])
            continue
        bed = [myentry[0], int(myentry[1]), int(myentry[2]), myentry[3], myentry[4], int(myentry[5])]   
        sup = tailsSearch(bam, bed, args)
        sup.extend(spotsSearch(bam, bed, args))
        #nulls
        if len(sup) == 0:
            sup = ['.']
        
        sys.stdout.write(line.strip() + "\t" + ",".join(sup) + "\n")
        i += 1
        if i % 250 == 0:
            sys.stdout.flush()
if __name__ == '__main__':
    run(sys.argv[1:])
