#!/usr/bin/env python

import sys, time, re, gc, argparse, logging
import pysam, numpy, h5py
from collections import defaultdict, Counter
from multiprocessing import Pool

from pbsuite.utils.setupLogging import *
from pbsuite.honey.SVMachines import *

USAGE = """\
Count number of errors and thier position releative to the reference. 
"""

TODO = """ 
Make multiprocessing.pool.map on a per region basis - The trick will be getting the
    numpy and h5 data put together
"""

columns = ["coverage", "matches", "mismatches", "insertions", \
           "insertionsize", "deletions", "tails", "avgmapq", "label"]
COV  = 0
MAT  = 1
MIS  = 2
INS  = 3
INSZ = 4
DEL  = 5
TAI  = 6
MAQ  = 7
LAB  = 8
           
"""
Proposed Metrics:
    I'm just going to use .m5 from now on
    1) TailMapping - Map reads, trim tails, map tails - document tailinfo
    break tails down to mappedTails and unmappedTails
    2) somehow need to document this primary and tails
"""
#Must not exceed 300,000 data points
CHUNKSHAPE = (7, 40000)

def markDups(bam):
    """
    Marks all reads that aren't the best read from a zmw as duplicate
    """
    names = {}
    numDups = 0
    for read in bam:
        n = "/".join(read.qname.split('/')[:2])
        try:
            cRead, cScore = names['/'.join(read.qname.split('/')[:2])]
            #myScore = sum([ord(y)-33 for y in read.qqual])/float(len(read.qqual))
            myScore = read.mapq
            if cScore > myScore:
                read.is_duplicate = True
            else:
                cRead.is_duplicate = True
            numDups += 1
        except KeyError:
            #myScore = sum([ord(y)-33 for y in read.qqual])/float(len(read.qqual))
            myScore = read.mapq
            names[n] = (read, myScore)
    logging.info("Marked %d ZMW duplicates" % (numDups))
    del(names)

def expandCigar(cigar):
    """
    Turns the abbreviated cigar into the full array
    0 = M; 1 = I; 2 = D
    """
    ret = []
    for t,s in cigar:
        if t < 3: #remove non mid (dangerous if blasr changes)
            ret.extend([t]*s)
    return ret

def expandMd(md):
    """
    Turns abbreviated MD into a full array
    """
    ret = []
    for i in re.findall("\d+|\^?[ATCGN]+", md):
        if i[0] == '^':
            d = list(i[1:])
        elif i[0] in ["A","T","C","G","N"]:
            d = list(i)
        else:
            d = xrange(int(i))
        #for j in d:
            #yield type(j) == int # is it a true Match
        ret.extend(d)
    return ret

def countErrors(reads, offset, size, mint, maxt, ignoreDups):
    """
    Essentially a pileup that much more organized between mid types
    """
    container = numpy.zeros( ( len(columns), size ) )
    numReads = 0
    
    for align in reads:
        numReads += 1
        if numReads % 1000 == 0:
            logging.info("parsed %d reads" % (numReads))
              
        seq = align.query
        cigar = expandCigar(align.cigar)
    
        mdTag = None
        for i in align.tags:
            if i[0] == "MD":
                mdTag = expandMd(i[1])
        
        if mdTag is None:# and alignment.target:
            logging.error(("MD tag is absent. Run samtools calmd"))
            exit(1)
         
        #get starts within region
        regionStart = max(0, align.pos - offset) 
        regionEnd = min(offset + size, align.aend)  # +1?
        
        #I'm always starting at the beginning of the read, but the beginning may hit before my regionStart
        #start = min(align.pos - offset, 0)
        start = align.pos
        #alignSpan = regionEnd - regionStart
        #myData = numpy.zeros( (4, alignSpan) )
        #check covering bases
        container[COV,  regionStart : regionEnd] += 1
        container[MAQ, regionStart : regionEnd] += align.mapq
        
        #previous base was an insert prevent multiple base ins
        pins = False
        curMd = 0
        
        for code in cigar:
            if start < regionStart or start >= regionEnd:
                if code != 1: 
                    curMd += 1
                    start += 1
                continue
            if curMd >= len(mdTag):#short circuit 
                break
            elif code == 0:
                if mdTag[curMd]: #mat
                    container[MAT, start] += 1
                else: #mis
                    container[MIS, start] += 1
                start += 1
                curMd += 1
                pins = False
            elif code == 1: #ins
                if not pins:
                    container[INS, start] += 1
                #insz
                container[INSZ, start] += 1
                pins = True
            elif code == 2: #del
                container[DEL, start] += 1
                start += 1
                curMd += 1
                pins = False
        
        #container[MAT:DEL+1, regionStart : regionEnd] += myData
        #check tails
        if ignoreDups and align.is_duplicate:
            continue
            
        if align.cigar[0][0] == 4:
            tailLen = align.cigar[0][1]
            
            #remove maxt because it's really dumb
            if tailLen >= mint:
                #get coordinates and prevent overstepping boundaries
                if align.is_reverse:
                    start = min(len(container[0], align.aend - offset)
                    end = min(len(container[0]), start + tailLen)
                else:
                    end = max(0, align.pos - offset)
                    start = max(0, end - tailLen) 
                    
                if  end - start > 0:
                    container[TAI, start:end] += 1
                    
        if align.cigar[-1][0] == 4:
            tailLen = align.cigar[-1][1]
            if tailLen >= mint:
                if not align.is_reverse:
                    start = min(len(container[0], align.aend - offset)
                    end = min(len(container[0]), start + tailLen)
                else:
                    end = max(0, align.pos - offset)
                    start = max(0, end - tailLen) 
                    
                if end - start > 0:
                    container[TAI, start:end] += 1
                 
    container[MAQ] = container[MAQ] / container[COV]
    
    logging.info("parsed %d reads" % numReads)
    
    
    return container, numReads
 
def signalTransform(data, cov, slopWindow, avgWindow):
    """
    Go through the signal processing changes on a 1d array
    """
    orig = numpy.convolve(data, avgWindow, "same") #smooth
    norm = orig / cov
    dat = numpy.convolve(norm, slopWindow, "same") #slope transform
    dat = numpy.abs(dat) # absolute value
    dat = numpy.convolve(dat, avgWindow, "same") #smooth
    dat = numpy.convolve(dat * orig, avgWindow, "same") # Take back to absolute space and smooth again
    
    return dat
    
def callHotSpots(data, threshPct, covThresh, binsize, offset):
    """
    """
    output = []
    
    sumWindow = numpy.ones(binsize)
    avgWindow = numpy.ones(binsize)/float(binsize)
    
    slopWindow = numpy.zeros(binsize)
    slopWindow[:binsize/3] = -1
    slopWindow[-binsize/3:] = 1

    data[LAB] = 0
    
    cov = numpy.convolve(data[COV], avgWindow, "same")
    
    mis = signalTransform(data[MIS], cov, slopWindow, avgWindow)
    truth = mis >= threshPct
    if numpy.any(truth):
        data[LAB, truth] += 2**MIS
    del(mis)
    
    ins = signalTransform(data[INS], cov, slopWindow, avgWindow)
    t = ins >= threshPct
    truth = numpy.any([truth, t], axis = 0)
    if numpy.any(t):
        data[LAB, t] += 2**INS
    
    insBas = signalTransform(data[INSZ]*ins, cov*binsize, slopWindow, avgWindow)
    t = insBas >= threshPct
    truth = numpy.any([truth, t], axis = 0)
    if numpy.any(t):
        data[LAB, t] += 2**INSZ
    del(insBas)
    del(ins)
    
    dele = signalTransform(data[DEL], cov, slopWindow, avgWindow)
    t = dele >= threshPct
    truth = numpy.any([truth, t], axis = 0)
    if numpy.any(t):
        data[LAB, t] += 2**DEL
    del(dele)
    
    tail = signalTransform(data[TAI], cov, slopWindow, avgWindow)
    t = tail >= threshPct
    truth = numpy.any([truth, t], axis = 0)
    if numpy.any(t):
        data[LAB, t] += 2**TAI
    del(tail)
    
    #maq = signalTransform(data[MAQ], cov, slopWindow, avgWindow)
    #truth = numpy.any([truth, maq >= threshPct], axis = 0)
    #del(maq)
    
    truth = numpy.all([truth, cov >= covThresh], axis = 0)
    
    points = makePoints(truth, binsize)
    output = summarizePoints(points, offset, data)
    
    return output
    
def makePoints(truth, binsize):
    """
    make the points for the truth set made from the data container
    truth = numpy.array() with boolean values
    """
    #prevent weirdness
    truth[-1] = False
    shift = numpy.roll(truth, 1)
    
    starts = truth & ~shift
    ends = ~truth & shift

    points = zip(numpy.nonzero(starts)[0], numpy.nonzero(ends)[0])
    npoints = []
    if len(points) == 0:
        return npoints
    curStart, curEnd = points[0]
    #compress:
    for start, end in points[1:]:
        if start - curEnd <= binsize:
            curEnd = end
        else:
            npoints.append((curStart, curEnd))
            curStart = start
            curEnd = end
    
    npoints.append((curStart, curEnd))
    
    return npoints
    
def summarizePoints(points, offset, container):
    """
    summarize the data in the points
    """
    output = []
    for start, end in points:
        p = container[:,start:end].mean(axis=1)
        label = Counter(container[LAB,start:end]).most_common()[0][0]
        output.append( (end-start, start + offset, end + offset, \
                        p[COV], \
                        p[MAT], \
                        p[MIS], \
                        p[INS], \
                        p[INSZ], \
                        p[DEL], \
                        p[TAI], \
                        p[MAQ], \
                        label) ) #I don't want an average here
    return output

def spotToString(chrom, spot):
    """
    """
    return "%s\t%s" % (chrom, 
          ("%d\t%d\t%d\t" + \
           "\t".join(["%.3f"] * len(columns))) % spot)

def parseArgs():
    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    
    parser.add_argument("-r", "--region", type=str, default=None,\
                        help="Only call spots in region.bed")
    
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Basename for output (BAM.hon)")
    parser.add_argument("-z", "--noZmwDedup", action="store_false", \
                        help="Allow tails from all subreads per ZMW (False)")
    
    parser.add_argument("-s", "--noCallSpots", action="store_true",\
                        help=("Don't call spots where error rates spike (False)"))

    parser.add_argument("-t", "--mintail", type=int, default=100, \
                        help=("Minimum number of bases in a tail before it's "
                              "considered discordant (100 bp)"))
    parser.add_argument("-T", "--maxtail", type=int, default=500, \
                        help=("Maximum number of bases in a tail before it no "
                              "longer contributes to the discordant count (500)"))
    
    parser.add_argument("-c", "--minCoverage", type=int, default=3, \
                        help="Minimum coverage for a spot call to be made (3)")
    parser.add_argument("-e", "--threshold",  type=float, default=2,
                        help="Spot Threshold (2))")
                        
    parser.add_argument("-b", "--binsize", type=int, default=50, \
                        help="binsize for window averaging (50)")
    
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    args = parser.parse_args()

    if args.output is None:
        args.output = args.bam[:-4]+".hon"
    
    return args

if __name__ == '__main__':
    numpy.seterr(all="ignore")
    args = parseArgs()
    setupLogging(args.debug)
    
    bam = pysam.Samfile(args.bam)
   
    if args.noZmwDedup:
        markDups(bam)
    
    MINTAIL = args.mintail
    MAXTAIL = args.maxtail
    EPTHRES = args.threshold
    CVTHRES = args.minCoverage
    BINSIZE = args.binsize
    
    if not args.noCallSpots:
        hotspots = open(args.output+".spots", 'w')
        hotspots.write(("chrom\tsize\tstart\tend\t%s\n" % "\t".join(columns)))
    
    regions = []
    if args.region:
        fh = open(args.region,'r')
        for line in fh.readlines():
            data = line.strip().split('\t')
            regions.append((data[3], data[0], int(data[1]), int(data[2])))
        fh.close()
    else:
        for chrom, size in zip(bam.references, bam.lengths):
            regions.append((chrom, chrom, 0, size))
    
    h5Name = args.output + ".h5"
    results = h5py.File(h5Name, 'w')
    results.attrs["columns"] = columns
    results.close()

    totReads = 0
    totSpots = 0
    for groupName, chrom, start, end in regions:
        size = end - start
        regName =  "%s:%d-%d" % (chrom, start, end)
        logging.info("making container for %s (%s %d bp)" % (groupName, regName, size))
        results = h5py.File(h5Name, 'a')
        out = results.create_group(groupName)
        
        out.attrs["reference"] = chrom
        out.attrs["start"] = start
        out.attrs["end"] = end
        
        logging.info("parsing bam" )
        
        #I'll need to extract this part
        reads = bam.fetch(chrom, start, end)
        #Might want to put this back
        #container = out.create_dataset("data", data=myData, chunks=CHUNKSHAPE, compression="gzip")
        myData, numReads = countErrors(reads, start, size, MINTAIL, \
                                  MAXTAIL, args.noZmwDedup)
        
        container = out.create_dataset("data", data = myData, \
                                chunks = CHUNKSHAPE, compression="gzip")
        totReads += numReads
        
        if numReads == 0:
            logging.warning("No reads found in %s" % groupName)
            continue
        
        if not args.noCallSpots:
            logging.info("calling spots")
            spots = callHotSpots(container, EPTHRES, CVTHRES, BINSIZE, start)
            #Down to here for multiprocessing
            logging.info("found %d spots" % (len(spots)))
            totSpots += len(spots)
            for i in spots:
                hotspots.write(spotToString(chrom, i)+"\n")
            hotspots.flush()
            
        results.flush()
        results.close()
        myData = None
        container = None
        gc.collect()
    
    logging.info("finished %d reads" % (totReads))
    if not args.noCallSpots:
        logging.info("finshed %d spots" % (totSpots))
        hotspots.close()
