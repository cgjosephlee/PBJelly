#!/usr/bin/env python

import sys, time, re, gc, argparse, logging, bisect
import pysam, numpy, h5py
from collections import defaultdict, Counter
from multiprocessing import Pool

from pbsuite.utils.setupLogging import *
#from pbsuite.honey.TGraf import * #future
VERSION = "13.10.15"

USAGE = """\
Detect 'smaller' SVs via measuring discordance between sample and reference in long reads.
"""

TODO = """ 
Make multiprocessing.pool.map on a per region basis - The trick will be getting the
    numpy and h5 data put together -- 
    Or we can split the input at the input .fastq level, merge the h5 files, then recallSpots
"""

columns = ["coverage", "matches", "mismatches", "insertions", \
           "insertionsize", "deletions", "avgmapq"]
COV  = 0
MAT  = 1
MIS  = 2  #4
INS  = 3  #8
INSZ = 4  #16  (which means 24 is best shown insertion)
DEL  = 5  #32
MAQ  = 6
           
"""
Proposed Metrics:
    I'm just going to use .m5 from now on
    1) TailMapping - Map reads, trim tails, map tails - document tailinfo
    break tails down to mappedTails and unmappedTails
    2) somehow need to document this primary and tails
"""
#Must not exceed 300,000 data points
CHUNKSHAPE = (7, 40000)

## {{{ http://code.activestate.com/recipes/511478/ (r1)
import math
import functools

def percentile(N, percent):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return N[int(k)]
    d0 = N[int(f)] * (c-k)
    d1 = N[int(c)] * (k-f)
    return d0+d1

# median is 50th percentile.

## end of http://code.activestate.com/recipes/511478/ }}}
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
    """ Slower method.. I hoped it'd be faster, but extend apparently isn't an expensive operation
        cigar = align.cigar
        #new iter method
        codPos = 0 #where I am in the cigar code
        codCnt = 0 #where I am in the specific cigar code's length
        code = cigar[codPos][0]
        while codPos < len(cigar):
            if cigar[codPos][0] >= 3:
                codPos += 1
                codCnt = 0
                try:
                    code = cigar[codPos][0]
                except IndexError:
                    pass
                continue
            elif codCnt >= cigar[codPos][1]:#Next cigar code
                codPos += 1
                codCnt = 0
                try:
                    code = cigar[codPos][0]
                except IndexError:
                    pass
                continue
            codCnt += 1 #Processing the Nth count of the code"""

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

def countErrors(reads, offset, size, args, readCount=None):
    """
    Sum the errors over any particular reference base
    """
    container = numpy.zeros( ( len(columns), size ) )
    numReads = 0.0 #number finished
    if readCount is not None:
        logging.info("%d reads to parse" % readCount)
    
    lastUpdate = 2
    for align in reads:
        numReads += 1
        if readCount is not None and int(numReads/readCount * 100) > lastUpdate:
            lastUpdate = int(numReads/readCount * 100)
            logging.info("parsed %d reads (%d %% complete)" % (numReads, lastUpdate))
              
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
        regionStart = 0 if align.pos < offset else align.pos - offset
        regionEnd = size if align.aend > (offset + size) else align.aend - offset
        
        #I'm always starting at the beginning of the read, but the beginning may hit before my regionStart
        start = align.pos - offset
        #check covering bases
        container[COV,  regionStart : regionEnd] += 1
        container[MAQ, regionStart : regionEnd] += align.mapq
        
        #previous base was an insert prevent multiple base ins
        pins = False
        pinsz = 0
        curMd = 0
        def pinsLoad(start, size):
            if not pins:
                return False, 0
            if size >= args.minInsz:
                container[INS, start-1] += 1
                container[INSZ, start-1] += size
            return False, 0
        
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
                pins, pinsz = pinsLoad(start, pinsz)
                #pins = False
            elif code == 1: #ins
                #if not pins:
                    #container[INS, start] += 1
                    #pins = True
                #insz
                #container[INSZ, start] += 1
                pins = True
                pinsz += 1
            elif code == 2: #del
                container[DEL, start] += 1
                start += 1
                curMd += 1
                #pins = False
                pins, pinsz = pinsLoad(start, pinsz)
    container[MAQ] = container[MAQ] / container[COV]
    
    logging.info("parsed %d reads" % numReads)
    
    return container, numReads

def checkTails():
    #check tails -- only those with a map
    if (ignoreDups and align.is_duplicate) or not (1 & align.flag):
        return
    if align.cigar[0][0] == 4:
        tailLen = align.cigar[0][1]
        #remove maxt because it's really dumb
        if tailLen >= mint:
            #get coordinates and prevent overstepping boundaries
            if align.is_reverse:
                start = min(len(container[0]), align.aend - offset)
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
                start = min(len(container[0]), align.aend - offset)
                end = min(len(container[0]), start + tailLen)
            else:
                end = max(0, align.pos - offset)
                start = max(0, end - tailLen) 
            if end - start > 0:
                container[TAI, start:end] += 1
     
def signalTransform(data, cov, slopWindow, avgWindow):
    """
    Go through the signal processing changes on a 1d array
    """
    orig = numpy.convolve(data, avgWindow, "same") #smooth
    #norm = orig / cov
    dat = numpy.convolve(orig, slopWindow, "same") #slope transform
    #dat = numpy.convolve(norm, slopWindow, "same") #slope transform
    #dat = numpy.abs(dat) # absolute value
    #dat = numpy.convolve(dat, avgWindow, "same") #smooth
    dat = numpy.convolve(dat * orig, avgWindow, "same") # / cov # Take back to absolute space and smooth again
    return numpy.nan_to_num(dat)
    #return dat
    
def callHotSpots(data, offset, args): #threshPct, covThresh, binsize, offset):
    """
    """
    sumWindow = numpy.ones(args.binsize)
    avgWindow = numpy.ones(args.binsize)/float(args.binsize)
    
    slopWindow = numpy.zeros(args.binsize)
    slopWindow[-args.binsize/2:] = 1
    slopWindow[:args.binsize/2] = -1
    
    slopWindow = slopWindow / float(args.binsize)
    
    ret = []
    
    #coverage
    cov = numpy.convolve(data[COV], avgWindow, "same")
    covTruth = numpy.all([cov >= args.minCoverage, cov <= args.maxCoverage], axis=0)
    logging.info("MaxCov:%d MeanCov:%d StdCov:%d MinCov:%d" \
            % (numpy.max(data[COV]), numpy.mean(data[COV]), numpy.std(data[COV]), numpy.min(data[COV])))
    
    #mis
    mis = signalTransform(data[MIS], cov, slopWindow, avgWindow)
    logging.info("MaxMis:%.3f MeanMis:%.3f StdMis:%.3f MinMis:%.3f" \
            % (numpy.max(mis), numpy.mean(mis), numpy.std(mis), numpy.min(mis)))
    ret.extend(makeSpotResults(mis, "MIS", cov, covTruth, args))
    
    #ins
    ins = signalTransform(data[INS], cov, slopWindow, avgWindow)
    logging.info("MaxIns:%.3f MeanIns:%.3f StdIns:%.3f MinIns:%.3f" \
            % (numpy.max(ins), numpy.mean(ins), numpy.std(ins), numpy.min(ins)))
    ret.extend(makeSpotResults(ins, "INS", cov, covTruth, args))
    del(ins)
    
    #insBas = signalTransform(data[INSZ]*ins, cov*binsize, slopWindow, avgWindow)
    insz = signalTransform(data[INSZ]*data[INS], cov, slopWindow, avgWindow)
    logging.info("MaxInsz:%.3f MeanInsz:%.3f StdInsz:%.3f MinInsz:%.3f" \
            % (numpy.max(insz), numpy.mean(insz), numpy.std(insz), numpy.min(insz)))
    ret.extend(makeSpotResults(insz, "INSZ", cov, covTruth, args))
    del(insz)
    
    #dele
    dele = signalTransform(data[DEL], cov, slopWindow, avgWindow)
    logging.info("MaxDel:%.3f MeanDel:%.3f StdDel:%.3f MinDel:%.3f" \
            % (numpy.max(dele), numpy.mean(dele), numpy.std(dele), numpy.min(dele)))
    ret.extend(makeSpotResults(dele, "DEL", cov, covTruth, args))
    del(dele)
    
    #maq = signalTransform(data[MAQ], cov, slopWindow, avgWindow)
    #logging.info("MaxMaq:%.3f MeanMaq:%.3f StdMaq:%.3f MinMaq:%.3f" \
            #% (numpy.max(maq), numpy.mean(maq), numpy.std(maq), numpy.min(maq)))
    #spots = makeSpotResults(maq, binsize, threshPct, covTruth, "MAQ")
    #del(maq)
    
    return ret
    

class SpotResult():
    def __init__(self, chrom=None, out_start=None, start=None, in_start=None, \
                                in_end=None, end=None, out_end=None, tags=None):
        self.chrom = chrom
        self.out_start = out_start
        self.start = start
        self.in_start = in_start
        self.in_end = out_end
        self.end = end
        self.out_end = out_end
        self.tags = tags if tags is not None else {}
        
    def offset(self, start):
        """
        moves the spot to an offset
        """
        if self.out_start is not None:
            self.out_start += start
        if self.start is not None:
            self.start += start
        if self.in_start is not None:
            self.in_start += start
        if self.in_end is not None:
            self.in_end += start
        if self.end is not None:
            self.end += start
        if self.out_end is not None:
            self.out_end += start
    
    def __str__(self):
        """
        changes a spot named tuple to a svp string
        """
        tag = []
        for key in self.tags:
            try:
                tag.append("%s=%0.3f" % (str(key), self.tags[key]))
            except TypeError:
                tag.append("%s=%s" % (str(key), str(self.tags[key])))
        tag = ";".join(tag)
        dat = [self.chrom, self.out_start, self.start, self.in_start, self.in_end, self.end, self.out_end, tag]
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(*dat).replace("None", ".")
    
def makeSpotResults(datpoints, label, cov, covTruth, args):
    """
    Find starts and ends ranges, then group the start/end pairs by closest proximity
    Need to explore the edge cases
    """
    entries = []
    startPoints = makePoints(numpy.all([datpoints <= -args.threshold, covTruth], axis = 0), args.binsize, 's')
    endPoints   = makePoints(numpy.all([datpoints >= args.threshold, covTruth], axis=0), args.binsize, 'e')
    sPos = 0
    ePos = 0
        
    points = []
    for i in startPoints:
        bisect.insort(points, i)
    for i in endPoints:
        bisect.insort(points, i)
    i = 0
    """
    need to fix
        V^V ^

    What if the procedure was to:
    1) turn starts into SpotResults and place into sorted list.
    2) for each endPoint, we find the nearest neighbor using insort
        if the start currently doesn't have a nearest neighbor
            if the end is contained within a start, the we throw the end away (hoping another end will be around)
            if the start is contained within an end, we throw away the start and keep searching for the end.
        if the end has a nearest neighbor
            we let the ends compete, the closest end will be put together with the start,
            the losing end will be put through the search 
    This doesnt work
    """
    while i < len(points):
        mySpot = SpotResult(tags={"label":label})
        if points[i][2] == 'e':
            mySpot.in_end  = points[i][0]
            mySpot.end     = datpoints[points[i][0]:points[i][1]].argmax() + points[i][0]
            mySpot.out_end = points[i][1]
            mySpot.tags["endCov"] = cov[points[i][0]:points[i][1]].mean()
            mySpot.tags["endSig"] = datpoints[points[i][0]:points[i][1]].mean()
            i += 1
        elif points[i][2] == 's':
            if i+1 < len(points) and points[i+1][2] == 'e':
                mySpot.out_start = points[i][0]
                mySpot.start     = datpoints[points[i][0]:points[i][1]].argmin()+points[i][0]
                mySpot.in_start  = points[i][1]
                mySpot.tags["startCov"] = cov[points[i][0]:points[i][1]].mean()
                mySpot.tags["startSig"] = datpoints[points[i][0]:points[i][1]].mean()
                mySpot.in_end    = points[i+1][0]
                mySpot.end     = datpoints[points[i+1][0]:points[i+1][1]].argmax()+points[i+1][0]
                mySpot.out_end   = points[i+1][1]
                mySpot.tags["endCov"] = cov[points[i+1][0]:points[i+1][1]].mean()
                mySpot.tags["endSig"] = datpoints[points[i+1][0]:points[i+1][1]].mean()
                i += 2
            else:
                mySpot.out_start = points[i][0]
                mySpot.start     = datpoints[points[i][0]:points[i][1]].argmin()+points[i][0]
                mySpot.in_start  = points[i][1]
                mySpot.tags["startCov"] = cov[points[i][0]:points[i][1]].mean()
                mySpot.tags["startSig"] = datpoints[points[i][0]:points[i][1]].mean()
                i += 1
        
        #Either allowing nonFull spots or mySpot is full
        if args.nonFull or (mySpot.start is not None and mySpot.end is not None):
            entries.append(mySpot)
        
    logging.info("%d %s entries" % (len(entries), label))
    return entries
            
def makePoints(truth, binsize, label):
    """
    make the points for the truth set made from the data container
    truth = numpy.array() with boolean values
    """
    #prevent weirdness
    truth[-1] = False
    shift = numpy.roll(truth, 1)
    
    starts = truth & ~shift
    #begins
    ends = ~truth & shift
    #ends
    #sort all of them together.a MEAKE YOUR NOTE HERE
    
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
            npoints.append((curStart, curEnd, label))
            curStart = start
            curEnd = end
    
    npoints.append((curStart, curEnd, label))
    
    return npoints

def filterINSZ(bam, spot, args):
    """
    Goes back through the insz called regions and checks to see
    if it's likely caused by a single large insertion

    filters on
        do at least 1/3rd of the reads have an insertion of MININSZ
        is the average insertion size more than 3x the median insertion size

    returns True if you should filter it out
    """
    ret = []
    
    #fetch left and right coords
    if spot.out_start is None:
        begin2 = max(0, spot.in_end - args.binsize)
        end2 = spot.out_end + args.binsize
        begin1 = begin2; end1 = end2
    elif spot.in_end is None:
        begin1 = max(0, spot.out_start - args.binsize)
        end1 = spot.in_start + args.binsize
        begin2 = begin1; end2 = end1
    else:
        begin1 = max(0, spot.out_start - args.binsize)
        end1 = spot.in_start + args.binsize
        begin2 = max(0, spot.in_end - args.binsize)
        end2 = spot.out_end + args.binsize
    
    #get every read over the region
    reads = bam.fetch(str(spot.chrom), min(begin1, begin2), max(end1, end2))
    totSizes = []
    coverage = 0
    nReadsIns = 0
    #count reads and insizes
        
    for i in reads:
        mySize = 0
        coverage += 1
        start = i.pos
        cigar = expandCigar(i.cigar)
        curSize = 0
        readHasIns = False
        
        for code in cigar: 
            #must be in region
            if ((start >= begin1 and start < end1) or \
                (start >= begin2 and start < end2)) and code == 1:
                curSize += 1
            else:
                if code != 1:
                    start += 1
                if curSize != 0:
                    if curSize >= args.minInsz:
                        readHasIns = True
                        mySize += curSize
                    curSize = 0
        if readHasIns:
            nReadsIns += 1
            totSizes.append(mySize)
        
    totSizes.sort()
    totSizes = numpy.array(totSizes)
    mean = totSizes.mean()
    median = numpy.percentile(totSizes, 50)
    firstQ = numpy.percentile(totSizes, 25)
    thirdQ = numpy.percentile(totSizes, 75)
    
    logging.debug(str(spot) )
    logging.debug("cov    %d" % coverage )
    logging.debug("sizes  %s" % str(totSizes) )
    logging.debug("mean   %d" % mean )
    logging.debug("median %d" % median)
    logging.debug("firstQ %d" % firstQ)
    logging.debug("thirdQ %d" % thirdQ)
    
    if len(totSizes) == 1:
        return True
    
    if len(totSizes) > (coverage * args.insPct):
        spot.tags["InszCount"]  = int(nReadsIns)
        spot.tags["InszMean"]   = int(mean)
        spot.tags["InszMedian"] = int(median)
        spot.tags["Insz1stQ"]   = int(firstQ)
        spot.tags["Insz3rdQ"]   = int(thirdQ)
        return False
    return True
           
def parseArgs(argv, established=False):
    parser = argparse.ArgumentParser(prog="spots", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    ioGroup = parser.add_argument_group("I/O Argument")
    ioGroup.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    if established:
        ioGroup.add_argument("hon", metavar="HON.H5", type=str, \
                        help="HON.h5 containing signal data")
        
    ioGroup.add_argument("-r", "--region", type=str, default=None,\
                        help="Only call spots in region.bed")
    
    ioGroup.add_argument("-o", "--output", type=str, default=None, \
                        help="Basename for output (BAM.hon)")
    
    ioGroup.add_argument("-s", "--noCallSpots", action="store_true",\
                        help=("Don't call spots where error rates spike (False)"))

    #parser.add_argument("-t", "--mintail", type=int, default=100, \
                        #help=("Minimum number of bases in a tail before it's "
                              #"considered discordant (100 bp)"))
    #parser.add_argument("-T", "--maxtail", type=int, default=500, \
                        #help=("Maximum number of bases in a tail before it no "
    procA = parser.add_argument_group("Processing Arguments")
    procA.add_argument("-b", "--binsize", type=int, default=50, \
                        help="binsize for window averaging (50)")
    
    pGroup = parser.add_argument_group("Spot-Calling Threshold Arguments")
    pGroup.add_argument("-e", "--threshold",  type=float, default=2,
                        help="Minimum Spot Threshold (2))")
    pGroup.add_argument("-c", "--minCoverage", type=int, default=3, \
                        help="Minimum coverage of a region (3)")
    pGroup.add_argument("-C", "--maxCoverage", type=int, default=200, \
                        help="Maximum coverage of a region (200)")
    pGroup.add_argument("-i", "--minInsz", type=int, default=20, \
                        help="Minimum insertion size (20)")
    pGroup.add_argument("-I", "--insPct", type=float, default=0.25, \
                        help="Minimum pct of spot coverage with insertion (0.25)")
    pGroup.add_argument("-f", "--nonFull", action="store_true", \
                        help="Allow calls with only putative starts xor ends")
    
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = args.bam[:-4]+".hon"
    
    return args

def run(argv):
    numpy.seterr(all="ignore")
    args = parseArgs(argv)
    setupLogging(args.debug)
    
    bam = pysam.Samfile(args.bam)
   
    if not args.noCallSpots:
        hotspots = open(args.output+".spots", 'w')
        hotspots.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tINFO\n")
    
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
    results.attrs["version"] = VERSION
    results.attrs["columns"] = columns
    results.attrs["parameters"] = str(args)
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
        readCount = bam.count(chrom, start, end)
        #Might want to put this back
        #container = out.create_dataset("data", data=myData, chunks=CHUNKSHAPE, compression="gzip")
        #myData, numReads = countErrors(reads, start, size, MINTAIL, \
                                  #MAXTAIL, args.noZmwDedup)
        myData, numReads = countErrors(reads, start, size, args, readCount)
        if size < CHUNKSHAPE[1]:
            chunk = (7, size-1)
        else:
            chunk = CHUNKSHAPE
        container = out.create_dataset("data", data = myData, \
                                chunks=chunk, compression="gzip")
        totReads += numReads
        
        if numReads == 0:
            logging.warning("No reads found in %s" % groupName)
            continue
        
        if not args.noCallSpots:
            logging.info("calling spots")
            spots = callHotSpots(container, start, args)
            fspots = 0
            logging.info("Filtering INSZ results")
            for spot in spots:
                spot.chrom = chrom
                #add tags for pval and container stats etc
                spot.offset(start)
                if spot.tags["label"] == "INSZ" and filterINSZ(bam, spot, args):
                    continue
                fspots += 1
                if groupName != chrom:
                    spot.tags["RegName"] = groupName
                hotspots.write(str(spot)+"\n")
            logging.info("found %d spots" % (fspots))
            totSpots += fspots
            hotspots.flush()
            #Down to here for multiprocessing
        logging.debug("flush")
        results.flush()
        logging.debug("close")
        results.close()
        myData = None
        container = None
        logging.debug("gc")
        gc.collect()
    
    logging.info("finished %d reads" % (totReads))
    if not args.noCallSpots:
        logging.info("finished %d spots" % (totSpots))
        hotspots.close()


if __name__ == '__main__':
    run(sys.argv[1:])
