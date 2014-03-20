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

#How big is an insertion error before we consider it something...
MINSINGLEINS = 5 #I don't want this to be a parameter

columns = ["coverage", "matches", "mismatches", "insertions", "deletions"]

### NUMPY ARRAY HDF5 COLUMNS AND SIZES
#Biggest integer I want to deal with
BIGINT  = 2000
BIGINTY = numpy.float32

COV  = 0
MAT  = 1
MIS  = 2  
INS  = 3  
DEL  = 4  

#Must not exceed 300,000 data points
CHUNKSHAPE = (5, 55000)

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
            if key == 'label':
                self.type = self.tags[key]
            else:
                try:
                    tag.append("%s=%0.3f" % (str(key), self.tags[key]))
                except TypeError:
                    tag.append("%s=%s" % (str(key), str(self.tags[key])))
        
        if self.type == 'INS':
            try:
                self.size = self.tags["InszMean"]
            except KeyError:
                self.size = "?"
        else:
            self.size = self.end - self.start
        
        tag = ";".join(tag)
        dat = [self.chrom, self.out_start, self.start, self.in_start, self.in_end, self.end, self.out_end, self.type, self.size, tag]
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(*dat).replace("None", ".")
    
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
    -- C translate candidate -- 
    """
    ret = []
    for t,s in cigar:
        if t < 3: #remove non mid (dangerous if blasr changes)
            ret.extend([t]*s)
    return ret

def expandMd(md):
    """
    Turns abbreviated MD into a full array
    --- C translate candidate --
    """
    ret = []
    for i in re.findall("\d+|\^?[ATCGN]+", md):
        if i[0] == '^':
            d = list(i[1:])
        elif i[0] in ["A","T","C","G","N"]:
            d = list(i)
        else:
            d = xrange(int(i))
        ret.extend(d)
    return ret

def countErrors(reads, offset, size, args, readCount=None):
    """
    Sum the errors over any particular reference base
    """
    container = numpy.zeros( ( len(columns), size ), dtype=BIGINTY )
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
        
        #I'm always starting at the beginning of the read,
        # but the beginning may hit before my regionStart
        start = align.pos - offset
        #check covering bases
        container[COV,  regionStart : regionEnd] += BIGINTY(1)
        #MAQ
        
        #previous base was an insert prevent multiple base ins
        pins = False
        pinsz = 0
        curMd = 0
        def pinsLoad(start, size):
            if not pins:
                return False, 0
            if size >= MINSINGLEINS:
                container[INS, start-1] += BIGINTY(1)
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
                    container[MAT, start] += BIGINTY(1)
                else: #mis
                    container[MIS, start] += BIGINTY(1)
                start += 1
                curMd += 1
                pins, pinsz = pinsLoad(start, pinsz)
            elif code == 1: #ins
                pins = True
                pinsz += 1
            elif code == 2: #del
                container[DEL, start] += BIGINTY(1)
                start += 1
                curMd += 1
                pins, pinsz = pinsLoad(start, pinsz)
    
    logging.info("parsed %d reads" % numReads)
    
    return container, numReads

def preprocessSignal(signal, coverage):
    """
    Normalize and print stats returning data and it's std
    """
    rate = numpy.convolve(signal/coverage, avgWindow, "same")
    rate[numpy.any([numpy.isinf(rate), numpy.isnan(rate)], axis=0)] = 0
    mu = numpy.mean(rate)
    sd = numpy.std(rate)
    logging.info("RateMean %f  -- RateStd  %f" % (mu, sd))
    return rate, mu, sd
    
def signalTransform(dat):
    """
    Go through the signal processing changes on a 1d array
    """
    return numpy.convolve(dat, slopWindow, "same") 

def postSigStats(sig):
    logging.info("MaxSig: %f MeanSig: %f StdSig %f MinSig: %f" \
                 % (numpy.max(sig), numpy.mean(sig), numpy.std(sig), \
                    numpy.min(sig)))
   
def makeKernals(binsize=100):
    """
    My Kernals for Convolution -- push global
    """
    global avgWindow
    global slopWindow
    avgWindow = numpy.ones(binsize, dtype=numpy.float16)/float(binsize)
    #slop - downstream average minus upstream average
    slopWindow = numpy.ones(binsize, dtype=numpy.float16) / (binsize/2)
    slopWindow[:binsize/2] = -slopWindow[:binsize/2]
    
def callHotSpots(data, offset, args): #threshPct, covThresh, binsize, offset):
    """
    """
    ret = []
    #coverage
    cov = numpy.convolve(data[COV], avgWindow, "same")
    covTruth = numpy.all([cov >= args.minCoverage, cov <= args.maxCoverage], axis=0)
    logging.info("MaxCov:%d MeanCov:%d StdCov:%d MinCov:%d" \
            % (numpy.max(data[COV]), numpy.mean(data[COV]), \
               numpy.std(data[COV]), numpy.min(data[COV])))
    
    #mis
    logging.info("MIS processing")
    mis, mu, sd = preprocessSignal(data[MIS], data[COV])
    mis = signalTransform(mis) 
    postSigStats(mis)
    ret.extend(makeSpotResults(mis, sd, "MIS", cov, covTruth, args))
    del(mis)
    
    #ins
    logging.info("INS processing")
    ins, mu, sd = preprocessSignal(data[INS], data[COV])
    ins = signalTransform(ins)
    postSigStats(ins)
    ret.extend(makeSpotResults(ins, sd, "INS", cov, covTruth, args))
    del(ins)
    
    #dele
    logging.info("DEL processing")
    dele, mu, sd = preprocessSignal(data[DEL], data[COV])
    dele = signalTransform(dele)
    postSigStats(dele)
    ret.extend(makeSpotResults(dele, sd, "DEL", cov, covTruth, args))
    del(dele)
    
    return ret
    
def makeSpotResults(datpoints, sd, label, cov, covTruth, args):
    """
    Find starts and ends ranges, then group the start/end pairs by closest proximity
    Need to explore the edge cases
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
    thresh = sd * args.threshold
    entries = []
    startPoints = makePoints(numpy.all([datpoints <= -thresh, covTruth], axis = 0), args.binsize, 's')
    endPoints   = makePoints(numpy.all([datpoints >=  thresh, covTruth], axis=0), args.binsize, 'e')
    sPos = 0
    ePos = 0
        
    points = []
    for i in startPoints:
        bisect.insort(points, i)
    for i in endPoints:
        bisect.insort(points, i)
    i = 0
    
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
    ends = ~truth & shift
    
    points = zip(numpy.nonzero(starts)[0], numpy.nonzero(ends)[0])
    npoints = []
    if len(points) == 0:
        return npoints
    curStart, curEnd = points[0]
    
    #compress: <-- Don't need anymore...?
    for start, end in points[1:]:
        #if start - curEnd <= binsize:
            #curEnd = end
        #else:
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
    
    # I think I could short circuit...?
    # are there enough reads with insertion
    
    #do the hard work
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
                    if curSize >= MINSINGLEINS:
                        readHasIns = True
                        mySize += curSize
                    curSize = 0
        if readHasIns and mySize >= args.minInsz:
            nReadsIns += 1
            totSizes.append(mySize)
        
    logging.debug("size %d %s" % (len(totSizes), str(totSizes)))
    if len(totSizes) == 0:
        logging.debug("no insertions found!? %s" % (str(spot)))
        return True # true you should filter
    
    totSizes.sort()
    totSizes = numpy.array(totSizes)
    mean = totSizes.mean()
    median = numpy.percentile(totSizes, 50)
    firstQ = numpy.percentile(totSizes, 25)
    thirdQ = numpy.percentile(totSizes, 75)
    
    logging.debug("cov    %d" % coverage )
    logging.debug("sizes  %s" % str(totSizes) )
    logging.debug("mean   %d" % mean )
    logging.debug("median %d" % median)
    logging.debug("firstQ %d" % firstQ)
    logging.debug("thirdQ %d" % thirdQ)
    
    if len(totSizes) == 1:
        logging.debug("single insertion found %s " % (str(spot)))
        return True
    
    if len(totSizes) > (coverage * args.insPct):
        spot.tags["InszCount"]  = int(nReadsIns)
        spot.tags["InszMean"]   = int(mean)
        spot.tags["InszMedian"] = int(median)
        spot.tags["Insz1stQ"]   = int(firstQ)
        spot.tags["Insz3rdQ"]   = int(thirdQ)
        return False
    logging.debug("not large pct insertion found %s " % (str(spot)))
    return True
           
def parseArgs(argv, established=False):
    parser = argparse.ArgumentParser(prog="Honey.py spots", description=USAGE, \
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

    procA = parser.add_argument_group("Processing Arguments")
    procA.add_argument("-b", "--binsize", type=int, default=100, \
                        help="binsize for window averaging (%(default)s)")
    
    pGroup = parser.add_argument_group("Spot-Calling Threshold Arguments")
    pGroup.add_argument("-e", "--threshold",  type=float, default=5,
                        help="Minimum Spot Threshold (%(default)s)")
    pGroup.add_argument("-c", "--minCoverage", type=int, default=5, \
                        help="Minimum coverage of a region (%(default)s)")
    pGroup.add_argument("-C", "--maxCoverage", type=int, default=BIGINT, \
                        help="Maximum coverage of a region (%(default)s)")
    pGroup.add_argument("-i", "--minInsz", type=int, default=50, \
                        help="Minimum insertion size (%(default)s)")
    pGroup.add_argument("-I", "--insPct", type=float, default=0.33, \
                        help="Minimum pct of spot coverage with insertion (%(default)s)")
    pGroup.add_argument("-f", "--nonFull", action="store_true", \
                        help="Allow calls with only putative starts xor ends")
    pGroup.add_argument("--sizeMin", type=int, default=150, \
                        help="Minimum Size of spot to be called (%(default)s)")
    pGroup.add_argument("--sizeMax", type=int, default=2000, \
                        help="Maximum Size of spot to be called (%(default)s)")
    
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    
    args = parser.parse_args(argv)
    if args.maxCoverage > BIGINT:
        logging.error("Max Coverge must be less than %d" % (BIGINT))
        exit(0)
    if args.output is None:
        args.output = args.bam[:-4]+".hon"
    
    return args

def run(argv):
    numpy.seterr(all="ignore")
    args = parseArgs(argv)
    setupLogging(args.debug)
    
    makeKernals(args.binsize)
    
    bam = pysam.Samfile(args.bam)
   
    if not args.noCallSpots:
        hotspots = open(args.output+".spots", 'w')
        hotspots.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
    
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
        #myData, numReads = countErrors(reads, start, size, MINTAIL, \
                                  #MAXTAIL, args.noZmwDedup)
        #container = out.create_dataset("data", data=myData, chunks=CHUNKSHAPE, compression="gzip")
        myData, numReads = countErrors(reads, start, size, args, readCount)
        if size < CHUNKSHAPE[1]:
            chunk = (5, size-1)
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
            logging.info("Filtering INS results")
            for spot in spots:
                spot.chrom = chrom
                #add tags for pval and container stats etc
                spot.offset(start)
                #filter on insertion stuff
                if spot.tags["label"] == "INS" and filterINSZ(bam, spot, args):
                    continue
                #Filter on size
                if spot.size < args.sizeMin or spot.size > args.sizeMax:
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
