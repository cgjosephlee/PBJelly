#!/usr/bin/env python

import sys, time, re, gc, argparse, logging
import pysam, numpy, h5py
from collections import defaultdict
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

#container columns
#I can't play with this anymore without changing countError
#COV = 0         #Coverage
#MAT = 1         #Matches
#MIS = 2         #Mismatches
#INS = 3         #Insertions
#INSZ = 4       #Number of bases inserted at position 
#DEL = 5         #Deletions
#TAI = 6         #Tails
#TER = 7         #Total Error #don't think I need
#SIM = 8         #Match percent --> MAT / ( COV + TAI ) #don't think I need
#ERR = 9         #Error Percent --> TER / COV #don't think I need
#RATIO = 10       #Error over Match Ratio --> TER / MAT #don't think I need
#DIFF = 11       #TotalError - Matches --> TER - MAT #don't think I need
#MAPQ = 12       #Average Map Quality 
#Maybe add a label field?

#columns = ["coverage", "matches", "mismatches", \
           #"insertions", "numinsbases", "deletions", "tails", \
           #"toterrors", "matchpct", "errorpct", \
           #"errmatratio", "difference", "avgmapq" ]


columns = ["coverage", "matches", "mismatches", \
           "insertions", "deletions", "tails", "avgmapq"]
COV = 0
MAT = 1
MIS = 2
INS = 3
DEL = 4
TAI = 5
MAQ = 6
#LAB = 7
           
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
        regionEnd = min(size, align.aend - offset)  # +1?
        alignSpan = regionEnd - regionStart
        #I'm always starting at the beginning of the read, but the beginning may hit before my regionStart
        start = min(align.pos - offset, 0)
        #alignSpan = regionEnd - regionStart
        myData = numpy.zeros( (4, alignSpan) )
        #check covering bases
        container[COV,  regionStart : regionEnd] += 1
        container[MAQ, regionStart : regionEnd] += align.mapq
 
        #previous base was an insert prevent multiple base ins
        pins = False
        curMd = 0
        
        for code in cigar:
            if start < 0 or start > alignSpan:
                if code != 1: 
                    curMd += 1
                    start += 1
            if curMd >= len(mdTag):#short circuit 
                break
            elif code == 0:
                if mdTag[curMd]: #mat
                    myData[0, start] += 1
                else: #mis
                    myData[1, start] += 1
                start += 1
                curMd += 1
                pins = False
            elif code == 1: #ins
                if not pins:
                    myData[2, start] += 1
                #insz
                #myData[3, start] += 1
                pins = True
            elif code == 2: #del
                myData[3, start] += 1
                start += 1
                curMd += 1
                pins = False
        
        container[MAT:DEL+1, regionStart : regionEnd] += myData
        
        #check tails
        if ignoreDups and align.is_duplicate:
            continue
            
        if align.cigar[0][0] == 4:
            tailLen = align.cigar[0][1]
            if tailLen >= mint and tailLen <= maxt:
                start = max(0, align.pos - tailLen) 
                end = align.pos 
                if start-end != 0:
                    container[TAI, start:end] += 1
        
        if align.cigar[-1][0] == 4:
            tailLen = align.cigar[-1][1]
            if tailLen >= mint and tailLen <= maxt:
                start = align.aend 
                end = min(len(container[2]), align.aend + tailLen)
                if start-end != 0:
                    container[TAI, start:end] += 1
    
    container[MAQ] = container[MAQ] / container[COV]
    
    logging.info("parsed %d reads" % numReads)
    
    
    return container, numReads
    #if numReads == 0:#Short circuiting. This doesn't usually happen, but still it's an easy time saver
    
    #logging.info("calculating results")
    #container[SIM] = (container[MAT]) / (container[COV] + container[TAI])
    #container[SIM][numpy.isinf(container[SIM])] = 0
    #container[SIM][numpy.isnan(container[SIM])] = 0
    
    #container[TER] = container[MIS] + container[DEL] + container[INS] + container[TAI]
            #
    #container[ERR] = container[TER] / container[COV]
    #container[ERR][numpy.isinf(container[ERR])] = 0
    #container[ERR][numpy.isnan(container[ERR])] = 0
    #container[RATIO] = container[TER] / container[MAT]

    #container[DIFF] = container[TER] - container[MAT]
    
    #return container, numReads


def expandAlign(alignment):
    """
    Takes a pysam Alignment and creates 
    (reference, query) alignments
    For example:
        query     =  ATCGC-GT
        reference =  AT-GCGGA
        Where C inserted, G deleted, A->T Sub
    """
    seq = alignment.query
    cigar = expandCigar(alignment.cigar)
    
    mdTag = None
    for i in alignment.tags:
        if i[0] == "MD":
            mdTag = expandMd(i[1])
    
    if mdTag is None:# and alignment.target:
        logging.error(("MD tag is absent. Run samtools calmd"))
        exit(1)
    
    qPos = 0
    tPos = 0
    tSeq = []
    qSeq = []
    for i in cigar:
        if i == 0:
            if mdTag[tPos] == '-':
                #mat
                tSeq.append(seq[qPos])
            else:
                #mis
                tSeq.append(mdTag[tPos])
            qSeq.append(seq[qPos])
            qPos += 1
            tPos += 1
        elif i == 1:
            #ins
            qSeq.append(seq[qPos])
            tSeq.append("-")
            qPos += 1
        elif i == 2:
            #del
            qSeq.append("-")
            tSeq.append(mdTag[tPos])
            tPos += 1
    #Expanding query seq and filling in target seq
    return (qSeq,tSeq)

def old_countErrors(reads, offset, size, mint, maxt, ignoreDups):
    """
    counts the base data in the read over the reference
    reads - the reads to process (sholud have been previously fetched)
    offset - the region's start
    size - size of the region
    mint - min tail len
    maxt - max tail contrib
    ignoreDups - don't count dups

    I can speed this up by not creating the expandAlign and instead
    doing direct comparisons between the expanded MD and Cigar
    """
    container = numpy.zeros( ( len(columns), size ) )
    numReads = 0
    
    for align in reads:
        numReads += 1
        if numReads % 1000 == 0:
            logging.info("parsed %d reads" % (numReads))
        
        container[COV, align.pos-offset : align.aend-offset] += 1
        container[MAPQ, align.pos-offset : align.aend-offset] += align.mapq
        #check covering bases
        start = align.pos - offset
        query, target = expandAlign(align)
        
        #previous base was an insert
        pins = False#prevent multiple base ins
        
        for qSeq,tSeq in zip(query, target):
            
            #filter out bases outside of our region
            if start < 0:#before our start
                if tSeq != '-': start += 1
                continue
            if start >= size:#after our end
                break
            
            if tSeq == '-':
                #insertion
                if not pins:
                    container[INS, start] += 1
                container[INSZ, start] += 1
                pins = True
            elif qSeq == '-':
                #deletion
                container[DEL, start] += 1
                start += 1
                pins = False
            elif qSeq == tSeq:
                #match
                container[MAT, start] += 1
                start += 1
                pins = False
            elif qSeq != tSeq:
                #mismatch
                container[MIS, start] += 1
                start += 1
                pins = False
            
        #check tails
        if ignoreDups and align.is_duplicate:
            continue
            
        if align.cigar[0][0] == 4:
            tailLen = align.cigar[0][1]
            if tailLen >= mint and tailLen <= maxt:
                start = max(0, align.pos - tailLen) 
                end = align.pos 
                container[TAI, start:end] += 1
        
        if align.cigar[-1][0] == 4:
            tailLen = align.cigar[-1][1]
            if tailLen >= mint and tailLen <= maxt:
                start = align.aend 
                end = min(len(container[2]), align.aend + tailLen)
                container[TAI, start:end] += 1
    
    logging.info("parsed %d reads" % numReads)
    
    if numReads == 0:#Short circuiting. This doesn't usually happen, but still it's an easy time saver
        return container, numReads
    
    logging.info("calculating results")
    container[SIM] = (container[MAT]) / (container[COV] + container[TAI])
    container[SIM][numpy.isinf(container[SIM])] = 0
    container[SIM][numpy.isnan(container[SIM])] = 0
    
    container[TER] = container[MIS] + container[DEL] + container[INS] + container[TAI]
            
    container[ERR] = container[TER] / container[COV]
    container[ERR][numpy.isinf(container[ERR])] = 0
    container[ERR][numpy.isnan(container[ERR])] = 0
    container[RATIO] = container[TER] / container[MAT]

    container[DIFF] = container[TER] - container[MAT]
    container[MAPQ] = container[MAPQ] / container[COV]
    
    return container, numReads

def makeSimpleSvmPt(points, binsize):
    """
    make the data or the spot
    #only doing three right now
    todo - all individually - super classification of the three types
    (4 if you include concordant)
    --- I'm going to need to create these in my numpy/h5 obj
    """
    avgWindow = numpy.ones(int(binsize))/float(binsize)
    data = numpy.zeros( (6, len(points[0])) )
    logging.debug("cov")
    #cov = numpy.convolve(data[COV], avgWindow, "same")
    data = points[MAT:TAI+1]/points[COV] #numpy.convolve(data[MAT], avgWindow, "same")/ cov
    #mis = data[MIS]/data[COV] #numpy.convolve(data[MIS], avgWindow, "same")/ cov
    #ins = data[INS]/data[COV] #numpy.convolve(data[INS], avgWindow, "same")/ cov
    #insz = data[INSZ]/data[COV] #numpy.convolve(data[INSZ], avgWindow, "same")/ cov
    #dele = data[DEL]/data[COV] #numpy.convolve(data[DEL], avgWindow, "same")/ cov
    #tail = data[TAI]/data[COV] #numpy.convolve(data[TAI], avgWindow, "same")/ cov
    
    return data

def makeSimpleSvmPt2(points, binsize):
    """
    make the data or the spot
    #only doing three right now
    todo - all individually - super classification of the three types
    (4 if you include concordant)
    --- I'm going to need to create these in my numpy/h5 obj
    """
    avgWindow = numpy.ones(int(binsize))/float(binsize)
    data = numpy.zeros( (3, len(points[0])) )
    logging.debug("cov")
    cov = numpy.convolve(points[COV], avgWindow, "same")
    logging.debug("mat")
    mat = numpy.convolve(points[MAT], avgWindow, "same")
    data[0] = mat/cov
    del(mat)
    
    logging.debug("ter")
    ter = numpy.convolve(points[TER], avgWindow, "same")
    data[1] = ter/cov
    del(ter)
    
    logging.debug("dif")
    dif = numpy.convolve(points[DIFF], avgWindow, "same")
    data[2] = dif/cov
    del(dif)
    
    del(cov)

def makeSimpleSvmPt_ecoli15(data, binsize):
    avgWindow = numpy.ones(binsize)/float(binsize)
    ret = numpy.zeros( (8, len(data[0])) )
    #cov = numpy.convolve(data[COV], avgWindow, "same")
    cov = data[COV]
    
    #errors = data[TER] / cov
    #errors = (numpy.convolve(data[TER], avgWindow, "same") )/ cov
    ret[0] = data[MAT]/data[COV]
    #mat = numpy.convolve(data[MAT], avgWindow, "same")/ cov
    ret[1] = data[MIS]/data[COV]
    #mis = numpy.convolve(data[MIS], avgWindow, "same")/ cov
    ret[2] = data[INS] / data[COV]
    #ins = data[INS]/data[COV] #numpy.convolve(data[INS], avgWindow, "same")/ cov
    ret[3] = data[INSZ] / data[COV]
    #insz = numpy.convolve(data[INSZ], avgWindow, "same")/ cov
    ret[4] = data[DEL]/data[COV]
    #dele = numpy.convolve(data[DEL], avgWindow, "same")/ cov
    ret[5] = data[TAI]/ cov
    #tail = numpy.convolve(data[TAI], avgWindow, "same")/ cov
    #diff = data[DIFF]/data[COV]
    ret[6] = (data[TER] - data[MAT]) /data[COV]
    #diff = numpy.convolve(data[DIFF], avgWindow, "same")/ cov
    
    insBas =  data[INS]/(data[INSZ]+1) 
    #insBas[numpy.where(insBas > binsize)] = binsize
    insBas = numpy.convolve(insBas, avgWindow, "same")
    ret[7] = insBas / cov
    
    
    return data


def makeSimpleSvmPtDict(points, binsize):
    data = makeSimpleSvmPt(points, binsize)
    ret = []
    for i in xrange(len(data[0])):
        d = {}
        for j in xrange(len(data)):
            d[j+1] = data[j][i]
        ret.append(d)
    return ret
            

def svmHotSpots(container, svmName, binsize, offset):
    """
    """
    machine = loadMachine(svmName)
    #this replaces writes to my container -- only memory efficient thing to do
    #avgWindow = numpy.ones(int(binsize))/float(binsize)
    #for i in xrange(len(container)):
        ##Don't need to do all...
        #container[i] = numpy.convolve(container[i], avgWindow, "same")
    """
    logging.info("Original")
    logging.debug(container[:,0])
    logging.info("Converting to SVM")
    points = makeSimpleSvmPtDict(container, binsize)
    logging.debug(points[0])
    logging.info("Normalizing Points")
    points = machine.normalizeInstances(points)
    logging.debug(points[0])
    logging.info("Running SVM Prediction")
    labels = numpy.array(machine.predict(points)[0])
    logging.info("Calling Ranges")
    logging.debug( labels )
    truth = labels != -1
    ranges = makePoints(truth, binsize)
    ret = summarizePoints(ranges, offset, container)
    logging.debug(ret[0])
    logging.info("found %d spots" % (len(ret)))
    """
    
    logging.info("Converting to SVM")
    points = makeSimpleSvmPt_ecoli15(container, binsize)
    logging.info("Normalizing Points")
    machine.normalize(points)
    logging.info("Running SVM Prediction")
    labels = machine.nppredict(points)
    logging.info("Calling Ranges")
    
    truth = labels != -1
    ranges = makePoints(truth, binsize)
    ret = summarizePoints(ranges, offset, container)
    logging.info("found %d spots" % (len(ret)))
    
    return ret
    
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
    
    avgWindow = numpy.ones(int(binsize))/float(binsize)
    
    slopWindow = numpy.zeros(binsize)
    slopWindow[:binsize/3] = -1
    slopWindow[-binsize/3:] = 1
    cov = numpy.convolve(data[COV], avgWindow, "same")
    
    mis = signalTransform(data[MIS], cov, slopWindow, avgWindow)
    truth = mis >= covThresh
    del(mis)
    
    #ins = numpy.convolve( data[INS] / data[COV], avgWindow, "same")
    #insBas = numpy.convolve(data[INSZ], avgWindow, "same")
    #insBas[numpy.where(insBas > binsize)] = binsize
    #insSz  = numpy.convolve(data[INSZ], avgWindow, "same")
    ins = signalTransform(data[INS], cov, slopWindow, avgWindow)
    truth = numpy.any([truth, ins >= covThresh], axis = 0)
    del(ins)
    
    dele = signalTransform(data[DEL], cov, slopWindow, avgWindow)
    truth = numpy.any([truth, dele >= covThresh], axis = 0)
    del(dele)
    
    tail = signalTransform(data[TAI], cov, slopWindow, avgWindow)
    truth = numpy.any([truth, tail >= threshPct], axis=0)
    del(tail)
    
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
        output.append( (end-start, start + offset, end + offset, \
                        p[COV], \
                        p[MAT], \
                        p[MIS], \
                        p[INS], \
                        p[DEL], \
                        p[TAI], \
                        p[MAQ]) )
                        #p[TER], \
                        #p[SIM], \
                        #p[ERR], \
                        #p[RATIO], \
                        #p[DIFF], \
                        #p[MAPQ], \
                        #p[INSZ]) )
    return output
      
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
    #parser.add_argument("-m", "--method", type=str, choices = \
    #                    ["SIM", "ERR", "RATIO", "DIFF", "MAPQ"], \
    #                    help="Method to call spots (RATIO) [not implemented]")

    parser.add_argument("-t", "--mintail", type=int, default=100, \
                        help=("Minimum number of bases in a tail before it's "
                              "considered discordant (100 bp)"))
    parser.add_argument("-T", "--maxtail", type=int, default=500, \
                        help=("Maximum number of bases in a tail before it no "
                              "longer contributes to the discordant count (500)"))
    parser.add_argument("-c", "--minCoverage", type=int, default=3, \
                        help="Minimum coverage for a spot call to be made (3)")
    parser.add_argument("-ep", "--errorPctThreshold",  type=float, default=1,
                        help="Error percent threshold (1))")
    #parser.add_argument("-em", "--errorMatThreshold", type=float, default=1.0,
                        #help="Error match ratio threshold (1.0)")
                        
    parser.add_argument("-b", "--binsize", type=int, default=50, \
                        help="binsize for window averaging (50)")
    #parser.add_argument("-B", "--binstep", type=int, default=1, \
                        #help="binstep for spot threshold calculations (1)")
    
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
    EPTHRES = args.errorPctThreshold
    #EMTHRES = args.errorMatThreshold
    CVTHRES = args.minCoverage
    BINSIZE = args.binsize
    #BINSTEP = args.binstep
    
    if not args.noCallSpots:
        hotspots = open(args.output+".spots", 'w')
        hotspots.write(("chrom\tsize\tstart\tend\t%s\n" % "\t".join(columns)))
    
    regions= []
    if args.region:
        fh = open(args.region,'r')
        for line in fh.readlines():
            data = line.strip().split('\t')
            regions.append((data[0], int(data[1]), int(data[2])))
        fh.close()
    else:
        for chrom, size in zip(bam.references, bam.lengths):
            regions.append((chrom, 0, size))
    
    h5Name = args.output + ".h5"
    results = h5py.File(h5Name, 'w')
    results.attrs["columns"] = columns
    results.close()

    totReads = 0
    for chrom, start, end in regions:
        size = end - start
        groupName =  "%s:%d-%d" % (chrom, start, end)
        logging.info("making container for %s (%d bp)" % (groupName, size))
        results = h5py.File(h5Name, 'a')
        out = results.create_group(groupName)
        
        out.attrs["reference"] = chrom
        out.attrs["start"] = start
        out.attrs["end"] = end
        out.attrs["size"] = size
        
        logging.info("parsing bam" )
        #I'll need to extract this part
        reads = bam.fetch(chrom, start, end)
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
            #svmachine = "hgsimsamp" 
            #svmachine = "ecoli_jul15"
            svmachine = "thresholding"
            logging.info("calling spots with svm %s" % (svmachine))
            #thresholding
            spots = callHotSpots(myData, EPTHRES, CVTHRES, BINSIZE, start)
            #spots = svmHotSpots(container, svmachine, BINSIZE, start)
            #Down to here for multiprocessing
        
            for i in spots:
                #size, start, end, cov, mat, mis, ins, del, tai, ter, sim, err
                hotspots.write(("%s\t%s\n" % (chrom, 
                                ("%d\t%d\t%d\t" + \
                                "\t".join(["%.3f"] * len(columns))) \
                                % i)) )
            hotspots.flush()
            
        results.flush()
        results.close()#my attempt to save memory
        myData = None
        container = None
        gc.collect()
    
    if not args.noCallSpots:
        hotspots.close()
    logging.info("finished %d reads" % (totReads))
    exit()
