#!/usr/bin/env python

import sys, time, re, gc, argparse, logging, bisect
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
           "insertionsize", "deletions", "avgmapq", "label"]
COV  = 0
MAT  = 1
MIS  = 2  #4
INS  = 3  #8
INSZ = 4  #16  (which means 24 is best insertion)
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

def countErrors(reads, offset, size):
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
            if size >= 5:
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
                if not pins:
                    container[INS, start] += 1
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
    dat = numpy.convolve(dat, avgWindow, "same") #smooth
    dat = numpy.convolve(dat * orig, avgWindow, "same") # Take back to absolute space and smooth again
    return numpy.nan_to_num(dat)
    
def callHotSpots(data, threshPct, covThresh, binsize, offset):
    """
    """
    
    sumWindow = numpy.ones(binsize)
    avgWindow = numpy.ones(binsize)/float(binsize)
    
    slopWindow = numpy.zeros(binsize)
    slopWindow[:binsize/2] = -1
    slopWindow[-binsize/2:] = 1
    
    slopWindow = slopWindow / float(binsize)
    
    #coverage
    cov = numpy.convolve(data[COV], avgWindow, "same")
    covTruth = cov >= covThresh
    logging.info("MaxCov:%d MeanCov:%d StdCov:%d MinCov:%d" \
            % (numpy.max(data[COV]), numpy.mean(data[COV]), numpy.std(data[COV]), numpy.min(data[COV])))
    #mis
    mis = signalTransform(data[MIS], cov, slopWindow, avgWindow)
    misSpots = makeSpotResults(mis, binsize, threshPct, covTruth, "MIS")
    logging.info("MaxMis:%.3f MeanMis:%.3f StdMis:%.3f MinMis:%.3f" \
            % (numpy.max(mis), numpy.mean(mis), numpy.std(mis), numpy.min(mis)))
    #ins
    ins = signalTransform(data[INS], cov, slopWindow, avgWindow)
    logging.info("MaxIns:%.3f MeanIns:%.3f StdIns:%.3f MinIns:%.3f" \
            % (numpy.max(ins), numpy.mean(ins), numpy.std(ins), numpy.min(ins)))
    insSpots = makeSpotResults(ins, binsize, threshPct, covTruth, "INS")
    
    #insBas = signalTransform(data[INSZ]*ins, cov*binsize, slopWindow, avgWindow)
    insz = signalTransform(data[INSZ], cov, slopWindow, avgWindow)# <-- poor
    logging.info("MaxInsz:%.3f MeanInsz:%.3f StdInsz:%.3f MinInsz:%.3f" \
            % (numpy.max(insz), numpy.mean(insz), numpy.std(insz), numpy.min(insz)))
    inszSpots = makeSpotResults(insz, binsize, threshPct, covTruth, "INSZ")
    del(insz)
    del(ins)
    
    #dele
    dele = signalTransform(data[DEL], cov, slopWindow, avgWindow)
    logging.info("MaxDel:%.3f MeanDel:%.3f StdDel:%.3f MinDel:%.3f" \
            % (numpy.max(dele), numpy.mean(dele), numpy.std(dele), numpy.min(dele)))
    deleSpots = makeSpotResults(dele, binsize, threshPct, covTruth, "DEL")
    del(dele)
    
    #maq = signalTransform(data[MAQ], cov, slopWindow, avgWindow)
    #logging.info("MaxMaq:%.3f MeanMaq:%.3f StdMaq:%.3f MinMaq:%.3f" \
            #% (numpy.max(maq), numpy.mean(maq), numpy.std(maq), numpy.min(maq)))
    #spots = makeSpotResults(maq, binsize, threshPct, covTruth, "MAQ")
    #del(maq)
    
    #points = makePoints(truth, binsize)
    #output = summarizePoints(points, offset, data)
    
    return misSpots, insSpots, inszSpots, deleSpots
    
def makeSpotResults(points, binsize, threshPct, covThresh, label):
    """
    Find starts and ends ranges, then group the start/end pairs by closest proximity
    Need to explore the edge cases
    """
    entries = []
    startPoints = makePoints(numpy.all([points <= -threshPct, covThresh], axis = 0), binsize, 's')
    #logging.critical('sps' + str(startPoints))
    endPoints   = makePoints(numpy.all([points >= threshPct, covThresh], axis=0), binsize, 'e')
    #logging.critical('eps' + str(endPoints))
    sPos = 0
    ePos = 0
        
    points = []
    for i in startPoints:
        bisect.insort(points, i)
    for i in endPoints:
        bisect.insort(points, i)
    i = 0
    while i < len(points):
        if points[i][2] == 'e':
            entries.append((".", ".", points[i][0], points[i][1]))
            i += 1
        elif points[i][2] == 's':
            if i+1 < len(points) and points[i+1][2] == 'e':
                entries.append((points[i][0], points[i][1], points[i+1][0], points[i+1][1]))
                i += 2
            else:
                entries.append((points[i][0], points[i][1], ".", "."))
                i += 1
                       
    logging.info("%d %s entries" % (len(entries), label))
    ret = []
    for ent in entries:
        ret.append(("%s\t.\t%s\t%s\t.\t%s\t"+label) % ent)
    logging.info("\n".join(ret))
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
    
def summarizePoints(points, offset, container):
    """
    summarize the data in the points
    I want to add a filter in column 4 (also points should be added in the .h5)
    Okay.... the filter should be (T) if the call is only tails, but no cluter.
             the filter should be (*) if the calls doesn't have any problems
             the filter should be INSZ if the call appears to be a single large insertion
             the filter should be 

    """
    output = []
    lab = ""
    for start, end in points:
        p = container[:,start:end].mean(axis=1)
        output.append( (end-start, start + offset, end + offset, \
                        p[COV], \
                        p[MAT], \
                        p[MIS], \
                        p[INS], \
                        p[INSZ],\
                        p[DEL], \
                        p[MAQ]) )
    return output

def spotToString(chrom, spot):
    """
    """
    return "%s\t%s" % (chrom, 
          ("%d\t%d\t%d\t" + \
           "\t".join(["%.3f"] * 7)) % spot)

def parseArgs():
    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    
    parser.add_argument("-r", "--region", type=str, default=None,\
                        help="Only call spots in region.bed")
    
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Basename for output (BAM.hon)")
    #parser.add_argument("-z", "--noZmwDedup", action="store_false", \
                        #help="Allow tails from all subreads per ZMW (False)")
    
    parser.add_argument("-s", "--noCallSpots", action="store_true",\
                        help=("Don't call spots where error rates spike (False)"))

    #parser.add_argument("-t", "--mintail", type=int, default=100, \
                        #help=("Minimum number of bases in a tail before it's "
                              #"considered discordant (100 bp)"))
    #parser.add_argument("-T", "--maxtail", type=int, default=500, \
                        #help=("Maximum number of bases in a tail before it no "
                              #"longer contributes to the discordant count (500)"))
    
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
   
    EPTHRES = args.threshold
    CVTHRES = args.minCoverage
    BINSIZE = args.binsize
    
    if not args.noCallSpots:
        hotspots = open(args.output+".spots", 'w')
        hotspots.write("chrom\tout_start\tstart\tin_start\tin_end\tend\tout_end\ttype\n")
    
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
        #myData, numReads = countErrors(reads, start, size, MINTAIL, \
                                  #MAXTAIL, args.noZmwDedup)
        myData, numReads = countErrors(reads, start, size)
        
        container = out.create_dataset("data", data = myData, \
                                chunks = CHUNKSHAPE, compression="gzip")
        totReads += numReads
        
        if numReads == 0:
            logging.warning("No reads found in %s" % groupName)
            continue
        
        if not args.noCallSpots:
            logging.info("calling spots")
            misSpots, insSpots, inszSpots, deleSpots = callHotSpots(container, EPTHRES, CVTHRES, BINSIZE, start)
            fspots = len(misSpots) + len(insSpots) + len(inszSpots) + len(deleSpots)
            #Down to here for multiprocessing
            logging.info("found %d spots" % (fspots))
            totSpots += fspots
            for s in misSpots:
                hotspots.write((chrom+"\t%s\t.\t%s\t%s\t.\t%s\tMIS\n") % s)
            for s in insSpots:
                hotspots.write((chrom+"\t%s\t.\t%s\t%s\t.\t%s\tINS\n") % s)
            for s in inszSpots:
                hotspots.write((chrom+"\t%s\t.\t%s\t%s\t.\t%s\tINSZ\n") % s)
            for s in deleSpots:
                hotspots.write((chrom+"\t%s\t.\t%s\t%s\t.\t%s\tDEL\n") % s)
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
