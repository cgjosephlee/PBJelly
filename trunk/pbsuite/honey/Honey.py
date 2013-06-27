#!/usr/bin/env python

import sys, time, re, pysam, numpy, argparse, logging, bisect, h5py
from collections import defaultdict
from multiprocessing import Pool

USAGE = """\
Count number of errors and thier position releative to the reference. 
"""

TODO = """ 
Make multiprocessing.pool.map on a per region basis - The trick will be getting the
    numpy and h5 data put together
"""

#container columns
COV = 0         #Coverage
MAT = 1         #Matches
MIS = 2         #Mismatches
INS = 3         #Insertions
DEL = 4         #Deletions
TAI = 5         #Tails
TER = 6         #Total Error
SIM = 7         #Match percent --> MAT / ( COV + TAI )
ERR = 8         #Error Percent --> TER / COV
RATIO = 9       #Error over Match Ratio --> TER / MAT
DIFF = 10       #TotalError - Matches --> TER - MAT
MAPQ = 11       #Average Map Quality 
INSZ = 12       #Number of bases inserted at position 

columns = ["coverage", "matches", "mismatches", \
           "insertions", "deletions", "tails", \
           "errors", "matchpct", "errorpct", \
           "errMatRatio", "difference", "avgMapQ", "numinsbases"]

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

def expandCigar(cigar):
    """
    Turns the abbreviated cigar into the full array
    
    0 = M
    1 = I
    2 = D
    """
    ret = []
    for t,s in cigar:
        #remove tails...
        if t not in [0,1,2]:
            continue
        ret.extend([t]*s)
    return ret

def expandMd(md):
    """
    Turns abbreviated MD into a full array
    """
    ret = []
    for i in re.findall("\d+|\^?[ATCGN]+", md):
        if i.startswith('^'):
            ret.extend(list(i[1:]))
        elif i[0] in ["A","T","C","G","N"]:
            ret.extend(list(i))
        else:
            ret.extend(['-']*int(i))
    return ret

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

def countErrors(reads, offset, size, mint, maxt, ignoreDups):
    """
    counts the base data in the read over the reference
    reads - the reads to process (sholud have been previously fetched)
    offset - the region's start
    size - size of the region
    mint - min tail len
    maxt - max tail contrib
    ignoreDups - don't count dups
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

def callHotSpots(container, errP_thresh, errM_thresh, cov_thresh, binsize, binstep, offset):
    """
    """
    def makePoints(truth, output):
        """
        make the points for this truth set
        """
        #prevent weirdness
        truth[-1] = False
        shift = numpy.roll(truth, 1)
        
        starts = truth & ~shift
        ends = ~truth & shift
        for start, end in zip(numpy.nonzero(starts)[0], numpy.nonzero(ends)[0]):
            p = container[:,start:end].mean(axis=1)
            output.append( (end-start, start + offset, end + offset, \
                            p[COV], \
                            p[MAT], \
                            p[MIS], \
                            p[INS], \
                            p[DEL], \
                            p[TAI], \
                            p[TER], \
                            p[SIM], \
                            p[ERR], \
                            p[RATIO], \
                            p[DIFF], \
                            p[MAPQ], \
                            p[INSZ]) )
    
    output = []
    window = numpy.ones(int(binsize))/float(binsize)
    #obs1 & obs2 & obs3
    #x$errorpct > 0.40 & x$errMatRatio >= 1 & x$difference >= coverage/2
    
    #for mismatches, deletions, and consistently placed insertions
    errpct = numpy.convolve(container[ERR], window, "same")
    errmatrat = numpy.convolve(container[RATIO], window, "same")
    errtruth = numpy.all([errpct >= errP_thresh, errmatrat >= errM_thresh], axis=0)
    
    #spurious insertions or (errmatch and errmatratio)
    #Right now I'm saying, are there more inserted bases in the window than what we'd expect
    #What I should say is if there are more inserted bases in the window per read
    #than what we'd expect -- that would involve remove the slideCov from expected
    #and dividing the slideSz j component by INS

    #for insertions that are spread across several base due to spurious matches
    #Nothing more than binsize -- 
    #this prevents single read insertion problems pretty well
    j = numpy.copy(container[INSZ])
    j[ numpy.where(container[INSZ] > binsize) ] = binsize 
    
    window2 = numpy.ones(binsize)#to get sum in window, not average
    #number of inserted bases per position per read
    slideInsCnt = numpy.convolve(container[INS], window2, "same")
    slideSz = numpy.convolve(j / slideInsCnt, window2, "same")
    
    slideCov = numpy.convolve(container[COV], window, "same")
    expected = (binsize * 0.15) * slideCov # more than the 15% err rate (should be a paramter)
    
    instruth = slideSz > expected
    
    truth = numpy.any( [ errtruth, instruth ], axis=0 )
    #truth = numpy.all( [ truth, container[COV] > cov_thresh ], axis=0 )
    
    makePoints(truth, output)
    
    return output
    
    #And I think what I should be doing is
    #slideSz = numpy.convolve(container[INSZ]/container[INS], window2, "same")
    #expected = (binsize * 0.20)
    #then, what I'm asking is when am 
           
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

    parser.add_argument("-t", "--mintail", type=int, default=50, \
                        help=("Minimum number of bases in a tail before it's "
                              "considered discordant (50bp)"))
    parser.add_argument("-T", "--maxtail", type=int, default=sys.maxint, \
                        help=("Maximum number of bases in a tail before it no "
                              "longer contributes to the discordant count (inf)"))
    parser.add_argument("-c", "--minCoverage", type=int, default=3, \
                        help="Minimum coverage for a spot call to be made (3)")
    parser.add_argument("-ep", "--errorPctThreshold",  type=float, default=0.40,
                        help="Error percent threshold (0.40))")
    parser.add_argument("-em", "--errorMatThreshold", type=float, default=1.0,
                        help="Error match ratio threshold (1.0)")
                        
    parser.add_argument("-b", "--binsize", type=int, default=10, \
                        help="binsize for spot threshold calculations (10)")
    parser.add_argument("-B", "--binstep", type=int, default=1, \
                        help="binstep for spot threshold calculations (1)")
    
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    args = parser.parse_args()

    if args.output is None:
        args.output = args.bam[:-4]+".hon"
    
    return args

def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))

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
    EMTHRES = args.errorMatThreshold
    CVTHRES = args.minCoverage
    BINSIZE = args.binsize
    BINSTEP = args.binstep
    
    results = h5py.File(args.output + ".h5", 'w')
    results.attrs["columns"] = columns
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
    
    totReads = 0
    for chrom, start, end in regions:
        size = end - start
        groupName =  "%s:%d-%d" % (chrom, start, end)
        logging.info("making container for %s (%d bp)" % (groupName, size))
        out = results.create_group(groupName)
        
        out.attrs["reference"] = chrom
        out.attrs["start"] = start
        out.attrs["end"] = end
        out.attrs["size"] = size
        
        logging.info("parsing bam" )
        #I'll need to extract this part
        reads = bam.fetch(chrom, start, end)
        logging.debug(start)
        myData, numReads = countErrors(reads, start, size, MINTAIL, MAXTAIL, \
                                       args.noZmwDedup)
        totReads += numReads
        container = out.create_dataset("data", data=myData )
        
        if numReads == 0:
            logging.warning("No reads found in %s" % groupName)
            continue
        if args.noCallSpots:
            continue
        
        logging.info("calling spots")
        spots = callHotSpots(myData, EPTHRES, EMTHRES, CVTHRES, BINSIZE, BINSTEP, start)
        #Down to here for multiprocessing
        logging.info("found %d spots" % (len(spots)))
        
        for i in spots:
            #size, start, end, cov, mat, mis, ins, del, tai, ter, sim, err
            hotspots.write(("%s\t%s\n" % (chrom, 
                            ("%d\t%d\t%d\t" + \
                             "\t".join(["%.3f"] * len(columns))) \
                             % i)) )
        
        hotspots.flush()
        
        
    hotspots.close()
    results.close()
    logging.info("finished %d reads" % (totReads))
    exit()
