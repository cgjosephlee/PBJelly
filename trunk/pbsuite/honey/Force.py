#!/usr/bin/env python
import sys, argparse, re
from collections import defaultdict, namedtuple, Counter
import pysam
import pbsuite.honey.TGraf as tails
import pbsuite.honey.HSpots as spots
from pbsuite.utils.setupLogging import *

USAGE = """
Checks if there are any reads in the given .bam that support predicted SVs

Takes a .bed with the first 6 columns being:
    chrom  start  end  name  svtype size

- svtype must be one of DEL, INS, MIS
- size is what the SV's size is estimated to be.  
- DEL and MIS size should equal the sv's span (end - start). INS is the 
    number of inserted bases

If you have single breakpoint events (such as translocations) specify --bedPE
Your input.bed's 9 columns become:
    chrom1  start1  end1 orient1  chrom2  start2  end2 orient2  name  svtype size 

- orient is the directionality of the sequence leading upto the breakpoint (+/-)

RegionBuffer is the +- space in which you consider reads for support 
    around predicted sv
Half of region buffer gives the area where reads must land to support 
    predicted sv.

Results are an extra column appended to the end of in the format REF[TAILS|SPOTS]
    REF:
        True if we found evidence of the reference over the region
        False if we had the opportunity to support the reference, but didn't.
        ? if we didn't have the opportunity to support the reference
    TAILS:
        A comma-separated list of chr:start-end(svtype)size coordinates for
        reads that have interrupted-mapping support of the SV
    SPOTS:
        A comma-separated list of chr:start-end(svtype)size coordinates for
        reads that have discordant-mapping support of the SV
"""



class Variant():
    def __init__(self, chrom, start, end, svtype, size, read=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.svtype = svtype
        self.size = size
        self.read = read
    
    def __str__(self):
        return "%s:%d-%d(%s)%d" % (self.chrom, self.start, self.end, self.svtype, self.size)
    
class BedEntry():
    def __init__(self, chrom, start, end, name, svtype, size):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.svtype = svtype
        self.size = int(size)
    
    def __str__(self):
        return "\t".join([str(x) for x in [self.chrom, self.start, self.end, self.name, self.svtype, self.size]])    

    def __repr__(self):
        return "<BedEntry '%s'>" % (str(self).replace('\t',' '))

class BedPEEntry():
    def __init__(self, chrom1, start1, end1, orient1, chrom2, start2, end2, orient2, name, svtype, size):
        self.chrom1 = chrom1
        self.start1 = int(start1)
        self.end1 = int(end1)
        self.orient1 = orient1
        self.chrom2 = chrom2
        self.start2 = int(start2)
        self.end2 = int(end2)
        self.orient2 = orient2
        self.name = name
        self.svtype = svtype
        self.size = int(size)

    def __str__(self):
        return "\t".join([str(x) for x in [self.chrom1, self.start1, self.end1, self.orient1, \
                                           self.chrom2, self.start2, self.end2, self.orient2, \
                                           self.name, self.svtype, self.size]])
    
    def __repr__(self):
        return "<BedPEEntry '%s'>" % (str(self).replace('\t',' '))

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
    parser.add_argument("-q", "--minMapq", type=int, default=100, \
                        help="Minimum mapping quality of a read and it's tail to consider (%(default)s)")
    parser.add_argument("-m", "--minErr", type=int, default=5, \
                        help="Minimum ins/del error size to consider (%(default)s)")
    parser.add_argument("-a", "--asm", action="store_true", \
                        help="Input reads are high-quality contigs")
    parser.add_argument("-p", "--bedPE", action="store_true", \
                        help="Input bed file is bedPE - only tails searching will be performed")
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    args = parser.parse_args(args)
    setupLogging(args.debug)
    
    return args

#MonkeyWrenching of sorts
FakeBread = namedtuple("FakeBread", "refKey, uBreak, dBreak, uDir, dDir, is_reverse, uTail, dTail, read")
Fr = namedtuple("fakeread", "is_reverse")
def tailsPESearch(bam, bed, args):
    """
    Populate the answer dictionary by looking for tails through the bam.
    I will need to setup orientation information, through.

    Based on paired end model
    """
    # These are all the element I'll be expected to make
    # refKey, uBreak, dBreak, uDir, dDir, is_reverse, uTail, dTail
    j = [bed.chrom1, bed.chrom2]; j.sort(); refKey = "_".join(j)
    
    #PE -> <-
    #+/- to pick the breakpoint locations
    bp1 = bed.start1 if bed.orient1 == '-' else bed.end1
    bp2 = bed.start2 if bed.orient2 == '-' else bed.end2
    #choose orientation of the first read and automatically flip the second
    orient1 = '5' if bed.orient1 == '-' else '3'
    orient2 = '5' if bed.orient2 == '+' else '3'
            
    j = [(bp1, orient1), (bp2, orient2)]; j.sort(); uBP1, dBP2 = j
    uBreak, uDir = uBP1
    dBreak, dDir = dBP2
    
    is_reverse = bed.orient1 == bed.orient2 
    fr = Fr(is_reverse)
    
    #My main read I'll fake cluster with everything
    fakey = FakeBread(refKey, uBreak, dBreak, uDir, dDir, is_reverse, 'i', 'e', fr)
    logging.info(fakey)
    reads = []
    
    #First breakpoint - 
    fetchS = max(0, bed.start1 - args.regionbuffer)
    fetchE = min(bed.end1 + args.regionbuffer, bam.lengths[bam.references.index(bed.chrom1)])
    reads.extend([x for x in bam.fetch(bed.chrom1, fetchS, fetchE)])
    
    fetchS = max(0, bed.start2 - args.regionbuffer)
    fetchE = min(bed.end2 + args.regionbuffer, bam.lengths[bam.references.index(bed.chrom2)])
    reads.extend([x for x in bam.fetch(bed.chrom2, fetchS, fetchE)])
    
    points = tails.makeBreakReads(reads, getrname = bam.getrname)
    nears = []
    for key in points:
        if key != refKey:
            continue
        for read in points[key]:
            #search for the original
            if read.near(fakey):
                anno = read.annotate()
                if anno in ['TLOC', 'INV']:
                    anno = 'MIS'
                nears.append(Variant("%s<->%s" % (read.uRef, read.dRef), read.uBreak, read.dBreak, anno, read.estsize))
    
    logging.info("Found %d penears" % (len(nears)))
    return len(reads) > 1, nears
    
def tailsSearch(bam, bed, args):
    """
    Populate the answer dictionary by looking for tails 
    through the bam

    Returns a list of pbsuite.honey.TGraf.Bnode that support
    """
    fetchS = max(0, bed.start - args.regionbuffer)
    fetchE = min(bed.end + args.regionbuffer, bam.lengths[bam.references.index(bed.chrom)])
    points = tails.makeBreakReads(bam.fetch(bed.chrom, fetchS, fetchE), getrname = bam.getrname)
    reads = []
    anyCoverage = False
    for key in points:
        anyCoverage = True
        #eventually will need tloc work
        if key.split('_')[0] != bed.chrom:
            continue
        #eventually will need a reference allele check for tails
        for read in points[key]:
            anno = read.annotate()
            if anno in ['TLOC', 'INV']:
                anno = 'MIS'
            
            #TLOCs...
            if bed.chrom != bam.getrname(read.read.tid):
                continue
            
            if anno != bed.svtype:
                #Not perfect..
                continue
            
            #within reciprocal ovl
            maxStart = max(bed.start, read.uBreak)
            minEnd   = min(bed.end, read.dBreak)
            if minEnd <= maxStart: #No overlap
                continue
            maxSpan  = max(bed.end-bed.start, read.dBreak - read.uBreak)
            recipOvl = abs(maxStart-minEnd) / float(maxSpan)
            logging.debug("predictVar [%d:%d] - tailRead [%d:%d]" \
                          % (bed.start, bed.end, read.uBreak, read.dBreak))
            if recipOvl < args.overlapbuffer:#not enough overlap
                continue
            
            anno  = read.annotate()    
            reads.append(Variant(bed.chrom, read.uBreak, read.dBreak, anno, read.estsize))
    
    #ret = ",".join(['t[%s]' % (str(x)) for x in reads])
    return anyCoverage, reads

def spotsSearch_asm(bam, bed, args):
    """
    find spots in high-accuracy contigs

    I'm going to have a problem with Insertion offsets before the variant
    if there are too many of them, I'm going to be effed
    """
    
    #MIS types can't be resolved from spots
    # EXCEPT, however, if they're actually INS/DEL 
    # except a single bp, they could
    if svtype == 'MIS':
        return False, '?',[]
    fetchS = max(0, bed.start - args.regionbuffer)
    fetchE = min(bed.end + args.regionbuffer, bam.lengths[bam.references.index(bed.chrom)])
    
    leeway = bed.size * args.sizebuffer
    ref = '?'
    vars = []
    anyCoverage = False
    for read in bam.fetch(bed.chrom, fetchS, fetchE):
        anyCoverage = True
        if read.pos > bed.start or read.aend < bed.end:
            #Not spanning our region
            logging.debug("%s doesn't span region" % (read.qname))
            continue
        ref = False
        if bed.size + leeway > 50000:
            #logging.warning("Variant is too long (%dbp), we are assuming reference" % (read.qname, len(read.seq)))
            continue
        #I'm going to need md if I get good a MIS
        cigar = spots.expandCigar(read.cigar)
        regionStart = max(read.pos, bed.start - (args.regionbuffer / 2))
        regionEnd = min(read.aend, bed.end + (args.regionbuffer / 2))
        readPosition = read.pos
        
        c = "".join([str(x) for x in cigar])
        logging.debug(c)
        logging.debug(c[regionStart-read.pos : regionEnd-read.pos])
        if svtype == 'INS':
            match = re.search("(^|[^1])1{%d,%d}([^1]|$)" % (bed.size-leeway, bed.size+leeway), c[regionStart-read.pos : regionEnd-read.pos])
        elif svtype == 'DEL':
            match = re.search("(^|[^2])2{%d,%d}([^2]|$)" % (bed.size-leeway, bed.size+leeway), c[regionStart-read.pos : regionEnd-read.pos])
        
        if match is None:
            ref = True
        else:
            #subtract insertion errors to correct the offset
            #Insertion offset subraction
            subtract = cigar[:(regionStart-read.pos) + match.start()].count(1)
            if bed.svtype == "INS":
                pos = match.start() + read.pos + (regionStart - read.pos) - subtract
                s,e = match.span(); size = e-s
                var = Variant(bed.chrom, pos, pos + 1, bed.svtype, size)
            if bed.svtype == "DEL":
                subtract = cigar[:(regionStart-read.pos)+ match.start()].count(1)
                spos = match.start() + read.pos + (regionStart -read.pos) - subtract
                epos = match.end() + read.pos + (regionStart - read.pos) - subtract
                s,e = match.span(); size = e-s
                var = Variant(bed.chrom, spos, epos, bed.svtype, size)

            vars.append(var)

    return anyCoverage, ref, vars
    
    
def spotsSearch(bam, bed, args):
    """
    take a pysam.Samfile and fetch reads in chrom/start/end region
    see if any reads support the call

    But this doesn't take into account that I have specific groupIds to use...
    """
    
    leeway = bed.size * args.sizebuffer
    
    fetchS = max(0, bed.start - args.regionbuffer)
    fetchE = min(bed.end + args.regionbuffer, bam.lengths[bam.references.index(chrom)])
    vars = []
    ref = '?'
    anyCoverage = False
    for read in bam.fetch(bed.chrom, fetchS, fetchE):
        if read.pos > bed.start or read.aend < bed.end:
            #Not spanning our region
            logging.debug("%s doesn't span region" % (read.qname))
            continue
        anyCoverage = True
        
        ref = False #we now have the opportunity to find the reference
         
        #I'm going to need md if I get good a MIS
        cigar = spots.expandCigar(read.cigar)
        regionStart = bed.start - (args.regionbuffer / 2)
        regionEnd = bed.end + (args.regionbuffer /2)
        readPosition = read.pos
        
       
        #classic finder
        numDel = 0
        numRef = 0
        before = False #need to see if we span
        after = False

        numIns = 0
        pins = False
        pinsz = 0
        curMd = 0
        
        pdel = False
        pdelz = 0
        #Find every insertion stretch
        #Find every deletion stretch
        #make a variant 
        #if any of the variants meet the criteria, we report it
        
        def pinsLoad(start, size):
            if not pins:
                return False, 0, 0
            rs = size if size >= args.minErr else 0
            return False, 0, rs

        def pdelLoad(start, size):
            if not pdel:
                return False, 0, 0
            rs = size if size >= args.minErr else 0
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
                #Did we just finish an ins/del
                pins, pinsz, t = pinsLoad(readPosition, pinsz)
                numIns += t
                pdel, pdelz, t = pdelLoad(readPosition, pdelz)
                numDel += t
            elif code == 1: #ins
                pins = True
                pinsz += 1
                #did we just finish a del?
                pdel, pdelz, t = pdelLoad(readPosition, pdelz)
                numDel += t
            elif code == 2: #del
                pdel = True
                pdelz += 1
                readPosition += 1
                #just finished an insertion
                pins, pinsz, t = pinsLoad(readPosition, pinsz)
                numIns += t
        pins, pinsz, t = pinsLoad(readPosition, pinsz)
        numIns += t
        pdel, pdelz, t = pdelLoad(readPosition, pdelz)
        numDel += t
        
        #what's the +- difference to validate the size
        leeway = bed.size * args.sizebuffer
        logging.debug(("regionStart, regionEnd, estSize, leeway, "
                       "estSize+leeway, estSize-leeway, numIns, numDel"))
        logging.debug("%d %d %d %d %d %d %d %d" % (regionStart, regionEnd, \
                        bed.size, leeway, size+leeway, size-leeway, numIns, numDel))
        
        if svtype == 'DEL':
            if bed.size + leeway >= numDel >= bed.size - leeway:
                vars.append(Variant(bed.chrom, bed.start, bed.end, "DEL", numDel))       
            else:
                ref = True
        elif svtype == 'INS':
            if bed.size + leeway >= numIns >= bed.size - leeway:
                vars.append(Variant(bed.chrom, bed.start, bed.end, "INS", numDel))       
            else:
                ref = True
        else:
            ref = True
    if len(vars) > 0:
        vars = [vars[0]]
    return anyCoverage, ref, vars

#Here are some helper methods for parsing force annotation results
forceRe = re.compile("(?P<ref>True|False|\?)\[(?P<tails>.*)\|(?P<spots>.*)\]")
def parseForce(data):
    """
    turns the force output into a dict
    """
    if data == 'no_cov' or data == '.':
        return None

    #will fail on malformed entries
    search = forceRe.search(data)
    if search is not None:
        d = search.groupdict()
    else:
        print "problem parsing", data
        return 'prob'

    if d["ref"] == '?':
        d["ref"] = None
    elif d["ref"] == 'True':
        d["ref"] = True
    elif d["ref"] == 'False':
        d["ref"] = False
    d["tails"] = [x for x in d["tails"].split(',') if x != '']
    d["spots"] = [x for x in d["spots"].split(',') if x != '']
    return d

def genoTyper(data):
    """
    """
    if data["ref"] is not None:
        if data["ref"]:
            if len(data["tails"]) > 0 or len(data["spots"]) > 0:
                genoType = "0/1"
            else:
                genoType = "0/0"
            
        elif not data["ref"]:
            if len(data["tails"]) > 0 or len(data["spots"]) > 0:
                genoType = "1/1"
            else:
                genoType = "./."
    else:
        if len(data["tails"]) > 0 or len(data["spots"]) > 0:
            genoType = "./1"
        else:
            genoType = "./."
    
    return genoType

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
        
        if not args.bedPE:
            myentry = BedEntry(*myentry)
        else:
            myentry = BedPEEntry(*myentry)
        
        if myentry.svtype not in vtypes:
            if myentry.svtype == 'UNK':
                logging.warning("Bed Entry %s is UNK and can't be forced... skipping" % (repr(myentry)))
                sys.stdout.write(line.strip() + "\t.\n")
                continue
            else:
                logging.error("Bed Entry %s svtype column isn't one of %s" % (repr(myentry), str(vtypes)))
                exit(1)
        
        if not args.bedPE:
            if myentry.chrom not in bam.references:
                logging.error("Invalid Chromosome %s" % myentry.chrom)
                continue
            
            anyCoverage1, tailVars = tailsSearch(bam, myentry, args)
            if args.asm:
                anyCoverage2, foundRef, spotVars = spotsSearch_asm(bam, myentry, args)
            else:
                anyCoverage2, foundRef, spotVars = spotsSearch(bam, myentry, args)
        
        else:
            anyCoverage1, tailVars = tailsPESearch(bam, myentry, args)
            foundRef = False
            anyCoverage2 = False
            spotVars = []
        
        #I'm outputing the variant reads and if we found the ref (True, False, ?) where ? means no evidence for or
        # against
        #
        if not anyCoverage1 and not anyCoverage2:
            annot = "no_cov"
        else:
            annot = "%s[%s|%s]" % (foundRef, ",".join([str(x) for x in tailVars]),\
                                             ",".join([str(x) for x in spotVars]))
        
        sys.stdout.write(line.strip() + "\t" + annot + "\n")
        i += 1
        if i % 250 == 0:
            sys.stdout.flush()

if __name__ == '__main__':
    run(sys.argv[1:])