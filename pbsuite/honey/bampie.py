#!/usr/bin/env python
import os, re, argparse, logging, tempfile
from collections import defaultdict
import pysam

from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import revComp
from pbsuite.utils.setupLogging import setupLogging


USAGE="""\
Extracts softclip bases from aligned reads and remaps them to the provided reference. 
Produces a unified bam with reads containing updated information about tail-mapping.
"""

def noSplitSubreads(readName):
    """
    Blasr won't give MapQ scores on alignments when noSplitSubreads 
    is specified -- so I gotta clean read names myself
    This is fixed
    """
    return "/".join(readName.split('/')[:-1])
    
def extractTails(bam, outFq, minLength=100):
    """
    0x1  -- template has multiple segments in sequencing
    0x40 -- first segment in template
    0x80 -- last segment in template

    Tail names will get _[pe][01].*:\d+
    on the end to hold metadata of:
    _    -- A delimieter
    [01] -- Strand of primary hit
    [pe] -- either p for prolog or e for epilog (strand ignorant)
    .*   -- refName of primary hit
    \d+  -- refPos of primary hit
    prolog, epilog, and the position of its primary (chr:pos)

    returns - nreads processed, ntails found, nreads with two ended tails
    """
    fout = open(outFq,'w')
    nreads      = 0
    ntails      = 0
    nmultitails = 0
    for read in bam:
        nreads += 1
        code, length = read.cigar[0]
        mateplace = bam.getrname(read.tid) 
        strand = 1 if read.is_reverse else 0
        hasTail = False
        if code == 4 and length >= minLength:
            hasTail = True
            ntails += 1
            if strand == 0:
                pos, tai = read.pos, 'p'
                seq = read.seq[:length]
                qal = read.qual[:length]
            else:
                pos, tai = read.pos, 'e'
                seq = read.seq[:length].translate(revComp)[::-1]
                qal = read.qual[:length][::-1]

            maq = int(read.mapq)
            loc = mateplace + ":" + str(pos)
            fout.write("@%s_%d%s%d%s\n%s\n+\n%s\n" % (read.qname, \
                       maq, tai, strand, loc, seq, qal))
                    
        code, length = read.cigar[-1]
        if code == 4 and length >= minLength:
            if hasTail:
                nmultitails += 1
            ntails += 1
            if strand == 0:
                pos, tai = read.aend, 'e'
                seq = read.seq[-length:]
                qal = read.qual[-length:]
            else:
                pos, tai = read.aend, 'p'
                seq = read.seq[-length:].translate(revComp)[::-1]
                qal = read.qual[-length:][::-1]
            maq = int(read.mapq)
            loc = mateplace + ":" + str(pos)
            fout.write("@%s_%d%s%d%s\n%s\n+\n%s\n" % (read.qname, \
                       maq, tai, strand, loc, seq, qal))
    fout.close()
    return nreads, ntails, nmultitails
    
    
def mapTails(fq, ref, nproc=1, out="tailmap.sam"):
    """
    automatically search for .sa
    """
    if os.path.exists(ref+".sa"):
        sa = "-sa " + ref + ".sa"
    else:
        sa = ""
    r,o,e = exe(("blasr %s %s %s -nproc %d -sam -bestn 1 -nCandidates 20 "
                 "-out %s -clipping soft -minPctIdentity 75 -sdpTupleSize 6"
                 " -noSplitSubreads") % (fq, ref, sa, nproc, out))
    if r != 0:
        logging.error("blasr mapping failed!")
        logging.error("RETCODE %d" % (r))
        logging.error("STDOUT %s" % (str(o)))
        logging.error("STDERR %s" % (str(e)))
        logging.error("Exiting")
        exit(r)
    
    logging.info(str([r, o, e]))

def uniteTails(origBam, tailSamFn, outBam="multi.bam"):
    """
    Put the tails and original reads into a single bam.
    Add tags uniting the pieces

    every read comprises upto three pieces
        X->Y->Z
    or
        prolog->primary->epilog
    
    each piece has 3 tags added: (R) ref - (P) pos - (S) strand
    
    prolog and eplog will only point to the primary and the primary will point to both
    
    """
    datGrab = re.compile("^(?P<rn>.*)_(?P<maq>\d+)(?P<log>[pe])(?P<strand>[01])(?P<ref>.*):(?P<pos>\d+)$")
    
    sam = pysam.Samfile(tailSamFn,'r')
    bout = pysam.Samfile(outBam, 'wb', template=origBam)
    checkout = defaultdict(list)
    nmapped = 0
    for read in sam:
        nmapped += 1
        readData = read.qname
        #trusting this doesn't fail
        data = datGrab.search(readData).groupdict()
        read.qname = data["rn"]
        read.tags += [("IR", data["ref"]), ("IP", int(data["pos"])), \
                      ("II", int(data["strand"])), ("IQ", int(data["maq"]))]
        
        ref = sam.getrname(read.tid)
        #primary or secondary
        if read.is_reverse:
            strand = 1
            if data["log"] == 'p':
                read.flag += 0x40
                pos = int(read.pos)
                code, len = read.cigar[-1]
                rmSeq = len if code == 4 else 0
            elif data["log"] == 'e':
                read.flag += 0x80
                pos = int(read.aend)
                code, len = read.cigar[0]
                rmSeq = len if code == 4 else 0
        else:
            strand = 0
            if data["log"] == 'p':
                read.flag += 0x40
                pos = int(read.aend)
                code, len = read.cigar[-1]
                rmSeq = len if code == 4 else 0
            elif data["log"] == 'e':
                read.flag += 0x80
                pos = int(read.pos)
                code, len = read.cigar[0]
                rmSeq = len if code == 4 else 0
            
        checkout[read.qname].append((data["log"], strand, ref, pos, int(read.mapq), rmSeq))
        bout.write(read)
    
    #add information to the primary
    for read in origBam:
        data = checkout[read.qname]
        if len(data) != 0:
            read.flag += 0x1
        for log, strand, ref, pos, maq, rmSeq in data:
            if log == 'p':
                read.tags += [("PR", ref), ("PP", pos), ("PI", strand), ("PQ", maq), ("PS", rmSeq)]
            if log == 'e':
                read.tags += [("ER", ref), ("EP", pos), ("EI", strand), ("EQ", maq), ("ES", rmSeq)]
        bout.write(read)
    
    bout.close()
    return nmapped

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="pie", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    parser.add_argument("ref", metavar="REFERENCE", type=str,\
                        help="REFERENCE to map tails to")
    
    parser.add_argument("-t", "--minTail", type=int, default=100,\
                        help="Minimum tail length to attempt remapping (100)")
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output Name (BAM.tails.bam)")
    parser.add_argument("-n", "--nproc", type=int, default=1,\
                        help="Number of processors to use (1)")
    parser.add_argument("--temp", type=str, default=tempfile.gettempdir(),
                        help="Where to save temporary files")
    parser.add_argument("--debug", action="store_true")
    
    args = parser.parse_args(argv)
    if args.output is None:
        args.output = args.bam[:-4] + ".tails.bam"
    setupLogging(args.debug)
    return args
    
def run(argv):
    args = parseArgs(argv)
    bam = pysam.Samfile(args.bam,'rb')
    
    logging.info("Extracting tails")
    tailfastq = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False, dir=args.temp)
    tailfastq.close(); tailfastq = tailfastq.name
    r, t, m = extractTails(bam, outFq=tailfastq, minLength=args.minTail)
    
    logging.info("Parsed %d reads" % (r))
    logging.info("Found %d tails" % (t))
    logging.info("%d reads had double tails" % (m))
    if t == 0:
        logging.info("No tails -- Exiting")
        exit(0)
    
    logging.info("Mapping Tails")
    tailmap = tempfile.NamedTemporaryFile(suffix=".sam", delete=False, dir=args.temp)
    tailmap.close(); tailmap = tailmap.name
    mapTails(tailfastq, args.ref, nproc=args.nproc, out=tailmap)
    bam.close() #reset
    
    logging.info("Consolidating alignments")
    bam = pysam.Samfile(args.bam,'rb')
    n = uniteTails(bam, tailmap, args.output)
    logging.info("%d tails mapped" % (n))
    
if __name__ == '__main__':
    run(sys.argv[:1])
