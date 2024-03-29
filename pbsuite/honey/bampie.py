#!/usr/bin/env python
import os
import re
import sys
import logging
import tempfile
import argparse
from collections import defaultdict

import pysam

from pbsuite.utils.CommandRunner import exe, partition
from pbsuite.utils.FileHandlers import revComp
from pbsuite.utils.setupLogging import setupLogging


#Edit this string to set which parameters blasr will use by default
#DO NOT! Set -nproc, -bestn, -clipping, or any output (e.g. -out -m 5)
#Remove -noSpotSubreads if your inputs are bax.h5 files [i think]
BLASRPARAMS = (" --affineAlign --noSplitSubreads --nCandidates 20 " \
               "--minPctIdentity 75 --sdpTupleSize 6")
#Parameters used in the eichler experiments
EEBLASRPARAMS = (" --maxAnchorsPerPosition 100 --advanceExactMatches 10 " \
               "--affineAlign --affineOpen 100 --affineExtend 0 " \
               "--insertion 5 --deletion 5 --extend --maxExtendDropoff 20 " \
               "--noSplitSubreads --nCandidates 20 ")
               #"-minPctIdentity 75 ") #didn't use this, but maybe should



VERSION = "17.x"

USAGE="""\
Maps Reads Using BLASRPARAMS to produce .sam file.

If input is a input.fofn, .fastq, or .fasta, we do the initial
mapping and tail mapping for you all at once.

If input is a .bam/sam, we extract only the softclipped
bases from aligned reads and remap them to the provided reference.

If your input is a .sam your output will be a .sam and if your
input is a .bam your output will be a .bam

Edit the $SWEETPATH/pbsuite/honey/bampie.py variable BLASRPARAMS
to change how reads are mapped.
"""
def checkBlasrParams(bp):
    """
    Ensure -bestn, -nproc, -clipping, -out are not specified
    """
    args = [" --bestn ", " --nproc ", " --clipping ", " --out ", " -m "]
    for i in args:
        if bp.count(i):
            logging.error("Do not specify %s through Honey.py pie", i)
            exit(1)
        #I have a problem here

def callBlasr(inFile, refFile, params, nproc=1, outFile="map.sam", stride=None):
    """
    fq = input file
    automatically search for .sa
    """
    if os.path.exists(refFile + ".sa"):
        sa = "-sa " + refFile + ".sa"
    else:
        sa = ""
    #Warning about large files
    sz = os.path.getsize(inFile)
    if sz > 4e9:
        if stride is None:
            logging.warning(("%.2f GB sequence file %s may cause "
                             "blasr memory problems. Consider using "
                             "--stride option") % (sz/1e9, inFile))
        elif sz/stride[1] > 4e9:
            logging.warning(("%.2f GB sequence file per stride %s may cause "
                             "blasr memory problems. Consider increasing "
                             "--stride option") % ((sz/stride[1])/1e9, inFile))
            
    logging.info("Running Blasr")
    cmd = ("blasr %s %s %s --nproc %d --bestn 1 "
           "--sam --clipping subread --out %s ") \
           % (inFile, refFile, sa, nproc, outFile)
    
    if stride != None:
        start, stride = stride
        cmd += "--start %d --stride %d " % (start, stride)
    #handle multiple versions of blasr
    #r, o, e = exe("blasr -version")
    #version = float(o.strip().split('\t')[1])
    #logging.critical(cmd)
    #if version >= 5:
        #cmd = cmd.replace(' -', ' --')
    #logging.critical(cmd)

    logging.debug(cmd)
    r, o, e = exe(cmd + params)

    #r,o,e = exe(("blasr %s %s %s --nproc %d --sam --bestn 1 --nCandidates 20 "
                 #"--out %s --clipping soft --minPctIdentity 75 "
                 #" --noSplitSubreads") % (fq, ref, sa, nproc, out))

    if r != 0:
        logging.error("blasr mapping failed!")
        logging.error("RETCODE %d", r)
        logging.error("STDOUT %s", str(o))
        logging.error("STDERR %s", str(e))
        logging.error("Exiting")
        exit(r)

    logging.info(str([r, o, e]))

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
    fout = open(outFq, 'w')
    nreads      = 0
    ntails      = 0
    nmultitails = 0
    for read in bam:
        nreads += 1
        #in case qualities are not present
        if read.qual is None:
            read.qual = b"!"*len(read.seq)
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
                qal = read.qual[:length].decode()
            else:
                pos, tai = read.pos, 'e'
                seq = read.seq[:length].translate(revComp)[::-1]
                qal = read.qual[:length][::-1].decode()

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
                qal = read.qual[-length:].decode()
            else:
                pos, tai = read.aend, 'p'
                seq = read.seq[-length:].translate(revComp)[::-1]
                qal = read.qual[-length:][::-1].decode()
            maq = int(read.mapq)
            loc = mateplace + ":" + str(pos)
            #fout.write("@%s_%d%s%d%s\n%s\n+\n%s\n" % (read.qname, \
                       #maq, tai, strand, loc, seq, qal))
            fout.write("@%d%s%d%s__%s\n%s\n+\n%s\n" % \
                (maq, tai, strand, loc, read.qname, seq, qal))
    fout.close()
    return nreads, ntails, nmultitails

def uniteTails(mappedFiles, outBam="multi.bam"):
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

    #datGrab = re.compile("^(?P<rn>.*)_(?P<maq>\d+)(?P<log>[pe])(?P<strand>[01])(?P<ref>.*):(?P<pos>\d+)$")
    datGrab = re.compile("^(?P<maq>\d+)(?P<log>[pe])(?P<strand>[01])(?P<ref>.*):(?P<pos>\d+)__(?P<rn>.*)$")

    bout = None
    
    #create your bout and update header from first of the files: 
    for i in mappedFiles:
        if i[0] is not None:
            orig = pysam.Samfile(i[0], 'r')
            break
    
    if orig is None:
        logging.critical("No mapped files!")
        exit(1)
    
    header = orig.header
    #should consider changing RG information also...
    header["PG"].append({"ID":"bampie.py", "VN":VERSION,"CL": " ".join(sys.argv)})
    if outBam.endswith('.bam'):
        bout = pysam.Samfile(outBam, 'wb', header=header)
    else:
        bout = pysam.Samfile(outBam, 'wh', header=header)
    
    for ibam, tbam in mappedFiles:
        if tbam is not None:
            sam = pysam.Samfile(tbam, 'r')
        else:
            sam = []
        
        #build lookup of tails
        checkout = defaultdict(list)
        nmapped = 0
        for read in sam:
            nmapped += 1
            readData = read.qname
            #trusting this doesn't fail
            data = datGrab.search(readData).groupdict()
            read.qname = data["rn"]
            read.set_tag("IR", data["ref"], 'Z')
            read.set_tag("IP", int(data["pos"]), 'i')
            read.set_tag("II", int(data["strand"]), 'i')
            read.set_tag("IQ", int(data["maq"]), 'i')
            #all = read.tags + [("IR", data["ref"], 's'), ("IP", int(data["pos"]), 'i'), \
                        #("II", int(data["strand"]), 'i'), ("IQ", int(data["maq"]), 'i')]
            #logging.debug(all)
            #read.set_tags(all)

            ref = sam.getrname(read.tid)
            #primary or secondary
            strand = 1 if read.is_reverse else 0
            if data["log"] == 'p':
                read.flag += 0x40
                if strand == 1:
                    pos = int(read.pos)
                else:
                    pos = int(read.aend)
                code, length = read.cigar[-1]
                rmSeq = length if code == 4 else 0
            elif data["log"] == 'e':
                read.flag += 0x80
                if strand == 1:
                    pos = int(read.aend)
                else:
                    pos = int(read.pos)
                code, length = read.cigar[0]
                rmSeq = length if code == 4 else 0

            checkout[read.qname].append((data["log"], strand, ref, int(pos), int(read.mapq), rmSeq))
            bout.write(read)

        #Open my ibam for merging into bout
        if ibam.endswith('.bam'):
            origBam = pysam.Samfile(ibam,'rb')
        elif ibam.endswith('.sam'):
            origBam = pysam.Samfile(ibam)
        else:
            logging.error("Cannot open input file! %s", ibam)
            exit(1)

        #add information to the initial alignment
        for read in origBam:
            data = checkout[read.qname]
            if len(data) != 0:
                read.flag += 0x1
            for log, strand, ref, pos, maq, rmSeq in data:
                logging.debug("%s has tail %s", read.qname, log)
                try:
                    if log == 'p':
                        adding = [("PR", ref, "Z"), ("PP", pos, "i"), ("PI", strand, "i"), \
                                    ("PQ", maq, "i"), ("PS", rmSeq, "i")]
                    elif log == 'e':
                        adding = [("ER", ref, "Z"), ("EP", pos, "i"), ("EI", strand, "i"), \
                                    ("EQ", maq, "i"), ("ES", rmSeq, "i")]
                    for i in adding:
                        read.set_tag(*i)
                    #read.set_tags(all)
                except IndexError:
                    logging.critical("Index Error at Tag Addition!?")
                    logging.critical("Dataset will be missing a %s tail on read %s", \
                                     log, read.qname)
                    logging.critical("This is one of %d tails", len(data))
                    logging.critical("Tag: %s", read.tags)
                    logging.critical("Adding: %s", str(adding))
                except OverflowError:
                    logging.critical("Overflow Error at Tag Addition!?")
                    logging.critical("Dataset will be missing a %s tail on read %s", log, read.qname)
                    logging.critical("This is one of %d tails", len(data))
                    logging.critical("Values: log - %s, strand - %s, ref - %s, pos - %s, mapq - %d, rmSeq - %d", \
                                    log, strand, ref, pos, maq, rmSeq)
                    logging.critical("Tag: %s", read.tags)
                    logging.critical("Adding: %s", str(adding))

            bout.write(read)

    bout.close()
    return nmapped

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="Honey.py pie", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("input", metavar="[SAM,BAM,FASTA,FASTQ,FOFN]", type=str, \
                        help="Input reads to be mapped")
    parser.add_argument("reference", metavar="REFERENCE", type=str,\
                        help="Reference to map tails")

    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output Name (BAM.tails.[sam|bam])")
    parser.add_argument("-t", "--minTail", type=int, default=100,\
                        help="Minimum tail length to attempt remapping (%(default)s)")
    parser.add_argument("-n", "--nproc", type=int, default=1,\
                        help="Number of processors to use (%(default)s)")
    parser.add_argument("-p", "--params", type=str, default=BLASRPARAMS, \
                        help="Specify custom blasr params. use -p=\"string\"")
    parser.add_argument("-s", "--stride", type=int, default=1, \
                        help=("Break the input fasta/q into stride pieces before "
                              "alignment to reduce the chances of blasr using "
                              "too much memory (%(default)s)"))
    parser.add_argument("--temp", type=str, default=tempfile.gettempdir(),
                        help="Where to save temporary files")

    parser.add_argument("--chunks", type=int, default=0, \
                        help=("Create N scripts containing commands to "
                              "each input of the fofn (%(default)s)"))
    parser.add_argument("--debug", action="store_true")

    args = parser.parse_args(argv)

    setupLogging(args.debug)
    checkBlasrParams(args.params)

    if args.output is None:
        ext = args.input[args.input.rindex('.'):]
        main = args.input[:args.input.rindex('.')]
        if ext in [".sam", ".bam"]:
            args.output = main + ".tails" + ext
        else:
            args.output = main + ".tails.sam"

    return args

def decipherInput(input, chunks=0, stride=1):
    """
    returns True if initial map needs to happen
    and list of inputFileNames
    in input.fofn and chunks, you'll have lists of lists
    """
    extension = input.split('.')[-1].lower()
    
    #Single Sam/Bam, only need to do tails
    if extension in ["bam", "sam"]:
        if chunks != 0:
            logging.error("chunks not applicable to %s files", extension)
            exit(1)
        return False, [input]
    
    if extension in ["fastq", "fasta", "fa", "fq", "h5"]:
        if chunks != 0:
            logging.error("chunks not applicable to %s files", extension)
            exit(1)
        
        return True, [input]*stride
            
    if extension == "fofn":
        inputs = [x.strip() for x in open(input)]
        if chunks != 0:
            return True, partition(inputs, chunks)
        else:
            return True, inputs
    logging.error("Do not recognize input. Exiting")
    exit(1)
    
def mapTails(bamFn, args):
    """
    Given a bamFn and argument, map the reads and 
    """
    if bamFn.endswith('.bam'):
        bam = pysam.Samfile(bamFn,'rb')
    elif bamFn.endswith('.sam'):
        bam = pysam.Samfile(bamFn)
    else:
        logging.error("Cannot open input file! %s", bamFn)
        exit(1)
    
    logging.info("Extracting tails")
    tailfastq = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False, dir=args.temp)
    tailfastq.close(); tailfastq = tailfastq.name
    r, t, m = extractTails(bam, outFq=tailfastq, minLength=args.minTail)
    
    logging.info("Parsed %d reads", r)
    logging.info("Found %d tails", t)
    logging.info("%d reads had double tails", m)
    if t == 0:
        logging.info("No tails -- short-circuiting")
        return None
    #I think I need to be able to chunk tails.
    tailmap = tempfile.NamedTemporaryFile(suffix=".sam", delete=False, dir=args.temp)
    tailmap.close(); tailmap = tailmap.name
    
    callBlasr(tailfastq, args.reference, args.params, args.nproc, tailmap)
    bam.close() #reset
    return tailmap

def run(argv):
    args = parseArgs(argv)

    steps, inputFiles = decipherInput(args.input, args.chunks, args.stride)
    try:
        index = argv.index("--chunks")
        argv.pop(index); argv.pop(index)
        argv.remove(args.input)
    except ValueError:
        pass

    if steps:
        #We need to do the full mapping
        mappedFiles = []
        for c, ifile in enumerate(inputFiles):
            if args.chunks != 0: #making commands
                fh = open("chunk%d.fofn" % (c), 'w')
                for indv in ifile:
                    fh.write(indv +'\n')
                fh.close()
                temp = list(argv)
                temp.insert(0, fh.name)
                print("Honey.py pie " + " ".join(temp))
            else:
                logging.debug("Mapping %s", ifile)
                #Need to put this in a tempFile
                outName = tempfile.NamedTemporaryFile(suffix="map%d.sam" % (c), \
                                                     delete=False, dir=args.temp)
                outName.close(); 
                outName = outName.name
                callBlasr(ifile, args.reference, args.params, args.nproc, outName, stride=(c, args.stride))
                mappedFiles.append(outName)
        if args.chunks != 0:#we've made the commands
            logging.info("Commands printed to STDOUT")
            exit(0)
    else:
        mappedFiles = inputFiles

    pairs = []
    logging.info("Mapping Tails")
    for file in mappedFiles:
        pairs.append((file, mapTails(file, args)))

    logging.info("Consolidating alignments")
    n = uniteTails(pairs, args.output)

    logging.info("%d tails mapped", n)
    logging.info("Finished!")

if __name__ == '__main__':
    run(sys.argv[:1])
