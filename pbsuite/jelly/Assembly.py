#!/usr/bin/env python
import os, sys, argparse, logging, json
from tempfile import NamedTemporaryFile

from collections import namedtuple

from pbsuite.utils.setupLogging import setupLogging
from pbsuite.utils.FileHandlers import FastqFile, M5File, revComp
from pbsuite.jelly.Support import AlignmentConnector, SUPPORTFLAGS
from pbsuite.banana.Polish import *


def blasr(query, target, bestn=200, nproc = 1, outname = "out.m5"):
    """
    Simple overlapper
    """
    c = ("blasr %s %s -m 5 -bestn %d -nCandidates %d -minMatch 6 "
                 "-nproc %d -noSplitSubreads -out %s -minPctIdentity 70 -minReadLength 5") % \
                 (query, target, bestn, bestn*5, nproc, outname)
    logging.debug(c)
    r,o,e = exe(c)
    logging.debug("blasr - %d - %s - %s" % (r, o, e))

def extractFlanks(reads, basedir="./"):
    """
    Takes FastqFile and separates the the reference reads (^ref)
    from the supporting reads
    returns queryFileName, targetFileName
    """
    query = NamedTemporaryFile(suffix=".fastq", delete=False, dir=basedir)
    target = NamedTemporaryFile(suffix=".fasta", delete=False, dir=basedir)
    for read in reads:
        if read.startswith("ref"):
            target.write(">%s\n%s\n" % (read, reads[read].seq))
        else:
            query.write(reads[read].toString())
    query.close()
    target.close()
    return query.name, target.name


def orderSeeds(seedNames):
    """
    Looks at the seed's names to figure out
    which one is upstream of the next and if alignments 
    should be on the same strand
    """
    if len(seedNames) == 1:
        seedNames.append(None)
    
    seed1, seed2 = seedNames
    
    logging.debug("Ordering %s and %s" % (seed1, seed2))
    if seed1 == None:
        logging.error("Seed1 must be non-None to AssessAssembly!")
        exit(5)
    
    #I will be returning a None, just need to know
    #if seed1 is trying to extend 5' or 3'
    if seed2 == None:
        sameStrand = True
        if seed1.endswith("e3"):
            ret = (None, seed1)
        elif seed1.endswith("e5"):
            ret = (seed1, None)
    elif seed1.endswith("e3") and seed2.endswith("e5"):
        sameStrand = True
        ret = (seed1, seed2)
    elif seed1.endswith("e5") and seed2.endswith("e3"):
        sameStrand = True
        ret = (seed2, seed1)
    else:
        #if seed1.endswith("e5") and seed2.endswith("e5"):
        #if seed1.endswith("e3") and seed2.endswith("e3"):
        #No way to know. Someone is reverse-compliment of
        #the other. -- One needs to be on the opposite Strand
        sameStrand = False
        ret = (seed1, seed2)
        
    logging.debug(("Seed Order %s - %s : strand -" % ret) + \
                    str(sameStrand))
    return sameStrand, ret

def createStats():
    """
    I just wanted to separate the stats so It is a little cleaner
    """
    #span, seed1, seed2
    return {"support":              [[], [], []], #keep all the flags I have  \
            "spanCount":            0,  
            "spanSeedName":         None, 
            "spanSeedScore":        0, 
            "spanSeedStart":        None,
            "spanSeedEnd":          None, 
            "spanSeedStrand1":      None, 
            "spanSeedStrand2":      None, 
            "avgSpanBases":         0, 
            "seed1":                None,
            "seed2":                None, 
            "predictedGapSize":     None,
            "sameStrand":           None,
            "extendF1Count":        0, 
            "avgExtF1Bases":        0, 
            "extendF1SeedName":     0, 
            "extendF1SeedScore":    0, 
            "extendF1SeedStart":    None, 
            "extendF1SeedEnd":      None, 
            "extendF1SeedStrand":   None, 
            "extendF2Count":        0, 
            "avgExtF2Bases":        0, 
            "extendF2SeedName":     0, 
            "extendF2SeedScore":    0, 
            "extendF2SeedStart":    None, 
            "extendF2SeedEnd":      None, 
            "extendF2SeedStrand":   None, 
            "extendSeq1":           None, 
            "extendSeq2":           None, 
            "fillSeq":              None, 
            "contribSeqs":          0, 
            "contribBases":         0, 
            "fillBases":            0,
            "isLowQ":               False }
   
def getSubSeqs(alignmentFile, readsFile, sameStrand, seeds, basedir="./"):
    """
    Finds the seqs that align to the flanks the best, creates a fastq of supporting reads
    and the seed
    
    Might have a problem with my best read no going off the edge fully
    so I put the maxFlank at 20
    I should do more strand correction here
    """
    connector = AlignmentConnector()
    aligns = connector.parseAlignments(M5File(alignmentFile))
    reads = FastqFile(readsFile)
    
    stats = createStats()
    stats["seed1"], stats["seed2"] = seeds
    stats["sameStrand"] = sameStrand
    
    sfout  = NamedTemporaryFile(suffix=".fastq", delete=False, dir=basedir)
    f1fout = NamedTemporaryFile(suffix=".fastq", delete=False, dir=basedir)
    f2fout = NamedTemporaryFile(suffix=".fastq", delete=False, dir=basedir)
    
    logging.debug("span subseqs writing to %s" % sfout.name)
    logging.debug("flank1 subseqs writing to %s" % f1fout.name)
    logging.debug("flank2 subseqs writing to %s" % f2fout.name)
    
    bestSpan = None
    bestF1E  = None
    bestF2E  = None
    tooShort = []
    for readGroup in aligns:
        #hit on each flank
        if len(readGroup) == 2:
            r1, r2 = readGroup
            if r1.tname == stats["seed2"]:
                r1, r2 = r2, r1
            #need to ensure that everything is on the right strand
            a = connector.extendsTarget(r1, maxFlank=10)
            b = connector.extendsTarget(r2, maxFlank=10)
            #it extends both flanks
            if a != SUPPORTFLAGS.none and b != SUPPORTFLAGS.none:
                #Need to ensure that it's extending in the correct orientation
                
                logging.debug("%s spans" % r1.qname)
                logging.debug("hit1- (%d, %d)" % (r1.qstart, r1.qend))
                logging.debug("hit2- (%d, %d)" % (r2.qstart, r2.qend))
                
                rStart = min(r1.qend, r2.qend)
                rEnd = max(r1.qstart, r2.qstart)
                sz = rEnd - rStart
                if sz < 50:
                    logging.info("fill seq is too short to call consensus")
                    tooShort.append(sz)
                    continue
                #check for negative gaps -- though this is impossible to detect with
                #how I do support and concordance
                
                stats["spanCount"] += 1
                stats["avgSpanBases"] += rEnd - rStart
                stats["support"][0].append(SUPPORTFLAGS.span)
                
                t = reads[r1.qname].subSeq(rStart, rEnd)
                sfout.write(str(t))
                
                #is it the best spanner
                score = r1.score + r2.score
                if score < stats["spanSeedScore"]:
                    stats["spanSeedScore"] = score
                    spanSeedName = r1.qname
                    stats["spanSeedStrand1"] = r1.tstrand
                    bestSpan = reads[r1.qname].subSeq(rStart, rEnd)
                    stats["spanSeedName"] = r1.qname
                    stats["spanSeedStart"] = rStart
                    stats["spanSeedEnd"] = rEnd
                    stats["spanSeedStrand2"] = r2.tstrand


        elif len(readGroup) == 1: # single extender
            a = readGroup[0]
            sup = connector.extendsTarget(a, maxFlank=10)
            
            if sup != SUPPORTFLAGS.none:
                #set the coordinates of the extending sequence
                if a.tname.endswith("e5"):
                    if sup not in [SUPPORTFLAGS.left, SUPPORTFLAGS.span]:
                        continue
                    if a.qstrand == '0':
                        mystart = 0
                        myend   = a.qstart
                    else:
                        mystart = a.qseqlength
                        myend   = a.qend
                elif a.tname.endswith("e3"):
                    if sup not in [SUPPORTFLAGS.right, SUPPORTFLAGS.span]:
                        continue
                    if a.qstrand == '0':
                        mystart = a.qend
                        myend   = a.qseqlength
                    else:
                        mystart = 0
                        myend   = a.qstart

                #what flank and is it the best
                if a.tname.replace('/','.') == stats["seed1"]:
                    stats["extendF1Count"] += 1
                    stats["avgExtF1Bases"] += a.qstart
                    stats["support"][1].append( sup )
                    
                    if a.score < stats["extendF1SeedScore"]:
                        stats["extendF1SeedScore"] = a.score
                        stats["extendF1SeedName"] = a.qname
                        stats["extendF1SeedStart"] = mystart
                        stats["extendF1SeedEnd"] = myend
                        stats["extendF1SeedStrand"] = a.tstrand
                        bestF1E = reads[a.qname].subSeq(mystart, myend)
                    myOut = f1fout    
                elif a.tname.replace('/','.') == stats["seed2"]:
                    stats["extendF2Count"] += 1
                    stats["avgExtF2Bases"] += a.qstart
                    stats["support"][2].append( sup )
                    
                    if a.score < stats["extendF2SeedScore"]:
                        stats["extendF2SeedScore"] = a.score
                        stats["extendF2SeedName"] = a.qname
                        stats["extendF2SeedStart"] = mystart
                        stats["extendF2SeedEnd"] = myend
                        stats["extendF2SeedStrand"] = a.tstrand
                        bestF2E = reads[a.qname].subSeq(mystart, myend)
                    myOut = f2fout
                else:
                    logging.info("Failed to make an acceptable assembly for %s" % (str(a)))
                    myOut = NullDevice()
                #write to the flank out
                myOut.write(str(reads[a.qname].subSeq(mystart, myend)))
        else:
            logging.warning("read %d gave too many alignments" % (readGroup[0].qname))
    
    sfout.close()
    sfout = sfout.name
    f1fout.close()
    f1fout = f1fout.name
    f2fout.close()
    f2fout = f2fout.name

    logging.info("%d reads span" % stats["spanCount"])
    logging.info("%d reads extend flank 1" % stats["extendF1Count"])
    logging.info("%d reads extend flank 2" % stats["extendF2Count"])
    
    nt = namedtuple("SubInfo", "stats spanReads flank1Reads flank2Reads spanSeed flank1Seed flank2Seed")
    #seeds out files
    ssfout  = None
    f1sfout = None
    f2sfout = None
    
    #replace too short with N's
    if stats["spanCount"] == 0 and len(tooShort) > (stats["extendF1Count"] + stats["extendF2Count"])/2:
        stats["avgSpanBases"] = sum(tooShort)/len(tooShort)
        stats["spanCount"] = len(tooShort)
        logging.info("estimated fill len %d" % (stats["avgSpanBases"]))
        stats["fillSeq"] = "N"*stats["avgSpanBases"]
        stats["contribSeqs"] = len(tooShort)
        stats["contribBases"] = sum(tooShort)
        stats["fillBases"] = stats["avgSpanBases"]
        stats["spanSeedScore"] = -500
        stats["spanSeedStrand1"] = '0'
        stats["spanSeedStart"] = 0
        stats["spanSeedName"] = "tooShortNs"
        stats["spanSeedStart"] = 0
        stats["spanSeedEnd"] = stats["avgSpanBases"]
        ret = nt(stats, None, None, None, None, None, None)
        return ret

    if stats["spanCount"] > 0:
        stats["avgSpanBases"] = stats["avgSpanBases"]/stats["spanCount"]
        logging.info("estimated fill len %d" % (stats["avgSpanBases"]))
        #write seed
        if len(bestSpan.seq) < 50:
            logging.warning("fill sequence is small (%dbp) can't call consensus" % (len(bestSpan.seq)))
            #I don't know what to return here
            
        ssfout = NamedTemporaryFile(prefix="span", suffix=".fasta", delete=False, dir=basedir)
        logging.debug("spanning with %s" % (bestSpan.name))
        ssfout.write(">%s\n%s\n" % (bestSpan.name, bestSpan.seq))
        ssfout.close()
        ssfout = ssfout.name

    if stats["extendF1Count"] > 0:
        stats["avgExtF1Bases"] = stats["avgExtF1Bases"]/stats["extendF1Count"]
        logging.info("estimated flank 1 extend len %d" % (stats["avgExtF1Bases"]))
        #write seed
        if len(bestF1E.seq) < 50:
            logging.warning("f1e sequence is small (%dbp) can't call consensus" % (len(bestF1E.seq)))
            #I don't know what to return here
        f1sfout = NamedTemporaryFile(prefix="flank1", suffix=".fasta", delete=False, dir=basedir)
        f1sfout.write(">%s\n%s\n" % (bestF1E.name, bestF1E.seq))
        f1sfout.close()
        f1sfout = f1sfout.name
        
    if stats["extendF2Count"] > 0:
        stats["avgExtF2Bases"] = stats["avgExtF2Bases"]/stats["extendF2Count"]
        logging.info("estimated flank 2 extend len %d" % (stats["avgExtF2Bases"]))
        #write seed
        if len(bestF2E.seq) < 50:
            logging.warning("f2e sequence is small (%dbp) can't call consensus" % (len(bestF2E.seq)))
            #I don't know what to return here
        f2sfout = NamedTemporaryFile(prefix="flank2", suffix=".fasta", delete=False, dir=basedir)
        f2sfout.write(">%s\n%s\n" % (bestF2E.name, bestF2E.seq))
        f2sfout.close()
        f2sfout = f2sfout.name
    
    #all of the info I need to return... refactor later and create useful objects
    ret = nt(stats, sfout, f1fout, f2fout, ssfout, f1sfout, f2sfout)
    #seeds writing
    return ret
     
def buildFillSeq(data, args):
    """
    Using all of the information in the namedtuple returned from getSubSeqs, 
    go through the process of building the filling sequence.

    load the filling sequence in to the data
    """
    #try to build span
    if SUPPORTFLAGS.span in data.stats["support"][0]:
        logging.debug("build span")
        alignFile = os.path.join(args.asmdir, "spancon.m5")
        blasr(data.spanReads, data.spanSeed, bestn = 1, nproc = args.nproc, outname=alignFile)
        aligns = M5File(alignFile)
        if len(aligns) > 0:
            con = consensus(aligns)
            #if successful we're done
            if con.contribBases > 0 and con.fillBases > 0:#must be
                sequence = con.sequence#strandCorrector(data.stats["spanSeedStrand1"], con.sequence)
                data.stats["fillSeq"] = sequence
                data.stats["contribSeqs"] = con.contribSeqs
                data.stats["contribBases"] = con.contribBases
                data.stats["fillBases"] = con.fillBases
                return 
    
    #no span -- we need to do flanks
    flank1Success = False
    flank2Success = False
    fl1Flag = SUPPORTFLAGS.left if data.stats["seed1"].endswith("e5") else SUPPORTFLAGS.right
    fl2Flag = SUPPORTFLAGS.left if data.stats["seed2"].endswith("e5") else SUPPORTFLAGS.right
    logging.debug((fl1Flag, fl2Flag))
    if fl1Flag in data.stats["support"][1]:
        logging.debug("build flank1 %d" % fl1Flag)
        alignFile = os.path.join(args.asmdir, "flank1con.m5")
        blasr(data.flank1Reads, data.flank1Seed, bestn=1, nproc=args.nproc, outname=alignFile)
        aligns = M5File(alignFile)
        if len(aligns) > 0:
            con = consensus(aligns)
            if con.contribBases > 0 and con.fillBases > 0:#must be
                sequence = con.sequence#strandCorrector(data.stats["extendF1SeedStrand"], con.sequence)
                data.stats["extendSeq1"] = sequence
                data.stats["contribSeqs"] += con.contribSeqs
                data.stats["contribBases"] += con.contribBases
                data.stats["fillBases"] += con.fillBases
                flank1Success = True
    
    if fl2Flag in data.stats["support"][2]:
        logging.debug("build flank2 %d" % fl2Flag)
        alignFile = os.path.join(args.asmdir, "flank2con.m5")
        blasr(data.flank2Reads, data.flank2Seed, bestn=1, nproc=args.nproc, outname=alignFile)
        aligns = M5File(alignFile)
        if len(aligns) > 0:
            con = consensus(aligns)
            if con.contribBases > 0 and con.fillBases > 0:#must be
                sequence = con.sequence#strandCorrector(data.stats["extendF2SeedStrand"], con.sequence)
                data.stats["extendSeq2"] = sequence
                data.stats["contribSeqs"] += con.contribSeqs
                data.stats["contribBases"] += con.contribBases
                data.stats["fillBases"] += con.fillBases
                flank2Success = True
    
    if flank1Success and flank2Success:
        logging.debug("mid unite")
        seq = singleOverlapAssembly(data, args)
        if seq is not None:
            data.stats["fillSeq"] = seq
    
    return

def strandCorrector(strand, sequence):
    """
    ensures that the sequence inside of data is from the same strand as the 
    first seed
    if -, flip it
    """
    logging.debug("Weird %s" % (strand))
    if strand == '1':
        sequence = sequence.translate(revComp)[::-1]
    return sequence
    
def singleOverlapAssembly(alldata, args):
    """
    
    """
    data = alldata.stats
    reads = NamedTemporaryFile(suffix=".fasta", delete=False, dir=args.asmdir)
    e1Seq = data["extendSeq1"]; e2Seq = data["extendSeq2"]
    reads.write(">%s\n%s\n>%s\n%s\n" % ("seq1", e1Seq, "seq2", e2Seq))
    reads.close()
    
    alignFn = NamedTemporaryFile(suffix=".m5", delete=False, dir=args.asmdir)
    blasr(reads.name, reads.name, nproc=args.nproc, outname=alignFn.name)
    aligns = M5File(alignFn)
    # find best hit between the two
    connector = AlignmentConnector()
    bestS = None
    bestA = 0
    for i in aligns:
        if i.qname != i.tname: 
            if connector.extendsTarget(i):
                if i.score < bestS: 
                    bestA = i
                    bestS = i.score
    if bestS is None:
        logging.info("no overlap between extenders")
        return
    
    #any of these steps could fail -- 
    #Ensure the hit is valid
    #(if + + and sameStrand we are okay, if - + and not sameStrand we are okay)
    if data["sameStrand"] == (bestA.qstrand == '0'):
        logging.info("bad overlap between extenders")
        return
    
    con = consensus([bestA])
    bestA = bestA[0]
    #strand correction...
    if bestA.qname == "seq1":
        if bestA.tstrand == '1':
            e2Seq = e2Seq[:bestA.tstart].translate(revComp)[::-1]
            seq = e1Seq[:bestA.qstart] + con.sequence.translate(revComp)[::-1] + e2Seq
        else:
            seq = e1Seq[:bestA.qstart] + con.sequence + e2Seq[bestA.tend:]
    else:
        if bestA.tstrand == '1':
            e2Seq = e2Seq[:bestA.qstart].translate(revComp)[::-1]
            seq = e1Seq[:bestA.tstart] + con.sequence + e2Seq
        else:
            seq = e1Seq[:bestA.qstart] + con.sequence + e2Seq[bestA.tstart:]
            
    return seq

    
    
def parseArgs():
    """
    input dir
    predicted gapsize
    if argument says that we need to extract the seeds we will have a single paramters
        extractFlanks
    """
    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("asmdir", metavar="DIR", type=str, \
                        help="Local assembly directory for a gap")
    parser.add_argument("-p", "--predictedGapSize", type=int, default=None)
    parser.add_argument("-n", "--nproc", type=int, default=1)
    parser.add_argument("--debug", action="store_true")
    
    args = parser.parse_args()

    if args.asmdir.endswith("/"):
        args.asmdir = args.asmdir[:-1]
    setupLogging(args.debug)

    return args

def run():
    args = parseArgs()
    
    dirName = os.path.basename(args.asmdir)
    sameStrand, seeds = orderSeeds(dirName.split('_'))
    
    inputReads = FastqFile(os.path.join(args.asmdir,"input.fastq"))
    supportFn, flankFn = extractFlanks(inputReads, basedir=args.asmdir)
    
    onFlank = os.path.join(args.asmdir, "onFlank.m5")
    blasr(supportFn, flankFn, bestn=2, nproc=args.nproc, outname=onFlank)
    
    data = getSubSeqs(onFlank, supportFn, sameStrand, seeds, basedir=args.asmdir)
    
    if data.stats["spanSeedName"] != "tooShortNs":
        buildFillSeq(data, args)
    #if data.stats["support"][0] == SUPPORTFLAGS.span:
        #logging.info("spanned gap")
    #else:
        #logging.info("seed1 extend %d - seed2 extend %d" % tuple(data.stats["support"][1:]))
    data.stats["predictedGapSize"] = args.predictedGapSize
    jOut = open(os.path.join(args.asmdir, "fillingMetrics.json"),'w')
    jOut.write(json.dumps(data.stats,indent=4))
    jOut.close()
    logging.info("Finished")
    
    
        
if __name__ == '__main__':
    run()
