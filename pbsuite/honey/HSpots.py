#!/usr/bin/env python

import sys, time, re, gc, argparse, logging, bisect, math, copy, shutil, multiprocessing
import pysam, numpy, h5py
from contextlib import contextmanager
from collections import defaultdict, Counter
from tempfile import NamedTemporaryFile

from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import M5File, revComp
from pbsuite.banana.Polish import consensus
from pbsuite.jelly.Assembly import blasr

VERSION = "14.10.14"
AUTHOR = "Adam English"

USAGE = """\
Detect 'smaller' SVs via measuring discordance between sample and reference in long reads.
"""

#########################
## --- Global Vars --- ##
#########################

#Biggest integer I want to deal with
BIGINT  = 2000
BIGINTY = numpy.float32
# NUMPY ARRAY HDF5 COLUMNS AND SIZES
COLUMNS = ["coverage", "matches", "insertions", "deletions"]
#Index of columns
COV  = 0
MAT  = 1
INS  = 2  
DEL  = 3  
#Must not exceed 300,000 data points
CHUNKSHAPE = (4, 70000)

##############################
## --- Helper Functions --- ##
##############################

#{{{ http://code.activestate.com/recipes/511478/ (r1)
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

class SpotResult():
    def __init__(self, chrom=None, out_start=None, start=None, in_start=None, \
                                in_end=None, end=None, out_end=None, \
                                svtype=None, size=None, tags=None):
        self.chrom = chrom
        self.out_start = out_start
        self.start = start
        self.in_start = in_start
        self.in_end = out_end
        self.end = end
        self.out_end = out_end
        self.svtype = "UNK" if svtype is None else svtype
        self.size = -1 if size is None else size
        self.tags = tags if tags is not None else {}
        self.varReads = []
   
    @classmethod
    def parseLine(cls, line):
        """
        Turn a line into a spot result
        """
        data = line.strip().split('\t')
        tags = {}
        for i in data[9].split(';'):
            k, v = i.split('=')
            try:
                v = int(v)
            except ValueError:
                try:
                    v = float(v)
                except ValueError:
                    pass
            tags[k] = v
        
        return cls(chrom=data[0], start=int(data[2]), end=int(data[5]), \
                   svtype=data[7], size=int(data[8]), tags=tags)
    
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
    
    def fetchbounds(self):
        """
        return (start,end) tuple of spot's boundaries
        """
        pnts = [x for x in [self.out_start, self.start, self.in_start, \
                            self.in_end, self.end, self.out_end] \
                            if x is not None]
        return min(pnts), max(pnts)
    
    def estimateSize(self):
        """
        estimate the sv's size by using the mean or bounds
        """
        if 'szMean' in self.tags:
            self.size = self.tags["szMean"]
        else:
            s,e = self.fetchbounds()
            self.size = e-s
    
    def qregstr(self):
        """
        returns quick region string chrom:start-end
        """
        return "%s:%d-%d" % (self.chrom, self.start, self.end)
    
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
        dat = [self.chrom, self.out_start, self.start, self.in_start, \
               self.in_end, self.end, self.out_end, self.svtype, self.size, \
               tag]

        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(*dat) \
                .replace("None", ".")

class HoneySpotH5():
    """
    Holds a HoneySport h5 file.
    keeps a semaphore to help prevent multiple access
    """
    
    def __init__(self, filename, version=None, columns=None, parameters=None, mode='w'):
        self.filename = filename
        if mode != 'r':
            self.results = h5py.File(self.filename, mode)
            self.results.attrs["version"] = version
            self.results.attrs["columns"] = columns
            self.results.attrs["parameters"] = parameters
            self.results.close()   
            
        self.lock = multiprocessing.Lock()

    @contextmanager
    def acquireH5(self, mode='r'):
        with self.lock:
            self.results = h5py.File(self.filename, mode)
            yield self.results
            self.results.close()
    
    def reopen(self, mode='r'):
        self.results = h5py.File(self.filename, mode)
        
    def close(self):
        self.results.close()
        
###################################
## --- Consumer/Task Objects --- ##
###################################

class Consumer(multiprocessing.Process):
    """
    Basic Consumer. Follow the two queues with your *args and **kwargs that should be sent
    to the task when __call__ 'd

    NOTE! args can't hold anything that isn't pickle-able for the subprocess
    """
    def __init__(self, task_queue, result_queue, bamName, referenceName, honH5): #*args, **kwargs):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.bam = pysam.Samfile(bamName)
        self.reference = pysam.Fastafile(referenceName)
        self.honH5 = honH5
        #self.args = args
        #self.kwargs = kwargs
        
    def run(self):
        try: 
            proc_name = self.name
            while True:
                next_task = self.task_queue.get()
                if next_task is None:
                    # Poison pill means shutdown
                    logging.info('Thread %s: Exiting\n' % proc_name)
                    self.task_queue.task_done()
                    break
                try:
                    next_task(self.bam, self.reference, self.honH5)
                except Exception as e:
                    logging.error("Exception raised in task %s - %s" % (next_task.name, str(e)))
                    next_task.failed = True
                    next_task.errMessage = str(e)
                
                self.task_queue.task_done()
                self.result_queue.put(next_task)
            
            return
        except Exception as e:
            logging.error("Consumer %s Died\nERROR: %s" % (self.name, e))
            return

#I can probably have a generic task class...
class ErrorCounter(): 
    """
    Counts all the errors in the region
    """
    def __init__(self, groupName, chrom, start, end, args):
        self.groupName = groupName
        self.chrom = chrom
        self.start = start
        self.end = end
        self.args = args
        self.name = groupName + ":ErrorCounter"
        self.failed = False
        self.errMessage = ""
    
    def countErrors(self, reads, offset, size, args):
        """
        Sum the errors over any particular reference base
        """
        container = numpy.zeros( ( len(COLUMNS), size ), dtype=BIGINTY )
        for align in reads:
            cigar = expandCigar(align.cigar)
            
            #get starts within region
            regionStart = 0 if align.pos < offset else align.pos - offset
            regionEnd = size if align.aend > (offset + size) else align.aend - offset
            
            #I'm always starting at the beginning of the read,
            # but the beginning may hit before my regionStart
            start = align.pos - offset
            
            #check covering bases
            container[COV, regionStart:regionEnd] += BIGINTY(1)
            #MAQ?
            
            #previous base was an insert prevent multiple base ins
            pins = False
            pinsz = 0
            pdel = False
            pdels = None
            
            def pinsLoad(start, size):
                if not pins:
                    return False, 0
                if size >= args.minIndelErr:
                    begin = max(0,start-(size/2))
                    end = min(start+(size/2), container.shape[1])
                    container[INS, begin:end] += BIGINTY(1)
                return False, 0
            
            def pdelLoad(startPos, curPos):
                if not pdel:
                    return False, None
                if curPos-startPos >= args.minIndelErr:
                    container[DEL, startPos:curPos] += BIGINTY(1)
                return False, None
            
            for code in cigar:
                if start < regionStart or start >= regionEnd:
                    if code != 1: 
                        start += 1
                    continue
                elif code == 0:
                    container[MAT, start] += BIGINTY(1)
                    start += 1
                    pins, pinsz = pinsLoad(start, pinsz)
                    pdel, pdels = pdelLoad(pdels, start)
                elif code == 1: #ins
                    pins = True
                    pinsz += 1
                    pdel, pdels = pdelLoad(pdels, start)
                elif code == 2: #del
                    if pdels is None:
                        pdel = True
                        pdels = start
                    start += 1
                    pins, pinsz = pinsLoad(start, pinsz)
        return container

    def __call__(self, bam, reference, honH5):
        """
        Takes a pysam.Samfile
        """
        logging.info("Starting %s" % (self.name))
        size = self.end - self.start
        regName =  "%s:%d-%d" % (self.chrom, self.start, self.end)
        
        logging.info("Making container for %s (%s %d bp)" % (self.groupName, regName, size))
                        
        logging.info("Parsing bam" )
        readCount = bam.count(self.chrom, self.start, self.end)
        reads = bam.fetch(self.chrom, self.start, self.end)
        if readCount == 0:
            logging.warning("No reads found in %s" % self.groupName)
            self.failed = True
            self.errMessage("No reads found in %s" % self.groupName)
            return
        else:
            logging.info("%d reads to parse in %s" % (readCount, self.groupName))
        
        myData = self.countErrors(reads, self.start, size, self.args)
        
        #request loc on honH5 and flush results
        with honH5.acquireH5('a') as h5dat:
            out = h5dat.create_group(self.groupName)
            
            out.attrs["reference"] = self.chrom
            out.attrs["start"] = self.start
            out.attrs["end"] = self.end
            
            if size < CHUNKSHAPE[1]:
                chunk = (CHUNKSHAPE[0], size-1)
            else:
                chunk = CHUNKSHAPE
            
            container = out.create_dataset("data", data = myData, \
                                    chunks=chunk, compression="gzip")
            h5dat.flush()

class SpotCaller():
    """
    Takes a full matrix from ErrorCounter and calls/filters spots
    """
    def __init__(self,  groupName, chrom, start, end, args):
        self.groupName = groupName
        self.chrom = chrom
        self.start = start
        self.end = end
        self.args = args
        self.name = groupName + ":SpotCaller"
        self.failed = False
        self.errMessage = ""

    def preprocessSignal(self, signal, coverage):
        """
        Normalize and print stats returning data and it's std
        """
        rate = numpy.convolve(signal/coverage, self.avgWindow, "same")
        rate[numpy.any([numpy.isinf(rate), numpy.isnan(rate)], axis=0)] = 0
        mu = numpy.mean(rate)
        sd = numpy.std(rate)
        logging.info("RateMean %f  -- RateStd  %f" % (mu, sd))
        return rate, mu, sd
        
    def callHotSpots(self, data, offset, bam, args): #threshPct, covThresh, binsize, offset):
        """
        """
        ret = []
        self.avgWindow = numpy.ones(args.binsize, dtype=numpy.float16)/float(args.binsize)
        #coverage
        cov = numpy.convolve(data[COV], self.avgWindow, "same")
        covTruth = numpy.all([cov >= args.minCoverage, cov <= args.maxCoverage], axis=0)
        logging.info("MaxCov:%d MeanCov:%d StdCov:%d MinCov:%d" \
                % (numpy.max(data[COV]), numpy.mean(data[COV]), \
                numpy.std(data[COV]), numpy.min(data[COV])))
        del(cov)
        
        #ins
        logging.info("INS processing")
        ins, mu, sd = self.preprocessSignal(data[INS], data[COV])
        #startPoints = self.makeSpots(numpy.all([ins*cov >= args.threshold, covTruth], axis=0), args.buffer)
        startPoints = self.makeSpots(numpy.all([ins >= mu+args.threshold*sd, covTruth], axis=0), args.buffer)
        for start, end in startPoints:
            mySpot = SpotResult(chrom=self.chrom, start=start, end=end, svtype="INS")
            if self.supportingReadsFilter(mySpot, bam, args):
                ret.append(mySpot)
        del(ins)
        
        #dele
        logging.info("DEL processing")
        dele, mu, sd = self.preprocessSignal(data[DEL], data[COV])
        #startPoints = self.makeSpots(numpy.all([dele*cov >= args.threshold, covTruth], axis=0), args.buffer)
        startPoints = self.makeSpots(numpy.all([dele >= mu+args.threshold*sd, covTruth], axis=0), args.buffer)
        for start, end in startPoints:
            mySpot = SpotResult(chrom=self.chrom, start=start, end=end, svtype="DEL")
            if self.supportingReadsFilter(mySpot, bam, args):
                ret.append(mySpot)
        del(dele)
        
        return ret
        
    def makeSpots(self, truth, buffer):
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
        #compress spots that are too close
        for start, end in points[1:]:
            if start - curEnd <= buffer:
                curEnd = end
            else:
                npoints.append((curStart, curEnd))
                curStart = start
                curEnd = end
        
        npoints.append((curStart, curEnd))
        return npoints
    
    def supportingReadsFilter(self, spot, bam, args):
        """
        filters insertions or deletion spots based on errors
        """
        if spot.svtype == "INS":
            errId = 1
            errLab = 'insertion'
        elif spot.svtype == "DEL":
            errId = 2
            errLab = 'deletion'
        else:#don't worry about other types
            return True
    
        begin, ending = spot.fetchbounds()
        buf = abs(begin-ending)*.5
        #begin -= abs(begin-ending)*.5
        #ending += abs(begin-ending)*.5
        
        reads = bam.fetch(str(spot.chrom), max(0, begin), ending)
        totSizes = []
        coverage = 0
        nReadsErr = 0
        #For tandem
        strandCnt = {True: 0, False: 0}
        
        #count reads and errSizes
        for read in reads:
            if not (read.pos < begin and read.aend > ending):
                continue
            
            mySize = 0
            coverage += 1
            start = read.pos - 1
            cigar = expandCigar(read.cigar)
            curSize = 0
            readHasErr = False
            
            for code in cigar: 
                if code != 1:
                    start += 1
                #must be in region
                if start < begin-buf:
                    continue
                if start >= ending+buf:
                    break
                
                if code == errId:
                    curSize += 1
                if curSize != 0 and code != errId:
                    if curSize >= args.minIndelErr:
                        readHasErr = True
                        mySize += curSize
                    curSize = 0
                
            if readHasErr and mySize >= args.minIndelSize:
                nReadsErr += 1
                totSizes.append(mySize)
                strandCnt[read.is_reverse] += 1
                spot.varReads.append(read.qname)
        
        spot.tags["strandCnt"] = "%d,%d" % (strandCnt[False], strandCnt[True])
        if len(totSizes) == 0:
            logging.debug("no %s found!? %s" % (errLab, str(spot)))
            return False # false - you should filter
        
        #if len(totSizes) < max(math.ceil(coverage * args.minIndelPct), args.minErrReads):
        if len(totSizes) < args.minErrReads:
            logging.debug("not large cnt %s found %s " % (errLab, str(spot)))
            return False
        
        totSizes.sort()
        totSizes = numpy.array(totSizes)
        mean = totSizes.mean()
        median = numpy.percentile(totSizes, 50)
        firstQ = numpy.percentile(totSizes, 25)
        thirdQ = numpy.percentile(totSizes, 75)
        
        logging.debug("PassFilt %s" % (str(spot)))   
        logging.debug("cov    %d" % coverage )
        logging.debug("size %d %s" % (len(totSizes), str(totSizes)))
        logging.debug("mean   %d" % mean )
        logging.debug("median %d" % median)
        logging.debug("firstQ %d" % firstQ)
        logging.debug("thirdQ %d" % thirdQ)
        
        spot.tags["coverage"] = coverage
        spot.tags["szCount"]  = int(nReadsErr)
        spot.tags["szMean"]   = int(mean)
        spot.tags["szMedian"] = int(median)
        spot.tags["sz1stQ"]   = int(firstQ)
        spot.tags["sz3rdQ"]   = int(thirdQ)
        return True
    
    def __call__(self, bam, reference, honH5):
        """
        """
        logging.info("Starting %s" % (self.name))
        with honH5.acquireH5('r') as h5dat:
            myData = numpy.array(h5dat[self.groupName]["data"])
        self.calledSpots = self.callHotSpots(myData, self.start, bam, self.args)
        
class ConsensusCaller():
    """
    For any particular spot, create a consensus
    """
    
    def __init__(self, spot, args):
        self.spot = spot
        self.args = args
        self.name = spot.qregstr() + ":ConsensusCaller"
        self.failed = False
        self.errMessage = ""
        
    def readTrim(self, read, start, end):
        """
        Trims a pysam.AlignedRead to only include the sequence that's aligned (or should be aligned)
        between start and end on reference
        returns the sequence and quality
        """
        score = 0
        if not read.is_unmapped:
            regTrim = 0
            upS = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
            dnS = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0
            
            trimS = None
            trimE = None
            if start > read.pos:
                for queryPos, targetPos in read.aligned_pairs:
                    if trimS is None and targetPos >= start:
                        trimS = queryPos
            else:
                score += abs(read.pos - start)
            if end < read.aend:
                for queryPos, targetPos in read.aligned_pairs[::-1]:
                    if trimE is None and targetPos <= end:
                        trimE = queryPos
            else:
                score += abs(read.aend-end)
            
            if trimS is not None:
                trimS = max(0, trimS) + upS
            else:
                trimS = 0
                    
            if trimE is not None:
                trimE = min(len(read.seq), trimE)  - dnS
            else:
                trimE = len(read.seq)
            seq = read.seq[trimS:trimE]
            qual = read.qual[trimS:trimE]
            if not read.is_reverse:
                seq = seq.translate(revComp)[::-1]
                qual = qual[::-1]
        
        return seq, qual
    
    def consensusCalling(self, spot, bam, reference, args):
        """
        Make a consensus of all the reads in the region and identify all of the SVs in the region
        """
        MAXNUMREADS = 100 #I don't think we'll need more than this many reads
        SPANBUFFER = 300#number of bases I want a read to span
        
        chrom, start, end = spot.chrom, spot.start, spot.end
        buffer = args.buffer
        #work
        supportReads = []
        spanReads = []
        #Fetch reads and trim
        totCnt = 0
        for read in bam.fetch(chrom, start-buffer, end+buffer):
            if read.qname not in spot.varReads:
                continue
            seq, qual = self.readTrim(read, start-buffer, end+buffer)
            if read.pos < start-SPANBUFFER and read.aend > end+SPANBUFFER:
                spanReads.append((len(seq), seq, qual))
            else:
                supportReads.append((seq, qual))
            totCnt += 1
            
        if len(spanReads) == 0:
            logging.info("noone spans - consensus aborted. %s" % (str(spot)))
            spot.tags["noSpan"] = True
            return [spot]
            
        spanReads.sort(reverse=True)
        refread = spanReads[0]
        logging.debug("%d reads %d support" % (totCnt, len(supportReads)))
        if len(spanReads) + len(supportReads) > MAXNUMREADS:
            if len(spanReads) > MAXNUMREADS:
                supportReads = [(x[1], x[2]) for x in spanReads[1:MAXNUMREADS]]
            else:
                supportReads = supportReads[:MAXNUMREADS-len(spanReads)]
                supportReads.extend([(x[1], x[2]) for x in spanReads[1:]])
        else:
            supportReads.extend([(x[1], x[2]) for x in spanReads[1:]])
        
        
        #read that spans most of the region goes first
        #use the rest for cleaning
        
        #building consensus sequence
        foutreads = NamedTemporaryFile(suffix=".fastq")
        for id, i in enumerate(supportReads):
            foutreads.write("@%d\n%s\n+\n%s\n" % (id, i[0], i[1]))
        foutreads.flush()
        
        foutref = NamedTemporaryFile(suffix=".fasta")
        foutref.write(">%s:%d-%d\n%s"%("ecoli", start, end, refread[1]))
        foutref.flush()
        
        alignOut = NamedTemporaryFile(suffix=".m5")
        
        blasr(foutreads.name, foutref.name, bestn=1, nproc=1, outname=alignOut.name)
        #shutil.copyfile(foutreads.name, "sup.fastq")
        #shutil.copyfile(foutref.name, "base.fasta")
        #shutil.copyfile(alignOut.name, "align.m5")
        if args.pbbanana:
            aligns = M5File(alignOut.name)
            con = ">con\n%s\n" % consensus(aligns).sequence
        else:
            logging.debug("pbdagcon")
            r, con, e = exe("pbdagcon -m 25 -c 1 -t 0 %s" % (alignOut.name), timeout=2)
            #logging.debug(str(r) + " - " + str(e))
            con = con[con.index("\n")+1:]
            #logging.debug("MySeq: " + con)
            #Check if con is blank
        
        conOut = NamedTemporaryFile(suffix=".fasta")
        conOut.write(con)
        conOut.flush()
        refOut = NamedTemporaryFile(suffix=".fasta")
        refOut.write(">%s:%d-%d\n%s\n" % (chrom, start, end, \
                    reference.fetch(chrom, start-buffer, end+buffer)))
        refOut.flush()
        
        #map consensus to refregion
        varSam = NamedTemporaryFile(suffix=".sam")
        cmd = "blasr %s %s -sam -bestn 1 -affineAlign -out %s" % (conOut.name, refOut.name, varSam.name)
        logging.debug(cmd)
        logging.debug(exe(cmd, timeout=5))
        
        foutreads.close()
        foutref.close()
        alignOut.close()
    
        #convert sam to bam
        input = pysam.Samfile(varSam.name)
        varBam = NamedTemporaryFile(suffix=".bam")
        output = pysam.Samfile(varBam.name, 'wb', template=input)
        nReads = 0
        for read in input:
            output.write(read)
            nReads += 1
        varSam.close()
        input.close()
        output.close()
        
        #do pileup for sequence
        pysam.sort(varBam.name, varBam.name[:-4])
        pysam.index(varBam.name)
        bam = pysam.Samfile(varBam.name, 'rb')
        
        mySpots = []
        for pos in bam.pileup():
            size = pos.pileups[0].indel
            if abs(size) < args.minIndelSize or size == 0:
                continue
            newspot = copy.deepcopy(spot)
            if size > 0 and spot.svtype == "INS":
                newspot.start = pos.pos + start - buffer
                newspot.end = pos.pos + start - buffer
                align = pos.pileups[0]
                newspot.tags["seq"] = align.alignment.seq[align.qpos : align.qpos + align.indel]
                newspot.size = size
                gt, gq = genotype(newspot)
                newspot.tags["GT"] = gt
                newspot.tags["GQ"] = gq
                mySpots.append(newspot)
            elif size < 0 and spot.svtype == "DEL":
                newspot.start = pos.pos + start - buffer
                newspot.end = pos.pos + abs(size) + start - buffer
                newspot.size = -size
                gt, gq = genotype(newspot)
                newspot.tags["GT"] = gt
                newspot.tags["GQ"] = gq
                mySpots.append(newspot)
        bam.close()
        varBam.close()
        #I need to filter the redundant spots...?
        #if len(mySpots) > 1:
            #filtSpots = []
            #knownSizes = Counter([newspot.size for newspot in mySpots])
            #for size in knownSizes:
                #if knownSizes[i] > 1:
            #for i in mySpots:
                
        logging.debug("%d consensus reads created %d spots" % (nReads, len(mySpots)))
        return mySpots
    
    def __call__(self, bam, reference, honH5):
        """
        """
        logging.info("Starting %s" % (self.name))
        self.newSpots = self.consensusCalling(self.spot, bam, reference, self.args)

def genotype(spot, priors=[['0/0',0.1], ['0/1',0.5], ['1/1',0.9]]):
    def log_choose(n, k):
        # swap for efficiency if k is more than half of n
        r = 0.0
        if k * 2 > n:
            k = n - k
            
        for  d in xrange(1,k+1):
            r += math.log(n, 10)
            r -= math.log(d, 10)
            n -= 1
        
        return r
    
    total = int(spot.tags["coverage"])
    alt = int(spot.tags["szCount"])
    ref = total - alt
    gtList = []
    for gt, p_alt in priors:
        gtList.append(log_choose(total, alt) + alt * math.log(p_alt, 10) + ref * math.log(1 - p_alt, 10))
    gt_idx = gtList.index(max(gtList))
    GL = gtList[gt_idx]
    gtList.remove(GL)
    GT = priors[gt_idx][0]
    GQ = -10 * (GL - sum(10**x for x in gtList))
    
    return GT, GQ

def parseArgs(argv, established=False):
    parser = argparse.ArgumentParser(prog="Honey.py spots", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    ioGroup = parser.add_argument_group("I/O Arguments")
    ioGroup.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    ioGroup.add_argument("--hon", metavar="HON.H5", type=str, default=None, \
                        help="HON.h5 containing Error data. Skips ErrorCouting.")
    ioGroup.add_argument("-r", "--region", type=str, default=None,\
                        help="Only call spots in region.bed")
    ioGroup.add_argument("--chrom", type=str, default=None, \
                        help="Only call spots on specified chromosomes (comma-separated) (%(default)s)")
    ioGroup.add_argument("-n", "--nproc", type=int, default=1, \
                        help="Number of processors to use (only for consensus) (%(default)s)")
    ioGroup.add_argument("-o", "--output", type=str, default=None, \
                        help="Basename for output (BAM.hon)")
    
    pGroup = parser.add_argument_group("Spot-Calling Threshold/Filtering Arguments")
    pGroup.add_argument("-b", "--binsize", type=int, default=100, \
                        help="binsize for window averaging (%(default)s)")
    pGroup.add_argument("-e", "--threshold",  type=float, default=3,
                        help="Minimum Spot Threshold (%(default)s)")
    pGroup.add_argument("-c", "--minCoverage", type=int, default=2, \
                        help="Minimum coverage of a region (%(default)s)")
    pGroup.add_argument("-C", "--maxCoverage", type=int, default=BIGINT, \
                        help="Maximum coverage of a region (%(default)s)")
    pGroup.add_argument("-m", "--minIndelErr", type=int, default=5,
                        help="Minimum size of an indel error to be counted (%(default)s)")
    pGroup.add_argument("-i", "--minIndelSize", type=int, default=50, \
                        help="Minimum indel SV size (%(default)s)")
    pGroup.add_argument("-E", "--minErrReads", type=int, default=5, \
                        help="Minimum number of reads with indel (max(minIndelPct,%(default)s))")
    pGroup.add_argument("--spanMax", type=int, default=2000, \
                        help="Maximum Size of spot to be called (%(default)s)")
    #pGroup.add_argument("-I", "--minIndelPct", type=float, default=0.20, \
                        #help="Minimum pct of reads with indel (max(%(default)s*cov,minErrReads)")
    
    aGroup = parser.add_argument_group("Consensus Arguments")
    aGroup.add_argument("--noConsensus", action="store_true", \
                        help="Turn off consensus calling, just report spots (False)")
    aGroup.add_argument("--buffer", default=1000, \
                        help="Buffer around SV to assemble (%(default)s)")
    aGroup.add_argument("--reference", default=None, type=str, \
                        help="Sample reference. Required with consensus calling (None)")
    aGroup.add_argument("--pbbanana", action="store_true", \
                        help="Use pbbanana for consensus. (default is pbdagcon)")
    aGroup.add_argument("--blasr", default="blasr", \
                        help="Path to blasr if it's not in the env")
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    
    args = parser.parse_args(argv)
    if args.maxCoverage > BIGINT:
        logging.error("Max Coverge must be less than %d" % (BIGINT))
        exit(0)
    
    #check bam is bamfile
    
    if args.output is None:
        #args.output = args.bam.filename[:-4]+".hon"
        args.output = args.bam[:-4]+".hon"
    
    if not args.noConsensus:
        if args.reference is None:
            logging.error("Reference is required with consensus calling")
            exit(0)
        #Check is fastafile
    if args.chrom is not None:
        args.chrom = args.chrom.split(',')
    return args

def run(argv):
    numpy.seterr(all="ignore")
    args = parseArgs(argv)
    setupLogging(args.debug)
    
    bam = pysam.Samfile(args.bam)
    try:
        if bam.header["HD"]["SO"] != "coordinate":
            logging.warning("BAM is not sorted by coordinates! Performance may be slower")
    except KeyError:
        logging.warning("Assuming BAM is sorted by coordinate. Be sure this is correct")
   
    hotSpots = open(args.output+".spots", 'w')
    hotSpots.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
    
    regions = {}
    if args.region:
        fh = open(args.region,'r')
        for line in fh.readlines():
            data = line.strip().split('\t')
            if args.chrom is None or data[0] in args.chrom:
                regions[data[3]] = (data[3], data[0], int(data[1]), int(data[2]))
                #regions.append((data[3], data[0], int(data[1]), int(data[2])))
        fh.close()
    else:
        for chrom, size in zip(bam.references, bam.lengths):
            if args.chrom is None or chrom in args.chrom:
                regions[chrom] = (chrom, chrom, 0, size)
    
        
    tasks = multiprocessing.JoinableQueue();
    results = multiprocessing.Queue();
    
    if args.hon is not None:
        honH5 = HoneySpotH5(args.hon, mode='r')
        gotoQueue = results
    else:
        honH5 = HoneySpotH5(args.output + ".h5", \
                            VERSION, \
                            COLUMNS, \
                            str(args))
        gotoQueue = tasks

    #ErrorCounting
    errConsumers = [ Consumer(tasks, results, args.bam, args.reference, honH5) for i in xrange(args.nproc) ]
    for w in errConsumers:
        w.start()
           
    num_jobs = 0
    for groupName, chrom, start, end in regions.values():
        gotoQueue.put(ErrorCounter( groupName, chrom, start, end, args))
        num_jobs += 1
               
    while num_jobs:
        result = results.get()
        num_jobs -= 1
        
        if result.failed:
            logging.error("Task %s Failed (%s)" % (result.name, result.errMessage))
        
        elif result.name.endswith("ErrorCounter"):#counted -> call
            logging.info("%s -> spot" % (result.name))
            groupName, chrom, start, end = regions[result.groupName]
            tasks.put(SpotCaller( groupName, chrom, start, end, args ))
            num_jobs += 1
            
        elif result.name.endswith("SpotCaller"):#call -> consensus
            nspot = 0
            for spot in result.calledSpots:
                spot.estimateSize()
                if spot.size < args.minIndelSize or spot.size > args.spanMax:
                    continue 
                if args.noConsensus:
                    hotSpots.write(result.spot + '\n')
                else:
                    nspot += 1
                    tasks.put(ConsensusCaller(spot, args))
                    num_jobs += 1
            if nspot > 0:
                logging.info("%s -> %d consensus" % (result.name, nspot))
            
        elif result.name.endswith("ConsensusCaller"):#consensus -> finish
            logging.info("Task %s -> finish" % (result.name))
            for x in result.newSpots:
                hotSpots.write(str(x)+'\n')
    
    #Poison the Consumers.. I'm done with them
    for i in xrange(args.nproc):
        tasks.put(None)
    #Wait to leave -- Ensuring everything is finished
    #logging.info("Joining")
    #tasks.join()
    logging.info("Finished")

if __name__ == '__main__':
    run(sys.argv[1:])

