#!/usr/bin/env python

import sys, logging, os, shutil, json, re, tempfile, copy
from optparse import OptionParser
from collections import defaultdict, namedtuple, deque

from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import *


USAGE="""%prog.py <inputDir> [--predictedGapSize] [--nproc]
Wraps Iteratively attempts to fill a gap by doing local assemblies.
"""

"""
Goals:
1) I want an object that can parse a level/ dir and know
    i) Layout of SeedContigs in Output Contigs
    ii) Number of reads and their names used per contig
    iii) Some qualiy score per contig
2) This object will replace the Remap/ReSupport Model  used now
3) Preserve the assessment performed in AssemblyAssessor.py - above object will work with it.
    Actually, mainly redo it.
4) Similarily, simplify the backtracking -- I believe these portions are tied in together
"""



class BankOutput():
    ContigLineRe = re.compile(("^#(?P<readName>.+?)\((?P<offset>\d+)\) " 
                               "\[(?P<rc>(RC)?)\] (?P<len>\d+) bases, "
                               "\d+ checksum\. {(?P<qstart>\d+) "
                               "(?P<qend>\d+)} <(?P<cstart>\d+) "
                               "(?P<cend>\d+)>$"))
    def __init__(self, bankPath, seedNames):
        """
        Parses all associated information with a bank... 
        """
        self.bankPath = bankPath
        self.bank = os.path.join(bankPath,"out.bank")
        self.inputFastq = os.path.join(bankPath,"inputReads.fastq")
        self.outputFasta = os.path.join(bankPath,"out.fasta")
        self.outputQual = os.path.join(bankPath,"out.qual")
        self.seed1Name, self.seed2Name = \
                        self.__orderSeeds__(seedNames)
        self.contigSeed1 = None
        self.contigSeed2 = None
        self.hasSeed1 = False
        self.hasSeed2 = False
        self.__parseBank()
    
    def __orderSeeds__(self, seedNames):
        """
        Looks at the seed's names to figure out
        which one is upstream of the next
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
            self.sameStrand = True
            if seed1.endswith("e3"):
                ret = (None, seed1)
            elif seed1.endswith("e5"):
                ret = (seed1, None)
        elif seed1.endswith("e3") and seed2.endswith("e5"):
            self.sameStrand = True
            ret = (seed1, seed2)
        elif seed1.endswith("e5") and seed2.endswith("e3"):
            self.sameStrand = True
            ret = (seed2, seed1)
        else:
            #if seed1.endswith("e5") and seed2.endswith("e5"):
            #if seed1.endswith("e3") and seed2.endswith("e3"):
            #No way to know. Someone is reverse-compliment of
            #the other. -- One needs to be on the opposite Strand
            self.sameStrand = False
            ret = (seed1, seed2)
            
        logging.debug(("Seed Order %s - %s : strand -" % ret) + \
                        str(self.sameStrand))
        return ret
     
    def __parseBank(self):
        c,o,e = exe("bank2contig %s -L" % (self.bank))
        self.contigs = []
        if c != 0:
            logging.error("Problem parsing bank %s" % (self.bank))
            return
        
        for line in o.split('\n')[:-1]:
            if line.startswith("##"):
                contigName = line.split(' ')[0][2:]
                curContig = ContigLayout(contigName)
                self.contigs.append(curContig)
            elif line.startswith("#"):
                match = self.ContigLineRe.match(line)
                if match is None:
                    logging.error("Problem parsing bank %s line %s" \
                                %(self.bank, line))
                    exit(1)
                curContig.addRead(match.groupdict())
            else:
                logging.debug("BankInfoLine (%s)" % line)
        
        for contig in self.contigs:
            contig.identifySeeds(self.seed1Name, self.seed2Name)
            if contig.readSeed1 is not None:
                self.contigSeed1 = contig
                logging.debug("Contig '" + contig.contigName + "' hasSeed1")
                self.hasSeed1 = True
            if contig.readSeed2 is not None:
                self.contigSeed2 = contig
                self.hasSeed2 = True
                logging.debug("Contig '" + contig.contigName + "' hasSeed2")
            
        logging.info("Made %d contigs " % (len(self.contigs)))
    
    def replaceOutput(self, reads):
        """
        Takes a bunch of FastqEntries and replaces my Bank
        """
        os.rename(self.outputFasta, os.path.join(self.bankPath,"_out.fasta"))
        os.rename(self.outputQual, os.path.join(self.bankPath,"_out.qual"))
        fa = open(self.outputFasta, 'w')
        fq = open(self.outputQual, 'w')
        for read in reads:
            fa.write(">%s\n%s\n" % (read, reads[read].seq))
            fq.write(">%s\n%s\n" % (read, " ".join(map(str, reads[read].translateQual()))))
        fa.close()
        fq.close()


    def extractSequencesToFile(self, contigs=None, reads=None, outName="output.fastq"):
        """
        takes out all of the contigs and reads in the lists
        Contigs are reads that are in self.outputFasta/Qual
        Reads are reads that are in self.inputFastq
        puts them into output.fastq
        """
        output = open(outName,'w')
        seqs = self.extractSequences(contigs, reads)
        for key in seqs:
            output.write(seqs[key][1].toString())
        output.close()
    
    def extractSequences(self, contigs=None, reads=None):
        """
        Takes out all of the contigs and reads and returns 
        a dictionary of
            {name : [ContigLayout, FastqEntry], ...}
        """
        ret = {}
        if contigs is not None and len(contigs) > 0:
            logging.info("Extracting %d contigs from %s %s" % \
                         (len(contigs), self.outputFasta, str(contigs)))
            inputReads = mergeFastaQual(self.outputFasta, self.outputQual)
            for i in contigs:
                for c in self.contigs:
                    if c.contigName == i:
                        layout = c;
                ret[i] = [layout, inputReads[i]]
                
        if reads is not None and len(reads) > 0:
            logging.info("Extracting %d reads from %s %s" % \
                         (len(reads), self.inputFastq, str(reads)))
            inputReads = FastqFile(self.inputFastq)
            for i in reads:
                ret[i] = [None, inputReads[i]]
                
        return ret

class Read():
    """
    Named tuples are immutable...
    Read = namedtuple("Read", "")
    """
    def __init__(self, readName, offset, isRC, expandLen, qstart, qend, \
                 queryLen, cstart, cend):
        self.readName = readName
        self.offset = offset
        self.isRC = isRC
        self.expandLen = expandLen
        self.qstart = qstart
        self.qend = qend
        self.queryLen = queryLen
        self.cstart = cstart
        self.cend = cend

    def translate(self, translate):
        """
        Updates this Read's Coordinates so that
        only the initial seed read is identified
        within the entire contig.
        For example:
            First Iteration made contig A using seed S 
            (where = represents seed S)
            
            A   --==--

            Second Iteration made contig B using contig A
            (where . represents contig A)

            B  ----......--

            Inside of B, we need to figure out the novel sequence by
            translating A in the newspace back down to the S sequence

            B  ----..==..--
            
            Now, if we wanted upstream sequence U, we'd pull out

            U ----..
        """
        logging.debug("Translating %s to %s" % (self.readName, translate.readName))
        logging.debug("Me" + str(self))
        logging.debug("Them" + str(translate))
        if (not self.isRC and not translate.isRC):
            self.cstart = self.cend - (self.qend - translate.cstart) - 1
            self.cend = self.cstart + (translate.cend - translate.cstart)
        
        if (self.isRC and translate.isRC):
            cstart = self.cstart
            self.cstart = cstart + translate.cstart
            self.cend = cstart + translate.cend
        
        if (self.isRC and not translate.isRC):
            cend = self.cend
            self.cstart = cend - translate.cend
            self.cend = cend - translate.cstart
        
        if (not self.isRC and translate.isRC):
            cend = self.cend
            self.cstart = cend - translate.cend
            self.cend = self.cend - translate.cstart 
        
        logging.debug("Me Again" + str(self))
            #What the dodo do I do about trims in 2nd iter?
        #else:
            #if self.isRC:
                #self is backwards and need to get back
            #else:
                #translate is backwards and needs to get back
        
    def __str__(self):
        #ref0000348e5(1308) [RC] 5393 bases, 00000000 checksum. {5211 1} <1309 6409>
        return "#%s(%d) [%s] %d expand bases, %d length. {%d %d} <%d %d>"\
                % (self.readName, self.offset, "RC" if self.isRC else "", \
                   self.expandLen, self.queryLen, self.qstart, self.qend, \
                   self.cstart, self.cend)
             
class ContigLayout():
    """
    Holds the layout of the reads in a contig
    HelperMethods:
        Get Fill Sequence Coordinates
    """
    def __init__(self, name):
        self.contigName = name
        self.reads = []
        self.readSeed1 = None
        self.readSeed2 = None
        self.isDegenerate = False
        #So dumb, but bank2contig doesn't give the actual 
        #contig length, but instead gives the inflated in
        #full by-base layout length.. so I gotta do work
        self.contigLength = 0

    def addRead(self, gDict):
        """
        Parses a groupdict to add a read to the layout
        """
        name = gDict["readName"]
        offset = int(gDict["offset"])
        isRC = True if gDict["rc"] == "RC" else False
        expandLen = int(gDict["len"])
        qstart = int(gDict["qstart"]) 
        qend = int(gDict["qend"])
        
        if isRC:
            qstart, qend = qend, qstart
        qstart -= 1
        queryLength = qend - qstart
        
        cstart = int(gDict["cstart"]) - 1
        cend = int(gDict["cend"])
        if cend > self.contigLength:
            self.contigLength = cend
            
        self.reads.append(Read(name, offset, isRC, expandLen, qstart, \
                               qend, queryLength, cstart, cend))
    
    def getRead(self, readName):
        """
        Returns the Read tuple associated with readName
        """
        logging.debug("Looking for read %s in contig %s" % (readName, self.contigName))
        for i in self.reads:
            if i.readName == readName:
                return i
                
        return None
    
    
    def identifySeeds(self, seed1Name, seed2Name):
        """
        Checks how many of the seeds in the provided list are inside 
        of this Contig's layout
        """
        #Number of seeds I have
        self.seedCount = 0
        
        read = self.getRead(seed1Name)
        if read is not None:
            self.readSeed1 = read
            self.seedCount += 1
        read = self.getRead(seed2Name)
        if read is not None:
            self.readSeed2 = read
            self.seedCount += 1
        
        if self.seedCount == 0:
            self.isDegenerate = True
            
        logging.debug("%d seeds found" % self.seedCount)
    
    def getFillSeqCoords(self, sameStrand=False):
        """
         --|-----|--fill-seq--|-----|-
           *seed1*            *seed2*
        """
        if not sameStrand:
            #Why? Learn to Document, Adam.
            logging.debug("COULD BE A PROBLEM")
            
        logging.debug("ReadSeed1 : %s" % (str(self.readSeed1)))
        logging.debug("ReadSeed2 : %s" % (str(self.readSeed2)))
        end   = max(self.readSeed1.cstart, self.readSeed2.cstart) 
        start = min(self.readSeed1.cend,   self.readSeed2.cend)
        logging.debug("Start %d - End %d" % (start, end))
        #Need to figure out what to do on Trims...
        return (start, end)
    
    def getExtendSeqCoords(self, seed):
        """
        downstream == True
            --|----|--extnd-seq---
              *seed*

        downstream == False
            --extnd-seq-|----|--
                        *seed*
        Note, must have only one seed in self.seeds
        
        #e3 means I want downstream sequence downstream = True
        #e5 means I want upstream sequence downstream = False
        """
        if seed.readName.endswith("e3"):
            downstream = True
        elif seed.readName.endswith("e5"):
            downstream = False
        else:
            logging.critical(("Safety Net for Read ending in "
                               "e[35]! " + seed))
            exit(7)
                
        if seed.isRC:
            downstream = not downstream
            
        if downstream:
            start = seed.cend
            end = None
        elif not downstream:
            start = 0
            end = seed.cstart
            
        return (start, end)
        
#Assembly Success codes
ASMSUC = enum(fill         = 4,
              trim         = 3,
              doubleExtend = 2,
              singleExtend = 1,
              none         = 0)

class AssemblyAssessor():
    """
    Takes in a bank and evaluates the kind of gap-fill it
    has produced

    AssemblyAssessor.success = ASMSUC

    Todo:
        how cool would AssemblyAssessor > AssemblyAsseor be!?!
    """
    def __init__(self, bankOutput, seedReadsFastq, predictedGapSize=None):
        self.bank = bankOutput
        self.seedReadsFastq = seedReadsFastq
        self.pGapSize = predictedGapSize
        self.success = ASMSUC.none
        self.contigSeed1 = None
        self.contigSeed2 = None
        # -- For the entire assembly (every contig) track what seeds
        #    we have so that a BankOperator can know if something 
        #    needs to be retreived
        self.metrics = {}
        self.__assess()
    
       
    def __assess(self):
        """
        fillType (extendL, extendR, span, trimL, trimR)
        fillSize
        sequence (eLSeq, eRSeq, spanSeq, tLAmt, tRamt)
        sequence strand (I won't automatically turn to + strand this time)
        seed1 strand
        seed2 strand
        """
        for contigLayout in self.bank.contigs:
            logging.debug("Looking at Contig %s" % (contigLayout.contigName))
            if contigLayout.seedCount == 2:
                #ensure that the strandedness of
                #the seeds is correct
                s1isRC = contigLayout.readSeed1.isRC
                s2isRC = contigLayout.readSeed2.isRC
                
                self.contigSeed1 = contigLayout
                self.contigSeed2 = contigLayout
                
                if not self.bank.sameStrand and s1isRC == s2isRC:
                    logging.info(("We're on the same strand, but "
                                  "shouldn't be! " + self.bank.bankPath))
                    self.success = ASMSUC.none
                    continue
                elif self.bank.sameStrand and s1isRC != s2isRC:
                    logging.info(("We're on different strands, but",
                                  "shouldn't be! " + self.bank.bankPath))
                    self.success = ASMSUC.none
                    continue
                    
                #Wee seem to have a span -- closed gap of some sort
                #start, end = contigLayout.getFillSeqCoords(self.bank.sameStrand)
                
                sequences = self.bank.extractSequences(contigs=[contigLayout.contigName])
                #Get the fastq entry I need
                entry = sequences.values()[0][1]

                readSeed1, readSeed2 = self.mapSeedsToFastqEntry(entry, \
                          contigLayout.readSeed1, contigLayout.readSeed2)
                
                end   = max(readSeed1.tstart, readSeed2.tstart) 
                start = min(readSeed1.tend,   readSeed2.tend)
                
                #Check seed overlaps:
                if readSeed1.tstart < readSeed2.tstart and readSeed1.tend > readSeed2.tstart:
                    logging.debug("Trim")
                elif readSeed2.tstart < readSeed1.tstart and readSeed2.tend > readSeed1.tstart:
                    logging.debug("Trim2")
                #if start > end ... trims... fuck  !!!!!!
                #self.success = ASMSUC.trim Need to get this
                 
                self.metrics["span"] = True
                self.metrics["contigName"] = contigLayout.contigName
                
                self.metrics[self.bank.seed1Name+"_Strand"] = "-" \
                            if contigLayout.readSeed1.isRC else "+"
                self.metrics[self.bank.seed2Name+"_Strand"] = "-" \
                            if contigLayout.readSeed2.isRC else "+"
                
                #Record Sequence, qualities, and it's coordinates
                self.metrics["fillStart"], self.metrics["fillEnd"] = start,end
                self.metrics["fillLength"] = end-start

                self.metrics["fillSeq"] = entry.seq[start:end]
                self.metrics["fillQual"] = entry.qual[start:end]
                
                self.metrics["numContigsProduced"] = len(self.bank.contigs)
                # -2 corrects for seeds in count
                self.metrics["numReadsInContig"] = len(contigLayout.reads) - 2
                self.success = ASMSUC.fill
                return 
                
            elif contigLayout.seedCount == 1:
                #We could have extended either seed 
                if contigLayout.readSeed1 is not None:
                    self.contigSeed1 = contigLayout
                    seedName = self.bank.seed1Name
                    seed = contigLayout.getRead(seedName)
                else:
                    self.contigSeed2 = contigLayout
                    seedName = self.bank.seed2Name
                    seed = contigLayout.getRead(seedName)
                
                self.metrics[seedName+"_Extend"] = True
                self.metrics[seedName+"_ContigName"] = contigLayout.contigName
                
                self.metrics[seedName+"_Strand"] = "-" \
                            if seed.isRC else "+"
                
                sequences = self.bank.extractSequences(contigs=[contigLayout.contigName])
                entry = sequences.values()[0][1]

                alignment = self.mapSeedsToFastqEntry(entry, seed)
                
                if seed.readName.endswith("e3"):
                    downstream = True
                elif seed.readName.endswith("e5"):
                    downstream = False
                else:
                    logging.critical(("Safety Net for Read ending in "
                                    "e[35]! " + seed))
                    exit(7)
                        
                if seed.isRC:
                    downstream = not downstream
                    
                if downstream:
                    start = alignment.tend
                    end = None
                elif not downstream:
                    start = 0
                    end = alignment.tstart
                    

                self.metrics[seed.readName+"_ExtendStart"] = start
                self.metrics[seed.readName+"_ExtendEnd"] = end
                self.metrics[seed.readName+"_ExtendSeq"] = entry.seq[start:end]
                self.metrics[seed.readName+"_ExtendQual"] = entry.qual[start:end]
                self.metrics[seed.readName+"_ExtendLen"] = len(entry.seq[start:end])

                self.metrics[seed.readName+"_NumContigsProduced"] = len(self.bank.contigs)
                self.metrics[seed.readName+"numReadsInContig"] = len(contigLayout.reads) - 1
                
                self.success += ASMSUC.singleExtend

    def mapSeedsToFastqEntry(self, fastqEntry, seed1, seed2=None):
        """
        Maps seeds back to read in order to get the fill seq coordinates
        """
        output = tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', \
                                            delete=False)
        output.write(">target\n%s\n" % (fastqEntry.seq))
        cmd = "blasr %s %s -bestn 1 -m 4 -noSplitSubreads" \
                    % (self.seedReadsFastq, self.bank.outputFasta)
        logging.debug("Running remap \"%s\"" % (cmd))
        r,o,e = exe(cmd)
         
        if r != 0:
            logging.error("Couldn't ReMap")
            logging.error(str(r))
            logging.error(str(o))
            logging.error(str(e))
            exit(10)
        
        ret = []
        for line in o.strip().split("\n"):
            logging.debug(line)
            ret.append(M4Line(line))

        if seed2 == None:#looking for single best
            for i in ret:
                if i.qname == seed1.readName:
                    return i
            logging.warning("Couldn't find %s when remapping" % (seed1.readName))
            exit(10)
        
        return ret

class BankOperator():
    """
    Object for holding information about what reads are used in each contig for an individual level's run
    And what operations (procedures) can be run on these levels.
    """
    def __init__(self, bankOutput, seed1, seed2):
        """
        """
        self.bank = bankOutput
        
        self.__seed1 = seed1
        self.__seed2 = seed2
        
        self.__availableOps = deque(OPERATIONORDER) 

    def updateOperations(self):
        """
        Looks at the bankOutput and figures out what procedures
        are available
        -- Need to have AssemblyAssessor the bank, first
        """
        if len(self.bank.contigs) < 2:
            #really can't do anything
            logging.debug("Bank can't push contigs")
            self.__availableOps.remove("PushContigs")
        
        if self.bank.hasSeed1 or self.bank.seed1Name is None:
            #Because we already have it
            logging.debug("Bank can't pull seed 1")
            #self.__availableOps.remove("PullSeed1")
        else:
            logging.debug("Pulling Seed1")
            self.pullSeed1()
        
        if self.bank.hasSeed2 or self.bank.seed2Name is None:
            #Because we already have it
            logging.debug("Bank can't pull seed 2")
            #self.__availableOps.remove("PullSeed2")
        else:
            logging.debug("Pulling Seed2")
            self.pullSeed2()
        
        hasDegen = False
        for contig in self.bank.contigs:
            if contig.isDegenerate:
                hasDegen = True
        if not hasDegen:
            logging.debug("Bank can't remove degenerate contigs")
            self.__availableOps.remove("RemoveDegenerateContigs")
            
    def clearOperations(self):
        """
        removes all operations
        """
        self.__availableOps = []
        
    def hasOperation(self):
        """
        Is there another operation available to perform?
        """
        return len(self.__availableOps) > 0

    def getAvailableOperations(self):
        """
        Returns a copy of availableOps
        """
        return list(self.__availableOps)
            
    def performNextOperation(self):
        """
        Looks at __availableOps and
            1
        """
        if not self.hasOperation():
            return False, None, None
        operName = self.__availableOps.popleft()
        
        logging.debug("Performing Operation %s on %s" % (operName,\
                                                self.bank.bankPath))
        operExe = OPERATIONS[operName]
        reads, contigs = operExe(self)
        
        return True, reads, contigs
    
    def renameContig(self, reads, contig, doOne): 
        """
        Helper name to proprogate seed names up through levels of the assembly

        reads is a dict of FastqEntries
        contig is the ContigLayout that has the seed we're renaming
        doOne boolean that says to rename it to seed1 if True
        """
        if doOne:
            reads[contig.readSeed1.readName] = reads[contig.contigName]
            reads[contig.readSeed1.readName].name = contig.readSeed1.readName
            del(reads[contig.contigName])
            contig.contigName = contig.readSeed1
        else:
            reads[contig.readSeed2.readName] = reads[contig.contigName]
            reads[contig.readSeed2.readName].name = contig.readSeed2.readName
            del(reads[contig.contigName])
            contig.contigName = contig.readSeed2
    
    #Begin Operations

    def pullSeed1(self):
        """
        find seed1 in inputs, pull it over.
        returns a dictonary of the FastqEntries to move forward with 
        """
        outputReads = mergeFastaQual(self.bank.outputFasta, self.bank.outputQual)
        inputReads = FastqFile(self.bank.inputFastq)
        allContigs = copy.deepcopy(self.bank.contigs)
        seedRead = inputReads[self.__seed1]
        outputReads[self.__seed1] = seedRead
        #Find the contig containing seed2 and rename it
        #if self.__seed2 is not None:#if there is a seed to get
            #for contig in allContigs:
                #if contig.readSeed2 is not None:
                    #if contig.readSeed2.isRC:
                        #outputReads[contig.contigName].reverseCompliment()
                    #self.renameContig(outputReads, contig, False)
                    #break
        #Replace output because we should always have seed2
        self.bank.replaceOutput(outputReads)
        #return outputReads, allContigs
    
    def pullSeed2(self):
        """
        find seed2 in inputs, pull it over
        returns a dictonary of the FastqEntries to move forward with 
        """
        outputReads = mergeFastaQual(self.bank.outputFasta, self.bank.outputQual)
        inputReads = FastqFile(self.bank.inputFastq)
        allContigs = copy.deepcopy(self.bank.contigs)
        seedRead = inputReads[self.__seed2]
        outputReads[self.__seed2] = seedRead
        #Find the contig containing seed1 and rename it
        #if self.__seed1 is not None:#if there is a seed to get
            #for contig in allContigs:
                #if contig.readSeed1 is not None:
                    #if contig.readSeed1.isRC:
                        #outputReads[contig.contigName].reverseCompliment()
                    #self.renameContig(outputReads, contig, True)
                    #break
        #Replace output because we should always have seed2
        self.bank.replaceOutput(outputReads)
        #return outputReads, allContigs

    def pushContigs(self):
        """
        Don't Touch Anything, Just push this forward
        returns a dictonary of the FastqEntries to move forward with 
        """
        outputReads = mergeFastaQual(self.bank.outputFasta, self.bank.outputQual)
        allContigs = copy.deepcopy(self.bank.contigs)
        for contig in allContigs:
            if not contig.isDegenerate:
                if contig.readSeed1 is not None:
                    if contig.readSeed1.isRC:
                        outputReads[contig.contigName].reverseCompliment()
                    self.renameContig(outputReads, contig, True)
                #both won't happen because?? I hope
                if contig.readSeed2 is not None:
                    if contig.readSeed2.isRC:
                        outputReads[contig.contigName].reverseCompliment()
                    self.renameContig(outputReads, contig, False)
        return outputReads, allContigs
           
    def removeDegenerateContigs(self):
        """
        removes all contigs in this level's output that
        do not contain the seed
        """
        logging.info("Removing Degenerate Contigs")
        outputReads = mergeFastaQual(self.bank.outputFasta, self.bank.outputQual)
        allContigs = copy.deepcopy(self.bank.contigs)
        for contig in allContigs:
            if contig.isDegenerate:
                logging.debug("%s is degenerate" % contig.contigName)
                del(outputReads[contig.contigName])
                continue
            if contig.readSeed1 is not None:
                if contig.readSeed1.isRC:
                    outputReads[contig.contigName].reverseCompliment()
                self.renameContig(outputReads, contig, True)
            #both won't happen because?? I hope
            if contig.readSeed2 is not None:
                if contig.readSeed2.isRC:
                    outputReads[contig.contigName].reverseCompliment()
                self.renameContig(outputReads, contig, False)
        return outputReads, allContigs
   

OPERATIONS = {"PushContigs"   : BankOperator.pushContigs ,\
              #"PullSeed1"  : BankOperator.pullSeed1 ,\
              #"PullSeed2" : BankOperator.pullSeed2 ,\
              "RemoveDegenerateContigs":BankOperator.removeDegenerateContigs}
            
OPERATIONORDER = ["PushContigs", "RemoveDegenerateContigs"]

class AssemblyIteration():
    
    def __init__(self, inputFastq, seedNames, translate=None, nproc=1, workDir=None, debug=False):
        """
        inputFastq is the input Fastq File
            or a list of [input.fasta, input.qual]
        Make my temporary directory
        Update OLCAssembly.py to work on input.fastq
        translate = a previous AssemblyIteration that can tell us how to 
            update this iteration's seed coordinates
        """
        self.inputFastq = inputFastq
        self.translate = translate
        self.nproc = nproc
        self.debug = debug
        if workDir is None:
            self.bankDir = tempfile.mkdtemp()
        else:
            self.bankDir = workDir
        if not os.path.exists(self.bankDir):
            raise OSError("Assembly Work Directory %s Does Not Exist" \
                            % self.bankDir)
        self.runAssembly(seedNames)
            
    def runAssembly(self, seedNames, threshold=800):
        """
        Runs OLCAssembly
        Makes self.BankOutput
        Makes self.BankOperations
        """
        if type(self.inputFastq) == list:
            fasta, qual = self.inputFastq
        elif type(self.inputFastq) == str:
            fasta = self.inputFastq
            qual = ""
        
        logging.info("Running OLCAssembly.py")
        dbg = "--debug" if self.debug else ""
        #cut the error rate in half if we're working with contigs and not raw reads
        errMod = "-e 0.07" if self.translate is not None else "-e 0.15"
        c, o, e = exe("OLCAssembly.py %s %s --nproc %s --workDir %s --threshold %d %s %s" \
                        % (fasta, qual, self.nproc, self.bankDir, threshold, errMod, dbg))
            
        if c != 0:
            logging.error("OLCAssembly Failed!")
            #logging.debug("RETCODE - %d\nSTDOUT - %s\nSTDERR - %s" % \
                            #(c, o, e))
            self.madeAssembly = False
        else:
            logging.info("Finished OLCAssembly.py")
            #logging.debug("RETCODE - %d\nSTDOUT - %s\nSTDERR - %s" % \
                        #(c, o, e))
            self.madeAssembly = True
        
        self.bank = BankOutput(self.bankDir, seedNames)
        self.operations = BankOperator(self.bank, self.bank.seed1Name, self.bank.seed2Name)
        #translate if necessary
        if self.translate is not None:
            self.translateSeedCoords(self.translate)

    
    def translateSeedCoords(self, translate):
        """
        Find the contig that has seed1 or seed2
        """
        logging.debug("Attempting To Translate Seed Coordinates Between Assemblies")
        if self.bank.hasSeed1 and translate.bank.hasSeed1:
            logging.debug("contig1 len orig %d new %d" % (self.bank.contigSeed1.contigLength,\
                            translate.bank.contigSeed1.contigLength))
            self.bank.contigSeed1.readSeed1.translate(translate.bank.contigSeed1.readSeed1)
        if self.bank.hasSeed2 and translate.bank.hasSeed2:
            logging.debug("contig2 len orig %d new %d" % (self.bank.contigSeed2.contigLength,\
                            translate.bank.contigSeed2.contigLength))
            self.bank.contigSeed2.readSeed2.translate(translate.bank.contigSeed2.readSeed2)


    def assessAssembly(self, readsFile):
        """
        Before hand I need to translate the coordinates?
        Runs and Makes self.AssemblyAssessor (metrics)
        readsFile holds the seeds so that I can map them to the output
        """
        if not self.madeAssembly:
            self.assessor = None
            self.operations.clearOperations()
        else:
            self.assessor = AssemblyAssessor(self.bank, readsFile)
            self.operations.updateOperations()
            if self.assessor.success == ASMSUC.none:
                self.operations.clearOperations()
    
    def hasNextOperation(self):
        """
        Checks to see if there is a next operation to try
        """
        return self.operations.hasOperation()

    def getNextOperation(self):
        """
        Makes self.operations perform the next operation,
        then it saves the resulting reads in a temporary file,
        then it returns the filename of the resulting reads we've created.
        """
        success, reads, contigs = self.operations.performNextOperation()
        if not success:
            return None
        
        output = tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', \
                                                            delete=False)
        logging.debug("Made File %s" % (output.name))
        for read in reads:
            output.write(reads[read].toString())
        output.close()
        return output.name
    
class Assembly():
    
    def __init__(self, args):
        self.parseArgs(args)
        setupLogging(self.options.debug)

    def parseArgs(self, argv):
        parser = OptionParser(USAGE)
        parser.add_option("--debug", action="store_true", \
                          help="Increases verbosity of logging" )
        parser.add_option("--nproc",type=int, default=1, \
                          help="Number of processors available for use [DEFAULT=1]")
        parser.add_option("--workDir", type=str, default=None,\
                          help="Working Directory. [DEFAULT=temporary folder]")
        parser.add_option("-p", "--predictedGapSize", type=int, default=None)
        
        self.options, args = parser.parse_args(argv)
        
        self.predictedGapSize = self.options.predictedGapSize
        
        if len(args) == 2:
            self.inputDir = args[1]
        else:
            parser.error("Incorrect number of arguments!")
            exit(1)
        
    def outputMetrics(self, assemblyIteration):
        """
        Save the metrics from a particular iteration
        """
        output = os.path.join(self.inputDir,"fillingMetrics.json")
        logging.info("Saving Filling Metrics to %s" % output)
        fout = open(output,'w')
        if self.predictedGapSize is not None:
            assemblyIteration.assessor.metrics["predictedGapSize"] = self.predictedGapSize
        json.dump(assemblyIteration.assessor.metrics, fout, indent=4)
        fout.close()
        
    def successfulAssembly(self, asIt):
        """
        Is this assembly as good as it can get?
        """
        if self.bestAssembly == ASMSUC.fill:
            if asIt.assessor.success == ASMSUC.fill:#or trim?
                if asIt.assessor.metrics["fillLength"] < 1:
                    return False
                return True
        if self.bestAssembly == ASMSUC.singleExtend:
            if asIt.assessor.success == ASMSUC.singleExtend:
                return True
        return False
    
    def stackTraceForNextOperation(self):
        """
        Gets the next non-None operation from the stack by popping
        if necessary
        """
        nextOperation = None
        #if len(self.assemblyStack) == 2|3: stack.pop() (don't want to get too deep)
        while nextOperation is None and len(self.assemblyStack) > 0:
            #nextOperation is a string to the next readInput
            nextOperation = self.assemblyStack[-1].getNextOperation()
            if nextOperation is None:
                logging.debug(("Assembly Iteration %d ran out of "
                               "Operations... popping") % \
                               ( len(self.assemblyStack) ))
                self.assemblyStack.pop() #<-- Done with this guy
                
        return nextOperation

    def chooseBestMetrics(self):
        """
        Tries to determine the best metrics we've found
        """
        #Obvi
        if len(self.assemblyCandidates) == 1:
            return self.assemblyCandidates[0]
        else:
            logging.warning("notimp")
            return self.assemblyCandidates[0]
        best = self.assemblyCandidates[0]
        for challenger in self.assemblyCandidates[1:]:
            #Predicted gap stuff would be good here
            #if challenger.has_key
            pass


    def run(self):
        """
        """
        #Stack of assemblies
        self.assemblyStack = []
        #A list of all the assemblies we've made
        self.assemblyCandidates = []
        readInput = os.path.join(self.inputDir, "input.fastq")
        
        self.seeds = [x.replace('.','/') for x in \
            os.path.abspath(self.inputDir).split('/')[-1].split('_')]
            
        self.bestAssembly = ASMSUC.fill
        if len(self.seeds) == 1:
            self.bestAssembly = ASMSUC.singleExtend
            self.seeds.append(None)
        
        #Make A temporary file for mapping seeds back to query
        fastq = FastqFile(readInput)
        tfile = tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', \
                                                    delete=False)
        SEEDSFASTQ = tfile.name
        
        logging.debug("Writing Seeds For Remap to %s" % SEEDSFASTQ)
        for i in self.seeds:
            if i is None:
                continue
            tfile.write(fastq[i].toString())
        tfile.close()
        
        logging.debug("Seeds! " + str(self.seeds))
        
        while True:
            prev = None if len(self.assemblyStack) == 0 else self.assemblyStack[-1]
            asIt = AssemblyIteration(readInput, self.seeds, translate=prev, \
                    nproc=self.options.nproc, workDir=self.options.workDir, \
                    debug=self.options.debug)
            
            logging.info("Assessing Assembly")
            asIt.assessAssembly(SEEDSFASTQ)
            logging.debug("Available Operations")
            logging.debug(asIt.operations.getAvailableOperations())
            
            #Did we even get an assembly out
            if not asIt.madeAssembly and not asIt.hasNextOperation():
                #No hope for this Assembly Path anyways,
                #Lower the threshold and see what we can do
                logging.info("Lowering Threshold")
                asIt.runAssembly(self.seeds, threshold=400 )
                asIt.assessAssembly(SEEDSFASTQ)
                
            
            #Still?
            if not asIt.madeAssembly:
                logging.debug(("Couldn't make an assembly... " \
                                "popping up the assembly stack"))
                if len(self.assemblyStack) < 1:
                    logging.info("Failed To Make A Single Assembly")
                    exit(1)
                
                nextOperation = self.stackTraceForNextOperation()
            else:
                self.assemblyStack.append(asIt)
                if self.successfulAssembly(asIt):
                    logging.info("Found Acceptable Assembly!")
                    self.outputMetrics(asIt)
                    break
                else:
                    #Let's try moving forward
                    logging.info("Assembly Success Code -> %d" % (asIt.assessor.success))
                    self.assemblyCandidates.append(asIt)
                    nextOperation = self.stackTraceForNextOperation()
            
            if nextOperation is None:
                logging.debug(("No Next Operation Found. Picking from "
                                "%d Assembly Candidates Found" % (len(self.assemblyCandidates))))
                
                best = self.chooseBestMetrics()
                self.outputMetrics(best)
                break
            else:
                readInput = nextOperation
                    
            
        
if __name__ == '__main__':
    runner = Assembly(sys.argv)
    runner.run()
    #print json.dumps(runner.assemblyStack[-1].assessor.metrics, indent=4)
"""
BankOutput(String bankPath)
    Holds Bank's Information
    
BankOperations(Static Init)
    Holds Operations that can be performed
 
AssemblyAssessor() 
    Similar to the GapSupporter, it compares what's expected (need ideas on how to make this) with what has been made (SingleAssembly)
   
AssemblyIteration(holds BankOutput, BankOperations, AssemblyAssessor)
    Holds A single assembly (a bank and a set of operations that could be / have been performed)

Assembly()
    Main class for running
    
    #predictedGapSize = sys.argv[2]
    while there is work to do:
        Make a temporary Directory?
        1) Build an AssemblyIteration
            --It makes it's own BankOperations
            --It makes it's own BankOutput
        
        2) If not first AssemblyIteration:
            Shift the SeedContigs

        3) Assess BankOutput 
        4) Decide what needs to be done Next -
            a) Success! Finished. (break)
            b) More to be done..
                Look at the current AssemblyIteraion's available BankOperations
                Make a new level based on the Operation
                
                NOTE! To keep the s1isRC s2isRC assumption safe,
                    through multiple iterations of OLCAssembly
                    I'll need to ReverseCompliment any RC seed
        
        Here is a little bit of what step 4 should look like
        Part of this work is going to be messing with BankOperations of 
        the current AssemblyIteration
        
        assessment = AssemblyAccessor(...)
        if assessment.success == ASMSUC.fill:
            finished
        if assessment.success == ASMSuc.doubleExtend:
            Check the extend lengths,
            see if we could potentially span with 
            a little overlap between the two
            or if there is no predicted gap size, try again anyway
            else break
        if assessment.success == ASMSUC.singleExtend:
            is this all that we need? If so finished
            if not, continue forward
        if assement.success == ASMSUC.none:
            I need to keep trying, manipulating the data

    Note: If I ever move backwards, I'm going to pop the stack into a list
        of best found candidates.
        If I get to the end, and haven't found a fantastic assembly,
        I'll just compare all of the assessments to get the best one.
        See the previous code for examples of 'best'
            -- Look into that (\d+) on the bank2contig -L.
            -- If it is perfectly preserved, is it 1000k?
"""
