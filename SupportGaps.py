#!/usr/bin/env python
import re, os, sys, logging, json
from bisect import bisect_left, bisect_right
from FileHandlers import *
from optparse import OptionParser
from collections import defaultdict
from Setup import refParser

USAGE = """SupportGaps.py <alignmentFile> <gapInfo> <outFile> [--spanOnly]

Parses an .m4 or .m5 alignment file and determines which gaps in gapInfo are
supported by reads. Results reported in outFile"""

"""
TODO

This Should be heavily refactored.
Consider putting start/ends of alignments in gapCanInfo so one can trim reads to what
supports the gap.

"""
#Max Distance a read's alignment can stop from a gap and still be considered 
#for support of the gap
MAXFLANK = 50
#Min number of bases a read needs to reach into a gap to be considered for 
#support of the gap
MINCOVERS = 25

#I need three classes here. SupportFlags is fine.
#I need to split Support Gaps into what is necessary for parsing an entire m4 file
#and something that handles a single group of alignments (for reuse in filling)

class SupportFlags(dict):
    """
    Keeps track of all the possible Flags for Reads 
    """  #Read's  Number Of Hits
    flags = ["UniqueMapping",
         "MultiMapping",
        #Read's Relation to any gap
        "LeftContig", 
        "RightContig", 
        "SpansGap", #Read Spans an entire gap
        #Read's ranking within all of it's hits
        "LQPA",
        "MostAccurate",
        "BestScore"]

    def __init__(self):
        super(dict)
        self.update({"UniqueMapping": 1,
                     "MultiMapping":  2,
                     "LeftContig":    4,
                     "RightContig":   8,
                     "SpansGap":      16,
                     "LQPA":          32,
                     "MostAccurate":  64,
                     "BestScore":     128})

    def translateFlag(self, flag, text=True):
        """
        Takes a flag and turns it into a string descriptor
        if not text, return an array instead. This will be useful for 
        populating Table Columns
        """
        ret = []#String
        for name in self.flags:
            if flag & self[name]:
                ret.append(name)
        if text:
            return " ".join(ret)
        return ret

class SupportClassifier():
    
    def __init__(self,  gapInfo, gapIndex = None):
        self.flagTable = SupportFlags()
        self.gapInfo = gapInfo
        self.gapIndex = gapIndex #self.gapInfo.getSortedGaps()


    def isConcordant(self, A, B):
        """
        Takes two reads (A,B) and checks their mapping concordancy
        """
        
        oscafName, aScaf, aContig = refParser.match(A.tname).groups()
        oscafName, bScaf, bContig = refParser.match(B.tname).groups()
        
        #We're looking at scaffolding
        if aContig == None and bContig == None: 
            if aScaf != bScaf:#Must align to the same scaffolding
                return 0
            #A check for negative gaps should be in here
            #Just see if the pieces are in order
            logging.debug("As, Ae, Ae, Bs, Bs, Be") 
            logging.debug(",".join(map(str,[A.qstart, A.qend, A.qend, \
                                            B.qstart, B.qstart, B.qend])))
            if A.tstrand != B.tstrand:
                return 0
            if A.tend < B.tstart and A.qend < B.qstart:
                return 1
            if B.tend < A.tstart and B.qend < A.qstart:
                return -1
            #if (A.qstart < A.qend and A.qend < B.qstart and B.qstart < B.qend):
                ##Also want to see if it's on the target correctly
                #if A.tstart < A.tend and A.tend < B.tstart and B.tstart < B.tend:
                    #return True
            return 0
            
        #Let's do the contig Work -- shouldn't happen anymore
        aContig, bContig = int(aContig), int(bContig)
        
        if aScaf != bScaf:#Put in across scaffolding support later
            return False
        
        if abs(aContig - bContig) == 1:#Adjacent
            if A.qstrand == B.qstrand and A.qstart < A.qend and A.qend < B.qstart:
                return True
        return False

    def mappingType(self, reads):
        """
        Sets Reads flag based on what type of Alignment Count it is.
        "UniqueMapping", "MultiMapping"
        """
        logging.debug("Determining Read's mapping type from hits")
        fInc = 0
        if len(reads) == 1:
            fInc = self.flagTable["UniqueMapping"]
        else:
            sorted = []
            for r in reads:
                oscafName, scaf, contig,  = refParser.match(r.tname).groups()
                if contig == None:
                    contig = 0
                sorted.append((int(scaf[3:]), int(contig), r.tstart, r.tend, \
                               r.qstart, r.qend, r))
            
            sorted = map(lambda x: x[-1], sorted)
            for i in range(0, len(sorted)-1):
                if self.isConcordant(sorted[i], sorted[i+1]):
                    fInc = self.flagTable["UniqueMapping"]
                else:
                    fInc = self.flagTable["MultiMapping"]
                    break
        
        for i in reads:
            i.flag += fInc
    
    def groupComparison(self, reads):
        """
        Find the best LQPA, ACC, and HISCORE
        """
        #The number behind the best; the list of bests (Cause there can be ties)
        logging.debug("Flagging hits based on group")
        LQPA = reads[0].queryPctAligned; lqpas = [reads[0]]
        ACC = reads[0].pctsimilarity; accs = [reads[0]]
        SCORE = reads[0].score; scores = [reads[0]]
        #Find the bests
        for read in reads[1:]:
            if read.queryPctAligned == LQPA:
                lqpas.append(read)
            elif read.queryPctAligned > LQPA:
                LQPA = read.queryPctAligned
                lqpas = [read]
                
            if read.pctsimilarity == ACC:
                accs.append(read)
            elif read.pctsimilarity > ACC:
                ACC = read.pctsimilarity
                accs = [read]
            
            if read.score == SCORE:
                scores.append(read)
            elif read.score < SCORE:
                SCORE = read.score
                scores = [read]
        #Set the apropriate flags
        for read in lqpas:
            read.flag += self.flagTable["LQPA"]
        for read in accs:
            read.flag += self.flagTable["MostAccurate"]
        for read in scores:
            read.flag += self.flagTable["BestScore"]

    def findSupport(self, read, scaffold, contig, spanOnly):
        """
        Find a read's support on the scaffolding
        """
        
        logging.debug("Finding read's support")
        if contig != None:
            return self.contigSupport(read, scaffold, contig)
        else:
            return self.scaffoldSupport(read, scaffold, spanOnly)
    
    def contigSupport(self, read, scaffold, contig):
        """
        Support for Contig Mapping
        """
        ret = defaultdict(GapCans)
        
        gapFromRight = Gap("ignore", -sys.maxint, 0, \
                           scaffold+"_"+str(int(contig)-1)+"_"+contig )
        gapFromLeft = Gap("ignore", read.tseqlength, sys.maxint, \
                          scaffold+"_"+contig+"_"+str(int(contig)+1) )
        
        if self.covers(read, gapFromLeft) != None:
            ret[gapFromLeft.name]["LeftContig"].append(read)
            read.flag += self.flagTable["LeftContig"]
        if self.covers(read, gapFromRight) != None:
            ret[gapFromRight.name]["RightContig"].append(read)
            read.flag += self.flagTable["RightContig"]
        
        return ret
        
    def scaffoldSupport(self, read, scaffold, spanOnly):
        """
        Support for Scaffold Mapping
        """
        if self.gapIndex != None:
            try:
                #Find range of gaps we can potentially support 
                #I play it safe and get upto 2 extra gaps we could support
                startIndex = max(0,bisect_left(self.gapIndex[scaffold], read.tstart) - 1)
                endIndex = bisect_right(self.gapIndex[scaffold], read.tend)+1
            except KeyError:
                return {}
            parse = self.gapIndex[scaffold][startIndex:endIndex]
        else:
            gaps = filter(lambda x: x.startswith(scaffold), self.gapInfo.keys())
            parse = []
            for key in gaps:
                parse.append(self.gapInfo[key])
        
        ret = defaultdict(GapCans)
        for gap in parse:#self.gapIndex[scaffold][startIndex:endIndex]: 
            #gap = self.gapInfo[gap.name]
            if self.spans(read, gap):
                logging.debug("Supports %s" % (str(gap)))
                logging.debug("SpansGap - gapLength %d" % (gap.length))
                ret[gap.name]["SpansGap"].append(read)
                read.flag += self.flagTable["SpansGap"]
            
            if not spanOnly:
                supportType = self.covers(read,gap)
                if supportType != None:
                    ret[gap.name][supportType].append(read)
                    read.flag +=  self.flagTable["SpansGap"]
        
        return dict(ret)
        
    def spans(self, read, gap):
        if read.tstart < gap.start and read.tend > gap.end:
            logging.debug("Supports %s" % (str(gap)))
            return True
        
        return False
    
    def covers(self, read, gap):
        #moving into gap
        if read.tstrand == "0":
            distanceFromContigEnd = gap.start - read.tend
            remainingReadSeq = read.qseqlength - read.qend - MINCOVERS
            if distanceFromContigEnd >= 0 and \
               distanceFromContigEnd < remainingReadSeq and \
               distanceFromContigEnd <= MAXFLANK:
                logging.debug("Supports %s" % (str(gap)))
                logging.debug("LeftContig - distanceFromEnd %d - remainingSeq %d " % \
                              (distanceFromContigEnd, remainingReadSeq))
                #Positive Strand Maps on Left Contig and enters gap     
                return "LeftContig"

            #moving outof gap
            distanceFromContigBeginning = read.tstart - gap.end 
            remainingReadSeq = read.qstart - MINCOVERS 
            
            if distanceFromContigBeginning >= 0 and \
               distanceFromContigBeginning < remainingReadSeq and \
               distanceFromContigBeginning <= MAXFLANK:
                #Positive Strand Maps on Right Contig and Exits Gap
                logging.debug("Supports %s" % (str(gap)))
                logging.debug("RightContig - distanceFromBeginning %d - remainingSeq %d " % \
                              (distanceFromContigBeginning, remainingReadSeq))
                return "RightContig"
        
        elif read.tstrand == "1":
            distanceFromContigEnd = gap.start - read.tend
            remainingReadSeq = read.qstart - MINCOVERS
            if distanceFromContigEnd >= 0 and \
               distanceFromContigEnd < remainingReadSeq and \
               distanceFromContigEnd <= MAXFLANK:
                logging.debug("Supports %s" % (str(gap)))
                logging.debug("LeftContig - distanceFromBeginning %d - remainingSeq %d " % \
                              (distanceFromContigEnd, remainingReadSeq))
                return "LeftContig"
            
            distanceFromContigBeginning = read.tstart - gap.end
            remainingReadSeq = read.qseqlength - read.qend - MINCOVERS
            if distanceFromContigBeginning >= 0 and \
               distanceFromContigBeginning < remainingReadSeq and \
               distanceFromContigBeginning <= MAXFLANK:
                logging.debug("Supports %s" % (str(gap)))
                logging.debug("RightContig - distanceFromBeginning %d - remainingSeq %d " % \
                              (distanceFromContigBeginning, remainingReadSeq))
                return "RightContig"
        
    
    def untangle(self, reads):
        """
        Given a group of subread's hits, see if we can eliminate the 
        spurious repeat matches, and create
        
        Find LQPA 
        recursively call(give me all the concordant neighbors to the left and to the right)
        continue until there isn't exactly one neighbor
        """
        if len(reads) == 1:
            return reads
        
        bestScore = []
        mostAccurate = []
        for read in reads:
            if read.flag & self.flagTable["BestScore"]:
                bestScore.append(read)
            if read.flag & self.flagTable["MostAccurate"]:
                mostAccurate.append(read)
        if len(bestScore) == 0:
            return []
        
        anchor = None
        if len(bestScore) != 1:#Need to find the next best
            if len(mostAccurate) == 1:
                for i in bestScore:
                    if i == mostAccurate[0]:
                        anchor = mostAccurate[0]
        #Just arbiturarily take one
        if anchor == None:
            anchor = bestScore[0]
        
        #now try to unchain to the right and left
        newReads = [anchor]
        newReads.extend(self.layout(anchor, reads, side=1))
        newReads.extend(self.layout(anchor, reads, side=-1))
        
        return newReads
        
    def layout(self, anchor, reads, side):
        """
        uses the anchor to find all unique concordant hits from the anchor 
        side = 1, check right
        side = -1, check left
        """
        reads = filter(lambda x: x != anchor, reads)
        neighbors = []
        for read in reads:
            if self.isConcordant(anchor, read) == side:
                neighbors.append(read)
        
        if len(neighbors) != 1:
            return []

        neighbors.extend(self.layout(neighbors[0], reads, side))

        return neighbors
    
    def removeDiscordantHits(self, readGroup):
        """
        Given a set of reads. Take out all of the reads that 
        have long unmapped tails that go over the reference
        
        DEPRICATED!
        """
        ret = []
        for read in readGroup:
            tailLenLeft = read.qstart
            tailLenRight = read.qseqlength - read.qend
            
            danger = False
            #Can this read save itself?
            #P.s. This is the similar code as entersExitsGap.. 
            # I'm sorry for my terrible programming. I'll clean
            #This next version.
            
            if tailLenLeft > 50:#you're in trouble
                if read.qstrand == '0' and read.tstart >= MAXFLANK:
                        danger = True
                elif read.qstrand =='1' and read.tseqlength-read.tend >= MAXFLANK:
                        danger = True
            if danger:#YOu're outta here, read!
                continue
            #Final chance for the read to prove itself
            if tailLenRight > 50:
                if read.qstrand == '0' and read.tseqlength - read.tend >= MAXFLANK:
                    danger = True
                elif read.qstrand == '1' and read.tstart >= MAXFLANK:
                    danger = True
            if danger:
                continue
            
            ret.append(read)
        return ret
    
    def classifyReads(self, reads, spanOnly=False):
        """
        Standard procedure to create how a set of reads supports a gap (if any)
        """
        logging.debug("Classifying %s" % (reads[0].qname))
        mySupport = defaultdict(GapCans)
        self.mappingType(reads)
        
        #self.reads[readGroup] = self.removeDiscordantHits(self.reads[readGroup])
        #if len(self.reads[readGroup]) == 0:
        #   continue
        
        self.groupComparison(reads)
        
        #Do "Mapping Untangling Here" also, filter out multimappers
        if reads[0].flag & self.flagTable["MultiMapping"]:
            reads = self.untangle(reads)
        
        for read in reads:
            oscafName, scaffold, contig = refParser.match(read.tname).groups()
            supportedGaps = self.findSupport(read, scaffold, contig, spanOnly)
            for i in supportedGaps.keys():
                mySupport[i].extend(supportedGaps[i])
            
        for gap in mySupport:
            mySupport[gap].consolidate()
        
        return dict(mySupport)

       
class SupportGaps():
    """
    Object for taking reads and classifying how map in conjunction with gaps
    """
    
    def __init__(self):
        self.parseArgs()
        self.__initLog()
    
    def __initLog(self):
        """Logging"""
        logLevel = logging.DEBUG if self.options.debug else logging.INFO
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
        logging.info("Running %s" % " ".join(sys.argv))
    
    def parseArgs(self):
        parser = OptionParser(USAGE)

        parser.add_option("--debug", action="store_true", \
                          help="Increases verbosity of logging" )
        parser.add_option("--spanOnly", action="store_true", \
                          help=("When using Scaffolding for mapping, span will"
                                " only allow support by reads that span an"
                                " entire gap. i.e. no contig extension."))
        self.options, args = parser.parse_args()
        
        if len(args) != 3:
            parser.error("Error! Incorrect number of arguments")
        
        if not os.path.isfile(args[0]):
            parser.error("Error! Alignment File Does Not Exist")
        self.alignmentFileName = args[0]
        
        if not os.path.isfile(args[1]):
            parser.error("Error! Gap Info File Does Not Exist") 
        self.gapFileName = args[1]
        
        if os.path.isfile(args[2]):
            parser.warning("Output File Being Overwritten!")
        self.outputFileName = args[2]
        
        self.gapInfo = GapInfoFile(self.gapFileName)
        if os.path.splitext(self.alignmentFileName)[1] == '.m4':
            #Temporarily set this as M5 since 
            self.alignmentFile = M4File(self.alignmentFileName)
        elif os.path.splitext(self.alignmentFileName)[1] == '.m5':
            self.alignmentFile = M5File(self.alignmentFileName)
        else:
            parser.error("Error! Alignment File Extension (%s) not recognized." \
                         % os.path.splitext(self.alignmentFileName)[1])
        
    
    def groupReadHits(self):
        logging.debug("Grouping Read Hits")
        reads = defaultdict(list)#readname: [hit hit hit]
        
        for line in self.alignmentFile:
            reads[line.qname].append(line)
        
        self.idAdapters(reads)
        
        return reads    
    
    def idAdapters(self, reads):
        """
        Using the alignments,
        Identify adapters and split those subreads into two alignments
        Consider a subread that overlaps with itself as a missed adapter, split it.
            example - for a single read
            hit1:
                qstart 0: qend 100: tstart = 400: tend=500: strand = 0
            hit2:
                qstart = 150: qend = 250: tstart = 400: tend = 500: strand = 1
        this is highly indicitive of a missed adapter

        right now I will just focus on one adapter being missed...
        new reads are added in place
        """
        update = defaultdict(list)
        for key in reads.keys():
            split = False
            for r1 in reads[key][1:]:
                for r2 in reads[key][:-1]:
                    if (r1.tname == r2.tname and r1.tstrand != r2.tstrand) \
                        and (abs(r1.tstart - r2.tstart) <= 50 or abs(r1.tend - r2.tend) <= 50): 
                        split = True
                        logging.warning("Missed adapter in read! %s" % r1.qname)
                        if (r1.qend < r2.qstart):
                            r1.trim = (0, r1.qend)
                            r1.qseqlength = r1.qend
                            
                            r2.trim = (r2.qstart, r2.qseqlength)
                            shift = r2.qend - r2.qstart
                            r2.qstart = 0
                            r2.qend -= shift
                            r2.qseqlength -= shift
                        else:
                            r2.trim = (0, r2.qend)
                            r2.qseqlength = r2.qend
                            
                            r1.trim = (r1.qstart, r1.qseqlength)
                            shift = r1.qend - r1.qstart
                            r1.qstart = 0
                            r1.qend -= shift
                            r1.qseqlength -= shift

            
            if split:#If there is any split, we're moving all of them apart by strandedness
                for r in reads[key]:
                    update[r.qname+'.'+r.tstrand].append(r)
                del(reads[key])
        reads.update(update)
                        
     
    def run(self):
        """
        Given a group of reads, put it through the paces, then assign it's flag
        """
        classifier = SupportClassifier(self.gapInfo, self.gapInfo.getSortedGaps())
        gapSupport = {}
        readGroups = self.groupReadHits()
        for group in readGroups:
            support = classifier.classifyReads(readGroups[group], \
                                               self.options.spanOnly)
            for key in support:
                if gapSupport.has_key(key):
                    gapSupport[key].extend(support[key])
                else:
                    gapSupport[key] = support[key]
        
                
        fout = open(self.outputFileName,'w')
        logging.debug("Dumping gap cans")
        fout.write(json.dumps(gapSupport, default=GapCansDecode, indent=4))
        logging.info("Finished Support!")


if __name__ == '__main__':
    runner = SupportGaps()
    runner.run()
