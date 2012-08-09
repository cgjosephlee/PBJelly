"""
TODO
remove loading the reference and making real_tseqlength - 
once I have extraction create the delta file, I'll need 
to load in the real_tseqlength

Need to add '--quick' option to WrapAssembly that gets
passed here

Complex fills aren't working at all.
I need to ensure that the gapCans returned
by resupport can have both hits of a read's
alignment reported when necessary
"""
import sys
import os
import logging
import re
import json
from CommandRunner import exe
from FileHandlers import *
from SupportGaps import SupportClassifier
from collections import defaultdict


MAX_TRIM = 150 #The most number of bases allowed to be trim from a contig

class AssemblyAssessor():
    """
    Assess an assembly and stores the best
    assembly it has seen in 
    pwd/[output.fasta, output.qual, remapping.m5, fillingMetrics.json]
    """
    def __init__(self, gapFn, gapInfo, nproc, supportType):
        """
        reference - the reference to map to
        gapFn - the GapInfo.bed holding all the gaps
        gapInfo - the single gap we're trying to assemble around
        nproc - number of processors to use
        supportType - the best possible result for this assembly
                Good metric for judging how much we have 
                assembled (like if we get only left or 
                right support)
        leftSeed rightSeed - The name of the reads that seed this assembly
        """
        #self.reference = reference
        self.gapFn = gapFn
        self.gapInfo = gapInfo

        self.myGap = self.gapInfo.values()[0]
        
        self.nproc = nproc
        self.supportType = supportType
        self.bestMetrics = None
        self.isSuccessful = False
        self.mySupporter = ReSupport(self.gapInfo, self.supportType)
    
    def assess(self, levelObj):
        """
        Returns True if we have an acceptable assembly.
        False if more work should try to be done.
        """
        
        logging.info("Entered Assessment")
        levelObj.beenAssessed = True
        #Short Circuiting. at least one of these needs to be true
        if self.shortCircuiting(levelObj):
            return False
        
        #gapCan = self.reMapAndSupport(levelObj)
        gapCan = self.quick_reMapAndSupport(levelObj)
        if gapCan.isEmpty():
            logging.info("Found No Resupport")
            return False
        #Need to remove Non-Seed reads
        logging.info("Found Gap Support")
        
        logging.debug("%s" % (json.dumps(gapCan, \
                        default=GapCansDecode, \
                        indent=4)))
        
        curMetrics = SupportMetrics(levelObj, self.myGap, \
                    gapCan, self.supportType)
        
        if not curMetrics["SpansGap"] \
                and not curMetrics["LeftContig"] \
                and not curMetrics["RightContig"]:
            
            logging.error(("Bad Support! We didn't actually " \
                           "get anything out. This happens" \
                           " when seeds are separate, but " \
                           "still 'span gap' OR the assembled "\
                           "contig trimmed too many bases."))
            return False
        
        if self.bestMetrics == None or self.bestMetrics.compare(curMetrics):
            logging.info("Found Better Assembly!")
            self.bestMetrics = curMetrics
            self.storeSuccess(levelObj)
        else:
            logging.info("Keeping Old Assembly")
        
        if self.supportType == "SpansGap" and curMetrics["SpansGap"]:
            return True
        if self.supportType == "LeftContig" and curMetrics["LeftContig"]:
            return True
        if self.supportType == "RightContig" and curMetrics["RightContig"]:
            return True
        
        return False
        
    def shortCircuiting(self, levelObj):
        """
        Try to find OBVIOUS shortcomins in the current 
        levelObj without a full (and quite expensive) 
        ReMap and ReSupport.
        """
        #We're not even looking at a good candidate!    
        if levelObj.leftSeedContig == None and levelObj.rightSeedContig == None:
            logging.debug("ShortCircuiting - No Seeds")
            return True
        #Being optimistic and saying we CAN get a span! Woo!
        if self.supportType == "SpansGap" and (levelObj.leftSeedContig==None \
                or levelObj.rightSeedContig==None):
            logging.debug("ShortCircuiting - Won't Span when could")
            return True 
        
        if self.supportType == "LeftContig" and levelObj.leftSeedContig == None:
            logging.debug("ShortCircuiting - Won't Support Left Correctly")
            return True
        
        if self.supportType == "RightContig" and levelObj.rightSeedContig == None:
            logging.debug("ShortCircuiting - Won't Support Right Correctly")
        
        if self.bestMetrics == None:
            return False
        #Best metrics work
        return False
        """
        -curmetrics could have a span and this guy doesn't have left 
         or right seedContig
        #This could be a problem when we're missing left or right 
        -support in the first place 
        #How do we check for that? fuck it.
        -curmetrics could have a good right and left support, and this 
         guy doesn't have left or right seedContig
        # Totally reasonable. This guy can't span and can't make a better 
         left/right
        -curMetrics used more reads than this guy? Can't implement this, yet.
        """
        
    def storeSuccess(self, levelObj):
        fasta = os.path.join(levelObj.bankPath,"out.fasta")
        qual = os.path.join(levelObj.bankPath,"out.qual")
        reMap = os.path.join(levelObj.bankPath, "remapping.m5")
        
        fout = open("fillingMetrics.json",'w')
        fout.write(json.dumps(dict(self.bestMetrics),indent=4))
        fout.close()
        
        logging.debug("Current Best Metrics -\n%s" % \
            json.dumps(dict(self.bestMetrics),indent=4))
        
        exe("cp %s %s" % (fasta, "output.fasta"))
        exe("cp %s %s" % (qual, "output.qual"))
        exe("cp %s ." % (reMap))
    
    def reMapAndSupport(self, levelObj):
        """
        Takes a level object and remaps, resupports.
        """
        logging.debug("ReMapping")
        remapFile = os.path.join(levelObj.bankPath, "remapping.m5")
        """New blasr -- only in custom build, not SMRTAnalysis 1.3.1
        logging.debug(exe(("blasr %s %s -noSplitSubreads -nproc %d -m 5 -out %s" \
             " -scoreMatrix \"-5 6 6 6 0 6 -5 6 6 0 6 6 -5 6 0 6 6 6 -5 0 0 0 0 0 0 \"" \
            ) % (levelObj.outputFasta,  self.reference, \
            self.nproc, remapFile)))
        
        """ #stable blasr
        logging.debug(exe(("blasr %s %s -sa %s.sa -noSplitSubreads -nproc %d -m 5 -out %s" \
            ) % (levelObj.outputFasta, self.reference, self.reference, \
            self.nproc, remapFile)))
        #"""
        
        #Resupport.
        #Should return or create metrics here?
        return self.mySupporter.resupport(remapFile, levelObj)

    def quick_reMapAndSupport(self, levelObj):
        """
        Takes a level object and remaps, resupports.
        This is quicker because it will map to a subset of the flanking contigs
        We will make Extraction provide the reference in <gapFolder>/contigRef.fasta
        #
        This needs to be severely refactored. 
        """
        logging.debug("Quick ReMapping")
        remapFile = os.path.join(levelObj.bankPath, "remapping.m5")
        contigRef = os.path.join(levelObj.bankPath, "../contigRef.fasta")
        
        logging.debug(exe(("blasr %s %s -noSplitSubreads -nproc %d -m 5 -out %s" \
             #" -scoreMatrix \"-5 6 6 6 0 6 -5 6 6 0 6 6 -5 6 0 6 6 6 -5 0 0 0 0 0 0 \"" \
            ) % (levelObj.outputFasta, contigRef, \
            self.nproc, remapFile)))
        
        #Let's shift the coordinates here.
        fh = open(os.path.join(levelObj.bankPath,"../delta"),'r')
        #We'll make this file during extraction. 
        #It will hold how to shift the alignment to their absolute position
        delta = int(fh.readline().strip())
        real_tseqlength = int(fh.readline().strip())
        fh.close()
        
        alignments = M5File(remapFile)
        fout = open(remapFile,'w')
        for align in alignments:
            align.tstart += delta
            align.tend += delta
            align.tseqlength = real_tseqlength
            #Need to make sure we output in the old coordinate system
            fout.write(str(align)+"\n")
        fout.close()
        #Resupport.
        #Should return or create metrics here?
        return self.mySupporter.resupport(remapFile, levelObj)



class ReSupport():
    """
    Take the realignment file, figure out best support and the
    Sequence/information that goes with it.
    """
    def __init__(self, gapInfo, supportType):
        """
        remappingFn - mapping results to evaluate
        gapInfo - Gap object of who we're looking for
        """
        self.gapInfo = gapInfo
        self.myGap = self.gapInfo.values()[0]
        logging.debug(self.myGap)
        self.supportType = supportType
        self.classifier = SupportClassifier(self.gapInfo)
    
    def resupport(self, remappingFn, levelObj):
        #Group same reads
        alignments = defaultdict(list) 
        mapping = M5File(remappingFn) 
        for hit in mapping: 
            alignments[hit.qname].append(hit) 
        
        alignments = dict(alignments)   
        
        gapCan = GapCans()
        #Support each one
        for group in alignments.keys():
            #First, force same scaffolding hits
            logging.debug("Running Resupport on read %s" % (group))
            filteredAlignments = []
            for hit in alignments[group]:
                if hit.tname == self.myGap.scaffold \
                        and (hit.qname == levelObj.leftSeedContig \
                        or hit.qname == levelObj.rightSeedContig):
                    filteredAlignments.append(hit)
            if len(filteredAlignments) == 0:
                logging.debug("Unfillable: " + alignments[group][0].qname +\
                    " didn't map the correct scaffolding")
                continue
            alignments[group] = filteredAlignments
            #Need to return a gapCans object here?
            #Then I can build metrics off of that.
            mySupport = self.classifier.classifyReads(alignments[group])
            for gap in mySupport.keys():
                gapCan.extend(mySupport[gap])
        
        return gapCan
        
FindExpandGap = re.compile("[^N\-]([N\-]*N[N\-]*)[^N\-]")   

class SupportMetrics(dict):
    
    def __init__(self, curLevel, myGap, gapCan, supportType):
        super(dict)
        self.curLevel = curLevel
        self.myGap = myGap
        self.gapCan = gapCan
        self.supportType = supportType
        self.buildMetrics()
    
    
    def buildMetrics(self):
        self["LeftContig"], self["RightContig"], self["SpansGap"] \
            = self.getSupportingReads()
        
        self["FillSequence"], self["FillQual"] = self.getFillSequence()
        self["FillLength"] = len(self["FillSequence"].replace('N',''))
        self["ContigsInOutput"] = len(self.curLevel.myContigs.values())
        self["Accuracy"] = self.getAccuracy()
        self["GapName"] = self.myGap.name
        self["GapPredictedSize"] = self.myGap.end - self.myGap.start
    
    def getSupportingReads(self):
        left = False
        right = False
        span = False
        
        for alignment in self.gapCan["LeftContig"]:
            if alignment.qname == self.curLevel.leftSeedContig:
                left = alignment.qname
        
        for alignment in self.gapCan["RightContig"]:
            if alignment.qname == self.curLevel.rightSeedContig:
                right = alignment.qname
        
        for alignment in self.gapCan["SpansGap"]:
            if alignment.qname == self.curLevel.rightSeedContig \
                    and alignment.qname == self.curLevel.leftSeedContig:
                span = alignment.qname
        
        return left, right, span
    
    def numReadsSupporting(self):
        """
        Current data structures don't allow this to be computed efficiently
        """
        count = 0
        if self["LeftContig"]:
            count += self.curLevel.leftSeedContigDepth
        if self["RightContig"]:
            count += self.curLevel.rightSeedContigDepth
        #this is inflating my numbers sometimes
        if self["SpansGap"]:
            #Don't want redundancy?
            if self.curLevel.rightSeedContig != self.curLevel.leftSeedContig:
                count += self.curLevel.leftSeedContigDepth
            count += self.curLevel.rightSeedContigDepth
        return count
     
    def getFillSequence(self):
        """
        Grabs the sequence that we've identified as going inside of the gap.
            --I hate how this code flows/looks.
        """
        if self["SpansGap"]:
            if len(self.gapCan["SpansGap"]) > 1:
                #Edge Case - There is no reason this can't happen, 
                #but it should never happen. So, I'm putting this 
                #in place as a safety net until I make it happen.
                logging.warning("Complex Span Extraction! Aborting.")
                sys.exit(1)
            
            m5 = self.gapCan["SpansGap"][0]
            
            #Find gap and realign around it
            for extGap in FindExpandGap.finditer(m5.targetSeq): 
                gapStart = extGap.start() + 1
                gapEnd = extGap.end() - 1
                count = m5.targetSeq.count('-',0, gapStart) 
                logging.debug("Found Gap at %d, %d" % (gapStart - count + m5.tstart, self.myGap.start))
                if (gapStart - count + m5.tstart) == self.myGap.start:
                    print m5.targetSeq
                    realign(m5, span=True, gapStart=gapStart, gapEnd=gapEnd)
                    print m5.targetSeq
                    break
            #Pulling sequence from m5 file offers better accuracy
            for extGap in FindExpandGap.finditer(m5.targetSeq): 
                #This has extGap.start gap.extGap.end 
                #If you have noticed how I use +1 -1 all around these
                #Gap Finding REs, it's because re.finditer will give
                #coordinates around the entire match, not the match's 
                #group. Therefore, there is always an extra base at the
                #beginning and end of the reMatch.span().

                gapStart = extGap.start() + 1
                gapEnd = extGap.end() - 1
                count = m5.targetSeq.count('-',0, gapStart) 
                logging.debug("where %d, %d" % (gapStart - count + m5.tstart, self.myGap.start))
                if (gapStart - count + m5.tstart) == self.myGap.start:
                    sequence = m5.querySeq[gapStart:gapEnd].replace('-','')
                    
                    start = m5.qstart + gapStart - count
                    end = start + len(sequence) 
                    if len(sequence) == 0:
                        #Negative gap handling.
                        #Find the coordinates where we think
                        #we actually have sequence again.
                        #Currently 2 consecutive bases need to be mapped to end trimming
                        left =  re.compile('.*[ATCGatcg]{2,}')
                        right = re.compile('[ATCGatcg]{2,}.*')
                        self["LeftTrim"] = left.search(m5.querySeq[:gapStart]).end()-gapStart
                        self["RightTrim"] = right.search(m5.querySeq[gapEnd:]).start()
                        logging.debug("LeftTrimAmt: %d" % self["LeftTrim"])
                        logging.debug("RightTrimAmt: %d" % self["RightTrim"])
                        if abs(self["LeftTrim"]) > MAX_TRIM or abs(self["RightTrim"]) > MAX_TRIM:
                            self["SpansGap"] = False

                    #short circuit
                    break
            
            logging.debug("SpanFill %d" % (end- start))
            qual = QualFile(self.curLevel.outputQual)
            self["SpanStart"] = start
            self["SpanEnd"] = end
            self["SpanStrand"] = '-' if m5.negStrand else '+'
            quality = " ".join(map(str,qual[m5.qname][start: end]))
        else:
            fasta = FastaFile(self.curLevel.outputFasta)
            qual = QualFile(self.curLevel.outputQual)
            leftSeq = rightSeq = ""
            leftQual = []; rightQual = []
            
            if self["LeftContig"]:
                logging.debug("LeftContig has %d supporting reads" % \
                    len(self.gapCan["LeftContig"]))
                
                m5 = self.gapCan["LeftContig"][0]
                realign(m5, left=True)
                if m5.tstrand == '0':
                    start = m5.qend + ( self.myGap.start-m5.tend )  
                    end = m5.qseqlength
                elif m5.tstrand == '1':
                    start = 0
                    end = m5.qstart
                
                logging.debug("LeftFill of %d" %(m5.qseqlength - start))
                self["LeftStart"] = start
                self["LeftEnd"] = end
                leftSeq = fasta[m5.qname][start:end]
                leftQual = qual[m5.qname][start:end]
                
                self["LeftStrand"] = '-' if m5.negStrand else '+'
                if m5.negStrand:
                    leftSeq = leftSeq.translate(revComp)
            
            if self["RightContig"]:
                logging.debug("RightContig has %d supporting reads" % \
                    len(self.gapCan["RightContig"]))
                
                m5 = self.gapCan["RightContig"][0]
                realign(m5, left=False)
                if m5.tstrand == '0':
                    start = 0
                    end = m5.qstart
                elif m5.tstrand == '1':
                    start = m5.qend
                    end = m5.qseqlength
                #start = 0
                #end = m5.qstart - ( m5.tstart-self.myGap.end )
                
                logging.debug("RightFill of %d" %(end))
                
                rightSeq = fasta[m5.qname][start:end]
                rightQual = qual[m5.qname][start:end]
                self["RightStart"] = start
                self["RightEnd"] = end
                self["RightStrand"] = '-' if m5.negStrand else '+'
                if m5.negStrand:
                    rightSeq = rightSeq.translate(revComp)
            
            #Consolidate Left and Right Seq
            gapSize = self.myGap.end - self.myGap.start
            remainingGap = gapSize - len(rightSeq) - len(leftSeq)
            
            logging.debug("GapSize - %d" % (gapSize))
            
            if remainingGap <= 0:
                
                logging.debug("Gap Size Underestimated")
                
                self["GapUnderestimated"] = -remainingGap
                remainingGap = 25
            
            sequence = leftSeq+("N"*remainingGap)+rightSeq
            leftQual.extend([0]*remainingGap)
            leftQual.extend(rightQual)
            quality = " ".join(map(str,leftQual))
        
        return sequence, quality
    
    def getAccuracy(self):
        """
        looks at matches, mismatches, insertions, deletions in alignment(s),
        excluding N's (not hurting us)
        """
        correct = 0
        total = 0.0
        gaps = 0
        for alignments in self.gapCan.values():
            for alignment in alignments:
                correct += alignment.nMatch
                total += alignment.nInsert \
                    + alignment.nDelete \
                    + alignment.nMismatch
                gaps += alignment.targetSeq.count('N')
        
        total += correct - min(gaps, self.myGap.end-self.myGap.start)
        return correct/total
    
    def compare(self, otherMetrics):
        """
        Returns True if otherMetrics are better than self.
            --Need more rigorous and intelligent approach   
        """
        """
        Keep This Simple For First Implementation
        first:
        span > left&right > left xor right > none
        second: (and if there is a tie in first)
        accuracyHigh > accuracyLow 
            Give a tolerance for this one. if they're close, check filling.
        amount of filling. -- more filling is better unless it's overfilling
        
        I need to stop those reads with Left and Right Seeds, but don't span
        from getting through
        
        I need to not have span automatically getting through over left and right
        when the span would put nothing into a big gap but left and right would.
        
        Finally I need to Save information inside of fillingMetrics about what kind
        of fill I have. This is super important.
        """
        logging.debug("Metics Manually - \nCurMetrics" + \
            json.dumps(dict(self),indent=4) + \
            "\n NewMetrics" + json.dumps(otherMetrics, indent=4))

        #Both SpanGap
        if self["SpansGap"] and otherMetrics["SpansGap"]:
            curAccuracy = self["Accuracy"]
            othAccuracy = otherMetrics["Accuracy"]
            curGapDev = abs(self["GapPredictedSize"] - self["FillLength"])
            othGapDev = abs(otherMetrics["GapPredictedSize"] \
                    - otherMetrics["FillLength"])
            #The new metrics will need to be more accurate and 
            #closer to the predicted gap size
            if othAccuracy > curAccuracy and othGapDev < curGapDev:
                return True
            #The new metrics will need to be MUCH more accurate, 
            #but still within 30% of the predicted gap size
            if (othAccuracy - .05) > curAccuracy \
                  and othGapDev < (otherMetrics["GapPredictedSize"]*.30):
                return True
            
            return False
        
        if self["SpansGap"] \
                and otherMetrics["LeftContig"] \
                and otherMetrics["RightContig"]:
            curAccuracy = self["Accuracy"]
            othAccuracy = otherMetrics["Accuracy"]
            curGapDev = abs(self["GapPredictedSize"] - self["FillLength"])
            othGapDev = abs(otherMetrics["GapPredictedSize"] \
                            - otherMetrics["FillLength"])
            #if current spans is outisde of 30% predicted gap size, 
            #we can't reasonably keep it
            if curGapDev > (self["GapPredictedSize"]*.60) \
                     and othGapDev < (otherMetrics["GapPredictedSize"]*.60):
                return True
            #otherwise, we'll keep it
            return False
        
        if otherMetrics["SpansGap"] \
               and self["LeftContig"] \
               and self["RightContig"]:
            
            curAccuracy = self["Accuracy"]
            othAccuracy = otherMetrics["Accuracy"]
            curGapDev = abs(self["GapPredictedSize"] - self["FillLength"])
            othGapDev = abs(otherMetrics["GapPredictedSize"] \
                            - otherMetrics["FillLength"])
            #if new span is within 30% of predicted gap size, take it.
            if othGapDev < (otherMetrics["GapPredictedSize"] * 0.60):
                return True
            #Or if span is closer than we are
            elif curGapDev >= othGapDev: 
                return True
            #elif self.has_key("RightTrim" or self.has_key("LeftTrim"):
            #otherwise, we'll keep it
            return False
        
        if self["LeftContig"] \
               and self["RightContig"] \
               and otherMetrics["LeftContig"] \
               and otherMetrics["RightContig"]:

            curAccuracy = self["Accuracy"]
            othAccuracy = otherMetrics["Accuracy"]
            curGapDev = abs(self["GapPredictedSize"] - self["FillLength"])
            othGapDev = abs(otherMetrics["GapPredictedSize"] \
                            - otherMetrics["FillLength"])
            #The new metrics will need to be more accurate \
            #and closer to the predicted gap size
            if othAccuracy > curAccuracy and othGapDev < curGapDev:
                return True
            #The new metrics will need to be MUCH more accurate, 
            #but still within 30% of the predicted gap size
            if (othAccuracy - .05) > curAccuracy \
                   and othGapDev < (otherMetrics["GapPredictedSize"]*.30):
                return True
                
            return False
        
        if self["SpansGap"] \
               and (otherMetrics["LeftContig"] \
               or otherMetrics["RightContig"]):

            return False

        if otherMetrics["SpansGap"] \
               and (self["LeftContig"] \
                or self["RightContig"]):

            return True
        
        if (self["LeftContig"] and self["RightContig"]) \
               and (otherMetrics["LeftContig"] or otherMetrics["RightContig"]):
            return False
        
        if (otherMetrics["LeftContig"] and otherMetrics["RightContig"]) \
               and (self["LeftContig"] or self["RightContig"]):
            return True
        
        return self["Accuracy"] > otherMetrics["Accuracy"]
        
#Need to clean this up, put it in an appropriate object
def pushSeq(align, stop):
    target = list(align.targetSeq)
    query = list(align.querySeq)
    compSeq = list(align.compSeq)
    
    for pos in range(len(align.targetSeq[:stop])):
        if align.targetSeq[pos] == '-':
            i = re.match('-+[ATCGatcg]?', align.targetSeq[pos:stop])
            if align.targetSeq[pos+i.end()-1] == align.querySeq[pos]:
                target[pos] = align.targetSeq[pos+i.end()-1]
                compSeq[pos] = '|'
                target[pos+i.end()-1] = '-'
                compSeq[pos+i.end()-1] = '*'
                align.targetSeq = "".join(target)
                align.compSeq = "".join(compSeq)
    

def realign(align, span=False, left=False, gapStart=None, gapEnd=None):
    orig = "".join(list(align.targetSeq))
    if span:#Push both sides of a span sequence
        pushSeq(align, gapStart)
        align.targetSeq = align.targetSeq[::-1]
        align.compSeq = align.compSeq[::-1]
        align.querySeq = align.querySeq[::-1]
        pushSeq(align, len(align.targetSeq) - gapEnd)
        align.targetSeq = align.targetSeq[::-1]
        align.compSeq = align.compSeq[::-1]
        align.querySeq = align.querySeq[::-1]
    elif left:#Push upto the end of the alignment
        pushSeq(align, len(align.targetSeq))
    else:#Push upto the start of the alignment (just use reverse of seqence)
        align.targetSeq = align.targetSeq[::-1]
        align.compSeq = align.compSeq[::-1]
        align.querySeq = align.querySeq[::-1]
        pushSeq(align,len(align.targetSeq))
        align.targetSeq = align.targetSeq[::-1]
        align.compSeq = align.compSeq[::-1]
        align.querySeq = align.querySeq[::-1]
    align.move = True


