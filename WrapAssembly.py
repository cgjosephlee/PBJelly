#!/usr/bin/env python
"""
TODO

need --quick option that allows mapping to 
contigRef.fasta 
AssemblyAssesor will do the work,
I'll just need to pass that the option is 
enabled to it.

Enable --tempDir.
    Working in on-node storage will speed things up and
    reduce cluster io stress.
    OLCAssembly already enables it, I just need levelObjects
    to track the folder
"""
import sys
import logging
import os
import shutil
import json
from optparse import OptionParser
from CommandRunner import exe
from collections import defaultdict
from FileHandlers import FastaFile, QualFile, GapInfoFile
from AssemblyAssessor import *

class WrapAssembly():
    
    def __init__( self ):
        self._parseOpts()
        self.__initLog()
        self.inputFasta = "input.fasta"
        self.inputQual = "input.qual"
        self.setupSeedKnowledge()
    
    def _parseOpts( self ): 
        parser = OptionParser( usage=USAGE ) 
        parser.add_option("--debug", action="store_true", \
                          help="Increases verbosity of logging" )
        parser.add_option("--nproc",type=int, default=1, \
                          help="Number of processors available for use")
        self.opts, args = parser.parse_args()
        if len(args) != 2:
            parser.error("Expected 2 Arguments!")
        self.inputDir, self.gapFn = args
        
        self.gapInfo = GapInfoFile(self.gapFn)
        if self.inputDir.endswith('/'):
            self.inputDir = self.inputDir[:-1]
        
        self.gapName = os.path.basename(self.inputDir)
        self.myGap = {self.gapName: self.gapInfo[self.gapName]}
        fh = open(os.path.join(self.inputDir,"support.gapCans"),'r')
        self.supportInfo = json.load(fh)
        fh.close()
        os.chdir(self.inputDir)


    def __initLog( self ):
        """Sets up logging based on command line arguments. """
        
        logLevel = logging.DEBUG if self.opts.debug else logging.INFO 
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        logging.basicConfig( stream=sys.stderr, \
                             level=logLevel, \
                             format=logFormat )
        logging.info("Running %s" % " ".join(sys.argv))
    
    def setupSeedKnowledge(self):
        """
        Using the flanking region information to figure out who 
        left and right, set 
        """
        #Pretty much like seqGrab?
        #I could shorten this...
        self.leftSeed = self.gapInfo[self.gapName].leftContig
        logging.debug("LeftSeed - %s" % self.leftSeed)
        self.rightSeed = self.gapInfo[self.gapName].rightContig
        logging.debug("RightSeed - %s" % self.rightSeed)
    
    def runAssembly(self, levelObj):
        """
        uses the previous level (levelObj)
        to setup and then run the new level (level)
        
        """
        name = "level%d" % (levelObj.levelId+1)
        os.mkdir(name)
        wDir = name
        
        levelObj.extractOutputs(wDir+"/input")
        
        t = "--threshold=200" if levelObj.threshold else ""
        c, o, e = exe(("OLCAssembly.py input.fasta input.qual " \
                        "--nproc=%d --rename %s --workDir=%s %s" \
                        % (self.opts.nproc, name, name, t)), timeout=30)
        logging.debug("OLCAssembly %s returned CODE %d \n STDOUT - %s \n STDERR - %s" \
                % (wDir, c, str(o), str(e)))
        if c == 214:#Timedout
            logging.error("OLCAssembly took more than 30 minutes. Something is likely wrong. Exiting %s" % self.gapName)
            sys.exit(1)
        return c, wDir
    
    def run(self):
        #Take information about support to setup expectations
        if len(self.supportInfo["SpansGap"]) > 0 \
              or (len(self.supportInfo["LeftContig"]) > 0 \
              and len(self.supportInfo["RightContig"]) > 0):
            self.supportType = "SpansGap"
        
        elif len(self.supportInfo["LeftContig"]) > 0:
            self.supportType = "LeftContig"
            PROCEDUREORDER.remove("PullRight")
        
        elif len(self.supportInfo["RightContig"]) > 0:
            self.supportType = "RightContig"
            PROCEDUREORDER.remove("PullLeft")
        else:
            #Should never happen.
            logging.warning("Couldn't find any support in support.gapCans")
            raise NotImplementedError("Word-up, dog.")
        
        self.myAssessor = AssemblyAssessor(self.gapFn, self.myGap, \
                    self.opts.nproc, self.supportType)
        
        self.myTrace = Trace(self.inputFasta, self.inputQual, self.myAssessor)
        
        while True:
            curLevel = self.myTrace.level + 1
            logging.info("Running Level%d Assembly" % (curLevel))
            
            retCode, bankPath = self.runAssembly(self.myTrace.curLevel) 
            
            logging.info("Evaluating Level%d" % (curLevel))
            #Need to handle retCode before trying to update the trace.
            threshold = self.myTrace.curLevel.threshold
            
            if curLevel == 0:
                leftSeed = self.leftSeed
                rightSeed = self.rightSeed
            else:
                leftSeed = self.myTrace.curLevel.leftSeedContig
                rightSeed = self.myTrace.curLevel.rightSeedContig
            
            self.myTrace.updateTrace(bankPath, threshold, rightSeed, leftSeed)

            #raw_input("Continue Evaluation of Level%d? " % (curLevel))
            
            if retCode == 0:
                #tell myTrace to figure out what is possible next
                isFinished = self.myTrace.evaluate()
            else:
                logging.error("OLCAssembly Returned Error Code %d" % (retCode))
                #tell trace to go back directly.
                isFinished = self.myTrace.backUp()
            
            if isFinished:
                break
            #if mytrace is successful, break, go do what's next.
            #if mytrace is failed, (eshausthed all options), break and report:
            #   Check for a best, if so break and go do what's next
            #   otherwise, quit, report failure
        
            #raw_input("Continue onto Level%d? " % (self.myTrace.level+1))
        logging.info("WrapAssembly found assembly == %s" % str(isFinished))
        #Extract Sequence? -- or give that to the next step to do?

#Maxum depth we can try to find a reasonable assembly at
MAXLEVEL = 5            
class Trace():
    """
    Keeps track of levels in order to easily trace back and 
    undo and assemblies performed.
    
    """
    
    def __init__(self, fasta, qual, assessor):
        self.level = -1
        self.traceDict = {self.level:StartLevel(fasta,qual)}
        self.curLevel = self.traceDict[self.level]
        self.assessor = assessor
    
    def updateTrace(self, bankPath, threshold, rightSeed, leftSeed):
        ##1 3 1953 bases, 00000000 checksum.
        #m110614_064618_00121_c100146272555500000315040906251190_s1_p0/1211/356_2148(0) [] **Trash(forthesepurposes)?**
        self.level += 1
        self.traceDict[self.level] = Level(self.level, bankPath, threshold)
        self.curLevel = self.traceDict[self.level]
        
        self.curLevel.setSeeds(rightSeed, leftSeed) 
       
    def evaluate(self):
        """
        Goes through all the information available at the current level
        and sets removes any procedures that are not viable options.
        """
        #Evaluate current level through Trace to see what is possible moving forward.
        #So Trace needs to keep logic of what is possible per level. Level can keep the actual
        #possibilities
        #if nothing, set self.myTrace.curLevel back one and update it's possibilities.
        #so, all this logic needs to go inside of Trace.
        # all of this logic needs to update the possibilites. 
        # if no possibilities, go back. Then go through again.
        
        if self.level == MAXLEVEL:
            self.backUp()
        
        #First thing I should do is evaluate the assembly,
        # if I have an acceptable (likely best) result, store it and stop.
        # or if I have a decent result, store it and continue.
        # otherwise, I should keep on keeping on. 
        #Check been assessed first to prevent double checking.
        if not self.curLevel.beenAssessed and self.assessor.assess(self.curLevel):
            return True
        
        self.curLevel.updateAvailableProcedures()
        
        if self.curLevel.runProcedure():
            return False
        else:
            return self.backUp()
        
    def backUp(self):
        #Tell trace to delete delete this current level 
        #Then go to evaluate so the previous level can check things out.
        if self.level == -1:
            logging.debug("Tried to backUp from StartNode -- Assembly Finished")
            return True
        
        logging.debug("BackingUp from %d" % (self.level))
        shutil.rmtree(self.curLevel.bankPath)
        del(self.traceDict[self.level])
        self.level -= 1
        self.curLevel = self.traceDict[self.level]
        
        return self.evaluate()
        

class StartLevel():
    """
    Place holder for original inputs to whole process
    
    #Something very unsafe about this, but python's 
    constructor overloading is not great (like Java's)
    also, Multiple Inheritance doesn't feel quite right.
    """
    def __init__(self, fasta, qual):
        self.levelId = -1
        self.outputFasta = fasta
        self.outputQual = qual
        self.threshold = False
        self.beenAssessed = False
        self.leftSeed = self.rightSeed = self.rightSeedContig = self.leftSeedContig = None
        self.myContigs = {}
        #push contigs automatically, so only threshold can be done
        self.availableProcedures = {"Threshold":StartLevel.threshold}#, \
                #   "Split":StartLevel.split}
        #Maybe remove inputs... we'll see how well I implement it.
    
    def extractOutputs(self, outName):
        """
        copies the reads left in this level to outName
        """
        exe("cp %s %s" % (self.outputFasta, outName+".fasta"))
        exe("cp %s %s" % (self.outputQual, outName+".qual"))
    
    def removeProcedure(self, proc):
        """
        procedure is a set of operations available to
        be run on this level.
        """
        if self.hasProcedure(proc):
            del(self.availableProcedures[proc])
    
    def hasProcedure(self, proc):
        """
        be sure procedure is still available
        """
        return self.availableProcedures.has_key(proc) 
    
    def runProcedure(self):
        """
        Remove, then run the procedure
        to be run on this level
        """
        #Something has been changed,
        #So no assessment, yet
        self.beenAssessed = False
        if len(self.availableProcedures) == 0:
            return False
        procName = self.availableProcedures.keys()[0]
        procExe = self.availableProcedures[procName]
        procExe(self)
        self.removeProcedure(procName)
        return True
        
    def threshold(self):
        logging.debug("Threshold in beginning")
        self.threshold = True
    
    def updateAvailableProcedures(self):
        if self.threshold:
            self.removeProcedure("Threshold")

class Level():
    """
    Object for holding information about what reads are used in each contig for an individual level's run
    And what operations (procedures) can be run on these levels.
    """
    def __init__(self, levelId, bankPath, threshold):
        self.levelId = levelId
        
        self.bankPath = bankPath
        self.bank = os.path.join(bankPath,"out.bank")
        self.inputFasta = os.path.join(bankPath,"input.fasta")
        self.inputQual = os.path.join(bankPath,"input.qual")
        self.outputFasta = os.path.join(bankPath,"out.fasta")
        self.outputQual = os.path.join(bankPath,"out.qual")
        
        #Was threshold used in this level?
        self.threshold = threshold
        
        self.rightSeed = None
        self.leftSeed = None
        self.rightSeedContig = None
        self.leftSeedContig = None
        self.removedContigs = {}
        self.pulledReads = []
        
        self.availableProcedures = list(PROCEDUREORDER) 
        self.beenAssessed = False
        self.__parseContigs()
        
    def __parseContigs(self):
        c,o,e = exe("bank2contig %s -L" % (self.bank))
        if c != 0:
            logging.error("Problem parsing bank %s" % (self.bank))
            self.myContigs = {}
            return
        
        self.myContigs = defaultdict(list)
        
        for line in o.split('\n')[:-1]:
            if line.startswith("##"):
                curContig = line.split(' ')[0][2:]
            elif line.startswith("#"):
                read = line.split(' ')[0][1:]
                read = read[:read.rindex('(')]
                self.myContigs[curContig].append(read)
            else:
                logging.debug("BankInfoLine (%s)" % line)
        
        self.myContigs = dict(self.myContigs)
        logging.debug("Produced Contigs:\n" + json.dumps(self.myContigs, indent=4))
    
    def setSeeds(self, right, left): 
        self.rightSeed = right
        self.leftSeed = left
        
        for contig in self.myContigs.keys():
            if right in self.myContigs[contig]:
                self.rightSeedContig = contig
            if left in self.myContigs[contig]:
                self.leftSeedContig = contig
        
        
        logging.debug("LeftSeedContig - %s" % (self.leftSeedContig))
        logging.debug("RightSeedContig - %s" % (self.rightSeedContig))
    
    def removeContig(self, contig):
        #Not checking is unsafe
        if self.myContigs.has_key(contig):
            self.removedContigs[contig] = list(self.myContigs[contig])
            del(self.myContigs[contig])
    
    def restoreContig(self, contig):
        #Not checking is unsafe
        if self.removedContigs.has_key(contig):
            self.myContigs[contig] = list(self.removedContigs[contig])
            del(self.myContigs[contig])
    
    def extractOutputs(self, outName):
        """
        takes out all of the contigs
        puts them into outName.fasta outName.qual
        """
        if len(self.removedContigs.keys()) == 0 and len(self.pulledReads) == 0:
            #Nothing has changed, just copy them
            self.copyOutputs(outName)
            return
        fout = open(outName+".fasta",'w')
        qout = open(outName+".qual",'w')
        
        #Keep what we've got (removed contigs won't be here)
        fin = FastaFile(self.outputFasta)
        qin = QualFile(self.outputQual, convert=False)
        for i in self.myContigs:
            fout.write(">"+i+"\n"+fin[i]+"\n")
            qout.write(">"+i+"\n"+qin[i]+"\n")
        
        #Are there also reads to pull?
        if len(self.pulledReads) != 0:
            fin = FastaFile(self.inputFasta)
            qin = QualFile(self.inputQual, convert = False)
            for i in self.pulledReads:
                fout.write(">"+i+"\n"+fin[i]+"\n")
                qout.write(">"+i+"\n"+qin[i]+"\n")  
        fout.close()
        qout.close()
        
    def pullRead(self, read):
        """
        Pulls a read from the into the output
        """
        #Not a save operation without knowing every input's name. 
        #it is is we make sure it's one of the seed's names!
        self.pulledReads.append(read)#you mess up in extractOutput
        
    def copyOutputs(self, outName):
        """
        copies outputs in this level to outName
        """
        exe("cp %s %s" % (self.outputFasta, outName+".fasta"))
        exe("cp %s %s" % (self.outputQual, outName+".qual"))
    
    def copyInputs(self, outName):
        """
        copies input in this level to outName
        should never be used
        """
        exe("cp %s %s" % (self.inputFasta, outName+".fasta"))
        exe("cp %s %s" % (self.inputQual, outName+".qual"))
    
    def removeProcedure(self, proc):
        """
        procedure is a set of operations available to
        be run on this level.
        """
        try:
            self.availableProcedures.remove(proc)
        except ValueError:
            #Totally acceptable that things
            #will be asked to be removed many times
            pass    
    
    def hasProcedure(self, proc):
        """
        be sure procedure is still available
        """
        return self.availableProcedures.count(proc) == 1
    
    def runProcedure(self):
        """
        Remove, then run the procedure
        to be run on this level
        """
        if len(self.availableProcedures) == 0:
            return False
        procName = self.availableProcedures[0]
        logging.debug("Performing %s on Level%d" % (procName, self.levelId))
        procExe = PROCEDURES[procName]
        procExe(self)
        self.removeProcedure(procName)
        return True
    
    def pushContigs(self):
        """
        Don't Touch Anything, Just push this forward
        """
        return  
    
    def pullLeft(self):
        """
        find the leftSeed in inputs, pull it over.
        """
        self.pulledReads.append(self.leftSeed)
    
    def pullRight(self):
        """
        find the rightSeed in inputs, pull it over
        """
        self.pulledReads.append(self.rightSeed)
        
    def removeExtraContigs(self):
        """
        removes all contigs in this level's output that
        do not contain the seed
        I'm not sure how I'm going to do this, but it IS! 
        the last resort (no going back from this)
        """
        for contig in list(self.myContigs.keys()):
            if contig != self.rightSeedContig and contig != self.leftSeedContig:
                self.removeContig(contig)
    
    def threshold(self):
        """
        lower the threshold
        """
        self.threshold = True
    
    def updateAvailableProcedures(self):
        if self.leftSeedContig == None and self.rightSeedContig == None:
            self.removeProcedure("PushContigs")
        
        if len(self.myContigs.keys()) <= 1:
            self.removeProcedure("PushContigs")
        
        if self.leftSeedContig != None or self.leftSeed == None:
            self.removeProcedure("PullLeft")
        
        if self.rightSeedContig != None or self.rightSeed == None:
            self.removeProcedure("PullRight")
        #Sorry to break the proramming style, but this logic different.
        remove = True
        for contig in self.myContigs:
            if contig == self.rightSeedContig or contig == self.leftSeedContig:
                continue
            remove = False; break#Just need a single non-seed contig and we can remove extras
        if remove:
            self.removeProcedure("RemoveExtraContigs")
        
        if len(self.myContigs.keys()) < 2 or self.threshold:
            self.removeProcedure("Threshold")
        
PROCEDURES = {"PushContigs":Level.pushContigs ,\
              "PullLeft":Level.pullLeft ,\
              "PullRight":Level.pullRight ,\
              "RemoveExtraContigs":Level.removeExtraContigs, \
              "Threshold":Level.threshold }
PROCEDUREORDER = ["PushContigs", "PullLeft", "PullRight", "RemoveExtraContigs", "Threshold"]


USAGE = """WrapAssembly.py <inputDir> <gapInfo>

    Also, look at OLCAssembly.py --help for more information 
    on parameters that can be passed into the assembler."""
if __name__ == "__main__":
    me = WrapAssembly()
    me.run()

