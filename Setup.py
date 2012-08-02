#!/usr/bin/env python
import sys, re, os, logging, subprocess
from optparse import OptionParser
from FileHandlers import FastaFile, QualFile, wrap, qwrap
from CommandRunner import exe

USAGE = """
USAGE: %prog <inputScaffolding.fasta> <inputScaffolding.qual> [outputContigs.fasta]

Take the input scaffolding and split it into contigs around all [Nn]+ 
"""

refParser = re.compile("(.*)\|(ref\d{7})/?(\d+)?$")

class Setup():
    def __init__(self):
        self.parseArgs()
        self.__initLog()
    
    def parseArgs(self):
        parser = OptionParser(USAGE)
        
        parser.add_option("-g", "--gapOutput", dest="gapOutput", \
            help="Create the table for gapInformation", default=None)
        parser.add_option("-i", "--index", dest="index", action="store_true", default=False, \
            help="Create the .sa index for faster blasr alignments")
        parser.add_option("--noRename", dest="rename", action="store_true", default=False, \
            help="Flag to indicate there will be no renaming of scaffolding.")
        parser.add_option("--debug",action="store_true", help="Increases verbosity of logging" )
        parser.add_option("--minGap", dest="minGap", type="int", default=25, help="Minimum number of consecutive Ns to be considered a gap. default=25")
        parser.add_option("--maxGap", dest="maxGap", type="string", default="", help="Maximum number of consecutive Ns to be considered a gap default=Inf")
        
        self.opts, args = parser.parse_args()
        if len(args) < 2:
            parser.error("Error! Incorrect number of arguments")
        if len(args) == 2:
            self.scaffInput = args[0]
            self.qualInput = args[1]
            if not self.scaffInput.endswith(".fasta"):
                parser.error("Reference must end in extension .fasta! Please rename it.")
        else:
            parser.error("Error! Incorrect number of arguments")
        
        qn = self.scaffInput[:self.scaffInput.rindex('.fasta')]
        self.contigsOutput = qn+".contigs.fasta"
        self.qualContigsOutput = qn+".contigs.qual"
        self.qualReferenceOutputName = qn+".qual"
        
        if not os.path.isfile(self.scaffInput):
            parser.error("Error! Scaffold File is not a file")

    def __initLog(self):
        """Logging"""
        logLevel = logging.DEBUG if self.opts.debug else logging.INFO
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
        logging.info("Running %s" % " ".join(sys.argv))
    
    def run(self):
        #Fasta Ref Output
        scaffTempName = self.scaffInput+".tempFasta"
        contigsOut = open(self.contigsOutput,'w')
        scaffOutput = open(scaffTempName, 'w')
        
        #Qual Ref Output
        qualTempName= self.qualInput+".tempQual"
        qualContigsOut = open(self.qualContigsOutput, 'w')
        qualOutput = open(qualTempName, 'w')
        
        #Gaps Output
        if self.opts.gapOutput != None:
            gapTableOut = open(self.opts.gapOutput,'w')
        else:
            gapTableOut = False
        
        logging.info("Creating reference sequence index names and contig split reference")
        
        refTemplate = "ref%07d"
        refId = 1
        
        reference = FastaFile(self.scaffInput)
        qualReference = QualFile(self.qualInput)    
        
        for key in reference:
            
            if self.opts.rename:
                scaffName, scaffIndex, null = refParser.match(key).groups()
            else:
                scaffIndex = refTemplate % refId
                scaffName = key.replace(' ','_')
            
            refId += 1
            
            scaffName = scaffName + "|" + scaffIndex
            scaffOutput.write(">"+scaffName+"\n"+wrap(reference[key])+"\n")
            
            qualOutput.write(">"+scaffName+"\n"+qwrap(qualReference[key])+"\n")
            
            gapCoords = []
            for gap in re.finditer("[^Nn]([Nn]{%d,%s})[^Nn]" % \
                    (self.opts.minGap, self.opts.maxGap), reference[key]):
                gapCoords.append([gap.start() + 1, gap.end() - 1])
            
            if len(gapCoords) == 0:#no Gaps
                contigName = scaffName+"/%i" % (0)
                contigsOut.write(">"+contigName+"\n"+wrap(reference[key])+"\n")
                qualContigsOut.write(">"+contigName+"\n"+qwrap(qualReference[key])+"\n")
                continue
            
            #Consolidate gaps that are too close -- indicating LQ regions.
            i = 0
            while i < len(gapCoords)-1:
                if gapCoords[i+1][0] - gapCoords[i][1] < 25:
                    gapCoords[i+1][0] = gapCoords[i][0]
                    del(gapCoords[i])
                else:
                    i += 1
            
            prevEnd = 0#Contig Start Tracking
            idx = 0
            #Make a contig and a gap for the first gap
            contigName = scaffName+"/%i" % (idx)
            newSeq = reference[key][prevEnd:gapCoords[0][0]]
            newQual = qualReference[key][prevEnd:gapCoords[0][0]]
            contigsOut.write(">"+contigName+"\n"+wrap(newSeq)+"\n")
            qualContigsOut.write(">"+contigName+"\n"+qwrap(newQual)+"\n")
            prevEnd = gapCoords[0][1]
            gapCoords[0][1]-gapCoords[0][0]
            
            if gapTableOut:
                gapTableOut.write("%s\t%i\t%i\t%s_%i_%i\n" \
                        % (scaffName, gapCoords[0][0], gapCoords[0][1], scaffIndex, idx, idx+1))

            #Now Go Through the rest of the gaps
            for i in range(1, len(gapCoords)):
                idx += 1
                contigName = scaffName+"/%i" % (idx)
                newSeq = reference[key][prevEnd:gapCoords[i][0]]
                newQual = qualReference[key][prevEnd:gapCoords[i][0]]
                contigsOut.write(">"+contigName+"\n"+wrap(newSeq)+"\n")
                qualContigsOut.write(">"+contigName+"\n"+qwrap(newQual)+"\n")
                prevEnd = gapCoords[i][1]
                gapCoords[i][1]-gapCoords[i][0]
            
                if gapTableOut:
                    gapTableOut.write("%s\t%i\t%i\t%s_%i_%i\n" \
                        % (scaffName, gapCoords[i][0], gapCoords[i][1], scaffIndex, idx, idx+1))
            
            #Finish up the last contig
            idx += 1
            contigName = scaffName+"/%i" % (idx)
            newSeq = reference[key][prevEnd:]
            newQual = qualReference[key][prevEnd:]
            contigsOut.write(">"+contigName+"\n"+wrap(newSeq)+"\n")
            qualContigsOut.write(">"+contigName+"\n"+qwrap(newQual)+"\n")
                   
        #Close shop
        scaffOutput.close()
        contigsOut.close()
        logging.debug(exe("mv %s %s" % (self.scaffInput, self.scaffInput+".original")))
        logging.debug(exe("mv %s %s" % (scaffTempName, self.scaffInput)))
        
        qualOutput.close()
        qualContigsOut.close()
        logging.debug(exe("mv %s %s" % (self.qualInput, self.qualInput+".original")))
        logging.debug(exe("mv %s %s" % (qualTempName, self.qualReferenceOutputName)))
        
        if gapTableOut:
            gapTableOut.close()
        
        if self.opts.index:
            logging.info("Creating .sa indexes for references")
            #logging.warning("There are unknown problems within creating indexes with sa writer. \
            #         If any stage that uses blasr and an index created here raises a \
            #        \"Segmentation Fault\" error, recreate the index with sawriter".replace('\t',''))
            
            #logging.debug(exe("sawriter %s.sa %s" % (self.contigsOutput, self.contigsOutput)))
            logging.debug(exe("sawriter %s.sa %s" % (self.scaffInput, self.scaffInput)))
        
        logging.info("Finished!")

if __name__ == '__main__':
    me = Setup()
    me.run()
