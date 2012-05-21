#!/usr/bin/env python
USAGE = """
Consolidates the reads over a single gap into a folder.

Extraction.py <reference.contigs.fasta> <reference.contigs.qual> \
<reference.gapInfo.bed> <inputDir> <outputDir> <inputFile> [<inputFile>...]
"""
"""
TODO
Change contigs.fasta to scaffolding.fasta
Then Extarct FLANKAMNT from gapInfoFile
second line of delta file should be scaffold length

instead of having delta as a file, it should be in support.gapCans
change support.gapCans to apriori.json

I'm going to use the marshall indexing method to randomly access the fasta files.
--This means I'm going to need two Fasta methods.
--You know, I only need this once, during extraction. So I'll just build a giant
index here using marshall. And then I'll use that to open whatever fasta file... right?
-- What if this results in too many open files? HOwW! grurr
"""
import sys, os, glob, logging, json, linecache
from optparse import OptionParser
from collections import defaultdict
from FileHandlers import FastaFile, QualFile, GapCans, GapInfoFile

FLANKAMT = 1000# Amount of flank to put in for the "Guided" Assembly
REFAMT = 5000# Amount of flank to remap to during assembly assessment
class Extraction():
    """
    Parses gapCans files and consolidates
    reads between files that map to the same gap.
    Also creates the qual file
    Final product is a set of folders for each
    gap that is supported so that it can be fed into 
    consensus
    ref0000001_0_1
        ref0000001_0_1.fasta
        ref0000001_0_1.qual
    """
    def __init__(self, args):
        self.parseArgs(args)
        self.__initLog()

    def __initLog(self):
        """Logging"""
        logLevel = logging.DEBUG if self.options.debug else logging.INFO
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
        logging.info("Running %s" % " ".join(sys.argv) )

    def parseArgs(self, args):
        parser = OptionParser(USAGE)
        parser.add_option("--debug",action="store_true",default=False)
        self.options, args = parser.parse_args(args)
        if len(args) < 5:
            parser.error("Invalid Number of Arguments.")
        self.fasta = args[0]
        self.qual = args[1]
        self.gapFile = args[2]
        self.inDir = args[3]
        self.outDir = args[4]
        self.jobDirs = args[5:] 
        
    def openGapCans(self):
        """
        Opens and consolidateds GapCans
        """
        gapSup = defaultdict(GapCans)
        #{gapSupportedBy: {gap: {read:flag}}}
        for file in glob.glob(os.path.join(self.inDir,"*.gapCans")):
            fh = open(file,'r')
            gapCans = json.load(fh)
            for gap in gapCans.keys():
                gapSup[gap].extend(gapCans[gap])
            fh.close()
        self.gapSup = dict(gapSup)
    
    def extractReads(self):
        for gapCan in self.gapSup.keys():
            logging.debug("Consolidating Gap %s" % gapCan)
            gapDir = os.path.join(self.outDir,gapCan)
            try:
                os.mkdir(gapDir)
            except OSError:
                logging.error("%s already has gap file. Skipping." % gapDir)
                continue
            
            fout = open(os.path.join(gapDir,"support.gapCans"),'w')
            fout.write(json.dumps(self.gapSup[gapCan], default=str, indent=4))
            fout.close()
            
            fastOut = open(os.path.join(gapDir,"input.fasta"),'w')
            qualOut = open(os.path.join(gapDir,"input.qual"),'w')
                
            for hit in self.gapSup[gapCan].getAllReads():
                fastOut.write(">" + hit + "\n" + \
                              self.allFasta[hit] + "\n")
                qualOut.write(">" + hit + "\n" + \
                              self.allQual[hit] + "\n")
            
            gap = self.gapInfo[gapCan]
            
            leftN = gap.leftContig
            fastOut.write(">" + leftN + "\n" + \
                          self.fastaRef[leftN][-FLANKAMT:] + "\n")
            qualOut.write(">" + leftN + "\n" + \
                          " ".join( map( str, \
                                    self.qualRef[leftN][-FLANKAMT:])) + "\n")
            
            rightN = gap.rightContig
            fastOut.write(">" + rightN + "\n" + \
                          self.fastaRef[rightN][:FLANKAMT] + "\n")
            qualOut.write(">" + rightN + "\n" + \
                          " ".join( map( str, \
                                    self.qualRef[rightN][:FLANKAMT])) + "\n")
            #Make quick assembly reference - used during remapping
            lContig = self.fastaRef[leftN][-REFAMT:]
            rContig = self.fastaRef[rightN][:REFAMT]
            
            contigRef = open(os.path.join(gapDir,"contigRef.fasta"),'w')
            contigRef.write(">"+gap.scaffold+"\n" + \
                            lContig + ("N"*gap.length) + rContig + "\n")
            contigRef.close()
            delta = gap.start - len(lContig)
            fout = open(os.path.join(gapDir,"delta"),'w')
            fout.write(str(delta))
            fout.close()
            
            fastOut.close()
            qualOut.close()

    def loadFastas(self):
        self.allFasta = cacheIndex()
        self.allQual = cacheIndex()
        for input in self.jobDirs:
            inQualName = input[:input.rindex('.fasta')]+".qual"
            self.allFasta.buildIndex(input)
            self.allQual.buildIndex(inQualName)


    def loadReferences(self):
        self.fastaRef = FastaFile(self.fasta)
        self.cleanReference()
        self.qualRef = QualFile(self.qual)
        self.gapInfo = GapInfoFile(self.gapFile)
    
    def cleanReference(self):
        """
        Remove IUB characters from the reference since
        we use contigs as part of the input reference for
        the local assembly and blasr doesn't like them
        """
        for entry in self.fastaRef.keys()
            self.fastaRef[entry] = self.fastaRef[entry].replace('M','C')
            self.fastaRef[entry] = self.fastaRef[entry].replace('R','A')
            self.fastaRef[entry] = self.fastaRef[entry].replace('W','T')
            self.fastaRef[entry] = self.fastaRef[entry].replace('S','G')
            self.fastaRef[entry] = self.fastaRef[entry].replace('Y','C')
            self.fastaRef[entry] = self.fastaRef[entry].replace('K','T')
            self.fastaRef[entry] = self.fastaRef[entry].replace('V','G')
            self.fastaRef[entry] = self.fastaRef[entry].replace('H','A')
            self.fastaRef[entry] = self.fastaRef[entry].replace('N','')

    def run(self):
        """
        Opens gaps. Loads Fastas. Extracts Reads.
        """
        logging.info("Opening Gap Objects")
        self.openGapCans()
        logging.info("Loading Fasta Files")
        self.loadFastas()
        logging.info("Loading References")
        self.loadReferences()
        logging.info("Extracting Reads")
        self.extractReads()
        logging.info("Finished Extraction")
        
class cacheIndex():
    def __init__(self):
        self.indices = {}
        
    def __getitem__(self, key):
        fn, start, numLines = self.indices[key]
        
        ret = ""
        for i in xrange(numLines):
            ret += linecache.getline(fn, i+start).strip()
        return ret

    def buildIndex(self, fn):
        fh = open(fn)
        fh.seek(0)
        prevHeader = None
        numLines = 0
        startLine = 0
        for line in fh.readlines():
            startLine += 1
            if line[0] == ">":
                header = line[1:].strip()
                self.indices[header] = [fn, startLine+1, None]
                if prevHeader != None:
                    self.indices[prevHeader][2] = numLines
                    numLines = 0
                prevHeader = header
            else:
                numLines += 1
        self.indices[prevHeader][2] = numLines
        fh.close()

if __name__ == '__main__':
    me = Extraction(sys.argv[1:])
    me.run()    

