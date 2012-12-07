#!/usr/bin/env python
USAGE = """
Consolidates the reads over a single gap into a folder.

Extraction.py <reference.fasta> <reference.qual> \
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
import sys, os, glob, logging, json, linecache, re
from optparse import OptionParser
from collections import defaultdict
from FileHandlers import FastaFile, QualFile, FastaIndex, GapCans, GapInfoFile

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
        if not os.path.exists(self.outDir):
            parser.error("%s does not exist!" % self.outDir)
        if not os.path.exists(self.inDir):
            parser.error("%s does not exist!" % self.inDir)
        if not os.path.exists(self.fasta):
            parser.error("%s does not exist!" % self.fasta)
        if not os.path.exists(self.qual):
            parser.error("%s does not exist!" % self.fasta)
        if not os.path.exists(self.gapFile):
            parser.error("%s does not exist!" % self.fasta)

    def openGapCans(self):
        """
        Opens and consolidateds GapCans
        """
        gapSup = defaultdict(GapCans)
        #{gap: [read, ...],...}
        needReads = {}
        #{readName:0, ...}
        for file in glob.glob(os.path.join(self.inDir,"*.gapCans")):
            fh = open(file,'r')
            gapCans = json.load(fh)
            for gap in gapCans.keys():
                gapSup[gap].extend(gapCans[gap])
                for read in gapCans[gap]:
                    needReads[read] = 0
            fh.close()
        
        self.gapSup = dict(gapSup)
        self.needReads = dict(needReads)
    
    def extractReads(self):
        getSub = re.compile("(.*)##(\d+)#(\d+)##$")
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
                if hit.endswith('#'):
                    name, start, end = getSub.match(hit).groups()
                    start = int(start)
                    end = int(end)
                else:
                    name = hit
                    start = 0
                    end = None
                """
                seq, qual = self.allSequence[name]
                qual = qual[start:end]
                fastOut.write(">" + hit + "\n"+ \
                                seq[start:end] + "\n")
                qualOut.write(">" + hit + "\n" + \
                                " ".join(map(lambda x: str(ord(x)-33), \
                                                    qual[start:end]+"\n")))
                """#FastaIndex
                fastOut.write(">" + hit + "\n" + \
                              self.allFasta.getFastaRead(name)[start:end] + "\n")
                qualOut.write(">" + hit + "\n" + \
                              " ".join(self.allQual.getQualRead(name)[start:end])+"\n")
                              #" ".join(map(str,self.allQual[name][start:end])) + "\n")
                #"""
            
            gap = self.gapInfo[gapCan]
            
            leftN = gap.leftContig
            leftStart = max(gap.start - FLANKAMT, 0)
            leftEnd = gap.start
            lContig = self.cleanSequence(self.fastaRef[gap.scaffold][leftStart:leftEnd])
            fastOut.write(">" + leftN + "\n" + lContig + "\n")
            qualOut.write(">" + leftN + "\n" + \
                          " ".join( map( str, \
                                    self.qualRef[gap.scaffold][leftStart:leftEnd])) + "\n")
            
            rightN = gap.rightContig
            rightStart = gap.end
            rightEnd = gap.end + FLANKAMT
            rContig = self.cleanSequence(self.fastaRef[gap.scaffold][rightStart:rightEnd])
            fastOut.write(">" + rightN + "\n" + rContig + "\n")
            qualOut.write(">" + rightN + "\n" + \
                          " ".join( map( str, \
                                    self.qualRef[gap.scaffold][rightStart:rightEnd])) + "\n")

            #Make quick assembly reference - used during remapping
            leftStart = max(gap.start - REFAMT, 0)
            leftEnd = gap.start
            lContig = self.fastaRef[gap.scaffold][leftStart:leftEnd]
            rightStart = gap.end
            rightEnd = gap.end + REFAMT
            rContig = self.fastaRef[gap.scaffold][rightStart:rightEnd]
            contigRef = open(os.path.join(gapDir,"contigRef.fasta"),'w')
            contigRef.write(">"+gap.scaffold+"\n" + \
                            lContig + ("N"*gap.length) + rContig + "\n")
            contigRef.close()
            qualOut.close()
            fastOut.close()
            
            delta = gap.start - len(lContig)
            fout = open(os.path.join(gapDir,"delta"),'w')
            fout.write(str(delta)+'\n')
            fout.write(str(len(self.fastaRef[gap.scaffold]))+'\n')
            fout.close()
            
    def loadFastas(self):
        self.allSequence = {}
        for input in self.jobDirs:
            inQualName = input[:input.rindex('.fasta')]+".qual"
            fasta = FastaFile(input)
            qual = QualFile(inQualName)
            for key in fasta:
                try:
                    x = self.needReads[key]
                    self.allSequence[key] = (fasta[key], \
                            "".join(map(lambda x: chr(x+33), qual[key])))
                except KeyError:
                    continue
   
    def loadFastas_Index(self):
        self.allFasta = FastaIndex()
        self.allQual = FastaIndex()
        for input in self.jobDirs:
            inQualName = input[:input.rindex('.fasta')]+".qual"
            self.allFasta.addFasta(input)
            self.allQual.addFasta(inQualName)
    
    def loadReferences(self):
        self.fastaRef = FastaFile(self.fasta)
        self.qualRef = QualFile(self.qual)
        self.gapInfo = GapInfoFile(self.gapFile)
    
    def cleanSequence(self,seq):
        """
        Remove IUB characters from the reference since
        we use contigs as part of the input reference for
        the local assembly and blasr doesn't like them
        """
        seq = seq.replace('M','C')
        seq = seq.replace('R','A')
        seq = seq.replace('W','T')
        seq = seq.replace('S','G')
        seq = seq.replace('Y','C')
        seq = seq.replace('K','T')
        seq = seq.replace('V','G')
        seq = seq.replace('H','A')
        seq = seq.replace('N','')
        return seq

    def run(self):
        """
        Opens gaps. Loads Fastas. Extracts Reads.
        """
        logging.info("Opening Gap Objects")
        self.openGapCans()
        logging.info("Loading Fasta Files")
        self.loadFastas_Index()
        logging.info("Loading References")
        self.loadReferences()
        logging.info("Extracting Reads")
        self.extractReads()
        logging.info("Finished Extraction")
        self.allFasta.close()
        self.allQual.close()

if __name__ == '__main__':
    me = Extraction(sys.argv[1:])
    me.run()    

