import tempfile
import logging
from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import *

"""
bamToFastq.py ~/english/StructuralVariation/CHM1/reads/CHM1.bam 1:72766317-72811839 > pacbio.fastq 
bamToFastq.py ~/english/StructuralVariation/CHM1/reads/CHM1.bam 1:72765317-72812839 > pacbio.fastq 
bamToFastq.py ~/english/StructuralVariation/CHM1/IlluminaReads/AllReads_300_350.bam 1:72765317-72812839 >
  illumina.fastq
bamToFastq.py ~/english/StructuralVariation/CHM1/IlluminaReads/AllReads_350_400.bam 1:72765317-72812839 >>
   illumina.fastq
bamToFastq.py ~/english/StructuralVariation/CHM1/IlluminaReads/AllReads_3kb.bam 1:72765317-72812839 >>
    illumina.fastq
bamToFastq.py ~/english/StructuralVariation/CHM1/IlluminaReads/AllReads_450-500.bam 1:72765317-72812839 >>
     illumina.fastq
./../SparseAssembler LD 0 NodeCovTh 1 EdgeCovTh 0 k 31 g 15 PathCovTh 100 GS 12000000 f illumina.fastq 
mkdir hybrid
./../../DBG2OLC k 17 KmerCovTh 2 MinOverlap 20 AdaptiveTh 0.002 LD1 0 Contigs ../Contigs.txt RemoveChimera 1
     f ../pacbio.fastq 
"""

class dbg2ovlAssembler(Assembler):
    
    def __init__(self, data, args):
        #buffer, tmpDir, timeout, *args, **kwargs):
        Assembler.__init__(self, data, args)
    
    def __assemble(self, sreads, lreads):
        """
        
        assembles
        reads results
        clears temp files
        returns results as a string
        Calls the assembler
        """
        self.myTmpFiles = []
        #Temporary Files
        sout = tempfile.NamedTemporaryFile(suffix=".fastq", mode="w", delete=False, dir=self.tmpDir)
        logging.debug(sout.name)
        self.myTmpFiles.append(sout.name)
        for name, seq, qual in reads:
            sout.write("@{0}\n{1}\n+\n{2}\n".format(name, seq, qual))
        sout.close()
        
        lout = tempfile.NamedTemporaryFile(suffix=".fastq", mode="w", delete=False, dir=self.tmpDir)
        logging.debug(fout.name)
        self.myTmpFiles.append(lout.name)
        for name, seq, qual in reads:
            lout.write("@{0}\n{1}\n+\n{2}\n".format(name, seq, qual))
        lout.close()

        r, o, e = exe("SparseAssembler LD 0 NodeCovTh 1 EdgeCovTh 0 k 31 g 15 " \
                      "PathCovTh 100 GS 12000000 f " + sout.name)
        r, o, e = exe("DBG2OLC k 17 KmerCovTh 2 MinOverlap 20 AdaptiveTh 0.002 "\
                      "LD1 0 Contigs ../Contigs.txt RemoveChimera 1
        r, o, e = exe("phrap %s -minmatch 6 -minscore 20" % (fout.name),\
                      timeout=self.timeout)
        self.myTmpFiles.extend([fout.name + ".contigs",  fout.name + ".contigs.qual", \
                           fout.name + ".problems", fout.name + ".problems.qual", \
                           fout.name + ".log",      fout.name + ".singlets"])
        if r == 214:
            super(PhrapAssembler, self).cleanupTmp()
            return "Failure - Assembly Timeout " + self.data.name
         
        results = mergeFastaQual(fout.name + ".contigs", fout.name + ".contigs.qual")
        
        #Try to push the problems through, too
        if os.stat(fout.name + '.problems').st_size != 0:
            pfile = fout.name + ".problems"
            r, o, e = exe("phrap %s -minmatch 6 -minscore 20" % (pfile), \
                        timeout=self.timeout)
                    
            self.myTmpFiles.extend([pfile + ".contigs",  pfile + ".contigs.qual", \
                            pfile + ".problems", pfile + ".problems.qual", \
                            pfile + ".log",      pfile + ".singlets"])
            if r == 214:
                super(PhrapAssembler, self).cleanupTmp()
                return "Failure - Assembly Timeout " + self.data.name
            
            results.update(mergeFastaQual(fout.name + ".problems.contigs", fout.name + ".problems.contigs.qual"))
        
        #save to file
        fout = tempfile.NamedTemporaryFile(prefix = "asm" + self.data.name, mode="w",\
                                    suffix=".fastq", delete=False, dir=self.tmpDir)
        for key in results:
            fout.write("@group" + self.data.name + "_" + key + "\n" + \
                        results[key].seq + '\n+\n' + \
                        results[key].qual + '\n')
        fout.close()
        self.results = fout.name
        
        #clean up
        super(PhrapAssembler, self).cleanupTmp()
        
        return self.results
    
    def __call__(self, nBams, tBams):
        #Fetch,
        logging.info("asm task groupid=%s start" % (self.data.name))
        sreads = {}
        lreads = {}
        chrom = self.data.chrom
        start = self.data.start - self.buffer
        start = max(0, start)
        end = self.data.end + self.buffer
        
        for bam in nBams:
            if self.data.start + self.buffer >= self.data.end - self.buffer:
                sreads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, start - self.buffer, end + self.buffer))
            else:
                sreads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.start - self.buffer), self.data.start + self.buffer))
                sreads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.end - self.buffer), self.data.end + self.buffer))
                
        if len(reads) > self.args.maxreads:
            logging.info("Downsampling %s" % (self.data.name))
            nreads = {}
            for i in random.sample(reads.keys(), self.args.maxreads):
                nreads[i] = reads[i]
            reads = nreads
        
        for bam in tBams:
            if self.data.start + self.buffer >= self.data.end - self.buffer:
                logging.critical("SmallSpan Problem")
                lreads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, start - self.buffer, end + self.buffer))
            else:
                lreads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                        max(0, self.data.start - self.buffer), self.data.start + \
                        self.buffer, trim=True))
                lreads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.end - self.buffer), self.data.end +
                             self.buffer, trim=True))
             
        sreads = reads.values() 
        lreads = reads.values()
        totReads = len(sreads) + len(lreads)
        
        #Assemble
        logging.info("assembling %d reads" % (totReads))
        self.result = self.__assemble(sreads, lreads)
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
 
