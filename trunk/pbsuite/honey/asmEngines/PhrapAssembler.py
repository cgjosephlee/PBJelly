import tempfile
import logging
from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import *
from pbsuite.honey.asmEngines.Assembler import Assembler

class PhrapAssembler(Assembler):
    
    def __init__(self, data, args):
        #buffer, tmpDir, timeout, *args, **kwargs):
        Assembler.__init__(self, data, args)
    
    def __assemble(self, reads):
        """
        writes temp files
        assembles
        reads results
        clears temp files
        returns results as a string
        Calls the assembler
        """
        self.myTmpFiles = []
        #Temporary Files
        fout = tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False, dir=self.tmpDir)
        logging.critical(fout.name)
        self.myTmpFiles.append(fout.name)
        qout = open(fout.name + '.qual', 'w')
        self.myTmpFiles.append(fout.name + '.qual')
        
        for name, seq, qual in reads:
            fout.write(">{0}\n{1}\n".format(name, seq))
            qout.write(">{0}\n{1}\n".format(name, qual))
        
        fout.close()
        qout.close()
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
        reads = {}
        chrom = self.data.chrom
        start = self.data.start - self.buffer
        start = max(0, start)
        end = self.data.end + self.buffer
        
        for bam in nBams:
            if self.data.start + self.buffer >= self.data.end - self.buffer:
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, start - self.buffer, end + self.buffer))
            else:
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.start - self.buffer), self.data.start + self.buffer))
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
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
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, start - self.buffer, end + self.buffer))
            else:
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                        max(0, self.data.start - self.buffer), self.data.start + \
                        self.buffer, trim=True))
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.end - self.buffer), self.data.end +
                             self.buffer, trim=True))
             
        reads = reads.values() 
        totReads = len(reads)
        
        #Assemble
        logging.info("assembling %d reads" % (len(reads)))
        self.result = self.__assemble(reads)
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
 
