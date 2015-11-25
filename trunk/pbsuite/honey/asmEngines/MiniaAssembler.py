import tempfile
import logging
from pbsuite.honey.asmEngines.Assembler import Assembler

class MiniaAssembler(Assembler):
    
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
        fout = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, dir=self.tmpDir, mode="w")
        self.myTmpFiles.append(fout.name)
        
        for name, seq, qual in reads:
            fout.write(">%s\n%s\n" % (name, seq))
        
        fout.close()
        resultOut = tempfile.NamedTemporaryFile(prefix="minia_", delete=False, dir=self.tmpDir, mode="w")
        
        estSize = self.buffer * 2
        if self.data.rest[0] != 'DEL':
            estSize += int(self.data.rest[1])
        
        r, o, e = exe("minia {reads} 20 3 {estSize} {outPrefix}".format(reads=fout.name, \
                      estSize=estSize, outPrefix=resultOut.name ), \
                      timeout=self.timeout)
        
        logging.debug("RET - %d\nOUT - %s\nERR- %s" % (r, o, e))
        
        self.myTmpFiles.extend([resultOut.name + ".contigs.fa",  resultOut.name + ".debloom", \
                           resultOut.name + ".debloom2", resultOut.name + ".false_positive_kmers", \
                           resultOut.name + ".reads_binary",      resultOut.name + ".solid_kmers_binary"])
        
        if r == 214:
            super(MiniaAssembler, self).cleanupTmp()
            return "Failure - Assembly Timeout " + self.data.name
         
        fasta = FastaFile(resultOut.name + ".contigs.fa")
        
        results = {}
        for key in fasta:
            results[key] = FastqEntry(key, fasta[key], '?' * len(fasta[key]))
        
        #save to file
        fout = tempfile.NamedTemporaryFile(prefix = "asm" + self.data.name, mode="w", \
                                    suffix=".fastq", delete=False, dir=self.tmpDir)
        for key in results:
            fout.write("@group" + self.data.name + "_" + key + "\n" + \
                        results[key].seq + '\n+\n' + \
                        results[key].qual + '\n')
            
        fout.close()
        self.results = fout.name
        
        #clean up
        super(MiniaAssembler, self).cleanupTmp()
        
        return self.results

    def __call__(self, nBams, tBams):
        #Fetch,
        logging.info("asm task groupid=%s start" % (self.data.name))
        reads = []
        chrom = self.data.chrom
        start = self.data.start - self.buffer
        start = max(0, start)
        #hope that fetching beyonde 3' boundary is okay
        end = self.data.end + self.buffer
        
        for bam in nBams:
            reads.extend(self.fetchReads(bam, chrom, start, end))
        if len(tBams) > 1:
            logging.warning("Minia isn't built to handle PacBio Reads.. results may be unoptimal")
        for bam in tBams:
            reads.extend(super(MiniaAssembler, self).fetchReads(bam, chrom, start, end, trim=True))
        
        if len(reads) > self.maxreads:
            return "Failure - Too Many Reads (%d) %s" % (len(reads), self.data.name)
        
        #Assemble
        logging.info("assembling %d reads" % (len(reads)))
        self.result = self.__assemble(reads)
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
 
