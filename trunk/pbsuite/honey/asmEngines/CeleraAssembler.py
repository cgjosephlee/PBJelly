
import tempfile

from pbsuite.honey.asmEngines.Assembler import *

class CeleraAssembler(Assembler):
    
    def __init__(self, data, buffer, tmpDir, timeout, *args, **kwargs):
        Assembler.__init__(self, data, buffer, tmpDir, timeout)
            
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
        #want to make a directory for temp files?
        
        """
        Celera procedure (current version does NOT use specs.. might be a bad idea)
        
        runCA -d directory -p prefix -s specfile <option=value> ... <input-files> ...

        input-files need to be output into FRG
        """
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
        super(CeleraAssembler, self).cleanupTmp()
        
        return self.results

    def __call__(self, nBams, tBams):
        #Fetch,
        logging.info("asm task groupid=%s start" % (self.data.name))
        reads = []
        chrom = self.data.chrom
        start = max(0, self.data.start - self.buffer)
        #hope that fetching beyond 3' boundary is okay
        end = self.data.end + self.buffer
        
        for bam in nBams:
            #I actually need to convert these into a FRG file..
            #make temp
            for name, seq, qual, in super(PhrapAssembler, self).fetchReads(bam, chrom, start, end):
                #write read to illumina file
                pass
        
        for bam in tBams:
            #I actually need to convert these into a FRG..
            #And to I need to trim anymore?
            for name, seq, qual, in super(PhrapAssembler, self).fetchReads(bam, chrom, start, end, trim=True):
                #write to pacbio file
                pass
        
        if totReads > self.maxreads:
            return "Failure - Too Many Reads (%d) %s" % (totReads, self.data.name)
        #fastqToCA -insertsize {muins} {stdins} -librayname ill -technology illumina -nonrandom
        #fastqToCA -librayname pac -technology pacbio-raw -nonrandom
        #Assemble
        logging.info("assembling %d reads" % (len(reads)))
        self.result = self.__assemble(illFn, pacFn)
        
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
 
