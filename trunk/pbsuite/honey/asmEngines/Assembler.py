import tempfile
import logging
 
class Assembler(object):
    
    def __init__(self, data, args):
        """
        Args is a namedtuple of every parameter needed
        required args are buffer, tmpDir, and timeout
        """
        #buffer, tmpDir, timeout):
        self.data = data
        self.args = args
        self.buffer = args.buffer
        self.tmpDir = args.temp
        self.timeout = args.timeout
    
    def __toQual(self, input):
        if type(input) == str:
            return " ".join([str(ord(x)-33) for x in input])
        elif type(input) == list and type(input[0]) == int:
            return "".join([chr(x+33) for x in input])
        raise TypeError("Expected string or list of ints for __toQual")
    
    def fetchReads(self, bam, chrom, start, end, trim=False):
        """
        Trim is for pacbio reads that have tails that were attempted to 
        be remapped -- helps reduce redundancy
        Also, I'm going to read tails that extend beyond my boundaries
        """
        logging.info("fetching %s from %s:%d-%d" % (bam.filename, chrom, start, end))
        ret = {}
        
        for id, read in enumerate(bam.fetch(reference=chrom, start=start, end=end)):
            name = read.qname 
            name += " DIRECTION: rev" if read.is_reverse else " DIRECTION: fwd"
            seq, qual = read.seq, read.qual 
            
            if trim and not read.is_unmapped:
                regTrim = 0
                
                trimS = None
                trimE = None
                if start > read.pos:
                    for queryPos, targetPos in read.aligned_pairs:
                        if trimS is None and targetPos is not None and targetPos >= start:
                            trimS = queryPos
                if end < read.aend:
                    for queryPos, targetPos in read.aligned_pairs[::-1]:
                        if trimE is None and targetPos is not None and targetPos <= end:
                            trimE = queryPos
                
                if trimS is not None:
                    upS = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
                    trimS = max(0, trimS) + upS
                else:
                    trimS = 0
                
                if trimE is not None:
                    dnS = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0
                    trimE = min(len(read.seq), trimE)  - dnS
                else:
                    trimE = len(read.seq)
                
                seq = seq[trimS:trimE]
                qual = qual[trimS:trimE]
            
            ret[name + seq[:10]] = [name, seq, self._toQual(qual)]
        logging.info("%d reads retreived" % len(ret))
        return ret
    
    def cleanupTmp(self):
        for i in self.myTmpFiles:
            if os.path.exists(i):
                if os.path.isfile(i):
                    os.remove(i)
                else:
                    shutil.rmtree(i)
    
    def __str__(self):
        return self.result


