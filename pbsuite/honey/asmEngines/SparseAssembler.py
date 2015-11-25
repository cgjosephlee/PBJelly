class SparseAssembler(Assembler):
    
    def __init__(self, data, args):
        #buffer, tmpDir, timeout, threads, maxreads):
        super(SpadesAssembler, self).__init__(data, args)
        self.threads = args.nproc
    
    def __fetchPEReads(self, bam, chrom, start, end):
        def write_read(read):
            """
            Write read to open FASTQ file.
            """
            COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}

            if read.is_reverse:
                ret = [read.qname + ("/%d" % (int(not read.is_read1) + 1)),
                       "".join((COMPLEMENT[b] for b in read.seq[::-1])),
                       read.qual[::-1]]
            else:
                ret = [read.qname + ("/%d" % (int(not read.is_read1) + 1)),
                       read.seq,
                       read.qual]
            return ret
        
        logging.info("fetching %s from %s:%d-%d" % (bam.filename, chrom, start, end))
        name = read_left = read_right = None
        
        sync_pairs = False
        for read in bam.fetch(reference=chrom, start = start, end = end):
            if name is not None and read.qname != name:
                if read_left and (not sync_pairs or read_right):
                    self.leftReads.append(write_read(read_left))
                if read_right and (not sync_pairs or read_left):
                    self.rightReads.append(write_read(read_right))
                read_left = read_right = None
            name = read.qname
            if read.is_read1:
                read_left = read
            else:
                read_right = read
        
        if read_left and (not sync_pairs or read_right):
            self.leftReads.append(write_read(read_left))
        if read_right and (not sync_pairs or read_left):
            self.rightReads.append(write_read(read_right))
        
    def __assemble(self):
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
        fout = tempfile.NamedTemporaryFile(prefix="spades_pe1", suffix=".fastq", delete=False, dir=self.tmpDir, mode="w")
        self.myTmpFiles.append(fout.name)
        for name, seq, qual in self.leftReads:
            fout.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        fout.close()
        
        fout2 = tempfile.NamedTemporaryFile(prefix="spades_pe2", suffix=".fastq", delete=False, dir=self.tmpDir, mode="w")
        self.myTmpFiles.append(fout2.name)
        for name, seq, qual in self.rightReads:
            fout2.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        fout2.close()
        
        foutp = tempfile.NamedTemporaryFile(prefix="spades_pb", suffix=".fastq", delete=False, dir=self.tmpDir, mode="w")
        self.myTmpFiles.append(foutp.name)
        for name, seq, qual in self.pbReads:
            foutp.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        foutp.close()
        
        #working here
        resultOut = tempfile.mkdtemp(prefix="spades", dir=self.tmpDir)
        
        estSize = self.buffer * 2
        if self.data.rest[0] != 'DEL':
            estSize += int(self.data.rest[1])
        
        #r, o, e = exe("dipspades.py -1 {pe1} -2 {pe2} --pacbio {pacbio} -o {output} "\
        r, o, e = exe("spades.py -1 {pe1} -2 {pe2} --pacbio {pacbio} -o {output} "\
                      .format(pe1=fout.name, pe2=fout2.name, pacbio=foutp.name, output=resultOut), \
                      timeout=self.timeout)
                    
        logging.debug("RET - %d\nOUT - %s\nERR- %s" % (r, o, e))
        #just the output dir, maybe?
        self.myTmpFiles.append(resultOut)
        if r == 214:
            super(SpadesAssembler, self).cleanupTmp()
            return "Failure - Assembly Timeout " + self.data.name
         
        outFsta = os.path.join(resultOut, "dipspades", "consensus_contigs.fasta")
        fasta = FastaFile(outFsta)
        
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
        super(SpadesAssembler, self).cleanupTmp()
        
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
        
        self.leftReads = []
        self.rightReads = []
        self.pbReads = []
        
        for bam in nBams:
            self.__fetchPEReads(bam, chrom, start, end)
        
        for bam in tBams:
            self.pbReads.extend(super(SpadesAssembler, self).fetchReads(bam, chrom, start, end, trim=True))
            #self.pbReads.extend(Assembler.fetchReads(self, bam, chrom, start, end, trim=True))
        
        #Assemble
        totReads = len(self.leftReads) + len(self.rightReads) + len(self.pbReads)
        if totReads > self.args.maxreads:
            return "Failure - Too Many Reads (%d) %s" % (totReads, self.data.name)
        logging.info("assembling %d reads" % (totReads))
        self.result = self.__assemble()
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
 
