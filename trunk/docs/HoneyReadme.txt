required stuff

    ?? means smtanalysis 2.0
    samtools ??
    blasr ??
    python2.7
    h5py ?? 
    pysam 0.7.4
    numpy ?? 

1) build bam from input reads mapped on reference
	>>  blasr reads.fastq reference.fasta -bestn 1 -noSplitSubreads -sam -clipping soft 
   you can edit the defaults if you're like, but the above arguments are
   required.
   next run sam2bam
   
    >>  sam2bam remap.sam

2) Run calmd on your bam
    >> samtools calmd -b remap.bam reference.fasta > remap_md.bam
    >> samtools index remap_md.bam

3) [Optional] Realign your bam 
    >> realign.py remap_md.bam
    # and calmd again (because I suck at recreating the MD tag)
    >> samtools calmd -b remap_md_realign.bam > remap_md.bam
    >> samtools index remap_md.bam

3) Run Honey.py to identify hotspots
    >> Honey.py -h
    

-- this is the annotator -- Comb.py
4) run Valid.py to assemble the region -- rename to Comb.py

   I can do this if I identify good seeds or if I perform super polish
  I can probably just use pbjpolish since I've already done all of this 
  realignment crap 
  Note! I should maybe use pbjPolish now! YAY!

