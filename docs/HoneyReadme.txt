required stuff
    
    samtools    0.1.17
    blasr       v1.3.1.121193
    python	2.7
    h5py        2.0.1
    pysam       0.7.4
    numpy       1.6

Process - 
1) build bam from input reads mapped on reference
	>>  blasr reads.fastq reference.fasta -bestn 1 -noSplitSubreads -sam -clipping soft 
   you can edit the defaults if you're like, but the above arguments are
   required.
   next run sam2bam
   
    >>  sam2bam remap.sam

2) Run calmd on your bam
    >> samtools calmd -b remap.bam reference.fasta > remap_md.bam
    >> samtools index remap_md.bam

3) Tails

4) Run Honey.py to identify hotspots
    >> Honey.py -h
    
5) Run TGraf.py to get the tails

-- todo --

unite 4 and 5 to create .svp file

and a .svp to .vcf creator

develop Comb.py -- Annotation (super specific break points) in Future--
   I can do this if I identify good seeds or if I perform super polish
   or use minimus or something -- maybe get dan's assembler

write the user documentation





