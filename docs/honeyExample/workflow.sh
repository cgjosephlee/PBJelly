#!/bin/bash

#Set the inputReads and reference -- change at your leisure
inputReads=filtered_subreads.fastq
reference=lambda_modified.fasta

echo "Mapping"
blasr $inputReads $reference -bestn 1 -nCandidates 15 -sdpTupleSize 6 \
      -minPctIdentity 75 -affineAlign -noSplitSubreads -nproc 1 -sam -clipping soft -out mapping.sam

echo "Extracting Tails"
Honey.py pie mapping.sam $reference

echo "Sam To Bam"
sam2bam $reference mapping.tails.sam

echo "Calling MD Tag"
samtools calmd -b mapping.tails.bam $reference > mappingFinal.bam
samtools index mappingFinal.bam

echo "Calling Tails"
Honey.py tails mappingFinal.bam

echo "Calling Spots"
Honey.py spots mappingFinal.bam
