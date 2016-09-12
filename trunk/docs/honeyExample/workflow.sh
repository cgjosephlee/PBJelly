#!/bin/bash

#Set the inputReads and reference -- change at your leisure
inputReads=filtered_subreads.fastq
reference=lambda_modified.fasta

echo "PIEMapping"
Honey.py pie $inputReads $reference -o mappingFinal.sam

echo "Sam To Bam"
sam2bam $reference mappingFinal.sam

echo "Calling Tails"
Honey.py tails mappingFinal.bam 

echo "Calling Spots"
Honey.py spots mappingFinal.bam --reference $reference

