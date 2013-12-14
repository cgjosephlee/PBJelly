Honey Documentation

== CONTENTS == 

I.   Using This README
II.  Requirements
III. Installation
IV.  Quick Start
V.   Details
VI.  Output Formats
VII. Extras

== I. Using This README ==

 Commands: 
    All commands are presented in the format
    > commandToExecute 
    where commandToExecute would be the actual command. 

== II. Requirements ==
    
 *  samtools    0.1.17 			*
 *  blasr       1.3.1.127046		
 *  python      2.7
 *  h5py        2.0.1
 *  pysam       0.7.4			*
 *  numpy       1.6

 Note: If you have PacBio's SMRTAnalysis suite v2.1, all of 
 these requirements will be met.

== III. Installation ==

 1) Edit setup.sh and change $SWEETPATH to the full directory where
    you've placed the package.
 
 2) To automatically place the package into your environment, add
    > source <path to>/setup.sh
    to your .bash_profile

  Be sure to source your .bash_profile (or just setup.sh) before 
  using Honey

== IV. Quick Start == 
 
 For more details on each step in the pipeline, see Section V
 below. 
 
 Beginning with a mapped bam, you'll use Honey.py to execute the stages 
 'pie', 'tails', 'spots'.

 1) Honey.py pie
    Extract the soft-clipped tails and attempt to remap them
 2) Honey.py tails
    Cluster the tail-mapping information to make genomic breakpoints
 3) Honey.py spots
    Look for genomic variants within the span of reads. 
    REQUIRED : The MD Tag MUST be present in the bam

== V. Details == 
 
 Below you'll find the full description of the procedure to create
 genomic variant calls from an input fastq
 
 1) Map your filtered subreads.
    Use blasr to map your filtered_subreads.fastq to your reference
    genome. At minimum, your command must have the following arguments:

    > blasr reads.fastq ref.fasta -bestn 1 -noSplitSubreads -sam -clipping soft 
    
    You can edit the other defaults however you like. I recommend using:
    * -nCandidates 15 
    * -sdpTupleSize 6
    * -minPctIdentity 75
 
 2) Turn the sam alignment into a bam.
    Honey comes with a utility called 'sam2bam'. Use this to quickly and 
    easily turn your sam results into a sorted and indexed bam.
    >  sam2bam reference.fasta mapresults.sam

 3) Extract the tails.
    Some number of your reads will have unmapped, soft-clipped tails. Use
    'pie' to extract, map, and consolidate those tails with your 
    results.
    > Honey.py pie mapresults.bam reference.fasta
    See Honey.py pie --help for details
 
 4) Call Honey Tails
    > Honey.py tails final.bam
 
 5) Call the Md tag.
    This step is optional if you're only interested in finding variants from
    tails. The command is:
    > samtools calmd -b mapresults.tails.bam reference.fasta > final.bam
    > samtools index final.bam
    
 6) Call Honey Spots
    > Honey.py spots final.bam
    
== VI. Output Format ==

 Honey spots columns:
 
    * CHROM       Reference entry where variant occurs
    * OUTERSTART  The 5' most boundary estimate of where variant begins
    * START       Best guess as the the exact start of the variant
    * INNERSTART  The 3' most boundary estimate of where variant begins
    * INNEREND    The 5' most boundary estimate of where variant ends
    * End         Best guess as the the exact end of the variant
    * OUTEREND    The 3' most boundary estimate of where the variant ends
    * TYPE        Variant type. One of MIS, INS, INSZ, or DEL
    * SIZE        Estimate size of the variant
    * INFO        More information associated with the calls

 
 Honey tails columns:

    * ID        Unique identifier of the call
    * CHROM     Reference entry where variant occurs
    * UTAIL     Upstream alignment piece of the seed is Prolog, Initial, or 
                Epilog (pie)
    * UMAPQ     Average mapqv of upstream reads
    * UDIR      Upstream seed alignment's direction
    * UBREAK    Average position of upstream breaks
    * ISINV     - if upstream and downstream alignment hit same strand, % if 
                inverted
    * DIR       Direction of the join between upstream and downstream hits
    * DBREAK    Average position of downstream breaks
    * DDIR      Downstream seed alignment's direction
    * DMAPQ     Average mapqv of downstream reads
    * DTAIL     Downstream alignment piece of the seed is Prolog, Initial, or 
                Epilog (pie)
    * EXSEQ     Average amount of sequence left between pie
    * NUMREADS  Number of unique reads in the cluster
    * NUMZMWS   Number of unique ZMWs in the cluster
    * EVIDENCE  Colon separated list of piece (p,i,e), direction (5,3), isInv
                (%,-),direction (<-,->), piece, direction(5,3) for each read
                in cluster
    
== VII. Extras ==

Future Directions:
develop Comb.py -- Annotation (super specific break points) in Future--
I can do this if I identify good seeds or if I perform super polish
or use minimus or something -- maybe get dan's assembler

write the user documentation





