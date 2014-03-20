Honey Documentation

== CONTENTS == 

I.   Using This README
II.  Requirements
III. Installation
IV.  Quick Start
V.   Details
VI.  Output Formats
VII. Interpreting Results
VIII.Extras

== I. Using This README ==

 Toy Data: 
    Provided with this distribution of Honey is a toy example 
    inside of docs/honeyExample directory. Use this once you've 
    setup Honey to test that everything is working as expected 
    and to become familiar with the software. 

 Commands: 
    All commands are presented in the format
    > commandToExecute 
    where commandToExecute would be the actual command. 

== II. Requirements ==
    
 *  samtools    0.1.17 			
 *  blasr       1.3.1.127046		
 *  python      2.7
 *  h5py        2.0.1
 *  pysam       0.7.4			
 *  numpy       1.6

 Note: If you have PacBio's SMRTAnalysis suite v2.1, all of 
 these requirements will be met.

== III. Installation ==
 
 1) Ensure all of the above requirements are in your environment.
 2) Edit setup.sh and change $SWEETPATH to the full directory where
    you've placed the package.
 
 3) To automatically place the package into your environment, add
    > source <path to>/setup.sh
    to your .bash_profile

  Be sure to source your .bash_profile (or just setup.sh) before 
  using Honey

== IV. Quick Start == 
 
 Here's a quick start guide to the Honey specific steps. For more 
 details on each step in the entire process, see Section V
 below. Otherwise, check out the --help for available parameters.
 
 Beginning with a bam with mapped PacBio reads, you'll use Honey.py to
 execute the stages 'pie', 'tails', 'spots'.

 1) Honey.py pie
    Extract the soft-clipped tails and attempt to remap them
 2) Honey.py tails
    Cluster the tail-mapping information to make genomic breakpoints
 3) Honey.py spots
    Look for genomic variants within the span of reads. 
    REQUIRED : The MD Tag MUST be present in the bam

== V. Details == 
 
 Below you'll find the full description of the procedure to create
 genomic variant calls from an input fastq. This README will be using
 the toy data in honeyExample for specific file names. Feel free to
 follow along using the toy data in honeyExample. See
 honeyExample/workflow.sh for the full set of commands.
 
 1) Map your filtered subreads.
    Use blasr to map your filtered_subreads.fastq to your reference
    genome. At minimum, your command must have the following arguments:

    > blasr reads.fastq ref.fasta -bestn 1 -noSplitSubreads -sam -clipping soft 
    
    You can edit the other defaults however you like. I recommend using:
    * -nCandidates 15 
    * -sdpTupleSize 6
    * -minPctIdentity 75
    * -affineAlign

 2) Extract the tails.
    Some number of your reads will have unmapped, soft-clipped tails. Use
    'pie' to extract, map, and consolidate those tails with your 
    results.
    > Honey.py pie mapping.sam reference.fasta
    For details on available parameters, see:
    > Honey.py pie --help 
   
 3) Turn the sam alignment into a bam.
    Honey comes with a utility called 'sam2bam'. Use this to quickly and 
    easily turn your sam results into a sorted and indexed bam.
    >  sam2bam reference.fasta mapping.tails.sam

 4) Call the MD tag and index your .bam
    Note: This step is optional if you only want the tails information
    > samtools calmd -b mappingSort.bam reference.fasta > mappingFinal.bam
    > samtools index mappingFinal.bam
 
 5) Call Honey Tails
    > Honey.py tails final.bam
    For details on available parameters, see:
    > Honey.py tails --help 
    
 6) Call Honey Spots
    > Honey.py spots final.bam
    For details on available parameters, see:
    > Honey.py spots --help 
    
== VI. Output Format ==

 Honey spots:
    .spots output -- Your variant calls with the format
    * CHROM       Reference entry where variant occurs
    * OUTERSTART  The 5' most boundary estimate of where variant begins
    * START       Best guess as the the exact start of the variant
    * INNERSTART  The 3' most boundary estimate of where variant begins
    * INNEREND    The 5' most boundary estimate of where variant ends
    * END         Best guess as the the exact end of the variant
    * OUTEREND    The 3' most boundary estimate of where the variant ends
    * TYPE        Variant type. One of MIS (missmatch), INS (insertion), or
    		  DEl (deletion)
    * SIZE        Estimated size of the variant
    * INFO        More information associated with the calls
    
    .h5 Output -- Contains each base-alignment type and coverage for each
    position in your reference. Use this if you wish to recall spots
    with different parameters (see recallSpots.py --help)

 Honey tails:
    .tails output -- Your variant calls with the format
    * id        Unique identifier of the call
    * chrKey    Key for what two chromosomes the breakpoint pair hits
    * uRef      First breakpoint location's reference
    * uBreak    First breakpoint location's coordinate
    * uMapq     Average mapping quality of reads that support the
                first breakpoint
    * dRef      First breakpoint location's reference
    * dBreak    First breakpoint location's coordinate
    * dMapq     Average mapping quality of reads that support the
                first breakpoint
    * remainSeq Average amount of sequence left unmapped between breakpoints
    * annot     Honey's predicted structural variant type. Note this
                isn't 100% accurate. Insertion (INS), deletion (DEL), or
		translocation (TLOC) are possible values.
    * numReads  Number of reads supporting structural variant
    * numZMWs   Number of ZMWs producing reads that support structrual
                variant
    * evidence  Semi-colon delimited list of read orientations around
                 breakpoint. See Section VII for details
    	
    Honey tails can produce two extra files with data that contains
    infomation which may help validation and annotation of tail
    structural variants. The two flags to create these files are:
    1) --fastq 
       Inside of a .tgz file, write a .fastq and .bam for each cluster
       identified as clu<id>.[bam|fastq] where <id> is column 1 of
       the tgraf output. These two files contain the full sequence and
       alignment information of each initial read that supported the
       structral variant.
    2) --verboseFile
       Write a .verbose with full details of every reads' alignment
       pair that supported the structural variant
       Each cluster begins with an identifier line
       
       ##Cluster <id> - <chrKey> 
	
       followed by 1+ lines fore each read's alignment. The coloumns
       are:
       
       #uRef uBreak uMapq dRef dBreak dMapq remainSeq break annot readName
       
       These are the same as what's reported above except on a
       per-read basis and the read's name is repoted as the final
       column.
           
== VII. Interpreting Results ==
  
  This section contains a tutorial on how to interpreate the data as well as
  some extra detail and advanced techniques for using Honey! To get the most
  out of this section, be sure you've run through Section V and created all of
  the results for the Toy Data
  
  = Spots Results =

  If you've followed along with the toy example, you should have found
  the following variants in the mappingFinal.hon.spots file
  
  #CHROM	OUTERSTART	START	INNERSTART	INNEREND	END	OUTEREND	TYPE	SIZE	INFO
  lambda_NEB3011	29917	29950	29982	30019	30050	30081	INS	106	startSig=-0.011;Insz3rdQ=112.000;InszMean=106.000;startCov=18.634;InszMedian=106.000;endSig=0.010;InszCount=16.000;endCov=18.132;Insz1stQ=101.000
  lambda_NEB3011	34956	35000	35046	35157	35200	35243	DEL	200	startCov=8.533;endSig=0.359;startSig=-0.381;endCov=8.875
  
  Honey has found an insertion somewhere between 29950-30050 of about 106bp
  and a deletion of 200bp between 35000-35200. For every spot-discovered
  structural variant, the info field has the following fields:
       * startSig, endSig - The average strength of the signal over the
                            start/end region
       * startCov, endCov - The average read coverage over the start/end
                            region
  For INS events, we have a few extra fields that describe the distribution of
  the insertion errors' size that Honey uses to predict the insert size:
       * InszCount - The number of reads that have an insertion in the region
       * Insz* - Description of the total insertion size per read over all
                 reads in the region
	
  = Visualizing Spots =
  
  If you wanted to create a graph showing the raw-errors and the transformed
  signal of our insertion region, use spotViz.py
      
      > spotViz.py mappingFinal.hon.h5 lambda_NEB3011 28917 31081
  
  Note: You have to provide a little bit of a buffer around the structural
  variant to get a good visualization (here I used 1kb).
  Two .png files are created rates.png and signals.png showing the original
  error rates in the region and the transformed signal.
  
  = Manual Inspection =
  
  Additionally, you can go into the bam itself to look at the alignments
  around the deletion breakpoints.
      > samtools tview mappedFinal.bam lambda_modified.fasta
     
  Press the '/' key to pull up the "Jump to Location" dialogue. Then type in
  "lambda_NEB3011:35000" then 'enter'. Move the cursor a little to the left (press 'j')
  and you should see something like this:
  
     34981          34991        35001	   35011      35021	35031	35041     35051   
  *C*ATTA*G*TG*AG*TTGA*TTG*AG*CTT*GGAATCAG*GAAGCTACGTT*CAACTCGA*****CTTATAAGGCGG*TGCCAGATG
   . .... . .. .. .... ... .. ... ..KK.KK. .KK.K.KT... KKKK.K.K     K..K.KK..K.. ..KKK.K..
  *,*,,,,*,*,,g,,*,,,,*,*,*,,*,,,*,,
  *.*....*.*..*..*....*...*..*...*..******************************************************
  *.*....*.*..*..*....T...*..*...*********************************************************
  *.*....*.**.*..*....*...*..C..**********************************************************
  ***....G.*..*..*....*...*..*...*
  *,*,,,,*,****,,*,*,,*,,,*,,*,,,*,,
  *,c,,,,*,*,,*,,g,,,,*,,,*,,*,,,t,*****,,t,**********************************************
  *,*,,,,*,g,,*,,*,,,,*,,,g,,*,,,*,*******************************************************
  *.*....*.*..*..*....*.C.*..**..*..************.GT...GA...*.*.AAAAC...********.T.T......*
  C.*....*.*..*..*....*...*..C...*********************************************************
  
  The top row is your position, second row is your reference sequence, third
  row is your consensus sequence, followed by a number of reads' alignments.
  As you can see, at about 35000, several reads either stop mapping or begin
  placing many deletion errors in the alignment. Because of the high error
  rate and alignment ambiguities, the exact point at which to begin the
  deletion varies on a per-read basis. This is straight forward in our
  simulated reference, but in real structural variants, there's usually
  microhomology around the breakpoints that ambiguity. Since our breakpoint
  here is so well definined, Honey predictes the exact breakpoint perfectly at
  position 35000. 
  
  = Recall Spots =
  
  The most time consuming part of spots calling is counting the errors. That's
  why Honey.py comes with a script called recallSpots.py that allows you to
  change the parameters you with to use for calling structural variants. For
  fun, run the following command:
  
      > recallSpots.py mappingFinal.bam mappingFinal.hon.h5 -e 2 > recall1.txt

  Inside of recall1, you can see that lowering the threshold from 5 to 2 gave
  us a couple more spots. These are both bad MIS results that span a 
  huge region. They are caused by a false-positive start signal and
  end signal that don't have a mate and Honey automatically tries to
  join the two points into a single call. Future versions of Honey
  will do a better job handling these issues, but for now you should
  continue to be aware of the issue.
  
  = Tails Results =
  
  Check out your tails results in mappingFinal.hon.tails and you'll
  see a total of 5 structural variants called (ids 0-3).

  ID=0 is simply evidence of lambda's circular genome. This is
  annotated as an Insertion because technically you can insert the
  sequence from the end of the reference at the beginning of the
  reference and still correctly represent your sample's genome.
  (However you wouldn't remove the INS SV call)

  ID=1 is also evidence of lambda's circular genome. However! this
  data indicates that a portion of the lambda molecules don't contain
  the first appx 900bp of the reference around their orgin. This is
  corectly called as a separate event.

  ID=2 Shows evidence of a 140bp (remainSeq) INS.
  ID=3 is an inversion. The simulated inversion was placed from
       9000-12000bp in the reference.
  ID=4 Shows a 1003bp deletion between 20000-21003, which is close to
       the exact coordinates of 20000-21003. See Manual Inspection
       above to see why the breakpoints may have been off by 3bp.
  
  = Tails Evidence =
  
  The final column in the Tails output shows how each read broke apart
  across the uBreak and dBreak. Every read starts with an initiail [i]
  alignment and any softclipped tails from the 5' of [i] that map
  create prologs [p] and from the 3' of [i] create epilogs [e]. The
  way to read a piece of evidence is to trace the read. In our sample
  genome - the sequence is always 
  
  	->p=i->i=e->
   
  Tracing down the 'pie' shows you the structure of the variant when
  represented in the reference space. For example, if you look at Tail
  ID=1, you'll see
  
        e->=->i
	i<-=<-e
	i->=->p
	p<-=
  This says that if you move on the direct strand to the dBreak (->i)
  in the reference, you'll pick back up ath the uBreak and move
  downstream on the direct strand (e->). The equal sign (=) in the middle
  represents that the two pieces of the alignment are on the same
  strand. If the pieces mapped to the compliment strand, the arrows
  would be backwards (<-) and if the pieces mapped on different
  strands, the equal sign would become a percent sign (%).
  
  = Important Tail Parameters =
 
  The best way to learn about the tail parameters are just to mess
  with them. Rerun Honey tails
  
      > Honey.py tails mappingFinal.bam -B 1000 -f -v -o tail2
     
  You'll see that we now make 17 structural variants! Using a buffer
  of 1000bp meant that our two structural variants that showed
  evidence of the circular genome has pooled because the reads all had
  the same orientation around the breakpoints and the breakpoints were
  within the specified distance. This is an important consideration
  when trying to find variants that border large repeats. If it's
  possible for the read to interrupt it's mapping over a large area, a
  larger buffer will collect all of the reads that support the event.
  However! You may incorrectly pool the breakpoints of two separate
  structural variants that happen to be in the same region within the
  buffer's size distance.
  
  = Missed Adapters =
  
  For the other events created with the modified parameters, these are
  evidence of a missed adapter. The way to tell is if the two
  breakpoints are near each other (usually within ~50bp or ~100bp) and
  the sample sequence switches strands and begins mapping back over
  itself (->%<- or <-%->). One way to possibly filter these out is to
  increase the number of reads and number of ZMWs required to support
  an event before it's reported. The below command will filter all
  but one of our missed adapter events.
  
      > Honey.py tails mappingFinal.bam -B 1000 -f -v -o tail3 -b 6 -z 6
  
  These can be identified as missed adapters by then looking at the
  verbose output provied in tail3.verbose. As you can see for the
  reads that support tail ID=3, they all look independent missed
  adapters. The problem is that they all occured close enough that
  they fell with the buffer. 

== VIII. Extras ==
  Change your calls into a .bed by using spotToBed.py and
  tailToBed.py

== IX. FAQ ==
  
  * Who can I report bugs to or ask questions?
    Please report your issues to the sourceforge ticketing system.

