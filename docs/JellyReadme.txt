	    PBJelly Documentation -- vTBD             

+-----------------+
#Using This README#
+-----------------+

 - Toy Data: 
 Provided with this distribution of PBJelly is a toy example 
 inside of the lambdaTest directory. Use this once you've 
 setup Jelly to test that everything is working as expected 
 and to become familiar with the software. 

 - Commands: 
 All commands are presented in the format

 >  [cmd]$: commandToExecute 

 where [cmd]$: is the tag that identifies a unix command and 
 commandToExecute would be the actual command. 

 - Alerts, Notes, and Protips 
 Throughout the document you will see tags for extra pieces 
 of information.
    + Alerts contain information that is crucial to running 
      PBJelly. 
    + Notes are extra tid-bits of knowledge about PBJelly 
      that may help you run PBJelly. 
    + Protips are considerations that one may keep in mind 
      to take full advantage of Jelly. 

+------------+
#Installation#
+------------+

 - Step 1 
 Download and Install the SMRTAnalysis software
	
 >  https://github.com/PacificBiosciences/SMRT-Analysis

 SMRTAnalysis v1.3.3 - v2.0.0 have all been vetted and should
 work with PBJelly. Use whichever distribution fits with you 
 system the best. 

 - Step 2 
 Edit exportPaths.sh and set two environment variables.
 >  set JELLYPATH to the path where you've placed PBJelly. 
 >  set SEYMOUR_HOME to the path where you've installed 
    SMRTAnalysis. 
 + Note: If SEYMOUR_HOME is already set, PBJelly assumes you 
   already sourced the SMRTAnalysis setup script. Be sure 
   this is the case. 

 - Step 3 
 To your .bash_profile , add the command.

 >  [cmd]$: source <path to PBJelly>/exportPaths.sh 

 Be sure to source your .bash_profile before beginning to 
 use PBJelly.
 + Alert: PBJelly runs using the python distributed with 
   SMRTAnalysis. Be sure that this python is your 
   environment's default after sourcing the .bash_profile 

+-----------+
#Quick Start#
+-----------+

 For more details on each step in the pipeline, see 
 "Running PBJelly" below. If, however, you'd like to just 
 run the program do the following.

 - Step 1 
 Create your Protocol.xml. To run the lambdaExample dataset 
 provided, edit the paths in the <reference> , <outputDir> 
 and the baseDir attribute in the inputs tag to the full 
 path that lambdaExample is sitting in. See 
 "Running PBJelly - Setup - Step 2" for more information 
 about creating Protocols. 
 
 - Step 2 
 Run each stage sequentially (meaning finish one stage 
 before continuing to the next). To run a stage, use the 
 command
 
 >  [cmd]$: Jelly.py <stage> yourProtocol.xml 

 The stages, in order, are

 >  setup 
 >  mapping 
 >  support 
 >  extraction 
 >  assembly 
 >  output 

 To get help with PBJelly, simply run Jelly.py --help.
 Or, for help with any stage, simply run 
 Jelly.py <stage> --help 

+---------------+
#Running PBJelly#
+---------------+

 - Pre Processing 
 If you would like to get some initial stats about your 
 reference (sizes and distributions of scaffolds, contigs, 
 and gaps) run: 

 >  [cmd]$: summarizeAssembly.py <reference.fasta> 
 >  [cmd]$: readSummary.py <Protocol.xml>
 
 see summarizeAssembly.py --help and readSummary.py for 
 details.
 
Setup Your Protocol
-------------------

 This is by far the longest and most involved step. Once you 
 get past this, PBJelly makes the rest of the workflow super 
 easy.
 
 - Step 1
 Gather the paths to your reference and all of your input 
 sequence files.
 
 + Alert: Every sequence file you use has the following 
   requirements:
 >  You should use the filtered_subreads produced by PacBio 
    SMRTAnalysis' protocols where SMRTBell adapters have 
    been removed. 
 >  Each set of filtered_subreads can be in a .fastq file or 
    in a .fasta file (no .fa, .fsta, etc extension.) that 
    has an associated qual file with the same name sitting 
    in the same directory beside the .fasta file. Qual files 
    should contain the Phred Scores of bases (0-93) and 
    should not be encoded (i.e. no Sanger/Solexa, only the 
    number for the score) If you do not have quals for your 
    sequences, see 
    [cmd]$: fakeQuals.py --help 
 > Each set of filtered_subreads need to have a unique file
   name. (e.g. filtered_subreads1.fastq, filtered_subreads2.fastq)
 > Sequence names should not have spaces.
 
 + ProTip: If you have a small number of very large sequence
   files and you want to speed up processing, split those 
   into several smaller files. PBJelly will submit one 
   mapping/support job per sequence file.
 
 - Step 2
 Create your Protocol.xml -- See TemplateProtocol.xml for an 
 idea of what a protocol should look like. You can name your 
 protocol whatever you'd like. Below are the elements needed 
 for a Protocol.

 <reference> : The reference tag contains the full path to 
 your reference.fasta. 
 + Note: All of the files PBJelly creates regarding your 
   reference will be placed in the same directory beside 
   the reference.fasta 
 + Note: You no longer need a .qual associated with your 
   reference! Yay! But if you do have one, just store it
   beside of your reference.fasta with the same base name
   and have it end with .qual

 <outputDir> : The output tag contains the full path to 
 where PBJelly will put the intermediate data for each stage 
 in the process. Protip: Placing your Protocol.xml into the 
 outputDir directory makes for excellent bookkeeping. Just 
 copy the provided TemplateProtocol.xml into the outputDir
 directory for each project you maintain 
 
 <cluster> [optional]: If you do not wish to submit jobs to 
 a cluster, do not include the <cluster> element in your 
 protocol. Otherwise, the cluster tag contains information 
 that PBJelly will use to submit jobs to your cluster. The 
 cluster tag contains two elements.

 <command> : This is a template that holds the structure of 
 how one submits jobs to his/her cluster. The example 
 provided in TemplateProtocol.xml is used on a MOAB job 
 management system. The command templatees has 4 REQUIRED 
 elements: 
 >  CMD     - The command one uses to execute on the cluster. 
 >  JOBNAME - The name to assign to the job being submitted. 
 >  STDOUT - The file that standard out will be directed to. 
 >  STDERR - The file that standard error will be directed to.
 
 + PROTIP: If you have a single, large machine and not a cluster,
   you can use the following command to submit all of your jobs
   to the background and parrallize operations. Just be careful
   about the number of jobs/resources you execute or you can 
   freeze your system
 >  ${CMD} ${JOBNAME} 2> ${STDERR} 1> ${STDOUT} &amp;

 <nJobs> : This is the maximum number of jobs PBJelly can 
 submit to the cluster. PBJelly tries to submit as many jobs 
 as possible. So, if you do not specify nJobs, PBJelly will 
 create a job for every single input file for any given stage. 
 For example, if you have 50 files containing read sequences 
 for mapping, and you specify nJobs = 5, PBJelly will submit 
 5 jobs, each will map 10 input files. If you do not specify 
 nJobs (i.e. don't include the nJobs tag or set nJobs to 0) 
 PBJelly will submit 50 jobs for mapping. 

 <blasr> : Place all of the mapping parameter to be sent to 
 blasr here. You can be comfortable with the default blasr 
 parameters that are provied in the TemplateProtocol.xml. 
 However, if you would like to customize the parameters, be 
 sure to read the blasr documentation (blasr --help). I also 
 highly recommend watching the video "Mapping and Alignment 
 using BLASR" located at http://www.smrtcommunity.com/Learn/Videos 
 
 + Alert : Always specify -noSplitSubreads in your blasr parameters. 

 <input> : Input contains information about where your input 
 data is located. First, there is the optional 'baseDir' 
 attribute. If all of your data has a common root path, 
 specifying baseDir=/That/Path will prevent redundancy in 
 the inputs. Inside of the <input> tag is the <job> tag. 
 This contains the path to each input file PBJelly will use. 

Setup Your Files
----------------
 Simply execute the command
 
 >  [cmd]$: Jelly.py setup Protocol.xm l 
 
 + Note : All of the files PBJelly creates regarding your 
 reference will be placed in the same directory beside the 
 reference.fasta 

Mapping Your Data
-----------------
 Run the mapping stage via the command

 > [cmd]$: Jelly.py mapping Protocol.xml 

 + Alert : Remember to wait until given stage is finished 
   before running the next stage.

 + Note: Standard Error and Standard Logs for each step are 
   placed next to the data PBJelly creates. For the mapping 
   step, this can be found in <outputDir>/mapping/ 

Support The Gaps
----------------
 Run the support stage 

 >  [cmd]$: Jelly.py support Protocol.xml 

Extract Useful Reads
--------------------
 Run the extraction stage

 >  [cmd]$: Jelly.py extraction Protocol.xml 

Assemble The Gaps
-----------------
 Run the assembly stage

 >  [cmd]$: Jelly.py assembly Protocol.xml 

 + Protip: If you have access to more than one core per gap 
 to be assembled, be sure to tell PBJelly to pass the nproc 
 parameter to the assembly stage via:
 >  [cmd]$: Jelly.py assembly Protocol.xml -x "--nproc=8"
 Where 8 can be replaced by the number of cores you're using.

Output Your Results
-------------------
 Run the output stage

 >  [cmd]$: Jelly.py output Protocol.xml 

 At the head of your log file, you can find information 
 about how many gaps were addressed, filled, etc. The output 
 stage collects all of your results into 3 files:

 <outputDir>/jellyOutput.fasta 
 -- The new reference 
 <outputDir>/jellyOutput.qual 
 -- the new reference qual. 
 <outputDir>/jellyLiftOverTable.json 
 -- For annotation lifting. 

 LiftOverTables are currently depricated.

+---------------+
#Post-Processing#
+---------------+
 removeOverFillGaps.py 
   This script will undo any of the changes PBJelly made in 
   overfilled gaps. For more information see 
   >  [cmd]$: removeOverFillGaps.py --help 

 summarizeAssembly.py 
   Use this script to get metrics describing the sizes and 
   distribution of your new reference's scaffolding, contigs,
   and gaps. Compare this to the original reference's summary
   to see the scale of the new reference's improvements.

+------+
#Extras#
+------+
 runBlasr.py: 
   If you have a large input.fasta file that you need to map,
   and access to a cluster, this script will do the heavy 
   lifting so that you can easily split the input into 
   several smaller pieces and all the jobs wil be submitted 
   to the cluster. This is essentially a neat wrapper around 
   blasr's -start and -stride commands.

 + Alert: Be sure to edit the clusterTemplate variable at 
 the head of the script.

 blasrToBed.py 
   This script will convert blasr's .m4 or .m5 format into a 
   BED Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format1 )
   If you would like to visualize the alignments, I 
   reccommend using IGB ( http://bioviz.org/igb/index.html ).

 bedToCoverageWig.py 
   Turn a bed file with alignments into a depth of coverage 
   WIG Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format6 ).

+--------------------+
#LiftOverTable Format#
+--------------------+
 The LiftOver Table is depricated until further notice.

+---+
#FAQ#
+---+

 Who can I report bugs to or ask questions?
   e-mail English@bcm.edu with any PBJelly related issue. 
   If, however your problem is with the SMRTAnalysis 
   software please consult http://www.smrtcommunity.com/

