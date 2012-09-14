PBJelly Documentation -- v12.9.14

Using This README 
Toy Data: 
Provided with this distribution of PBJelly is a toy example inside of the lambdaTest directory. Use this once you've setup Jelly to test that everything is working as expected and to become familiar with the software. 
Commands: 
All commands are presented in the format

[cmd]$: commandToExecute 
where [cmd]$: is the tag that identifies a unix command and commandToExecute would be the actual command. 
Alerts, Notes, and Protips 
Throughout the document you will see tags for extra pieces of information.

Alerts contain information that is crucial to running PBJelly. 
Notes are extra tid-bits of knowledge about PBJelly that may help you run PBJelly. 
Protips are considerations that one may keep in mind to take full advantage of Jelly. 

Installation 
Step 1 
Download and Install the SMRTAnalysis software v1.3.1 from http://www.smrtcommunity.com/SMRT-Analysis/Software/SMRT-Analysis 

Use whichever distribution fits with you system the best. 
Step 2 
Edit exportPaths.sh

set JELLYPATH to the path where you've placed PBJelly. 
set SEYMOUR_HOME to the path where you've installed SMRTAnalysis. 

Note: If SEYMOUR_HOME is already set, PBJelly assumes you already sourced the SMRTAnalysis setup script. Be sure this is the case. 
Step 3 
To your .bash_profile , add the command.

[cmd]$: source <path to PBJelly>/exportPaths.sh 

Be sure to source your .bash_profile before beginning to use PBJelly.

Alert: PBJelly runs using the python distributed with SMRTAnalysis. Be sure that this python is your environment's default after sourcing the .bash_profile 

Quick Start 
For more details on each step in the pipeline, see "Running PBJelly" below. If, however, you'd like to just run the program do the following.

Step 1 
Create your Protocol.xml. To run the lambdaExample dataset provided, edit the paths in the <reference> , <outputDir> and the baseDir attribute in the inputs tag to the full path that lambdaExample is sitting in. See "Running PBJelly - Setup - Step 2" for more information about creating Protocols. 
Step 2 
Run each stage sequentially (meaning finish one stage before continuing to the next). 

To run a stage, use the command

    [cmd]$: Jelly.py <stage> yourProtocol.xml 

The stages, in order, are

setup 
mapping 
support 
extraction 
assembly 
output 

To get help with PBJelly, simply run Jelly.py --help.

To get help with any stage, simply run Jelly.py <stage> --help 

Running PBJelly 
Pre Processing 
If you would like to get some initial stats about your reference (sizes and distributions of scaffolds, contigs, and gaps) run: 

    [cmd]$: summarizeAssembly.py <reference.fasta> 
see summarizeAssembly.py --help for details.

Alert: Every sequence file you use has the following requirements

It must be in fasta format with the full .fasta file extension (no .fa or .fsta etc) 
Every fasta file requires an associated qual sitting in the same directory. 
 Qual files should contain the Phred Scores of bases (0-93) and should not be encoded (i.e. no Sanger/Solexa, only the number for the score) If you do not have quals for your sequences, see [cmd]$: fakeQuals.py --help 

If you have a fastq file as your beginning reference, you can use fastqSplit.py to divide the file into a fasta and qual.

Sequence names will have all spaces replaced with underscores. 

Setup 
This is by far the longest and most involved step. Once you get past this, PBJelly makes the rest of the workflow super easy.

 Gather the paths to your reference and all of your fasta files. 
 Create your Protocol.xml -- See TemplateProtocol.xml for an idea of what a protocol should look like. You can name your protocol whatever you'd like. Below are the elements needed for a Protocol.

<reference> : The reference tag contains the full path to your reference.fasta. Note: All of the files PBJelly creates regarding your reference will be placed in the same directory beside the reference.fasta 
<outputDir> : The output tag contains the full path to where PBJelly will put the intermediate data for each stage in the process. Protip: Placing your Protocol.xml into the <outputDir> directory makes for excellent bookkeeping. Just copy the provided TemplateProtocol.xml into the <outputDir> directory for each project you maintain 
<cluster> [optional]: If you do not wish to submit jobs to a cluster, do not include the <cluster> element in your protocol. Otherwise, the cluster tag contains information that PBJelly will use to submit jobs to your cluster. The cluster tag contains two elements.

<command> : This is a template that holds the structure of how one submits jobs to his/her cluster. The example provided in TemplateProtocol.xml is used on a MOAB job management system. The command templatees has 4 REQUIRED elements: CMD - The command one expects to execute on the cluster; JOBNAME - The name to assign to the job being submitted; STDOUT - The name to assign to the job being submitted; STDERR - The file that the standard error will be directed to. 
<nJobs> : This is the maximum number of jobs PBJelly can submit to the cluster. PBJelly tries to submit as many jobs as possible. So, if you do not specify nJobs, PBJelly will create a job for every single input file for any given stage. For example, if you have 50 files containing read sequences for mapping, and you specify nJobs = 5, PBJelly will submit 5 jobs, each will map 10 input files. If you do not specify nJobs (i.e. don't include the nJobs tag or set nJobs to 0) PBJelly will submit 50 jobs for mapping. 

<blasr> : Place all of the mapping parameter to be sent to blasr here. You can be comfortable with the default blasr parameters that are provied in the TemplateProtocol.xml. However, if you would like to customize the parameters, be sure to read the blasr documentation (blasr --help). I also highly recommend watching the video "Mapping and Alignment using BLASR" located at http://www.smrtcommunity.com/Learn/Videos . Alert : Always specify -noSplitSubreads in your blasr parameters. 
<input> : Input contains information about where your input data is located. First, there is the optional 'baseDir' attribute. If all of your data has a common root path, specifying baseDir=/That/Path will prevent redundancy in the inputs. Inside of the <input> tag is the <job> tag. This contains the path to each input file PBJelly will use. 

Run the Setup stage by executing [cmd]$: Jelly.py setup Protocol.xm l Note : All of the files PBJelly creates regarding your reference will be placed in the same directory beside the reference.fasta 

Mapping 
Run the mapping stage via the command

    [cmd]$: Jelly.py mapping Protocol.xml 

(Told you things would get easier)

Alert : Remember to wait until given stage is finished before running the next stage.

Note: Standard Error and Standard Logs for each step are placed next to the data PBJelly creates. For the mapping step, this can be found in <outputDir>/mapping/ 

Support 
Run the support stage 

[cmd]$: Jelly.py support Protocol.xml 

Extraction 
Run the extraction stage

[cmd]$: Jelly.py extraction Protocol.xml 

Alert: PBJelly currently doesn't handle the data in a memory efficent mannor. It simply loads the reference and all of the input files into memory. This means you need to be sure the machine you're running this stage on has enough memory to handle all of your data.

Assembly 
Run the assembly stage

[cmd]$: Jelly.py assembly Protocol.xml 

Protip: If you have access to more than one core per gap to be assembled, be sure to tell PBJelly to pass the nproc parameter to the assembly stage via:
[cmd]$: Jelly.py assembly Protocol.xml -x "--nproc=8"

Where 8 can be replaced by the number of cores you're using.

Output 
Run the output stage

[cmd]$: Jelly.py output Protocol.xml 

At the head of your log file, you can find information about how many gaps were addressed, filled, etc. The output stage collects all of your results into 3 files:

<outputDir>/assembly/jellyOutput.fasta -- The new reference 
<outputDir>/assembly/jellyOutput.qual -- the new reference qual. 
<outputDir>/assembly/jellyLiftOverTable.txt -- For annotation lifting. 

See the LiftOverTable documentation below for more information.

Post-Processing 
removeOverFillGaps.py This script will undo any of the changes PBJelly made in overfilled gaps. For more information see [cmd]$: removeOverFillGaps.py --help 

summarizeAssembly.py 

Use this script to get metrics describing the sizes and distribution of your new reference's scaffolding, contigs, and gaps. Compare this to the original reference's summary to see the scale of the new reference's improvements.

Extras 
runBlasr.py: 

If you have a large input.fasta file that you need to map, and access to a cluster, this script will do the heavy lifting so that you can easily split the input into several smaller pieces and all the jobs wil be submitted to the cluster. This is essentially a neat wrapper around blasr's -start and -stride commands.

Alert: Be sure to edit the clusterTemplate variable at the head of the script.

blasrToBed.py 

This script will convert blasr's .m4 or .m5 format into a BED Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format1 ). If you would like to visualize the alignments, I reccommend using IGB ( http://bioviz.org/igb/index.html ).

bedToCoverageWig.py 

Turn a bed file with alignments into a depth of coverage WIG Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format6 ).

LiftOverTable Format 
The LiftOverTable is a tab-delimited file that contains all the information about changes Jelly made to the input reference in order to produce the upgraded reference genome. A header exists on the first line of the LiftOverTable and the first character is a '#'.

The LiftOverTable uses a 0-based coordinate system. The best definition of this coordinate system is found in the SAM Format Documentation ( http://samtools.sourceforge.net/SAM1.pdf )
0-based coordinate system: A coordinate system where the first base of a 
sequence is zero. In this coordinate system, a region is specified by a half-closed-half-open interval. For example, the region between the 3rd and the 7th bases inclusive is [2, 7). The BAM, BED, and PSL formats are using the 0-based coordinate system. 
Each line of the LiftOverTable represents one featureType of a reference.

The column definitions (& data type) are:

scaffoldName (string) -- The scaffold on which this feature is located. 
oStart (int) -- The start coordinate for where this feature was located in the original reference. 
oEnd (int) -- The end coordinate for where this feature was located in the original reference. 
nStart (int) -- The start coordinate for where this feature is located in the new reference. May be 'na' for trim features 
nEnd (int) -- The end coordinate for where this feature is located in the new reference. May be 'na' for trim features 
featureType (string) -- The type of feature represented in this range. featureType is one of 8 values.

contig: A contig sequence. 
new_sequence: Sequence added into the assembly by PBJelly. This feature's oStart and oEnd are always 'na' 
trim: Sequence that was trimmed from the boundary of a gap. This features's nStart and nEnd are always 'na' 
gap_reduced: A gap that PBJelly reduced in size by extending one or both of the flanking contigs. 
gap_filled: A gap that PBJelly closed and removed from the new reference. This feature's nStart and nEnd are always 'na' 
gap_overfilled: A gap that PBJelly reduced, but put in more sequence than the original predicted gap size. For example, consider a gap predicted to be 100 bp. PBJelly identifies 75 bp of new_sequence that extends the 5' flanking contig and 75 bp of new_sequence that extends the 3' contig, but cannot identify an overlap between the extenstions, and therefore cannot close the gap. This means PBJelly overfilled the gap by putting 150 bp into a 100 bp predicted gap. 
gap_failedAssembly: A gap that was supported by reads, but PBJelly was unable to assemble a contig that improved the gap. 
gap_unaddressed: A gap that had no reads mapping in a manner to 

FAQ 

Who can I report bugs to or ask questions?

e-mail English@bcm.edu with any PBJelly related issue. If, however your problem is with the SMRTAnalysis software please consult http://www.smrtcommunity.com/
