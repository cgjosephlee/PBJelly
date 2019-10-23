# PBJelly
***This is not official repo. All updates here are for in-house using and no warranty.***

Official repo: https://sourceforge.net/projects/pb-jelly/

## Requirements
Here listed the tested version.

- python 3.7.3
- networkx 1.7
- blasr 5.3.3

## Installation
Using a conda enviroment is encouraged.

```bash
conda create -n PBJelly python blasr
conda activate PBJelly
git clone https://github.com/cgjosephlee/PBJelly
cd PBJelly
python setup.py install
```

## Quick start
More details in `docs/`.

Example `Protocol.xml`:

```xml
<jellyProtocol>
    <reference>/FULL/PATH/TO__/PBJelly/data/reference/lambda.fasta</reference>
    <outputDir>/FULL/PATH/TO__/PBJelly/lambdaExample/</outputDir>
    <cluster>
        <command notes="For single node, multi-core machines" >${CMD} ${JOBNAME} 2> ${STDERR} 1> ${STDOUT} &amp;</command>
        <command notes="For PBS/Moab">echo '${CMD}' | msub -N "${JOBNAME}" -o ${STDOUT} -e ${STDERR} -l nodes=1:ppn=8,mem=48000mb</command>
        <nJobs>1</nJobs>
    </cluster>
    <blasr>--minMatch 8 --sdpTupleSize 8 --minPctIdentity 75 --bestn 1 --nCandidates 10 --maxScore -500 --nproc 8 --noSplitSubreads</blasr>
    <input baseDir="/FULL/PATH/TO__/PBJelly/lambdaExample/data/reads/">
        <job>pacbioReads.fasta</job>
    </input>
</jellyProtocol>
```

Example commands:

```bash
cd PBJelly/docs/jellyExample
# edit Protocol.xml according to your path
Jelly.py setup Protocol.xml
Jelly.py mapping Protocol.xml
Jelly.py support Protocol.xml
Jelly.py extraction Protocol.xml
Jelly.py assembly Protocol.xml -x "--nproc=8"
Jelly.py output Protocol.xml
```

## Reminder
- Only PBJelly is tested.
- Sometimes it fails without returning error code, go check each `xxx.err`.