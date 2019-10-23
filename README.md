# PBJelly
***This is not official repo. All updates here are for in-house using and no warranty.***

Official repo: https://sourceforge.net/projects/pb-jelly/

## Requirements
Here listed the tested version.

- python 3.7.3
- networkx 1.7
- blasr 5.3.3

## Installation
Using a conda environment is encouraged.

```bash
conda create -n PBJelly python blasr
conda activate PBJelly
git clone https://github.com/cgjosephlee/PBJelly
cd PBJelly
python setup.py install
```

## Usage
`<stage>` is one of
- setup
- mapping
- support
- extraction
- assembly
- output

```
Jelly.py -h
Jelly.py <stage> -h
```

## Quick start
- Input sequence files must be end with `.fasta` or `.fastq`.
- If you are inputing fasta files, use `fakeQuals.py` to generate dummy quality scores.
- More details in `docs/`.

Example `Protocol.xml`:
```xml
<jellyProtocol>
    <reference>/FULL/PATH/TO__/PBJelly/data/reference/lambda.fasta</reference>
    <outputDir>/FULL/PATH/TO__/PBJelly/lambdaExample/</outputDir>
    <blasr>--minMatch 8 --sdpTupleSize 8 --minPctIdentity 75 --bestn 1 --nCandidates 10 --maxScore -500 --nproc 8 --noSplitSubreads</blasr>
    <input baseDir="/FULL/PATH/TO__/PBJelly/lambdaExample/data/reads/">
        <job>pacbioReads.fasta</job>
    </input>
</jellyProtocol>
```
See also `docs/TemplateProtocol.xml` and `docs/jellyExample/Protocol.xml`.

Example data:
```bash
cd PBJelly/docs/jellyExample
# Edit Protocol.xml according to your path
fakeQuals.py data/reference/lambda.fasta data/reference/lambda.qual
Jelly.py setup Protocol.xml
Jelly.py mapping Protocol.xml
Jelly.py support Protocol.xml
Jelly.py extraction Protocol.xml
Jelly.py assembly Protocol.xml -x "--nproc=8"
Jelly.py output Protocol.xml
```

## Reminder
- Only PBJelly works.
- PBJelly will edit `input.fasta` and make a backup `input.fasta.original`, use alias if you feel it's annoying.
- Sometimes it fails without returning error code, go check each `xxx.err`.
