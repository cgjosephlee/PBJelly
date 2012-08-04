#!/usr/bin/env python

import sys, os, glob, subprocess
from string import Template
from optparse import OptionParser

#Edit the template below to customize the submission for your cluster.
clusterTemplate = Template("echo '${CMD}' | msub -N \"${JOBNAME}\" -o ${STDOUT} -e ${STDERR} -l nodes=1:ppn=8,mem=48000")


#These are your default blasr parameters. Adjust at will.
parameters = "-bestn 3 -nproc 8 -minSubreadLength 200 -nCandidates 20 -minMatch 8"


command = Template("blasr ${FAS} ${REF} ${SA} -m 4  -out ${OUT} -start ${START} -stride ${STRIDE} ${EXTRAPARAMS}")

USAGE="""runBlasr.py <reads.fasta> <reference.fasta> --output <outName> [--sa <reference.fasta.sa> --stride]
This script builds and submits cluster commands that will perform blasr mapping."""

def _exe(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdoutVal, stderrVal =  proc.communicate()
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal

def parseArgs():
    
    parser = OptionParser(USAGE) 
    parser.add_option("--sa", dest="sa", action="store_true", help="Use reference index. MUST EXIST BESIDE OF REFERENCE SPECIFIED. [Optional]", default=False)
    parser.add_option("--output",dest="output", help="Output results name", default=None)
    parser.add_option("--stride", dest="stride", type="int", help="See -stride and -start in blasr -h", default = 1)
    parser.add_option("-p","--params",type="str", help="Parameters to pass to blasr. Surround string of params with \"'s", default=parameters)
    opts, args = parser.parse_args()
        
    if len(args) != 2:
        parser.error("Error! Expceted 2 arguments")
    reads = os.path.abspath(args[0])
    reference = os.path.abspath(args[1])
    if opts.sa == True:
        refIndex = "-sa " + reference+".sa" #hackey shack
    else:
        refIndex = ""
    
    if opts.output == "None":
        parser.error("Error! Must specify output")
    else:
        output = os.path.abspath(opts.output)
    
    stride = opts.stride
    params = opts.params
    return reads, reference, refIndex, output, stride, params
    

if __name__ == '__main__':
    reads, reference, refIndex, outFile, stride, params = parseArgs()
    
    
    for i in range(stride):
        #Build the stuff for the fasta Command
        myOutFile = outFile + ".chunk_%d.m4" % i
        myParams = {"REF":reference, "SA": refIndex, "FAS":reads, "OUT":myOutFile, "START":i, "STRIDE":stride, "EXTRAPARAMS":params}
        #Build the stuff for the cluster Command
        myCommand = {"CMD":command.substitute(myParams), \
            "JOBNAME":myOutFile, \
            "STDOUT":myOutFile+".out", \
            "STDERR":myOutFile+".err" } 
        #Submit the cluster command
        print _exe(clusterTemplate.substitute(myCommand))
