#!/usr/bin/env python
"""
Run the normal Jelly pipeline through Setup, Mapping, and Support.
For Setup, set minGap=1 & maxGap=24
For mapping, set blasr parameter for output format to -m 5 #Not easy, yet.
             set bestn parameter to 1
For Support, specify --spanOnly

Then run this script with the parameters 
    <reference.fasta> <reference.qual> <gapInfo.bed> <inputDir>
where inputDir is where you're running the Jelly pipeline.

This does a quick and dirty consensus of all the reads over lowQuality
N's to make a reference base.

This IS NOT an assembly procedure. The accuracy here is lower than normal.
"""
import sys
import os
import json
import re
from glob import glob
from collections import defaultdict
from FileHandlers import M5File, GapInfoFile
from numpy import *

def sortM5(file):
    ret = {}
    for x in file:
        ret[x.qname] = x
    return ret

FindExpandGap = re.compile("[^N\-]([N\-]*N[N\-]*)[^N\-]")

if __name__ == '__main__':
    sys.stderr.write("Warning! This code is currently unsupported.\n")
    lookup = {'-':0, 'A':1, 'T':2, 'C':3, 'G':4, 'N':5}
    revLook = {0:'', 1:'A', 2:'T', 3:'C', 4:'G', 5:'N'}

    gapInfoFile = GapInfoFile(sys.argv[1])
    supportFiles = glob("support/*.gapCans")
    nCorrectingArray = {}
    for supportFileName in supportFiles:
        mappingFile = sortM5(M5File(os.path.join("mapping", \
                                   os.path.basename(supportFileName) \
                                   .replace("gapCans","m4"))))
        supportFile = json.load(open(supportFileName,'r'))
        for gapName in supportFile.keys():
            gap = gapInfoFile[gapName]
            if not nCorrectingArray.has_key(gapName):
                nCorrectingArray[gapName] = zeros( (6, gap.end - gap.start) )

            for readName in supportFile[gapName]["SpansGap"]:
                read = mappingFile[readName]
                start = read.tstart
                gapBase = 0
                for pos,base in enumerate(zip(read.querySeq, read.targetSeq)):
                    qbase, tbase = base
                    if tbase == '-':
                        start -= 1
                        continue
                    if pos+start >= gap.start and pos+start < gap.end:
                        if qbase == '-':#Better
                            continue
                        nCorrectingArray[gapName][lookup[qbase]][gapBase] += 1
                        gapBase += 1

    for gapName in nCorrectingArray.keys():
        out = ""
        for i in nCorrectingArray[gapName].argmax(axis=0):
            out += revLook[i]
        print str(gapInfoFile[gapName]) + "\t" + out
        #print nCorrectingArray[gapName]

def allReplacements(maybe):
    file[entry] = file[entry].replace('M','C')
    file[entry] = file[entry].replace('R','A')
    file[entry] = file[entry].replace('W','T')
    file[entry] = file[entry].replace('S','G')
    file[entry] = file[entry].replace('Y','C')
    file[entry] = file[entry].replace('K','T')
    file[entry] = file[entry].replace('V','G')
    file[entry] = file[entry].replace('H','A')
    file[entry] = file[entry].replace('N','')   
