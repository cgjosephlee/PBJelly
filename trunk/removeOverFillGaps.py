#!/usr/bin/env python
import sys
from collections import namedtuple
from optparse import OptionParser
from math import sqrt, ceil
from FileHandlers import FastaFile, QualFile, wrap, LiftOverTable, LiftOverEntry

USAGE = """USAGE: %prog <input.fasta> <input.qual> <liftOverTable.txt> [--options]
Removes sequence put into overfilled gaps. 
Renames the Feature Type to gap_overfilled_undo"""

LiftEntry = namedtuple("LiftOverEntry", "scaffold oStart oEnd nStart nEnd gType")

def getStdv(x):
    n, mean, std = len(x), 0, 0 
    if n == 0:
        sys.stderr.write("Error! - No filled gaps to estimate standard deviation from")
        exit(1)
    obs = []
    for i in x: 
        a = (i.oStart - i.oEnd) - (i.nStart - i.nEnd)
        obs.append(a)
        mean = mean + a 
    mean = mean / float(n) 
    
    for a in obs: 
        std = std + (a - mean)**2 
    std = sqrt(std / float(n-1)) 
    
    #se = std/sqrt(len(x))
    return int(mean+ceil(std))
    #return mean, std, se

if __name__ == '__main__':
    parser = OptionParser(USAGE)
    parser.add_option("-m", "--max", default=None, type="int",
            help=("Maximum overfill amount that isn't removed\n"\
                  "Default is one standard deviation in the distribution \n"\
                  "of amount of sequence put into filled gaps \n"\
                  "minus predicted gap size."))
    parser.add_option("-o","--output",default="reference",
            help=("Name of file to output (DEFAULT=reference)"))
    
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("Expected exactly 3 arguments.")

    fastaName, qualName, liftTableName = args
    
    fasta = FastaFile(fastaName)
    #Since we're changing individual bases
    for entry in fasta.keys():
        fasta[entry] = list(fasta[entry])
    
    qual = QualFile(qualName)

    liftTable = LiftOverTable(liftTableName)
    
    fGaps = []#Filled
    oGaps = []#Overfilled
    for entry in liftTable:
        if entry.gType == "gap_filled":
            fGaps.append(entry)
        elif entry.gType == "gap_overfilled":
            oGaps.append(entry)

    if opts.max == None:
        opts.max = getStdv(fGaps)

    sys.stderr.write("Removing Overfills greater than %d bp\n" % (opts.max))
    sys.stderr.write("length %d\n" % (len(oGaps)))
    oGaps.sort(cmp = lambda x,y: x.nStart - y.nStart, reverse=True)
    nCleaned = 0
    for gap in oGaps:
        if ((gap.nEnd - gap.nStart) - (gap.oEnd - gap.oStart)) > opts.max:
            nCleaned += 1
            
            fasta[gap.scaffold][gap.nStart:gap.nEnd] = 'N' * (gap.oEnd - gap.oStart)
            qual[gap.scaffold][gap.nStart:gap.nEnd] = [0] * (gap.oEnd - gap.oStart)
            
            gap.gType += "_undo"
            #Amount put in - amount putting back
            newEnd = gap.nStart + (gap.oEnd - gap.oStart)
            #shift = (gap.nStart - gap.nEnd) - (gap.oEnd - gap.oStart) + nEnd
            shift = gap.nEnd - newEnd
            gap.nEnd = newEnd
            liftTable.updateScaffold(gap, -shift)
    
    sys.stderr.write("Removed %d overfilled gaps\n" % (nCleaned))
    
    liftOut = open(opts.output+".liftOver.txt",'w')
    liftOut.write(str(liftTable))
    liftOut.close()
    
    fout = open(opts.output+".fasta",'w')
    qout = open(opts.output+".qual",'w')
    for entry in fasta.keys():
        fout.write(">"+entry+"\n"+wrap("".join(fasta[entry]))+"\n")
        qout.write(">"+entry+"\n"+" ".join(map(str,qual[entry]))+"\n")

    fout.close()
    qout.close()
    


