#!/usr/bin/env python
import sys
from collections import namedtuple
from optparse import OptionParser
from math import sqrt, ceil
from FileHandlers import FastaFile, QualFile, wrap, LiftOverTable

USAGE = """USAGE: %prog <input.fasta> <input.qual> <liftOverTable.txt> [--options]
Removes sequence put into overfilled gaps. 
Renames the Feature Type to gap_overfilled_undo"""


def getStdv(x):
    n, mean, std = len(x), 0, 0 
    if n == 0:
        sys.stderr.write("Error! - No filled gaps to estimate standard deviation from")
        exit(1)
    obs = []
    for i in x: 
        mean = mean + i
    mean = mean / float(n) 
    
    for a in x: 
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
            help=("Name of file to output (DEFAULT=reference)\nWarning! Output files are overwriten!"))
    
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
    outT = []
    oGaps = []#Overfilled
    for entry in liftTable:
        if entry.gType == "gap_closed":
            if entry.next.gType == "new_sequence":
                newSeq = entry.getNext("new_sequence")
                fGaps.append((entry.oEnd - entry.oStart) - (newSeq.nEnd - newSeq.nStart))
        elif entry.gType == "gap_overfilled":
            oGaps.append(entry)

    if opts.max == None:
        opts.max = getStdv(fGaps)

    sys.stderr.write("Removing Overfills greater than %d bp\n" % (opts.max))
    sys.stderr.write("length %d\n" % (len(oGaps)))
    
    oGaps.sort(cmp = lambda x,y: x.nStart - y.nStart, reverse=True)
    nCleaned = 0
    for gap in oGaps:
        fillAmt = 0
        if gap.prev.gType == 'new_sequence':
            fillAmt += gap.prev.nEnd - gap.prev.nStart
        if gap.next.gType == 'new_sequence':
            fillAmt += gap.next.nEnd - gap.next.nStart
            
        if fillAmt > opts.max:
            nCleaned += 1
            oGapSize = gap.oEnd - gap.oStart
            
            pSize = nSize = 0
            if gap.prev.gType == "new_sequence":
                pSize = gap.prev.nEnd - gap.prev.nStart
                pSeq = gap.prev
                liftTable.removeEntry(gap.prev)
            else:
                pSeq = gap
            if gap.next.gType == "new_sequence":
                nSize = gap.next.nEnd - gap.next.nStart
                nSeq = gap.next
                liftTable.removeEntry(gap.next)
            else:
                nSeq = gap
            
            #Amount of sequence removed minus amount we're putting back in
            shift = pSize + nSize + (gap.nEnd - gap.nStart) - oGapSize
            
            gap.gType += "_flagged"
            gap.nStart = pSeq.nStart
            gap.nEnd = gap.nStart + oGapSize
            
            fasta[gap.scaffold][pSeq.nStart:nSeq.nEnd] = 'N' * (oGapSize)
            qual[gap.scaffold][pSeq.nStart:nSeq.nEnd] = [0] * (oGapSize)
           
            liftTable.updateScaffold(gap, -shift)
    
    sys.stderr.write("Removed %d overfilled gaps\n" % (nCleaned))
    
    liftOut = open(opts.output+".liftOver.txt",'w')
    liftOut.write("#scaffoldName\toStart\toEnd\tnStart\tnEnd\tfeatureType\n")
    liftOut.write(str(liftTable))
    liftOut.close()
    
    fout = open(opts.output+".fasta",'w')
    qout = open(opts.output+".qual",'w')
    for entry in fasta.keys():
        fout.write(">"+entry+"\n"+wrap("".join(fasta[entry]))+"\n")
        qout.write(">"+entry+"\n"+" ".join(map(str,qual[entry]))+"\n")

    fout.close()
    qout.close()
    


