#!/usr/bin/env python

import sys, os
from collections import namedtuple
from optparse import OptionParser

sys.path.append("/users/english/english/Jelly/Jelly")

from FileHandlers import FastaFile, QualFile, LiftOverTable, LiftOverEntry, wrap, qwrap
"""
Given a fasta file, qual file, and liftover table, 
We undo any changes that result in new_sequences that are < 10bp long (which happens)

Todo:
    make user friendly.
    Not in this script - but we should never make 0bp new_sequences in Collection.py
"""

USAGE = """USAGE: %prog <input.fasta> <input.qual> <liftOverTable.txt> [--options]
Removes small sequence put into gaps. 
Repalces gap_closed feature gap_closed_flagged
Renames gap_reduced to gap_reduced_u5_u3, 
    Where u5 means 5' new_sequence was removed.
    Where u3 means 3' new_sequence was removed."""

if __name__ == '__main__':
    #I have to keep all changes within
    #entries of the liftOver table and can't use update scaffold until
    #everyone has their information
    parser = OptionParser(USAGE)
    parser.add_option("-m", "--min", default=10, type="int",
            help=("Minimum size of new_sequences that isn't removed\n" \
                  "Default is 10bp"))
    parser.add_option("-o","--output",default="reference",
            help=("Name of file to output (DEFAULT=reference)\nWarning! Output files are overwriten!"))
    
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("Expected exactly 3 arguments.")

    fastaName, qualName, liftTableName = args
    
    Change = namedtuple("change", "entry shift start end length size")

    fasta = FastaFile(fastaName)
    qual = QualFile(qualName)
    table = LiftOverTable(liftTableName)
    
    for entry in fasta:
        fasta[entry] = list(fasta[entry])
    
    changes = []
    for entry in table:
        
        if entry.gType == 'gap_reduced' or entry.gType == 'gap_overfilled':
            changed = False
            size = 0
            if entry.prev.gType == 'new_sequence' and entry.prev.nEnd - entry.prev.nStart < opts.min:
                changed = True
                size += entry.prev.nEnd - entry.prev.nStart
                entry.gType += "_u5"
                entry.nStart = entry.prev.nStart
                table.removeEntry(entry.prev)
            
            if entry.next.gType == 'new_sequence' and entry.next.nEnd - entry.next.nStart < opts.min:
                changed = True
                size += entry.next.nEnd - entry.next.nStart
                entry.gType += "_u3"
                entry.nEnd = entry.next.nEnd
                table.removeEntry(entry.next)
                
            if changed:
                changes.append(Change(entry, 0, entry.nStart, entry.nEnd, entry.nEnd - entry.nStart, size))
                        
        elif entry.gType == 'gap_closed' and entry.next.gType == 'new_sequence' and entry.next.nEnd - entry.next.nStart < opts.min:
            entry.gType = 'gap_closed_flagged'
            #replace the gap
            size = entry.next.nEnd - entry.next.nStart
            gapLen = entry.oEnd - entry.oStart
            
            entry.nStart = entry.next.nStart
            entry.nEnd = entry.nStart + gapLen
            
            shift = gapLen - (entry.next.nEnd - entry.next.nStart) 
            start = entry.next.nStart
            end = entry.next.nEnd
            table.removeEntry(entry.next)
            changes.append(Change(entry, shift, start, end, gapLen, size))
    
    changes.sort(cmp = lambda x,y: x.start - y.start, reverse=True)
    
    sys.stderr.write("Removing %d new_sequence less than %d bp\n" % (len(changes), opts.min))
    sys.stderr.write("Total number of bases removed: %d\n" % (sum(map(lambda x: x.size, changes))))
    
    for c in changes:
        fasta[c.entry.scaffold][c.start:c.end] = 'N' * c.length
        qual[c.entry.scaffold][c.start:c.end] = [0] * c.length
        table.updateScaffold(c.entry, c.shift)
    

    liftOut = open("%s.liftOver.txt" % (opts.output),'w')
    liftOut.write("#scaffoldName\toStart\toEnd\tnStart\tnEnd\tfeatureType\n")
    liftOut.write(str(table))
    liftOut.close()
    
    fout = open("%s.fasta" % (opts.output),'w')
    qout = open("%s.qual" % (opts.output),'w')
    for entry in fasta.keys():
        fout.write(">"+entry+"\n"+wrap("".join(fasta[entry]))+"\n")
        qout.write(">"+entry+"\n"+" ".join(map(str,qual[entry]))+"\n")
    fout.close()
    qout.close()

            