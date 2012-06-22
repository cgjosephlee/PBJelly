import sys

sys.path.append("/users/english/english/Jelly/Jelly")
from collections import namedtuple
from FileHandlers import FastaFile, QualFile, LiftOverTable, LiftOverEntry, wrap, qwrap
"""
Given a fasta file, qual file, and liftover table, 
We undo any changes that result in new_sequences that are < 10bp long (which happens)

Todo:
    make user friendly.
    Not in this script - but we should never make 0bp new_sequences in Collection.py
"""

if __name__ == '__main__':
    #I have to keep all changes within
    #entries of the liftOver table and can't use update scaffold until
    #everyone has their information
    Change = namedtuple("change", "entry shift start end length")

    fasta = FastaFile(sys.argv[1])
    qual = QualFile(sys.argv[2])
    table = LiftOverTable(sys.argv[3])
    
    for entry in fasta:
        fasta[entry] = list(fasta[entry])
    changes = []
    for entry in table:
        
        if entry.gType == 'gap_reduced':
            changed = False
            if entry.prev.gType == 'new_sequence' and entry.prev.nEnd - entry.prev.nStart < 10:
                changed = True
                entry.gType += "_u5"
                entry.nStart = entry.prev.nStart
                table.removeEntry(entry.prev)
            
            if entry.next.gType == 'new_sequence' and entry.next.nEnd - entry.next.nStart < 10:
                changed = True
                entry.gType += "_u3"
                entry.nEnd = entry.next.nEnd
                table.removeEntry(entry.next)
                
            if changed:
                changes.append(Change(entry, 0, entry.nStart, entry.nEnd, entry.nEnd - entry.nStart))
                        
        elif entry.gType == 'gap_closed' and entry.next.gType == 'new_sequence' and entry.next.nEnd - entry.next.nStart < 10:
            entry.gType = 'gap_undoC'
            #replace the gap
            gapLen = entry.oEnd - entry.oStart
            
            entry.nStart = entry.next.nStart
            entry.nEnd = entry.nStart + gapLen
            
            shift = gapLen - (entry.next.nEnd - entry.next.nStart) 
            start = entry.next.nStart
            end = entry.next.nEnd
            table.removeEntry(entry.next)
            changes.append(Change(entry, shift, start, end, gapLen))
    
    changes.sort(cmp = lambda x,y: x.start - y.start, reverse=True)
    
    for c in changes:
        fasta[c.entry.scaffold][c.start:c.end] = 'N' * c.length
        qual[c.entry.scaffold][c.start:c.end] = [0] * c.length
        table.updateScaffold(c.entry, c.shift)
    
    liftOut = open("smallRemove.liftOver.txt",'w')
    liftOut.write("#scaffoldName\toStart\toEnd\tnStart\tnEnd\tfeatureType\n")
    liftOut.write(str(table))
    liftOut.close()
    
    fout = open("smallRemove.fasta",'w')
    qout = open("smallRemove.qual",'w')
    for entry in fasta.keys():
        fout.write(">"+entry+"\n"+wrap("".join(fasta[entry]))+"\n")
        qout.write(">"+entry+"\n"+" ".join(map(str,qual[entry]))+"\n")
    fout.close()
    qout.close()

            
