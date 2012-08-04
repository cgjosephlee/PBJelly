#!/usr/bin/env python
import sys, os
from collections import defaultdict
sys.path.append("/users/english/english/Jelly/Jelly")
from FileHandlers import FastaFile, QualFile, wrap, qwrap

"""
After using agpLiftMerge.py to combine an agp file and a liftover table,

This will create an AGP 2.0

Todo:
    Currently don't preserve linkage evidence, gap type, or linkage from
    the original agp (mainily because agpLiftMerge doesn't work with 2.0, yet.)
    
    make user friendly
"""
keys = ["name", "oStart", "oEnd", "nStart", "nEnd", "gType", "name2", \
        "t1", "t2", "part", "w", "conName", "start", "end", "strand"] 
lookup = {}
for pos,i in enumerate(keys):
    lookup[i] = pos

def makeWLine(entry, partNumber, id=None, start=None, end=None):
    """
    Create an agp W line either soley from a liftoverentry
    or with extra information
    """
    #object object_beg  object_end  part_numer  component_type  component_id component_beg   component_end
    #   orientation
    s = []
    s.append(entry.scaffold)
    s.append(str(entry.nStart+1))
    s.append(str(entry.nEnd))
    s.append(str(partNumber))
    s.append('W')
    if id == None:
        s.append(entry.agpInfo.conName)
        s.append(entry.agpInfo.subStart)
        s.append(entry.agpInfo.subEnd)
        s.append(entry.agpInfo.strand)
    else:
        s.append(id)
        s.append(str(start))
        s.append(str(end))
        s.append('+')
    
    return "\t".join(map(str,s))
    
def makeNLine(entry, partNumber):
    """
    Make an N line from a liftoverentry
    """
    #object object_beg object_end part_number component_type gap_length component_beg component_end 
    #gap_type ='scaffold' linkage ='yes' linkage_evidence ='paired_ends'
    s = []
    s.append(entry.scaffold)
    start = int(entry.nStart) + 1
    end = int(entry.nEnd)
    s.append(str(entry.subStart))
    s.append(str(entry.subEnd))
    s.append(str(partNumber))
    s.append('N')
    s.append(str(entry.subEnd - entry.subStart + 1))
    s.append('scaffold')
    s.append('yes')
    s.append('paired-ends')
    return "\t".join(s)
    
def createSubmissionFiles(agpLftMrg):
    """
    I'm building this around the agpLftMrg that was is build and carried
    from the previous two scripts
    setGapNamesInLift.py (poorly named)
    agpLiftMerge.py (aptly named)
    
    Previously, these scripts had to be run sequenctially, but this is my 
    effort to do it all at once.
    """
    partNumber = 1
    prevScaf = None
    uid = 1
    fout = open("newAssemblySeqs.fasta", 'w')
    qout = open("newAssemblySeqs.qual", 'w')
    aout = open("newAssembly.agp", 'w')
    for entry in agpLftMrg:
        if prevScaf != entry.name:
            partNumber = 1
            prevScaf = entry.name

        if entry.gType == 'contig':#Write W line
            ptrim = 0
            if entry.prev.gType == 'trim':
                ptrim = entry.prev.oEnd - entry.prev.oStart
            
            ntrim = 0
            if entry.next.gtype == 'trim':
                ntrim = entry.next.oEnd - entry.next.oStart
            
            if ntrim or ptrim:
                if entry.agpInfo.strand == '-':
                    entry.agpInfo.subStart += ntrim
                    entry.agpInfo.subEnd -= ptrim
                else:
                    entry.agpInfo.subStart += ptrim
                    entry.agpInfo.subEnd -= ntrim

            aout.write(makeWLine(entry, partNumber)+'\n')
            #Finished Making W Line
            partNumber += 1

        elif entry.gType == 'new_sequence':
            if entry.nEnd - entry.nStart == 0:
                continue
            seqId = "PBJ%07d" % uid
            #I think from here I can use the liftTable new_seqInfo to get this info
            seq = entry.new_seqInfo.contigSeq 
            qual = entry.new_seqInfo.contigQual
            fout.write(">"+seqId+"\n"+wrap(seq)+"\n")
            qout.write(">"+seqId+"\n"+qwrap(qual)+"\n")
            
            aout.write(makeWLine(data, partNumber, id=seqId, \
                start=entry.new_seqInfo.fillStart, end=entry.new_seqInfo.fillEnd)+'\n')
            
            uid += 1
            partNumber += 1
            
        elif entry.gType.startswith("gap") and not entry.gType == 'gap_closed':
            aout.write(makeNLine(entry, partNumber)+'\n')
            partNumber += 1
            
    fout.close()
    qout.close()
    aout.close()
    agpLftMrg.close()

if __name__ == '__main__':
    #agpLiftMerge.txt
    agpLiftMerge = open(sys.argv[1],'r')
    fasta = FastaFile(sys.argv[2])
    qual = QualFile(sys.argv[3])
    createSubmissionFiles(agpLiftMerge, fasta, qual)
    
