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

def makeWLine(data, partNumber, id=None, length = None):
    #object object_beg  object_end  part_numer  component_type  component_id component_beg   component_end
    #   orientation
    s = []
    s.append(data[lookup["name"]])
    s.append(str(int(data[lookup["nStart"]])+1))
    s.append(str(data[lookup["nEnd"]]))
    s.append(str(partNumber))
    s.append('W')
    if id == None:
        s.append(data[lookup["conName"]])
        s.append(data[lookup["start"]])
        s.append(data[lookup["end"]])
        s.append(data[lookup["strand"]])
    else:
        s.append(id)
        s.append('1')
        s.append(str(length))
        s.append('+')
    return "\t".join(map(str,s))
    
def makeNLine(data, partNumber):
    #object object_beg object_end part_number component_type gap_length component_beg component_end 
    #gap_type ='scaffold' linkage ='yes' linkage_evidence ='paired_ends'
    s = []
    s.append(data[lookup["name"]])
    start = int(data[lookup["nStart"]]) + 1
    end = int(data[lookup["nEnd"]])
    s.append(str(start))
    s.append(str(end))
    s.append(str(partNumber))
    s.append('N')
    s.append(str(end-start+1))
    s.append('scaffold')
    s.append('yes')
    s.append('paired-ends')
    return "\t".join(s)
    
def customGapInfoFile(fn):
    """
    Return a dictionary with structure:
    [seqName]:[start:end]=gapName
    """
    ret = defaultdict(dict)
    fh = open(fn,'r')
    for line in fh.readlines():
        seq,start,end,name = line.strip().split('\t')
        ret[seq][start+':'+end] = name
    fh.close()
    return dict(ret)
        

if __name__ == '__main__':
    #agpLiftMert.txt
    fh = open(sys.argv[1],'r')
    fasta = FastaFile(sys.argv[2])
    qual = QualFile(sys.argv[3])
    fh.readline()#header
    #Gap Info File
    #gapInfo = customGapInfoFile(sys.argv[3])
    #Assembly Folder
    #assemblyFolder = sys.argv[4]
    
    partNumber = 1
    prevScaf = None
    uid = 1
    fout = open("dpse_sequences.fasta",'w')
    qout = open("dpse_sequences.qual",'w')
    aout = open("dpse_Agp.agp",'w')
    lines = fh.readlines()
    for pos,line in enumerate(lines):
        data = line.strip().split('\t')
        if prevScaf != data[lookup["name"]]:
            partNumber = 1
            prevScaf = data[lookup["name"]]

        if data[lookup["gType"]] == 'contig':#Write W line
            ptrim = 0
            if pos != 0:
                pdata = lines[pos-1].strip().split('\t')
                if pdata[lookup["gType"]] == 'trim':
                    ptrim = int(pdata[lookup["oEnd"]]) - int(pdata[lookup["oStart"]])
            ntrim = 0
            if pos != len(lines)-1:
                ndata = lines[pos+1].strip().split('\t');
                if ndata[lookup["gType"]] == 'trim':
                    ntrim = int(ndata[lookup["oEnd"]]) - int(ndata[lookup["oStart"]])
            
            if ntrim or ptrim:
                if data[lookup["strand"]] == '-':
                    data[lookup["start"]] = int(data[lookup["start"]]) + ntrim
                    data[lookup["end"]] = int(data[lookup["end"]]) - ptrim
                else:
                    data[lookup["start"]] = int(data[lookup["start"]]) + ptrim
                    data[lookup["end"]] = int(data[lookup["end"]]) - ntrim
            
            aout.write(makeWLine(data, partNumber)+"\n")
            partNumber += 1

        elif data[lookup["gType"]] == 'new_sequence':
            start = int(data[lookup["nStart"]])
            end = int(data[lookup["nEnd"]])
            if end - start == 0:
                continue
            seqId = "PBJ%07d" % uid
            
            #search for the chromosome name directly 
            #if not found, we'll have to parse for it
            #once found, set it's real name to what it should be
            #note that _seqName_ will help except for chromosome 3 & 2
            #gapInfo[
            seq  =  fasta[data[lookup["name"]]][start:end]
            q =  qual[data[lookup["name"]]][start:end]
            fout.write(">"+seqId+"\n"+wrap(seq)+"\n")
            qout.write(">"+seqId+"\n"+qwrap(q)+"\n")
            aout.write(makeWLine(data, partNumber, seqId, len(seq))+'\n')
            
            uid += 1
            partNumber += 1
            
        elif data[lookup["gType"]].startswith("gap") and not data[lookup["gType"]] == 'gap_closed':
            aout.write(makeNLine(data, partNumber)+'\n')
            partNumber += 1
            
    fout.close()
    qout.close()
    aout.close()
    fh.close()
            
