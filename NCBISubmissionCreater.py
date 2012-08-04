#!/usr/bin/env python

import sys, json, os
from collections import defaultdict, namedtuple
from FileHandlers import *

NewSeqInfo = namedtuple("NewSeqInfo", "contigSeq contigQual fillStart fillEnd")

class AgpInfo():
    def __init__(self, scaffName, start, end, part, featType, conName, subStart, subEnd, strand):
        self.scaffName = scaffName
        self.start = int(start)-1
        self.end = int(end)
        self.part = part
        self.featType = featType
        self.conName = conName
        self.subStart = int(subStart)
        self.subEnd = int(subEnd)
        self.strand = strand
    def __str__(self):
        """
        Returns a string of the agp info.
        Also returns target coordinates to 1 based.
        """
        return "\t".join([self.scaffName, str(self.start+1), str(self.end), self.part, \
                          self.featType, self.conName, str(self.subStart), str(self.subEnd),\
                          self.strand])


def customGapInfoFile(fn):
    """
    Return a dictionary with structure:
    [seqName]:[start:end]=gapName
    -- This helps find the associated gap assembly folder
    """
    ret = defaultdict(dict)
    fh = open(fn,'r')
    for line in fh.readlines():
        seq,start,end,name = line.strip().split('\t')
        ret[seq][start+':'+end] = name
    fh.close()
    return dict(ret)
    
def findGapName(gapInfo, entry):
    """
    figure out a gap's name from the LiftOverEntry's
    information -- 
    This needs to to be refactored. Currently, it 
    is assuming it is working with Dpse reference,
    but that won't always be the case.
    """
    if gapInfo.has_key(entry.scaffold):
        return gapInfo[entry.scaffold]["%d:%d" % (entry.oStart, entry.oEnd)]
    elif entry.scaffold == 'Ch2':
        for key in gapInfo:
            if key.count('chromosome_2') == 1:
                gapInfo[entry.scaffold] = gapInfo[key]
                return gapInfo[key]["%d:%d" % (entry.oStart, entry.oEnd)]
    elif entry.scaffold == 'Ch3':
        for key in gapInfo:
            if key.count('chromosome_3') == 1:
                gapInfo[entry.scaffold] = gapInfo[key]
                return gapInfo[key]["%d:%d" % (entry.oStart, entry.oEnd)]
    else:
        for key in gapInfo.keys():
            try:
                index = key.index('_'+entry.scaffold+'_')
                gapInfo[entry.scaffold] = gapInfo[key]
                return gapInfo[key]["%d:%d" % (entry.oStart, entry.oEnd)]
            except ValueError:
                pass
    print "couldn't find", str(entry)
    exit(1)
    
def getNewSequenceInfo(assemblyDir, entry, gapName, supportType):
    """
    Builds a NewSeqInfo
    for each gap that is supported.
    """
    if supportType == "Span":
        SEQ, STRAND, START, END = ("SpansGap", "SpanStrand", "SpanStart", "SpanEnd")
    elif supportType == "Left":
        SEQ, STRAND, START, END = ("LeftContig", "LeftStrand", "LeftStart", "LeftEnd")
    elif supportType == "Right":
        SEQ, STRAND, START, END = ("RightContig", "RightStrand", "RightStart", "RightEnd")
        
    fh = open(os.path.join(assemblyDir,gapName,"fillingMetrics.json"),'r')
    fillMetrics = json.load(fh)
    fh.close()
    
    fasta = FastaFile(os.path.join(assemblyDir,gapName,"output.fasta"))
    qual = QualFile(os.path.join(assemblyDir,gapName,"output.qual"))
    seqName = fillMetrics[SEQ]
    try:
        contigFasta = fasta[seqName]
    except KeyError:
        print str(entry), gapName
    contigQual = qual[seqName]
    if False and fillMetrics[STRAND] == '-':
        contigFasta = contigFasta.translate(revComp)[::-1]
        contigQual = contigQual[::-1]
    
    return NewSeqInfo(contigFasta, contigQual, fillMetrics[START], fillMetrics[END])

def addNewSeqInfo(gapInfo, liftOverTable, assemblyDir):
    """
    Goes through a liftOverTable and for each closed, reduced, or overfilled gap
    it will create a NewSeqInfo and put it in  LiftOverEntry.new_seqInfo
    """
    for entry in liftOverTable:
        if entry.gType == 'gap_closed':
            gapName = findGapName(gapInfo, entry)
            seqEntry = entry.next
            seqEntry.new_seqInfo = getNewSequenceInfo(assemblyDir, seqEntry, gapName, "Span")
                        
        elif entry.gType == 'gap_overfilled' or entry.gType == 'gap_reduced':
            gapName = findGapName(gapInfo, entry)
            prevEntry = entry.prev
            if prevEntry.gType == "new_sequence":
                prevEntry.new_seqInfo = getNewSequenceInfo(assemblyDir, prevEntry, gapName, "Left")
            
            nextEntry = entry.next
            if nextEntry.gType == "new_sequence":
                nextEntry.new_seqInfo = getNewSequenceInfo(assemblyDir, nextEntry, gapName, "Right")
    
    return liftOverTable

def mergeAgpWithLift(table, agpfile):
    """
    LiftOverTable
    and a file handler holding an agpfile
    """
    curScaffold = None
    line =  agpfile.readline()
    while line != "":
        data = line.strip().split('\t')
        if data[4] != 'W':
            line = agpfile.readline()
            continue
        
        data = AgpInfo(*data)

        if curScaffold != data.scaffName:
            curScaffold = data.scaffName
            curContig = table.scaffoldRoots[curScaffold]
        else:
            curContig = curContig.getNext("contig")
            if curContig == None:
                curContig = table.scaffoldRoots[curScaffold]
        
        #lft ----
        #agp ----
        if curContig.oStart == data.start and curContig.oEnd == data.end:
            curContig.agpInfo = data
            line = agpfile.readline()
        
        #lft    ------
        #agp  -----
        #or agp -----
        elif data.start <= curContig.oStart and curContig.oStart < data.end and curContig.oEnd > data.end:
            #change this lft to the overlap,
            #create a new lft that is the 3' of lft
            #trim the curContig back
            #Need to track the trim to the first contig so that
            #downstream contigs don't pay the price
            threeTrim = 0
            if curContig.prev != None and curContig.prev.gType == 'trim':
                threeTrim = curContig.prev.oEnd - curContig.prev.oStart
            
            newContig = LiftOverEntry(data.scaffName, data.end, curContig.oEnd, \
                                      curContig.nStart + (data.end-curContig.oStart) - threeTrim,\
                                      curContig.nEnd, \
                                      'contig')
            curContig.oEnd = data.end
            curContig.nEnd = newContig.nStart
            data.subStart = str(data.subEnd - (curContig.oEnd - curContig.oStart) + 1) 
            curContig.agpInfo = data
            table.insertEntry(curContig, newContig)
            line = agpfile.readline()

        #lft  ------
        #agp--------
        elif data.start < curContig.oStart and curContig.oEnd == data.end:
            if data.strand == '-':
                newEnd = data.subStart + (curContig.oEnd - curContig.oStart) - 1
                data.subEnd = str(newEnd)
            else:
                data.subStart = str(data.subEnd - (curContig.oEnd - curContig.oStart) + 1) 

            curContig.agpInfo = data
            line = agpfile.readline()
        #lft -----
        #agp -------
        elif data.start == curContig.oStart and curContig.oEnd < data.end:
            if data.strand == '-':
                newStart = data.subEnd - (curContig.oEnd - curContig.oStart) + 1
                data.subStart = str(newStart)
            else:
                data.subStart = str(data.subStart + (curContig.oStart - data.start))
                data.subEnd = str(data.subEnd - (data.end - curContig.oEnd))
            curContig.agpInfo = data
        #lft  ------
        #agp----------
        elif data.start < curContig.oStart and curContig.oEnd < data.end:
            #change the sStart sEnd for this lft's agpInfo,
            if data.strand == '-':
                newStart = data.subStart + (data.end - curContig.oEnd)#(curContig.oStart - start) 
                newEnd = data.subEnd - (curContig.oStart - data.start)#(end - curContig.oEnd)
                data.subStart = str(newStart)
                data.subEnd = str(newEnd)
            else:
                data.subStart = str(data.subStart + (curContig.oStart - data.start))
                data.subEnd = str(data.subEnd - (data.end - curContig.oEnd))
            curContig.agpInfo = data
            #keep the same agpline
        
        #lft  ------
        #agp    ------ 
        elif curContig.oStart < data.start and curContig.oEnd > data.start and curContig.oEnd < data.end:
            #change this lft to the overlap
            #use same agp and nextLft automatically)
            curContig.oStart += start-curContig.oStart
            curContig.nStart += start-curContig.oStart
            data.subEnd = str(data.subEnd - (end-curContig.oEnd))
            curContig.agpInfo = data
        
        #lft --------        
        #agp   ----
        elif data.start > curContig.oStart and data.end < curContig.oEnd:
            print "I don't think this happens but it did -- so refactor your code"
            newContig = LiftOverEntry(data.scaffName, data.end, curContig.oEnd, \
                                      curContig.nStart + (end-curContig.oStart),\
                                      curContig.nEnd, \
                                      'contig')
            curContig.nStart -= curContig.oStart - data.start
            curContig.nEnd -= curContig.oEnd - data.end
            curContig.oStart = data.start
            curContig.oEnd = data.end
            curContig.agpInfo = data
            table.insertEntry(curContig, newContig)
            line = agpfile.readline()
 
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
    #Gap Info
    gapInfo = customGapInfoFile(sys.argv[1])
    #Lift Over Table
    liftOverTable = LiftOverTable(sys.argv[2])
    #Assembly Directory
    assemblyDir = sys.argv[3]
    #Agp File for Original Assembly
    agpFileHandle = open(sys.argv[4],'r')
    
    addNewSeqInfo(gapInfo, liftOverTable, assemblyDir)
    mergeAgpWithLift(liftOverTable, agpFileHandle)
    #outputAgpMergeLiftTable(table)
    createSubmissionFiles(liftOverTable)
