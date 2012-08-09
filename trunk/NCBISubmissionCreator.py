#!/usr/bin/env python

import sys, json, os
from collections import defaultdict, namedtuple
from FileHandlers import *
from optparse import OptionParser

USAGE = """USAGE: %prog <parameters> 
Using information about your original assembly and PBJelly gap-filled assembly, 
this program helps create:
    newAssemblySeqs.fasta -- The new contig sequence PBJelly Created
    newAssemblySeqs.qual -- The new contig qualities PBJelly Created
    newAssembly.agp -- An agp description of the upgraded assembly.

Required Parameters are:
    gapInfo, liftOverTable, assemblyDir, agpFile
"""

"""
Todo: 
    make user friendly (optparse) -- k
    make documentation -- eh
    make a --new option, for new assemblies that doesn't have an existing .agp 
        this will mean I need to create an .agp .fasta .qual for everything.
    
    Work towards preserving gap_type evidence from original .agp
        This will require differentiating between 1.1 and 2.0 agp versions.
"""

def __parseOpts(myArgs):
    parser = OptionParser(USAGE)
    parser.add_option("-g", "--gapInfo", type="str", default=False,
                        help = ("The gapInfo.bed file created during the " \
                                "setup stage for your original reference"))
    parser.add_option("-l", "--liftOverTable", type="str", default=False,
                        help = ("The liftOverTable.txt file created during the "\
                                "output stage"))
    parser.add_option("-d", "--assemblyDir", type="str", default=False, \
                        help = ("The assembly/ directory created during the " \
                                "assembly stage"))
    parser.add_option("-a", "--agpFile", type="str", default=False, \
                        help = ("The AGP for the original reference"))
    opts, args = parser.parse_args(myArgs)
    
    if not opts.gapInfo:
        parser.error("--gapInfo is required")
    if not opts.liftOverTable:
        parser.error("--liftOverTable is required")
    if not opts.assemblyDir:
        parser.error("--assemblyDir is required")
    if not opts.agpFile:
        parser.error("--agpFile is required")
    
    return opts.gapInfo, opts.liftOverTable, opts.assemblyDir, opts.agpFile

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
    contigFasta = fasta[seqName]
    
    contigQual = qual[seqName]
    
    if fillMetrics[STRAND] == '-':
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
                        
        elif entry.gType.startswith('gap_overfilled') or entry.gType.startswith('gap_reduced'):
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
        
        #lft    --------
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
            data.subStart = data.subEnd - (curContig.oEnd - curContig.oStart) + 1
            curContig.agpInfo = data
            table.insertEntry(curContig, newContig)
            line = agpfile.readline()

        #lft  ------
        #agp--------
        elif data.start < curContig.oStart and curContig.oEnd == data.end:
            if data.strand == '-':
                data.subEnd = data.subStart + (curContig.oEnd - curContig.oStart) - 1
            else:
                data.subStart = data.subEnd - (curContig.oEnd - curContig.oStart) + 1

            curContig.agpInfo = data
            line = agpfile.readline()
        #lft -----
        #agp -------
        elif data.start == curContig.oStart and curContig.oEnd < data.end:
            if data.strand == '-':
                data.subStart = data.subEnd - (curContig.oEnd - curContig.oStart) + 1
            else:
                data.subStart = data.subStart + (curContig.oStart - data.start)
                data.subEnd = data.subEnd - (data.end - curContig.oEnd)
            curContig.agpInfo = data
        #lft  ------
        #agp----------
        elif data.start < curContig.oStart and curContig.oEnd < data.end:
            #change the sStart sEnd for this lft's agpInfo,
            if data.strand == '-':
                data.subStart = data.subStart + (data.end - curContig.oEnd)#(curContig.oStart - start) 
                data.subEnd = data.subEnd - (curContig.oStart - data.start)#(end - curContig.oEnd)
            else:
                data.subStart = data.subStart + (curContig.oStart - data.start)
                data.subEnd = data.subEnd - (data.end - curContig.oEnd)
            curContig.agpInfo = data
            #keep the same agpline
        
        #lft  ------
        #agp    ------ 
        elif curContig.oStart < data.start and curContig.oEnd > data.start and curContig.oEnd < data.end:
            #change this lft to the overlap
            #use same agp and nextLft automatically)
            curContig.oStart += start-curContig.oStart
            curContig.nStart += start-curContig.oStart
            data.subEnd = data.subEnd - (end-curContig.oEnd)
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
    s.append(str(start))
    s.append(str(end))
    s.append(str(partNumber))
    s.append('N')
    s.append(str(end - start + 1))
    s.append('scaffold')
    s.append('yes')
    s.append('paired-ends')
    return "\t".join(s)
    
def createSubmissionFiles(agpLftMrg):
    """
    I'm building this around the agpLftMrg that was is built and carried
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
        if prevScaf != entry.scaffold:
            partNumber = 1
            prevScaf = entry.scaffold

        if entry.gType == 'contig':#Write W line
            ptrim = 0
            if entry.prev != None and entry.prev.gType == 'trim':
                ptrim = entry.prev.oEnd - entry.prev.oStart
            
            ntrim = 0
            if entry.next != None and entry.next.gType == 'trim':
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
            try:
                seq = entry.new_seqInfo.contigSeq 
            except AttributeError:
                print 'no seqInfo!', str(entry)
                exit(1)
            qual = entry.new_seqInfo.contigQual
            fout.write(">"+seqId+"\n"+wrap(seq)+"\n")
            qout.write(">"+seqId+"\n"+qwrap(qual)+"\n")
            
            aout.write(makeWLine(entry, partNumber, id=seqId, \
                start=entry.new_seqInfo.fillStart + 1, end=entry.new_seqInfo.fillEnd)+'\n')
            
            uid += 1
            partNumber += 1
            
        elif entry.gType.startswith("gap") and not entry.gType == 'gap_closed':
            aout.write(makeNLine(entry, partNumber)+'\n')
            partNumber += 1
            
    fout.close()
    qout.close()
    aout.close()

    
if __name__ == '__main__':
    gapInfo, liftOverTable, assemblyDir, agpFile = __parseOpts(sys.argv)
    #Gap Info
    #gapInfo = customGapInfoFile(sys.argv[1])
    #Lift Over Table
    #liftOverTable = LiftOverTable(sys.argv[2])
    #Assembly Directory
    #assemblyDir = sys.argv[3]
    #Agp File for Original Assembly
    #agpFileHandle = open(sys.argv[4],'r')
    
    #Make new_seqInfo for new_sequence entries
    addNewSeqInfo(gapInfo, liftOverTable, assemblyDir)
    #Makes agpInfo to all contig entries
    mergeAgpWithLift(liftOverTable, agpFileHandle)
    #outputAgpMergeLiftTable(table)
    #output newSeqs.fasta and newSeqs.qual and newAssembly.agp
    createSubmissionFiles(liftOverTable)

