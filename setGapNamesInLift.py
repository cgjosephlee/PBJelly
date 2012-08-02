import sys, json
from collections import defaultdict, namedtuple

sys.path.append("/users/english/english/Jelly/Jelly")
from FileHandlers import FastaFile, QualFile, LiftOverTable, revComp

"""
Figure out the gap name for every gap-entry in a liftOverTable

Let's work on renamed stuff first.

Then we'll rename whatever new_seq we got to a name that is easily identifiable to createAGP
for when we need to create those stupid coordinates

This will  likely need to be merged with agpLiftMerge creation so that I can update the coordinates used in the pacbio
sequence.
"""
NewSeqInfo = namedtuple("NewSeqInfo", "contigSeq contigQual fillStart fillEnd")
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
    
def findGapName(gapInfo, entry):
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
    
    if supportType == "Span":
        SEQ, STRAND, START, END = ("SpansGap", "SpanStrand", "SpanStart", "SpanEnd")
    elif supportType == "Left":
        SEQ, STRAND, START, END = ("LeftContig", "LeftStrand", "LeftStart", "LeftEnd")
    elif supportType == "Right":
        SEQ, STRAND, START, END = ("RightContig", "RightStrand", "RightStart", "RightEnd")
        
    fh = open(os.path.join(assemblyDir,gapName,"fillingMetrics.json"),'r')
    fillMetrics = json.load(fh)
    fh.close()
    
    fasta = FastaFile(os.path.join(assemblyDir,gapName,"ouptut.fasta"))
    qual = QualFile(os.path.join(assemblyDir,gapName,"ouptut.qual"))
    newName = "%s,%d,%d" % (entry.scaffold, entry.oStart, entry.oEnd)
    seqName = fillMetrics[SEQ]
    contigFasta = fasta[seqName]
    contigQual = qual[seqName]
    if fillMetrics[STRAND] == '-':
        contigFasta = contigFasta.translate(revComp)[::-1]
        contigQual = contigQual[::-1]
    
    return NewSeqInfo(contigFasta, contigQual, fillMetrics[START], fillMetrics[END])

def addNewSeqInfo(gapInfo, liftOverTable, assemblyDir)
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
                nextEntry.new_seqInfo = getNewSequenceInfo(assemblyDir, prevEntry, gapName, "Right")
    
    return liftOverTable

if __name__ == '__main__':
    gapInfo = customGapInfoFile(sys.argv[1])
    liftOverTable = LiftOverTable(sys.argv[2])
    assemblyDir = sys.argv[3]
    addNewSeqInfo(gapInfo, liftOverTable, assemblyDir)
       
