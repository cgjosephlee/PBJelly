#!/usr/bin/env python

import sys, re
from collections import defaultdict, namedtuple
from FileHandlers import M4File, LiftOverTable, LiftOverEntry

"""
Merges a liftOverTable with an associated agp file.

In Development.
ToDo:
    Make compatible with agp v2.0 (just make agp FileHandlers)
    Write docs / make user friendly with optparse
WARNING!
    Say there is a contig Jelly identifies because of gaps.
    (represented by 'X' = sequence and '-' = gaps)
    And this contig is actually composed of multiple contigs
    in the agp (represented by A and B).A
    And the first base of the Jelly contig (X) is trimmed
    (represented by a [:space:])

agp  ---AAABBB---
orig ---XXXXXX----
jelly--- XXXXX----
    
    This trim affects all indexing of bases downstream within
    the original contig.
        
    So our united liftTable and Agp will look like:

    ...
    
    scaffold    3   4   na  na   trim
    scaffold    3   6   3   5   contig  contigA
    scaffold    6   9   5   9
    
"""
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
    
def outputAgpMergeLiftTable(table):
    print "#scaffoldName\toStart\toEnd\tnStart\tnEnd\tfeatureType\tagpInfo"
    
    for entry in table:
        if entry.gType != 'contig':
            print str(entry)
        else:
            try:
                ai = entry.agpInfo
                try:
                    print str(entry)+"\t"+str(ai)
                except TypeError:
                    print str(ai)
                    exit(1)
            except AttributeError:
                print str(entry)+"\t"+"noExtras?"

if __name__ == '__main__':
    table = LiftOverTable(sys.argv[1])
    agpFileHandle = open(sys.argv[2],'r')
    mergeAgpWithLift(table, agpFileHandle)
    outputAgpMergeLiftTable(table)
