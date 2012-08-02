#!/usr/bin/env python

import sys, re
from collections import defaultdict 
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
if __name__ == '__main__':
    table = LiftOverTable(sys.argv[1])
    #AGP File
    fh = open(sys.argv[2],'r')
    curScaffold = None
    line = fh.readline()
    while line != "":
        data = line.strip().split('\t')
        if data[4] != 'W':
            line = fh.readline()
            continue

        if curScaffold != data[0]:
            curScaffold = data[0]
            curContig = table.scaffoldRoots[curScaffold]
        else:
            curContig = curContig.getNext("contig")
            if curContig == None:
                curContig = table.scaffoldRoots[curScaffold]
        
        scaf, start, end, part, w, name, sStart, sEnd, strand = data
        start = int(start)-1
        end = int(end)
        sStart = int(sStart)
        sEnd = int(sEnd)
        
        #lft ----
        #agp ----
        if curContig.oStart == start and curContig.oEnd == end:
            curContig.extras = data
            line = fh.readline()
        
        #lft    ------
        #agp  -----
        #or agp -----
        elif start <= curContig.oStart and curContig.oStart < end and curContig.oEnd > end:
            #change this lft to the overlap,
            #create a new lft that is the 3' of lft
            #trim the curContig back
            #Need to track the trim to the first contig so that
            #downstream contigs don't pay the price
            threeTrim = 0
            if curContig.prev != None and curContig.prev.gType == 'trim':
                threeTrim = curContig.prev.oEnd - curContig.prev.oStart
            
            newContig = LiftOverEntry(scaf, end, curContig.oEnd, \
                                      curContig.nStart + (end-curContig.oStart) - threeTrim,\
                                      curContig.nEnd, \
                                      'contig')
            curContig.oEnd = end
            curContig.nEnd = newContig.nStart
            data[6] = str(sEnd - (curContig.oEnd - curContig.oStart) + 1) 
            curContig.extras = data
            table.insertEntry(curContig, newContig)
            line = fh.readline()

        #lft  ------
        #agp--------
        elif start < curContig.oStart and curContig.oEnd == end:
            if data[8] == '-':
                newEnd = sStart + (curContig.oEnd - curContig.oStart) - 1
                data[7] = str(newEnd)
            else:
                data[6] = str(sEnd - (curContig.oEnd - curContig.oStart) + 1) 

            curContig.extras = data
            line = fh.readline()
        #lft -----
        #agp -------
        elif start == curContig.oStart and curContig.oEnd < end:
            if data[8] == '-':
                newStart = sEnd - (curContig.oEnd - curContig.oStart) + 1
                data[6] = str(newStart)
                #data[7] = sEnd
                #data[7] = str(sEnd -  #str(newStart + (sEnd - sStart))
            else:
                data[6] = str(sStart + (curContig.oStart - start))
                data[7] = str(sEnd - (end - curContig.oEnd))
            curContig.extras = data
        #lft  ------
        #agp----------
        elif start < curContig.oStart and curContig.oEnd < end:
            #change the sStart sEnd for this lft's extras,
            if data[8] == '-':
                newStart = sStart + (end - curContig.oEnd)#(curContig.oStart - start) 
                newEnd = sEnd - (curContig.oStart - start)#(end - curContig.oEnd)
                data[6] = str(newStart)

                data[7] = str(newEnd)
            else:
                data[6] = str(sStart + (curContig.oStart - start))
                data[7] = str(sEnd - (end - curContig.oEnd))
            curContig.extras = data
            #keep the same agpline
        
        #lft  ------
        #agp    ------ 
        elif curContig.oStart < start and curContig.oEnd > start and curContig.oEnd < end:
            #change this lft to the overlap
            #use same agp and nextLft automatically)
            curContig.oStart += start-curContig.oStart
            curContig.nStart += start-curContig.oStart
            data[7] = str(sEnd - (end-curContig.oEnd))
            #data[7] = str(sStart + (end-curContig.oStart))
            curContig.extras = data
        
        #lft --------        
        #agp   ----
        elif start > curContig.oStart and end < curContig.oEnd:
            print "I don't think this happens but it did -- refactor"
            newContig = LiftOverEntry(scaf, end, curContig.oEnd, \
                                      curContig.nStart + (end-curContig.oStart),\
                                      curContig.nEnd, \
                                      'contig')
            curContig.nStart -= curContig.oStart - start
            curContig.nEnd -= curContig.oEnd - end
            curContig.oStart = start
            curContig.oEnd = end
            curContig.extras = data
            table.insertEntry(curContig, newContig)
            line = fh.readline()

    
    print "#scaffoldName\toStart\toEnd\tnStart\tnEnd\tfeatureType\tagpFields"
    
    for entry in table:
        if entry.gType != 'contig':
            print str(entry)
        else:
            try:
                x = entry.extras
                try:
                    print str(entry)+"\t"+"\t".join(x)
                except TypeError:
                    print x
                    exit(1)
            except AttributeError:
                print str(entry)+"\t"+"noExtras?"
