#!/usr/bin/env python
import sys, bisect
import pysam
import networkx as nx

BUFFER = 500
class Bnode():
    """
    A node for a sorted list where a cluster of similar breakpoint
    are pooled 
    """
    def __init__(self, bread=None):
        """
        start with a point
        """
        self.breads = [bread]
        self.upoint = bread.uBreak #furthest upstream breakpoint
        self.dpoint = bread.dBreak #furthest downstream breakpoint
        
    def addBread(self, bread):
        self.breads.append(bread)
        if bread.uBreak < self.upoint:
            self.upoint = bread.uBreak
        if bread.dBreak > self.dpoint:
            self.dpoint = bread.dBreak
        
    def __eq__(self,other):
        return True if self.tid == other.tid and self.overlap(other) else False

    def __ne__(self,other):
        return False if self == other else True
        
    def __gt__(self,other):
        if self != other:
            if self.tid != other.tid:
                return self.tid > other.tid
            return self.RIPstart > other.RIPstart
        return False
    
    def __lt__(self,other):
        if self != other:
            if self.tid != other.tid: 
                return self.tid < other.tid
            return self.RIPstart < other.RIPstart
        return False
    
    def __ge__(self,other):
        return True if self == other or self > other else False
    
    def __le__(self,other):
        return True if self < other or self == other else False



class Bread():
    """
    Holds a read that has a break in it
    and all relevant information for clustering
    """
    def __init__(self, read):
        """
        extract information from pysam.AlignedRead 
        """
        self.read = read
        
        if read.is_reverse:
            begin = read.aend
            end = read.pos
            strand = 1
        else:
            begin = read.pos
            end = read.aend
            strand = 0
            
        one = False
        self.proref = getTag(read, "XR")
        self.prostr = getTag(read, "XI")
        self.propos = getTag(read, "XP")
        if self.propos is not None:
            one = True
            if self.propos <= begin:
                s, e, a, b = (self.propos, begin, "p", "i")
                ud = 3 if self.prostr == 0 else 5
                dd = 5 if not self.read.is_reverse else 5
            else:
                dd = 3 if self.prostr == 0 else 5
                ud = 5 if not self.read.is_reverse else 5
                s, e, a, b = (begin, self.propos, "i", "p")
            inv = False if self.prostr == strand else True    
            
        self.epiref = getTag(read, "ZR")
        self.epistr = getTag(read, "ZI")
        self.epipos = getTag(read, "ZP")
        if self.epipos is not None:
            if self.epipos <= end:
                s, e, a, b = (self.epipos, end, "e", "i")
                ud = 5 if self.epistr == 0 else 3
                dd = 5 if not self.read.is_reverse else 3
            else:
                s, e, a, b = (end, self.epipos, "i", "e")
                dd = 5 if self.epistr == 0 else 3
                ud = 5 if not self.read.is_reverse else 3
            inv = False if self.epistr == strand else True
        
        if self.proref is not None and self.epiref is not None:
            sys.stderr.write("FUCK!... we got a two tail %s\n" % (self.read.qname))
            #I can move this up to just after where I make the epi stuff
            #and decision fork reassigning the ud belo
        self.uBreak = s
        self.dBreak = e
        self.uTail = a
        self.dTail = b
        self.uDir = ud
        self.dDir = dd
        self.isInverted = inv
            
    def getEpilog(self):
        pass
        
    def getProlog(self):
        """
        returns the Bread that is paired with this guy
        """
        pass
    
    def nearRead(self, other):
        """
        Is this Bread and it's mate near the other Bread
        """
        #Same target
        if self.read.tid != other.read.tid:
            return False
        # are our components within buffer bp of each other
        if abs(self.uBreak - other.uBreak) > BUFFER:
            return False
        if abs(self.dBreak - other.dBreak) > BUFFER:
            return False
        # are we moving in the same direction
        if self.uDir != other.dDir or self.dDir != other.dDir:
            return False
        
        return True
        
    def nearNode(self, node):
        """
        Is this Bread near a node -- tricky
        """
        pass
   
    def __str__(self):
        return "%s %d %s %s %d %s\t%s" % (self.uTail, self.uBreak, \
                                  "%" if self.isInverted else "-", \
                           "<-" if self.read.is_reverse else "->", \
                            self.dBreak, self.dTail, self.read.qname)
    
    def __lt__(self, other):
        return self.uBreak < other.uBreak
    
    def __gt__(self, other):
        return self.dBreak > other.uBreak
        
def makeBreakReads(bam):
    """
    Extracts all of the tail-mapped reads from a bam and crates break reads (Bread)
    """
    ret = []
    for read in bam:
        if read.flag & 0x1:
            if read.flag & 0x40:   #pro
                continue
                matepos = getTag(read, "YP")
            elif read.flag & 0x80: #epi
                continue
                matepos = getTag(read, "YP")
            else: #just primary
                ret.append( Bread(read) )
                continue
                
    return ret
    
def getTag(read, tagId):
    """
    Returns the specified tag or None from an AlignmentRecord
    """
    for i in read.tags:
        if i[0] == tagId:
            return i[1]
    return None
    
if __name__ == '__main__':
    bam = pysam.Samfile(sys.argv[1],'rb')
    points = makeBreakReads(bam)
    genomeBPs = []
    for i in points:
        bisect.insort(genomeBPs, i)
    
    for i in genomeBPs:
        print i
    
    #clusterBreads(points)
