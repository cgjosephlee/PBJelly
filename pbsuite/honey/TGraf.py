#!/usr/bin/env python
import sys, bisect, argparse
from collections import Counter
import pysam
from pbsuite.utils.setupLogging import setupLogging

USAGE="""Parse a BAM with tails mapped to create clusters of tails
Use parameters to filter reported clusters."""

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
        self.promaq = getTag(read, "XQ")
        if self.propos is not None:
            one = True
            if self.propos <= begin:
                s, e, a, b, uq, dq = (self.propos, begin, "p", "i", self.promaq, self.read.mapq)
                ud = 3 if self.prostr == 0 else 5
                dd = 5 if self.read.is_reverse else 3
            else:
                dd = 3 if self.prostr == 0 else 5
                ud = 5 if self.read.is_reverse else 3
                s, e, a, b, uq, dq = (begin, self.propos, "i", "p", self.read.mapq, self.promaq)
            inv = False if self.prostr == strand else True    
            
        self.epiref = getTag(read, "ZR")
        self.epistr = getTag(read, "ZI")
        self.epipos = getTag(read, "ZP")
        self.epimaq = getTag(read, "ZQ")
        if self.epipos is not None:
            if self.epipos <= end:
                s, e, a, b, uq, dq = (self.epipos, end, "e", "i", self.epimaq, self.read.mapq)
                ud = 3 if self.epistr == 0 else 5
                dd = 5 if self.read.is_reverse else 3
            else:
                s, e, a, b, uq, dq = (end, self.epipos, "i", "e", self.read.mapq, self.epimaq)
                dd = 3 if self.epistr == 0 else 5
                ud = 5 if self.read.is_reverse else 3
            inv = False if self.epistr == strand else True
        
        if self.proref is not None and self.epiref is not None:
            pass
        
        self.uBreak = s
        self.dBreak = e
        self.uTail = a
        self.dTail = b
        self.uDir = ud
        self.dDir = dd
        self.uMapq = uq
        self.dMapq = dq
        self.isInverted = inv
            
    def getEpilog(self):
        pass
        
    def getProlog(self):
        """
        returns the Bread that is paired with this guy
        """
        pass
    
    def near(self, other):
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
        # this creates 2 cluters - one per strand
        if self.uDir != other.uDir or self.dDir != other.dDir:
            #I need a compliment function
            if self.read.is_reverse != other.read.is_reverse:
                if self.uDir != other.uDir and self.dDir != other.dDir:
                    return True
            return False
        
        return True
        
        
    def __str__(self):
        return "%s %d %d' %d %s %s %d %d' %d %s\t%s" % (self.uTail, self.uMapq, \
                    self.uDir, self.uBreak, "%" if self.isInverted else "-", \
                    "<-" if self.read.is_reverse else "->", \
                     self.dBreak, self.dDir, self.dMapq, self.dTail, self.read.qname) 
    
    def __lt__(self, other):
        return self.uBreak < other.uBreak
    
    def __gt__(self, other):
        return self.uBreak > other.uBreak
 
class Bnode(Bread):
    """
    A node for a sorted list where a cluster of similar breakpoint
    are pooled 
    """
    def __init__(self, bread):
        """
        start with a point
        """
        Bread.__init__(self, bread.read)
        self.breads = [bread]
        
    def breadMatch(self, bread):
        """
        Checks if the bread is near any of the bread in this node
        """
        for i in self.breads:
            if bread.near(i):
                return True
        return False
    
    def addBread(self, bread):
        """
        put more bread in the basket
        """       
        self.breads.append(bread)
    
    def numReads(self):
        """
        returns the number of reads supporting in the node
        """
        return len(self.breads)
        
    def numUniqueReads(self):
        """
        returns the number of unique reads in the node
        """
        count = Counter([x.read.qname for x in self.breads])
        return len(count)

    def numUniqueZMWs(self):
        x = []
        for y in self.breads:
            x.append('/'.join(y.read.qname.split('/')[:2]))
        return len(Counter(x))
    
    def avgMapq(self):
        """
        return average mapping quality
        """
        x = []
        for y in self.breads:
            x.append(y.uMapq); x.append(y.dMapq)
        return sum(x)/float(len(x))
        
    def __str__(self):
        ret = "Bnode w/ %d Breads %d unique sub %d unique zmws\n" % \
              (len(self.breads), self.numUniqueReads(), self.numUniqueZMWs())
        for i in self.breads:
            ret += str(i)+"\n"
        return ret
       
def makeBreakReads(bam, buffer=500):
    """
    Extracts all of the tail-mapped reads from a bam and crates break reads (Bread)
    that are then bisect placed inside of
    """
    ret = {}
    for read in bam:
        refName = bam.getrname(read.tid)
        if refName not in ret.keys():
            ret[refName] = []
        clist = ret[refName]
        if read.flag & 0x1:
            if read.flag & 0x40 or read.flag & 0x80: 
                continue; matepos = getTag(read, "YP")
                
            else: #just primary
                pan = Bread(read)
                point = bisect.bisect_left(clist, pan)
                
                unear = False; dnear = False
                #while moving upstream and I'm within buffer, 
                #see if I've got someone to merge with
                lpoint = point
                while lpoint > 0 and abs(pan.uBreak - clist[lpoint - 1].uBreak) <= BUFFER:
                    if clist[lpoint-1].breadMatch(pan):
                        unear = True
                        break
                    lpoint -= 1
                    
                dpoint = point
                while dpoint < len(clist) and abs(pan.uBreak - clist[dpoint].uBreak) <= BUFFER:
                    if clist[dpoint].breadMatch(pan):
                        dnear = True
                        break
                    dpoint += 1
                
                if not (unear or dnear):
                    bisect.insort( clist, Bnode(pan) )
                elif unear and not dnear:
                    clist[lpoint-1].addBread( pan )
                elif not unear and dnear:
                    clist[dpoint].addBread(pan)
                elif unear and dnear:
                    node = Bnode( pan )
                    for i in clist[lpoint-1].breads:
                        node.addBread( i )
                    for i in clist[dpoint].breads:
                        node.addBread( i )
                    del(clist[dpoint])
                    del(clist[lpoint-1])
                    bisect.insort( clist, node )
    return ret
    
def bNodeMerge(node1, node2):
    ret = Bnode(node1.read)
    
def getTag(read, tagId):
    """
    Returns the specified tag or None from an AlignmentRecord
    """
    for i in read.tags:
        if i[0] == tagId:
            return i[1]
    return None

def parseArgs():
    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    parser.add_argument("-B", "--buffer", type=int, default=500, \
                        help=("Buffer around breaks reads must fall "
                              "within to become clustered (500)"))
    parser.add_argument("-b", "--minBreads", type=int, default=4,\
                        help="Minimum number of reads (4)")
    parser.add_argument("-z", "--minZMWs", type=int, default=1, \
                        help="Minimum number of unique ZMWs (1)")
    parser.add_argument("-q", "--minMapq", type=int, default=200, \
                        help="Minimum average map quality score (200)")
    args = parser.parse_args()
    global BUFFER
    BUFFER = args.buffer
    #if args.output is None:
        #args.output = args.bam[:-4]+".hon"
    
    return args


if __name__ == '__main__':
    args = parseArgs()
    bam = pysam.Samfile(args.bam,'rb')
    points = makeBreakReads(bam)
    
    for i in points:
        print "chrom", i, "-", len(points[i]),"clusters"
        clu = 0
        for j in points[i]:
            if j.numUniqueReads() >= args.minBreads \
               and j.numUniqueZMWs >= args.minZMWs \
               and j.avgMapq >= args.minMapq:
                print "-"*10
                print "cluster %d" % clu, j
                fout = open("chrom%s_clu%d.fastq" % (i, clu),'w')
                for r in j.breads:
                    r = r.read
                    fout.write("@%s\n%s\n+\n%s\n" % (r.qname, r.seq, r.qual))
                fout.close()
                clu += 1
        print "#"*10

    
    #clusterBreads(points)
