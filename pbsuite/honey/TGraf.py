#!/usr/bin/env python
import sys, bisect, argparse, tarfile, StringIO, os, pwd, grp, logging, time
from tempfile import NamedTemporaryFile
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
        self.proref = getTag(read, "PR")
        self.prostr = getTag(read, "PI")
        self.propos = getTag(read, "PP")
        self.promaq = getTag(read, "PQ")
        if self.propos is not None:
            one = True
            if self.propos <= begin:
                s, e, a, b, uq, dq = (self.propos, begin, "p", "i", self.promaq, self.read.mapq)
                ud = '3' if self.prostr == 0 else '5'
                dd = '5' if self.read.is_reverse else '3'
            else:
                dd = '3' if self.prostr == 0 else '5'
                ud = '5' if self.read.is_reverse else '3'
                s, e, a, b, uq, dq = (begin, self.propos, "i", "p", self.read.mapq, self.promaq)
            inv = False if self.prostr == strand else True    
            
        self.epiref = getTag(read, "ER")
        self.epistr = getTag(read, "EI")
        self.epipos = getTag(read, "EP")
        self.epimaq = getTag(read, "EQ")
        
        if self.epipos is not None and self.epimaq > self.promaq:
            if self.epipos <= end:
                s, e, a, b, uq, dq = (self.epipos, end, "e", "i", self.epimaq, self.read.mapq)
                ud = '3' if self.epistr == 0 else '5'
                dd = '5' if self.read.is_reverse else '3'
            else:
                s, e, a, b, uq, dq = (end, self.epipos, "i", "e", self.read.mapq, self.epimaq)
                dd = '3' if self.epistr == 0 else '5'
                ud = '5' if self.read.is_reverse else '3'
            inv = False if self.epistr == strand else True
        
        if self.proref is not None and self.epiref is not None:
            #Doubles
            pass
        
        #Points
        self.uBreak = s 
        self.dBreak = e
        #P,I,E
        self.uTail = a
        self.dTail = b
        #Strands
        self.uDir = ud
        self.dDir = dd
        #Mapqs
        self.uMapq = uq if uq is not None else 255
        self.dMapq = dq if dq is not None else 255
        #Inv
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
        
        
    def getInvStr(self):
        return "%" if self.isInverted else "-"
    
    def getRevStr(self):
        return "<-" if self.read.is_reverse else "->"
    
    #def header(self):
        #return "uTail uMapq uDir uBreak isInv dir dBreak dDir dMapq dTail readName"
        
    def anyNone(self):
        if self.uTail is None:
            logging.debug("uTail none %s" % str(self.read.qname))
        elif self.uMapq is None:
            logging.debug("uMapq none %s" % str(self.read.qname))
        elif self.uDir is None:
            logging.debug("uDir  none %s" % str(self.read.qname))
        elif self.uBreak is None:
            logging.debug("uBrea none %s" % str(self.read.qname))
        elif self.getInvStr() is None:
            logging.debug("getIn none %s" % str(self.read.qname))
        elif self.getRevStr() is None:
            logging.debug("getRe none %s" % str(self.read.qname))
        elif self.dBreak is None:
            logging.debug("dBrea none %s" % str(self.read.qname))
        elif self.dDir is None:
            logging.debug("dDir  none %s" % str(self.read.qname))
        elif self.dMapq is None:
            logging.debug("dMapq none %s" % str(self.read.qname))
        elif self.dTail is None:
            logging.debug("dTail none %s" % str(self.read.qname))
    
    def __str__(self):
        
        return "%s %d %s %d %s %s %d %s %d %s\t%s" % (self.uTail, self.uMapq, \
                    self.uDir, self.uBreak, self.getInvStr(), self.getRevStr(), \
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
        
    def estimateEnds(self):
        """
        Go through all of the breaksand try to make the best estimate
        of where they are.

        I'm going to create a mode average median of the breaks
        """
        upstream = map(lambda x: x.uBreak, self.breads)
        self.upPos = Counter(upstream).most_common()
        self.avgUpPos = sum(upstream)/len(upstream)
        
        dnstream = map(lambda x: x.dBreak, self.breads)
        self.dnPos = Counter(dnstream).most_common()
        self.avgDnPos = sum(dnstream)/len(dnstream)
        
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
        
    #def prettyHeader(self):
        #data = self.breads[0].header
        #data = data[:data.rindex(' ')]
        #return data + " numReads numZmws"
    
    def toPrettyStr(self):
        """
        Need to make this average mapq and starts and stuff
        """
        self.estimateEnds()
        
        self.uBreak = self.avgUpPos
        self.dBreak = self.avgDnPos
        
        #This should be in a method...
        
        self.uMapq = sum(map(lambda x: x.uMapq, self.breads)) / len(self.breads)
        self.dMapq = sum(map(lambda x: x.dMapq, self.breads)) / len(self.breads)
        
        readData = []
        for i in self.breads:
            readData.append(i.uTail + i.uDir + i.getInvStr() + i.getRevStr() +  i.dTail + i.dDir)
        
        #uTails = set(map(lambda x: x.uTail, self.breads))
        #self.uTail = ""
        #for i in ['p','i','e']:
            #if i in uTails:
                #self.uTail += i
        
        #dTails = set(map(lambda x: x.dTail, self.breads))
        #self.dTail = "" 
        #for i in ['p','i','e']:
            #if i in dTails:
                #self.dTail += i
        
        #uDirs = set(map(lambda x: x.uDir, self.breads))
        #self.uDir = ""
        #for i in ['5', '3']:
            #if i in uDirs:
                #self.uDir += str(i)+"'"
        
        #dDirs = set(map(lambda x: x.dDir, self.breads))
        #self.dDir = ""
        #for i in ['5', '3']:
            #if i in dDirs:
                #self.dDir += str(i)+"'"
        
        data = Bread.__str__(self)
        data, read = data.split('\t')
        #us = ",".join(map(lambda x: "%d:%d" % x, self.upPos))
        #ds = ",".join(map(lambda x: "%d:%d" % x, self.dnPos))
        data += " %d %d %s" % (self.numUniqueReads(), self.numUniqueZMWs(), ":".join(readData))
                                  #us, ds)
        return data.replace(' ','\t')
        
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
                continue; matepos = getTag(read, "IP")
                
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
    parser.add_argument("-B", "--buffer", type=int, default=200, \
                        help=("Buffer around breaks reads must fall "
                              "within to become clustered (200)"))
    parser.add_argument("-b", "--minBreads", type=int, default=2,\
                        help="Minimum number of reads (2)")
    parser.add_argument("-z", "--minZMWs", type=int, default=2, \
                        help="Minimum number of unique ZMWs (2)")
    parser.add_argument("-q", "--minMapq", type=int, default=100, \
                        help="Minimum average map quality score (100)")
    parser.add_argument("-f", "--fastq", action="store_true", \
                        help="Write fastq for each cluster into a .tgz archive (false)")
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output file to write results (BAM.tgraf)")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true", \
                        help="Print each read inside of a cluster to stdout.")
    args = parser.parse_args()
    global BUFFER
    BUFFER = args.buffer
    if args.output is None:
        args.output = args.bam[:-4]+".tgraf"
    setupLogging(args.debug)
    return args


if __name__ == '__main__':
    args = parseArgs()
    bam = pysam.Samfile(args.bam,'rb')
    logging.info("Parsing Reads")
    points = makeBreakReads(bam)
    
    if args.fastq:
        tarOut = tarfile.open(args.output+".tgz", 'w:gz')
    
    fout = open(args.output,'w')
    fout.write("#Args: %s\n" % str(args))
    #uChrom dChrom
    fout.write("#id chrom uTails uMapq uDirs uBreak dir isInv dir dBreak dDirs dMapq dTails numReads numZMWs\n")
    logging.info("Writing Results")
    clu = 0
    for chrom in points:
        #print "chrom", i, "-", len(points[i]),"clusters"
        logging.info("Chrom %s made %d pre-filter clusters" % (chrom, len(points[chrom])))
        postCnt = 0
        for j in points[chrom]:
            if j.numUniqueReads() >= args.minBreads \
               and j.numUniqueZMWs() >= args.minZMWs \
               and j.avgMapq() >= args.minMapq:
                postCnt += 1
                fout.write(str(clu) + " " + chrom + " " + j.toPrettyStr()+"\n")
                if args.fastq:
                    fastq = StringIO.StringIO()
                    tfn = NamedTemporaryFile(suffix=".bam", delete=False).name
                    align = pysam.Samfile(tfn, 'wb', template=bam)
                    for r in j.breads:
                        read = r.read
                        fastq.write("@%s\n%s\n+\n%s\n" % (read.qname, read.seq, read.qual))
                        align.write(read)
                    info = tarOut.tarinfo()
                    info.name  = "clu%d.fastq" % (clu)
                    info.uname = pwd.getpwuid(os.getuid())[0]
                    info.gname = grp.getgrgid(os.getgid())[0]
                    info.size  = fastq.len
                    info.mtime = time.time()
                    #print dir(info)
                    fastq.seek(0)
                    align.close()
                    tarOut.addfile(info, fastq)
                    tarOut.add(tfn, arcname="clu%d.bam" % clu)
                    os.remove(tfn)
                if args.verbose:
                    print "##Cluster %d - %s" % (clu, chrom)
                    for r in j.breads:
                        print r

                    
                clu += 1
        logging.info("Chrom %s made %d post-filter clusters" % (chrom, postCnt))
    
    fout.close()
    if args.fastq:
        tarOut.close()
    
    logging.info("Finished")
