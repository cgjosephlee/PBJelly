#!/usr/bin/env python
import sys, bisect, argparse, tarfile, StringIO, os, pwd, grp, logging, time
from tempfile import NamedTemporaryFile
from collections import Counter
from string import maketrans
import pysam
from pbsuite.utils.setupLogging import setupLogging

USAGE="""\
Parse and cluster mapped tails from a bam to create breakpoint candidates."""

class Bread():
    """
    Holds a read that has a break in it
    and all relevant information for clustering
    """
    def __init__(self, read, readRef):
        """
        extract information from pysam.AlignedRead 
        """
        self.read = read
        self.readRef = readRef
        
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
        self.prorem = getTag(read, "PS")
        if self.propos is not None:
            one = True
            if self.propos <= begin:
                s, e, a, b, uq, dq, uR, dR = (self.propos, begin, \
                                              "p", "i", self.promaq, \
                                              self.read.mapq, self.proref, \
                                              readRef)
                ud = '3' if self.prostr == 0 else '5'
                dd = '5' if self.read.is_reverse else '3'
            else:
                s, e, a, b, uq, dq, uR, dR = (begin, self.propos, \
                                              "i", "p", self.read.mapq, \
                                              self.promaq, readRef, \
                                              self.proref)
                dd = '3' if self.prostr == 0 else '5'
                ud = '5' if self.read.is_reverse else '3'

            rmSeq = self.prorem if self.prorem is not None else 0
            inv = False if self.prostr == strand else True    
            
        self.epiref = getTag(read, "ER")
        self.epistr = getTag(read, "EI")
        self.epipos = getTag(read, "EP")
        self.epimaq = getTag(read, "EQ")
        self.epirem = getTag(read, "ES")
        #Choose higher quality
        if self.epipos is not None and self.epimaq > self.promaq:
            if self.epipos <= end:
                s, e, a, b, uq, dq, uR, dR = (self.epipos, end, \
                                              "e", "i", self.epimaq, \
                                              self.read.mapq, self.epiref, \
                                              readRef)
                ud = '3' if self.epistr == 0 else '5'
                dd = '5' if self.read.is_reverse else '3'
            else:
                s, e, a, b, uq, dq, uR, dR = (end, self.epipos, \
                                              "i", "e", self.read.mapq, \
                                              self.epimaq, readRef, \
                                              self.epiref)
                dd = '3' if self.epistr == 0 else '5'
                ud = '5' if self.read.is_reverse else '3'
            rmSeq = self.epirem if self.epirem is not None else 0
            inv = False if self.epistr == strand else True
        
        self.uRef = uR
        self.dRef = dR
        #reference key; for sorting
        j = [uR, dR]; j.sort()
        self.refKey = "_".join(j)
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
        #Remain
        self.remainSeq = rmSeq
    
    def near(self, other):
        """
        Is this Bread and it's mate near the other Bread
        """
        #Same target
        if self.refKey != other.refKey:
            return False
        
        # are our components within buffer bp of each other
        if abs(self.uBreak - other.uBreak) > BUFFER:
            return False
        if abs(self.dBreak - other.dBreak) > BUFFER:
            return False
        """
            del
            ->p| |i->       ud=3,dd=3
            ->i| |e->       ud=3,dd=3
            <-i| |p<-       ud=5,dd=5
            <-e| |i<-       ud=5,dd=5
            
            ins gain same as del but close can estimate size from remaining

            ins sequence
            |i<-  <-e|      ud=5,dd=5   
            |p<-  <-i|      ud=5,dd=5
            |e->  ->i|      ud=3,dd=3
            |i->  p->|      ud=3,dd=3
            
            inv sequence
            ->p|    <-i|    ud=3,dd=5
            ->i|    <-e|    ud=3,dd=5
               |i->    |<-p ud=3,dd=5
               |e->    |<-i ud=3,dd=5

               |p<-    |i-> ud=5,dd=3
               |i<-    |e-> ud=5,dd=3
            <-i|    ->p|    ud=5,dd=3
            <-e|    ->i|    ud=5,dd=3

            
            missed adapter evidence is any one of these. (a || b)
            both would suggest some kind of duplication inversion
            but that requires non-local information.
            These are the same as the inversion/but the sequence is 
            on top of itself
            ->i|        a       ud=3
            <-e|        a       dd=5

            ->p|        a       ud=3
            <-i|        a       dd=5
            
               |e->     b       ud=3
               |i<-     b       dd=5
            
               |i->     b       ud=3
               |p<-     b       dd=5
        """
        # are we moving in the same direction
        # this creates 2 cluters - one per strand
        if self.uDir == other.uDir and self.dDir == other.dDir:
            return True
        elif self.read.is_reverse != other.read.is_reverse:
            #If we're on opposite strands, 
            #but our pieces are pointing together
            if self.uDir != other.uDir and self.dDir != other.dDir:
                    return True
        else:
            if self.uTail == 'i' and other.uTail in ['p', 'e'] or \
               self.dTail == 'i' and other.dTail in ['p', 'e']:
                if self.uDir == other.dDir and self.dDir == other.uDir \
                   and self.read.is_reverse == other.read.is_reverse:
                    return True
        return False
        
    def getInvStr(self):
        return "%" if self.isInverted else "="
    
    def getRevStr(self):
        return "<-" if self.read.is_reverse else "->"
    
    def annotate(self):
        """
        based on the properties of orientation, create annotation
        of what possible variant is here
        """
        if self.uRef == self.dRef:
            if self.uDir != self.dDir:
                return "INV"
            if self.dBreak - self.uBreak < 100 and self.remainSeq >= 100:
                return "INS"
            elif self.uDir == self.dDir:
                ut, dt = self.__tailtoint__()
                if (self.uDir == '3' and ut > dt) or (self.dDir == '5' and ut < dt):
                    return "INS"
                else:
                    return "DEL"
        else:
            return "TLOC"
            
    def __tailtoint__(self):
        """
        returns the uTail and dTail as ints
        """
        trans = {"p":1,
                 "i":2,
                 "e":3}
        x = trans[self.uTail]
        y = trans[self.dTail]
        return x, y
    
    def bpStr(self):
        def swap(a, b):
            trans = maketrans("<>","><")
            if a in ['p','i','e']:
                b = b.translate(trans)[::-1]
            elif b in ['p','i','e']:
                a = a.translate(trans)[::-1]
            return (b, a)
        
        x, y = self.__tailtoint__()
        if (x < y):
            uC = ("->", self.uTail)
            dC = (self.dTail, "->")
        elif (x > y):
            uC = (self.uTail, "->")
            dC = ("->", self.dTail)

        if self.uDir == '5':
            uC = swap(*uC)
        if self.dDir == '5':
            dC = swap(*dC)
        
        return "".join(uC) + self.getInvStr() + "".join(dC)
        
    def anyNone(self):
        """
        This is an old debugging method
        """
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
    
    def getBriefData(self):
        """
        gets the relevant data in a list
        """
        return [self.uRef, self.uBreak, self.uMapq, \
                self.dRef, self.dBreak, self.dMapq, \
                self.remainSeq, self.read.qname]
    
    def __str__(self):
        """
        Need to update this...
        Outputs the string representation of the bread

        chrom   uTails  uMapq   uDirs   uBreak  isInv   dBreak  dDirs   dMapq   dTails  remainSeq
        chrom_chrom uRef    uTail   uMapq   uDir uBreak isInv   remainSeq   dBreak  dDir    dMapq   dTail   dRef
        
        #id chrKey 
        
        uRef uBreak uMapq dRef dBreak dMapq remainSeq numReads numZMWs evidence
        """
        data = self.getBriefData()
        data.insert(-1, self.bpStr())
        data.insert(-1, self.annotate())
        #       ur ub uq dr db dq rs bp an rn
        return "%s %d %d %s %d %d %d %s %s %s" % tuple(data)
        
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
        Bread.__init__(self, bread.read, bread.readRef)
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
    
    def avgMapq(self, threshold=None):
        """
        return average mapping quality
        """
        if threshold is None:
            x = []
            for y in self.breads:
                x.append(y.uMapq); x.append(y.dMapq)
            return sum(x)/float(len(x))
        else:
            us = []
            ds = []
            for y in self.breads:
                us.append(y.uMapq)
                ds.append(y.dMapq)
            us = sum(us)/float(len(us))
            ds = sum(ds)/float(len(ds))
            return (us >= threshold) and (ds >= threshold)
    
    def avgRemainSeq(self):
        """
        return the average length of the remaining sequence
        """
        x = 0
        for y in self.breads:
            x += y.remainSeq
        return x/len(self.breads)

    def toPrettyStr(self):
        """
        Gets the average positions and mapqs, places that into
        the parent Bread, and returns a fully formatted string
        """
        self.estimateEnds()
        
        self.uBreak = self.avgUpPos
        self.dBreak = self.avgDnPos
        
        #This should be in a method...
        
        self.uMapq = sum(map(lambda x: x.uMapq, self.breads)) / len(self.breads)
        self.dMapq = sum(map(lambda x: x.dMapq, self.breads)) / len(self.breads)
        
        readData = set() #getBPstr
        for i in self.breads:
            readData.add(i.bpStr())
        readData = list(readData)
        readData.sort()
        self.remainSeq = self.avgRemainSeq()
              
        data = Bread.getBriefData(self)[:-1]
        data.append(Bread.annotate(self))
        data.append(self.numUniqueReads())
        data.append(self.numUniqueZMWs())
        data.append(";".join(readData))
        
        return "\t".join([str(x) for x in data])
        
            
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
        if not (read.flag & 0x1) or (read.flag & 0x40 or read.flag & 0x80):
            continue; 
        
        #just primary
        pan = Bread(read, refName)
        refKey = pan.refKey
        if refKey not in ret.keys():
            ret[refKey] = []
        clist = ret[refKey]
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

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="tails", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    parser.add_argument("-B", "--buffer", type=int, default=200, \
                        help=("Buffer around breaks reads must fall "
                              "within to become clustered (%(default)s)"))
    parser.add_argument("-b", "--minBreads", type=int, default=3,\
                        help="Minimum number of reads (%(default)s)")
    parser.add_argument("-z", "--minZMWs", type=int, default=3, \
                        help="Minimum number of unique ZMWs (%(default)s)")
    parser.add_argument("-q", "--minMapq", type=int, default=100, \
                        help="Minimum average map quality score per breakpoint (%(default)s)")
    parser.add_argument("-f", "--fastq", action="store_true", \
                        help="Write fastq for each cluster into a .tgz archive (%(default)s)")
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output file to write results (BAM.hon.tails)")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("-v", "--verboseFile", action="store_true", \
                        help="Print each read inside of a cluster to <output>.verbose (%(default)s)")
    args = parser.parse_args(argv)
    global BUFFER
    BUFFER = args.buffer
    if args.output is None:
        args.output = args.bam[:-4] + ".hon.tails"
    if args.verboseFile:
        args.verboseFile = args.output + '.verbose'
        
    setupLogging(args.debug)
    return args


def run(argv):
    args = parseArgs(argv)
    bam = pysam.Samfile(args.bam,'rb')
    logging.info("Parsing Reads")
    points = makeBreakReads(bam)
    
    if args.fastq:
        tarOut = tarfile.open(args.output+".tgz", 'w:gz')
    
    fout = open(args.output,'w')
    fout.write("#Args: %s\n" % str(args))
    #uChrom dChrom
    fout.write(("#id\tchrKey\tuRef\tuBreak\tuMapq\tdRef\tdBreak\tdMapq"
                "\tremainSeq\tannot\tnumReads\tnumZMWs\tevidence\n"))
    logging.info("Writing Results")

    if args.verboseFile:
        vOut = open(args.verboseFile, 'w')
        vOut.write("#uRef uBreak uMapq dRef dBreak dMapq remainSeq break annot readName\n")
    clu = 0
    for chrom in points:
        #print "chrom", i, "-", len(points[i]),"clusters"
        logging.info("Chrom %s made %d pre-filter clusters" % (chrom, len(points[chrom])))
        postCnt = 0
        for j in points[chrom]:
            if j.numUniqueReads() >= args.minBreads \
               and j.numUniqueZMWs() >= args.minZMWs \
               and j.avgMapq(args.minMapq):
                postCnt += 1
                fout.write(str(clu) + "\t" + chrom + "\t" + j.toPrettyStr()+"\n")
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
                if args.verboseFile:
                    vOut.write("##Cluster %d - %s\n" % (clu, chrom))
                    for r in j.breads:
                        vOut.write(str(r)+'\n')
                    
                clu += 1
        logging.info("Chrom %s made %d post-filter clusters" % (chrom, postCnt))
    
    fout.close()
    if args.fastq:
        tarOut.close()
    if args.verboseFile:
        vOut.close()
    
    logging.info("Finished")

if __name__ == '__main__':
    run(sys.argv[1:])
