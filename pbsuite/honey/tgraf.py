import sys, bisect

import pysam

import networkx as nx

def getPoints(bam):
    """
    for every flag & 0x1:
        if primary:
            extract pro tail
            make tuple (propos, read.pos)
            
            extract epi tail
            make tuple (read.aend, epipos)
            
        if pro:
            extract pro tail
            make tuple (read.aend, pripos)
        if epi:
            extract pro tail
            make tuple (read.pripos, read.pos)
            
    """
    ret = []
    for read in bam:
        if read.flag & 0x1:
            if read.flag & 0x40:   #pro
                continue
                matepos = getTag(read, "YP")
                
                end, direction = (read.pos, "<-") if read.is_reverse else (read.aend, "->")
                
                s, e = (matepos, end) if matepos <= end else (end, matepos)
                ret.append( (s, direction, e) )
                
            elif read.flag & 0x80: #epi
                continue
                matepos = getTag(read, "YP")
                
                begin, direction = (read.aend, "<-") if read.is_reverse else (read.pos, "->")
                
                s, e = (matepos, begin) if matepos <= begin else (begin, matepos)
                
                ret.append( (s, direction, e) )
                
            else: #just primary
                if read.is_reverse:
                    begin = read.aend
                    end = read.pos
                    direction = "<-"
                    strand = 1
                else:
                    begin = read.pos
                    end = read.aend
                    direction = "->"
                    strand = 0
                    
                propos = getTag(read, "XP")
                prostr = getTag(read, "XI")
                if propos is not None:
                    s, e = (propos, begin) if propos <= begin else (begin, propos)
                    inversion = "-" if prostr == strand else "%"
                    ret.append( (s, direction, inversion, e, read.qname) )
                
                epipos = getTag(read, "ZP")
                epistr = getTag(read, "ZI")
                if epipos is not None:
                    s, e = (epipos, end) if epipos <= end else (end, epipos)
                    inversion = "-" if epistr == strand else "%"
                    ret.append( (s, direction, inversion, e, read.qname) )
                
    return ret
    
def getTag(read, tagId):
    """
    Returns tag or None
    """
    for i in read.tags:
        if i[0] == tagId:
            return i[1]
    return None
    
if __name__ == '__main__':
    BUFFER = 500
    bam = pysam.Samfile(sys.argv[1],'rb')
    points = getPoints(bam)
    points.sort()
    for i in points:
        print "%d %s %s %d\t%s" % i
    
