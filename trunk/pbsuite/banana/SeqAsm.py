import sys

revComp = str.maketrans("ATCGNatcgn","TAGCNtagcn")

"""
Given reads aligned to a range,

I want to assemble them into
the highest quality contig possible

Then I want to be able to put that contig back 
on the reference and show clear differences


"""

import pysam
chrom = "1"
start = 1048122
end = 1049122
buffer = 1000

bam = "/users/english/english/StructuralVariation/LupskiBreakpoints/pacbio/HS1011.full.bam"
bam = pysam.Samfile(bam)

ref = "/users/english/pacbio/data/references/hs37d5/sequence/hs37d5.fasta"
reference = pysam.Fastafile(ref)
print(">ref\n" + reference.fetch(chrom, start , end ))


for read in bam.fetch(chrom, start - buffer, end + buffer):
    score = 0
    trimS = None
    trimE = None
    if start > read.pos:
        for queryPos, targetPos in read.aligned_pairs:
            #print(queryPos, targetPos, trimS, start)
            if trimS is None and targetPos is not None and targetPos >= start:
                trimS = queryPos
            else:
                score += abs(read.pos - start)
    
    if end < read.aend:
        for queryPos, targetPos in read.aligned_pairs[::-1]:
            if trimE is None and targetPos is not None and targetPos <= end:
                trimE = queryPos
            else:
                score += abs(read.aend-end)

            if trimS is not None:
                #trimS = max(0, trimS) + upS
                trimS = max(0, trimS)
            else:
                trimS = 0

            if trimE is not None:
                #trimE = min(len(read.seq), trimE)  - dnS
                trimE = min(len(read.seq), trimE)
            else:
                trimE = len(read.seq)
    #print(trimS, trimE)
    #print(read.query)
    seq = read.query[trimS:trimE]

    qual = read.qqual[trimS:trimE]
    #if read.is_reverse:
        #seq = seq.translate(revComp)[::-1]
        #qual = qual[::-1]
    print(">" + read.qname + '\n' + seq)
    #print("alignLen trimStart, trimEnd, readLen, trimLen, seq")
    #print(start, end)
    #print(read.pos, read.aend, read.aend - read.pos, trimS, trimE, len(read.seq), len(seq), seq)
    #print(qual)
