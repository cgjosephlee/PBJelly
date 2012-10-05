#!/usr/bin/env python
import sys, re, math
from optparse import OptionParser
from string import Template
from FileHandlers import FastaFile

USAGE="""%prog <file.fasta> [options]
Returns basic statistics (like N50s) about an assembly"""

def parseArgs():
    parser = OptionParser(USAGE)
    parser.add_option("-b", "--binsize", dest="binsize", type="int", default=0,\
                        help="Bin size for creating gap frequency data. (Default is to not print the frequency)")
    parser.add_option("-m", "--min", dest="min", type="int", default=25,\
                        help="Minimum gap size to be considered. DEFAULT=25")
    parser.add_option("-M", "--max", dest="max", type="str", default="",\
                        help="Maximum gap size to be considered. DEFAULT=inf")
    parser.add_option("-c", "--consolidate",dest="consolidate", type=int, default=25,\
                        help=("Concolidate gaps within XXbp as a single gap since this is more" \
                              "indicative of a single LowQuality Region than multiple gaps. DEFAULT=25 [0 means off]"))
    #Put this in. I think it would be nice.
    #parser.add_option("-o","--output", dest="output", type="str", default=None,
                        #help="File name to output a Bed File with chromosome, start, end of gaps features. Default=Off")
    opts, args = parser.parse_args()
    
    if len(args) != 1:
        parser.error("No Fasta Specified!")
    if opts.min < 1:
        parser.error("Minimum gap size must be at least 1.")
    return opts, args[0]
   
def getStats(seqLengths):
    data = {}

    seqLengths.sort(reverse=True)
    
    data["numSeqs"] = len(seqLengths)
    data["totalLength"] = sum(seqLengths)
    tl = data["totalLength"]
    n50_mark = data["totalLength"] * .5
    n90_mark = data["totalLength"] * .90
    n95_mark = data["totalLength"] * .95
    
    data["n50"] = None
    data["n90"] = None
    data["n95"] = None
    basesSeen = 0
    
    for n in seqLengths:
        basesSeen += n
        if data["n50"] is None and basesSeen > n50_mark:
            data["n50"] = n
        if data["n90"] is None and basesSeen > n90_mark:
            data["n90"] = n
        if data["n95"] is None and basesSeen > n95_mark:
            data["n95"] = n
            break
    #may not have gaps
    if data["numSeqs"] == 0:
        return data
    data["min"] = seqLengths[-1]
    data["FstQu"] = seqLengths[ int(math.floor(data["numSeqs"]*.75)) ]
    median = data["numSeqs"]*.50
    data["median"] = int( (seqLengths[ int(math.floor(median)) ] + \
                           seqLengths[ int(math.floor(median)) ]) / 2)
    data["mean"] = data["totalLength"]/data["numSeqs"]
    data["TrdQu"] = seqLengths[ int(math.floor(data["numSeqs"]*.25)) ] 
    data["max"] = seqLengths[0]

    return data
            
def printBins(seq, binsize):
    """
    Print Bin Sizes.
    """
    if binsize < 1:
        exit(0)

    seq.sort()
    
    BINSIZE = binsize
    bin_mark = BINSIZE
    bincount = 0
    i = 0
    while i < len(seq):
        if seq[i] <= bin_mark:
            bincount += 1
            i += 1
        else:
            if bincount != 0:
                print str(bin_mark-BINSIZE+1)+"bp : "+str(bin_mark)+"bp\t"+str(bincount)
            bincount = 0
            bin_mark += BINSIZE

    if bincount != 0:
        print str(bin_mark-BINSIZE+1)+"bp : "+str(bin_mark)+"bp\t"+str(bincount)

if __name__ == '__main__':
    opts, ref = parseArgs()
    
    reference = FastaFile(ref)
    
    gapLengths = []
    contigLengths = []
    scaffoldLengths = []
    gapRE = re.compile("[^Nn]([Nn]{%d,%s})[^Nn]" % (opts.min, opts.max))
    for entry in reference:
        seq = reference[entry]
        scaffoldLengths.append( len(seq) )
        gapCoords = []
        
        for gap in gapRE.finditer( seq ):
            #Finditer gives the full span of the match.
            #The first and last characters of the match are not N's
            #Therefore they are not part of the gap
            gapCoords.append([gap.start() + 1, gap.end() - 1])
        
        if len(gapCoords) == 0:
            continue
        
        #Consolidate gaps that are too close
        i = 0
        while i < len(gapCoords)-1:
            if gapCoords[i+1][0] - gapCoords[i][1] < opts.consolidate:
                gapCoords[i+1][0] = gapCoords[i][0]
                del(gapCoords[i])
            else:
                i += 1
        
        contigLengths.append(gapCoords[0][0])
        gapLengths.append(gapCoords[0][1]-gapCoords[0][0])
        
        for i in range(1, len(gapCoords)):
            contigLengths.append(gapCoords[i][0] - gapCoords[i-1][1])
            gapLengths.append(gapCoords[i][1] - gapCoords[i][0])
        contigLengths.append(len(seq) - gapCoords[-1][1])
            
        #prevStart = 0 # previous contig start
        #contigLengths.append(gap.start() - prevStart - 1)
        #prevStart = gap.end() - 1
        #contigLengths.append(len(seq) - prevStart)
    
    scafStats = getStats(scaffoldLengths)
    contStats = getStats(contigLengths)
    gapStats = getStats(gapLengths)
    
    report = Template("#Seqs\t$numSeqs\n" + \
                      "Min\t$min\n" + \
                      "1st Qu.\t$FstQu\n" + \
                      "Median\t$median\n" + \
                      "Mean\t$mean\n" + \
                      "3rd Qu.\t$TrdQu\n" + \
                      "Max\t$max\n" + \
                      "Total\t$totalLength\n" + \
                      "n50\t$n50\n" + \
                      "n90\t$n90\n" + \
                      "n95\t$n95\n")
    
    if scafStats["numSeqs"] == 0:
        print "="*20
        print "No Scaffolding!"
        print "="*20
    else:
        print "="*20
        print "Scaffold Stats"
        print "="*20
        print report.substitute(scafStats)
        print "="*20
    
    if contStats["numSeqs"] == 0:
        print "No Contigs! (or gaps betwen them)"
        print "="*20
    else:
        print "Contig Stats"
        print "="*20
        print report.substitute(contStats)
        print "="*20
    
    if gapStats["numSeqs"] == 0:
        print "No Gaps!"
        print "="*20
    else:
        print "Gap Stats"
        print "="*20
        print report.substitute(gapStats)
        print "="*20

    if opts.binsize != 0:
        printBins(gapLengths, opts.binsize)
