#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from StringIO import StringIO
from collections import defaultdict
from FileHandlers import FastaFile, QualFile, revComp, wrap, qwrap

USAGE = """Usage: %prog --agp <input.agp> --fasta <contigs.fasta> [--qual <contigs.qual> --output <outputName>]
Creates a scaffolding from contigs.fasta/qual using input.agp
"""

def __parseArgs():
    parser = OptionParser(usage=USAGE)
    parser.add_option("-a", "--agp", type="string", default=None, \
                      help="Input AGP File")
    parser.add_option("-f", "--fasta", type="string", default=None, \
                      help="Input Fasta File")
    parser.add_option("-q", "--qual", type="string", default=None, \
                      help="Input Qual File. default is None")
    parser.add_option("-o", "--output", type="string", default="output", \
                      help="Basename for output files. default is output")
    
    opts, args = parser.parse_args(sys.argv)
    
    if opts.agp is None or not os.path.exists(opts.agp):
        parser.error("Agp file specified is invalid. (%s)" % (opts.agp))
    if opts.fasta is None or not os.path.exists(opts.fasta):
        parser.error("Fasta file specified is invalid. (%s)" % (opts.fasta))
    if opts.qual is not None and not os.path.exists(opts.qual):
        parser.error("Qual file specified is invalid. (%s)" % (opts.qual))
    
    return opts
        
if __name__ == '__main__':
    opts = __parseArgs()
    sys.stderr.write("Loading Sequences\n")
    fastaOut = defaultdict(StringIO)
    contigsFasta = FastaFile(opts.fasta)
    if opts.qual is not None:
        contigsQual = QualFile(opts.qual)
        qualOut = defaultdict(list)
    
    #For each agp line
    sys.stderr.write("Parsing Agp\n")
    fh = open(opts.agp, 'r')
    for line in fh.readlines():
        #Skip Header
        if line.startswith('#'):
            continue
        
        data = line.strip().split('\t')
        #We're looking at a gap
        if data[4] in ['N', 'U']:
            sys.stderr.write("Building gap (len: %s, scaf: %s)\n" % (data[5], data[0]))
            fastaOut[data[0]].write('N'* int(data[5]))
            if opts.qual is not None:
                qualOut[data[0]].extend([0] * int(data[5]))
            continue
        
        sys.stderr.write("Building contig (len: %s, scaf: %s)\n" % (data[5], data[0]))
        #Grab the sequence from the contigs file
        seq = contigsFasta[data[5]][int(data[6])-1:int(data[7])]
        if opts.qual is not None:
            qual = contigsQual[data[5]][int(data[6])-1:int(data[7])]
        
        #Reverse Compliment if necessary
        if data[8] == '-':
            seq = seq.translate(revComp)[::-1]
            if opts.qual is not None:
                #Something breaks here
                try:
                    qual = qual[::-1]
                except IndexError:
                    sys.stderr.write("Error!\n\t%s\n\t%s\n" % (qual, line))
                    exit(1)
        
        #Append it
        fastaOut[data[0]].write(seq)
        if opts.qual is not None:
            qualOut[data[0]].extend(qual)
        
    #Write It
    fout = open(opts.output+".fasta", 'w')
    if opts.qual is not None:
        qout = open(opts.output+".qual", 'w')
        
    for entry in fastaOut:
        fout.write(">%s\n%s\n" % (entry, wrap(fastaOut[entry].getvalue())))
        if opts.qual is not None:
            qout.write(">%s\n%s\n" % (entry, qwrap(qualOut[entry])))
    sys.stderr.write("Finished!\n")
    #Done
    fout.close()
    if opts.qual is not None: 
        qout.close()
    
