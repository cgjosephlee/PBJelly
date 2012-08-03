#!/usr/bin/env python

import sys
from optparse import OptionParser
from collections import namedtuple
from FileHandlers import wrap, qwrap


USAGE = """Usage: %prog <input.fastq> [--output baseName]
Splits a fastq into <baseName>.fasta and <baseName>.qual
Assumes Sanger Encoded Phred Scores in fastq
"""

def __parseArgs():
    parser = OptionParser(usage=USAGE)
    parser.add_option("-o", "--output", type="string", default="output", \
                      help="Basename for .fasta and .qual output files")
    
    opts, args = parser.parse_args(sys.argv)
    if len(args) != 2: parser.error('Expected 1 arguments')
    
    return opts.output, args[1]

def fastqIter( fn ):
    fh = open(fn, 'r')
    FastQEntry = namedtuple("FastQEntry", "name seq qual")
    while True:
        name = fh.readline().strip()[1:]
        seq = fh.readline().strip()
        plus = fh.readline().strip()
        qual = fh.readline().strip()
        
        if qual == "": break

        yield FastQEntry(name, seq, qual)

def phredToQual( qual ):
    """
    Take a qual string that is phred/sanger encoded
    turn it into a list of quals
    """
    return map(lambda x: ord(x)-33, list(qual))
    
if __name__ == '__main__':
    output, fastq = __parseArgs()
    
    fout = open(output+".fasta", 'w')
    qout = open(output+".qual", 'w')

    for entry in fastqIter(fastq):
        fout.write(">%s\n%s\n" % (entry.name, wrap(entry.seq)))
        qout.write(">%s\n%s\n" % (entry.name, qwrap(phredToQual(entry.qual))))
    
    fout.close()
    qout.close()
