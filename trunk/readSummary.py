#!/usr/bin/env python
import sys, shutil, os, json
from optparse import OptionParser
from xml.etree import ElementTree

from CommandRunner import exe
from summarizeAssembly import getStats


def countBases(file):
    r,o,e = exe('grep -v "^>" %s | wc -c' % (file))
    if r != 0:
        print r
        print o
        print e
        exit(r)


    return int(o.split(' ')[0])

def countReads(file):
    r,o,e = exe('grep -c "^>" %s' % (file))
    if r != 0:
        print r
        print o
        print e
        exit(r)

    return int(o.split(' ')[0])
    
def readLengths(file):
    r,o,e = exe("fastalength %s | cut -f1 -d\  " % (file))
    if r != 0:
        print r
        print o
        print e
        exit(r)

    data = [int(x) for x in o.strip().split('\n')]
    return data
    
def renameFiles(path, file, name):
    r,o,e = exe('fastqSplit.py {0} -o {1}'.format(\
                os.path.join(path, 'filtered_subreads.fastq'),\
                os.path.join(path, file)))
    if r != 0:
        print r
        print o
        print e
        exit(r)

if __name__ == '__main__':
    parser = OptionParser("Get Summary of all reads pointed to inside of a Protocol.xml")
    opts, args = parser.parse()
    if len(args) != 1:
        parser.error("Expected one argument... the Protocol.xml")
    x = ElementTree.parse(args[0])
    root = x.getroot()
    j = root.find("input")
    os.chdir(j.attrib["baseDir"])
    bases = 0
    reads = 0
    lengths = []
    for dir in j.getchildren():
        f = dir.text.replace('.fasta','.fastq')
        path, file = os.path.split(dir.text)
        name = f.split('/')[-3]
        #print f, path, file, name
        lengths.extend(readLengths(dir.text))

    #print "bases:",bases
    #print "reads:",reads
    print json.dumps(getStats([x for x in lengths if x >= 100]),indent=4)
    
