#!/usr/bin/env python

import sys, json, os
from glob import glob
from optparse import OptionParser

USAGE = """USAGE: %prog <supportFolderOrig> <supportFolderNew> [--options]
Removes all reads that support gaps in supportFolderNew 
that also support gaps supportFolderOrig.
"""

if __name__ == '__main__':
    parser = OptionParser(usage=USAGE)
    parser.add_option("-i","--inplace", type="bool", default=False, \
            help="Overwrite the .gapCan files in <supportFolderNew>")
    opts, args = parser.parse_args()
    
    if len(args) != 2:
        parser.error("Expected 2 arguments")
        
    setAfolder, setBfolder = args
    
    setA = {}
    for file in glob(os.path.join(setAfolder, "*.gapCans")):
        fh = open(file)
        data = json.load(fh)
        for key in data:
            for supType in data[key]:
                for item in data[key][supType]:
                    setA[item] = 1
        fh.close()
    
    print "%d reads used in set A" % (len(setA.keys()))
    removed = 0
    total = 0
    for file in glob(os.path.join(setAfolder, "*.gapCans")):
        fh = open(file)
        output = {}
        data = json.load(fh)
        fh.close()
        for key in data:
            output[key] = {}
            for supType in data[key]:
                output[supType] = []
                for item in data[key][supType]:
                    total += 1
                    try:
                        x = setA[item]
                    except KeyError:
                        output[supType].append(item)
                        removed += 1
        
        if opts.inplace:
            fout = open(file,'w')
        else:
            fout = open(file+".removed",'w')
        
        fout.write(json.dumps(output,indent=4))
        fout.close()
    
    print "Removed %d of %d reads (%.3f%%)" % (removed, total, float(removed)/total)
