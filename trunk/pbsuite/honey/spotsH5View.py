#!/usr/bin/env python
import sys
import h5py

#view with chr:start-end

if __name__ == '__main__':
    h5 = h5py.File(sys.argv[1])
    #reg = sys.argv[2]
    #d = reg.split(':')
    #if len(d) == 1:
        #chrom = d[0]
        #start = 0
        #end = None
    #elif len(d) == 2:
        #chrom = d[0]
        #e = d[1].split('-')
        #if len(e) != 2:
            #print "incorrect region... expected chrom:start-end or just chrom"
            #exit(0)
        #start, end = e
        #start = int(start); end = int(end)
    #print chrom, start, end
    print "chrom\tposition\t"+"\t".join(h5.attrs["columns"])
    for chrom in h5.keys():
        for pos,i in enumerate(h5["/%s/data" % chrom].value.transpose().tolist()):
            print chrom + "\t" + str(pos) + '\t' + "\t".join(map(str, i))
