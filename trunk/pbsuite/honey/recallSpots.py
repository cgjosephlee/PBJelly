#!/usr/bin/env python

import sys, logging, argparse
import pysam
from pbsuite.honey.HSpots import *
from pbsuite.utils.setupLogging import setupLogging
USAGE = "Recall spots in a hon.h5 file"


args = parseArgs(sys.argv[1:], established=True)

setupLogging(args.debug)

f = h5py.File(args.hon,'a')
bam = pysam.Samfile(args.bam)
tsp = 0
makeKernals(args.binsize)
print "#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO"
for chrom in f.keys():
    logging.info("Calling %s" % (chrom))
    container = f[chrom]["data"]
    spots = callHotSpots(container, f[chrom].attrs["start"], args)
    logging.info("Filtering INSZ Spots")
    fspot = 0
    for spot in spots:
        spot.chrom = chrom
        spot.offset(start)
        if spot.tags["label"] == "INS" and filterINSZ(bam, spot, args):
            continue
        if (spot.size < args.sizeMin or spot.size == -1) or spot.size > args.sizeMax:
            continue
        fspot += 1 
        if groupName != chrom:
            spot.tags["RegName"] = groupName
        print spot
    tsp += fspot
    logging.info("Found %d spots" % (fspot))
logging.info("Finished %d spots" % (tsp))
