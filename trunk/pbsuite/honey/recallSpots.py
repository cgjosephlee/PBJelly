#!/usr/bin/env python

import sys, logging, argparse
import pysam
from pbsuite.honey.Honey import *
from pbsuite.utils.setupLogging import setupLogging
USAGE = "Recall spots in a hon.h5 file"

parser = argparse.ArgumentParser(description=USAGE, \
        formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("bam", metavar="BAM", type=str, \
                    help="Bam containing original reads")
parser.add_argument("hon", metavar="HON", type=str, \
                    help="Hon.h5 containing error counts")
parser.add_argument("-c", "--minCoverage", type=int, default=3, \
                    help="Minimum coverage for a spot call to be made (3)")
                    
parser.add_argument("-e", "--threshold",  type=float, default=2,
                    help="Spot Threshold (2))")
                    
parser.add_argument("-b", "--binsize", type=int, default=50, \
                    help="binsize for window averaging (50)")
    
parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
                    
args = parser.parse_args()

setupLogging()

f = h5py.File(args.hon,'a')
bam = pysam.Samfile(args.bam)
tsp = 0
for chrom in f.keys():
    logging.info("Calling %s" % (chrom))
    container = f[chrom]["data"]
    spots = callHotSpots(container, args.threshold, args.minCoverage, args.binsize, f[chrom].attrs["start"])
    print "\n".join(map(str, spots))
    logging.info("filtering insz spots")
    fspot = 0
    for spot in spots:
        spot.chrom = chrom
        if spot.tags["label"] == "INSZ" and not filterINSZ(bam, spot, 50, 10, 0.33):
            continue
        fspot += 1 
        print spot
    tsp += fspot
    logging.info("found %d spots" % (fspot))
logging.info("Finished %d spots" % (tsp))
