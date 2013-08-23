#!/usr/bin/env python

import sys, logging, argparse
from pbsuite.honey.Honey import *
from pbsuite.utils.setupLogging import setupLogging
USAGE = "Recall spots in a hon.h5 file"

parser = argparse.ArgumentParser(description=USAGE, \
        formatter_class=argparse.RawDescriptionHelpFormatter)

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
tsp = 0
for chrom in f.keys():
    logging.info("Calling %s" % (chrom))
    container = f[chrom]["data"]
    points = callHotSpots(container, args.threshold, args.minCoverage, args.binsize, f[chrom].attrs["start"])
    cnt = len(points)
    logging.info("%d spots found" % (cnt))
    tsp += cnt
    for i in points:
        print spotToString(chrom,i)
logging.info("Finished %d spots" % (tsp))
