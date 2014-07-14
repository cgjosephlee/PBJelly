#!/usr/bin/env python

import argparse
from pbsuite.utils.setupLogging import *
from pbsuite.honey import bampie, TGraf, HSpots, Force

STAGES = {"pie":   bampie.run, \
          "tails": TGraf.run, \
          "spots": HSpots.run, \
          "force": Force.run }

USAGE = """\
   Honey - genomic variant calling with long sequencing reads

   STAGE is one of
     pie        Extract and map soft-clipped Tails from a bam.
     tails      Cluster mapped tails to make break-points of larger events.
     spots      Find genomic variants within reads' spans.
     force      Given a BedFile of predicted variants, force search for matching
   
   See HoneyReadme.txt for documentation or --help for details\
"""

def parseArgs():
    parser = argparse.ArgumentParser(prog="Honey.py", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("stage", metavar="STAGE", choices=STAGES.keys(), type=str, \
                        help="Stage to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,\
                         help="Options to pass to the stage")

    args = parser.parse_args()

    STAGES[args.stage](args.options)
    
if __name__ == '__main__':
    parseArgs()