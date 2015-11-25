#!/usr/bin/env python
import os
import sys
import shutil
import random
import argparse
import tempfile
import multiprocessing
from collections import namedtuple

import pysam

from pbsuite.utils.BedIO import *
from pbsuite.utils.FileHandlers import *
from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import *
from pbsuite.honey.asmEngines import Assembler, MiniaAssembler, PhrapAssembler, SpadesAssembler

USAGE = """\
Takes a list of putative SVs and a set of BAMs builds
as many of your sites as possible.
Remap Files are stored in
<--output>.sam
<--output>.tails.sam
"""

CRITICAL = """\
Need to separate arguments
Need to make it explicit that Celera is ran blandly
Only the first bam's insert size is considered
"""

class Consumer(multiprocessing.Process):

    def __init__(self, task_queue, result_queue, args):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.args = args

    def run(self):
        try:
            proc_name = self.name
            #Open the bams
            nBams = [] # nonTrim Bams
            tBams = [] # trim Bams
            for i in self.args.bam:
                nBams.append(pysam.Samfile(i))
            for i in self.args.pacBam:
                tBams.append(pysam.Samfile(i))

            while True:
                next_task = self.task_queue.get()
                if next_task is None:
                    # Poison pill means shutdown
                    logging.info('Thread %s: Exiting\n' % proc_name)
                    self.task_queue.task_done()
                    break
                try:
                    answer = next_task(nBams, tBams)
                except Exception as e:
                    logging.error("Exception raised in task %s" % (str(e)))
                    self.task_queue.task_done()
                    self.result_queue.put("Failure - UNK - %s" % str(e))
                    logging.info("fail in groupid=%s" % next_task.data.name)
                    continue
                self.task_queue.task_done()
                self.result_queue.put(answer)
            return
        except Exception as e:
            logging.error("Consumer %s Died\nERROR: %s" % (self.name, e))
            return

def insertDist(bam):
    """
    Samples reads to get mean insert size and standard deviation
    """
    num_samp = 1000000
    counter = 0
    skip = 5000000
    skip_counter = 0
    ins_list = []
    for read in bam.fetch():
        if skip_counter < skip:
            skip_counter += 1
            continue
        if read.is_proper_pair and not read.is_reverse and not read.is_secondary:
            ins_list.append(read.tlen)
            counter += 1
        if counter == num_samp:
            break
    mean = sum(ins_list)/float(len(ins_list))
    v = 0
    for i in ins_list:
        v += (i-mean)**2
    variance = v/float(len(ins_list))
    stdev = variance**(0.5)
    return (mean, stdev)

def parseArgs(argv):
    parser = argparse.ArgumentParser(description=USAGE, \
                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("putative", metavar="BED", type=str, \
                        help="Bed of regions to assemble")
    parser.add_argument("-b", "--bam", type=str, nargs="*", \
                        help="Input Bam (NonTrim)")
    parser.add_argument("-p", "--pacBam", type=str, nargs="*", \
                        help="PacBio Bam")
    parser.add_argument("-a", "--assembler", type=str, default='phrap', choices=["phrap", "minia", "spades"],
                        help="Assembly program to use (%(default)s)")
    parser.add_argument("-B", "--buffer", type=int, default=1000, \
                        help="Amount of buffer sequence around the variant to use (%(default)s)")
    parser.add_argument("-n", "--nproc", type=int, default=1, \
                        help="Number of processors to use (%(default)s)")
    parser.add_argument("-o", "--output", default="asm.fastq",\
                        help="Where to write the resultant assemblies (%(default)s)")
    parser.add_argument("-r", "--reference", default=None, \
                        help="Reference to map to (optional if --noRemap)")
    parser.add_argument("--noRemap", action="store_false", \
                        help="Do not remap assembly")
    parser.add_argument("--noSplitMap", action="store_false", \
                        help="Do not map tails from remapped assembly (off if --noRemap)")
    parser.add_argument("--timeout", type=int, default=30, \
                        help="Timeout assembly after N minutes (%(default)s)")
    parser.add_argument("--maxspan", type=int, default=100000, \
                        help="Maximum Span of SV to attempt assembling (%(default)s)")
    parser.add_argument("--maxreads", type=int, default=500, \
                        help="Maximum number of Illumina reads used to attempt assembling (%(default)s)")
    parser.add_argument("--temp", type=str, default=tempfile.gettempdir(),
                            help="Where to save temporary files")
    parser.add_argument("--start", type=int, default=0,
                        help="Index of the first variant to begin assembling. (%(default)s)")
    parser.add_argument("--stride", type=int, default=1,
                        help="Assemble one every N reads (%(default)s)")
    parser.add_argument("--debug", action="store_true",\
                        help="Verbose Logging")

    #parser.add_argument("--insertsize", type=int, default=None, \
                        #help=("Celera - insert size for PE Illumina reads (auto_detect)"))
    #parser.add_argument("--insertstd", type=float, default=None, \
                        #help=("Celera - insert std for PE Illumina reads (auto_detect)"))

    args = parser.parse_args(argv)
    setupLogging(args.debug)

    # Parameter checks
    if args.bam is None and args.pacBam is None:
        logging.error("Expected at least one BAM argument")
        exit(1)

    if not args.output.endswith(".fastq"):
        logging.error("Output needs to end with .fastq")
        exit(1)

    if not os.path.exists(args.putative):
        logging.error("Input {inp} does not exist".format(inp=args.putative))
        exit(1)

    if args.noRemap and args.reference == None:
        logging.error("Cannot remap without --reference")
        exit(1)

    if args.reference and not os.path.exists(args.reference):
        logging.error("Reference {ref} does not exist".format(ref=args.reference))
        exit(1)

    if args.bam is None:
        args.bam = []
        #if args.insertsize is None and args.bam is not None:
            #j = pysam.Samfile(args.bam[0])
            #mu,std = insertDist(j)
            #j.close()
            #args.insertsize = mu
            #args.insertstd = std if args.insertstd is None else args.insertstd

    if args.pacBam is None:
        args.pacBam = []

    return args


def run(argv):
    args = parseArgs(argv)
    # Establish communication queues
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    # Start consumers
    num_consumers = args.nproc
    consumers = [ Consumer(tasks, results, args)
                    for i in range(num_consumers) ]
    for w in consumers:
        w.start()

    # Enqueue jobs
    num_jobs = 0

    #open the putative file
    bed = BedFile.fromFile(args.putative)

    entryPos = args.start

    ### Set the Assembler
    if args.assembler == 'phrap':
        MyAssembler = PhrapAssembler.PhrapAssembler
    elif args.assembler == 'minia':
        MyAssembler = MiniaAssembler.MiniaAssembler
    elif args.assembler == 'spades':
        MyAssembler = SpadesAssembler.SpadesAssembler
    elif self.assembler == 'dipspades':
        raise NotImplementedError("Working on it")
    else:
        logging.error("%s - Not a valid assembler. Seek --help" % (assembler))
        exit(1)

    eout = open(args.output + ".err", 'w')
    while True:
        if entryPos >= len(bed):
            break

        entry = bed[entryPos]

        if abs(entry.start - entry.end) > args.maxspan:
            eout.write("Too Big %s %d\n" % (entry.name, abs(entry.start-entry.end)))
        else:
            tasks.put(MyAssembler(entry, args))
            num_jobs += 1

        #stride over
        #for i in range(args.stride):
            #entryPos += 1
        entryPos += args.stride


    # Add a poison pill for each consumer
    for i in range(num_consumers):
        tasks.put(None)

    logging.info("%d Tasks" % (num_jobs))
    # Wait for all of the tasks to finish -- I might not need to
    #logging.info("Joining")
    #tasks.join()

    fout = open(args.output, 'w')
    #tracking who failed
    asmFails = 0
    # Consolidate results
    logging.info("Consolidating")
    while num_jobs:
        result = results.get()
        num_jobs -= 1

        if result.startswith("Failure "):
            logging.error(result)
            eout.write(result + '\n')
            asmFails += 1
            continue

        fh = open(result,'r')
        fout.write(fh.read())
        fh.close()

        os.remove(result)

    eout.close()
    logging.error("%d assemblies failed" % (asmFails))
    fout.close()
    #for each, I'm going to need to
    #Grabbing my reads
    if args.noRemap:
        logging.info("PIE mapping")
        r, o, e = exe(("Honey.py pie --temp {tempDir} --nproc {np} "
                       "--minTail 50 {reads} {ref}" \
                       .format(tempDir=args.temp, reads=args.output, \
                               ref=args.reference, np=args.nproc)))
        if r != 0:
            logging.error("Honey pie quit (%d)\n%s\n%s" % (r, o, str(e)))
            exit(r)
        logging.info("Honey Log:\n" + o)

    logging.info("Finished")



if __name__ == '__main__':
    args = parseArgs(sys.argv[1:])
    takeMassivePhrap(args)
