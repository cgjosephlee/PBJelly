#!/usr/bin/env python
import sys, os, logging, pickle, re, json
from glob import glob
from optparse import OptionParser
from collections import defaultdict
from FileHandlers import FastaFile, QualFile, GapInfoFile, wrap
from Setup import refParser
from StringIO import StringIO

USAGE= ("Collection.py <inputDir> <reference.contigs.fasta> " \
        "<reference.contigs.qual> <gapInfo.bed>")
DEBUGON = False

"""
TODO:

Need to Refactor all this code - 
especially since I've got LiftTable Objects.
"""
def outputNewScaffold(allFilling, contigsFasta, \
                      contigsQual, gapInfoFile, inputDir):
    logging.info("Loading Fasta/Qual Files")
    inputFasta = FastaFile( contigsFasta )
    inputQual = QualFile( contigsQual)
    gapInfo = GapInfoFile(gapInfoFile)
    
    allGapMetrics = {}
    #Clean off absolute path in gap name
    for n in allFilling:
        allGapMetrics[n.split('/')[-1]] = allFilling[n]
    
    logging.info("Finished Loading Files")
    
    originalNames = defaultdict()
    groupedFasta= defaultdict(dict)
    groupedQual = defaultdict(dict)
        
    #Group scaffolds
    for key in inputFasta.keys():
        origScaf, scafIndex, contig = refParser.match(key).groups()
        groupedFasta[scafIndex][contig] = inputFasta[key]
        groupedQual[scafIndex][contig] = inputQual[key]
        originalNames[scafIndex] = origScaf+"|"+scafIndex
    
    groupedFasta = dict(groupedFasta)
    groupedQual = dict(groupedQual)
    
    outputFasta = open(os.path.join(inputDir,"jellyOutput.fasta"),'w')
    outputQual  = open(os.path.join(inputDir,"jellyOutput.qual"),'w')
    shiftTable  = open(os.path.join(inputDir,"jellyLiftOverTable.txt"),'w')
    shiftTable.write("#scaffoldName\toStart\toEnd\tnStart\tnEnd\tFeatureType\n")

    for key in groupedFasta:
        outputFasta.write(">"+originalNames[key]+"\n")
        outputQual.write(">"+originalNames[key]+"\n")
        
        curScafFasta = StringIO()
        curScafQual = StringIO()
        
        if len(groupedFasta[key]) == 1:
            outputFasta.write(wrap(groupedFasta[key]['1'])+"\n")
            outputQual.write(" ".join(map(str,groupedQual[key]['1']))+"\n")
            shiftTable.write("%s\t0\t%d\t0\t%d\t%s\n" % (originalNames[key], \
                             len(groupedFasta[key]['1']), \
                             len(groupedFasta[key]['1']), \
                             "contig") )
            continue
        
        oStart = 0
        oEnd = None
        nStart = 0
        nEnd = None
        allShift = 0 #Where things have moved
        for i in xrange(len(groupedFasta[key])-1):
            gapName = key + "_" + str(i) + "_" + str(i+1)
            
            oEnd = oStart + len(groupedFasta[key][str(i)])
            #Shift based on what we've seen up to here
            nStart = oStart + allShift
            nEnd = oEnd + allShift
            
            leftOfGapTrim = ""
            rightOfGapTrim = ""
            
            fivePrimeTrim = 0#how much is this contig's 5' trimmed
            threePrimeTrim = 0#how much is this contig's 3' trimmed
            
            prevGap = key + "_" + str(i-1) + "_" + str(i)
            
            if allGapMetrics.has_key(prevGap) and allGapMetrics[prevGap].has_key("RightTrim"):
                #Trim this contig's 5' using previous gap's rightTrim information
                fivePrimeTrim = allGapMetrics[prevGap]["RightTrim"]
            
            if allGapMetrics.has_key(gapName) and allGapMetrics[gapName].has_key("LeftTrim"):
                threePrimeTrim = allGapMetrics[gapName]["LeftTrim"]

            gapSize = gapInfo[gapName].length
            
            if not allGapMetrics.has_key(gapName):
                logging.info("Didn't Address %s" % (gapName))
                fillSequence = "N" * gapSize
                fillQual = "0 " * gapSize
                nGapType = "gap_unaddressed"
                fillShift = 0
            elif allGapMetrics[gapName].has_key("NoFillingMetrics"):
                logging.info("Didn't Pass Assembly %s" % (gapName))
                fillSequence = "N" * gapSize
                fillQual = "0 " * gapSize
                nGapType = "gap_unimproved"
                fillShift = 0
            #Same Safety Net From Assembly -- Shouldn't be necessary
            elif allGapMetrics[gapName]["FillSequence"] == "complex":
                logging.info("Complex Fill on %s" % (gapName))
                fillSequence = "N" * gapSize
                fillQual = "0 " * gapSize
                nGapType = "gap_unimprovedC"
                fillShift = 0
            else:
                logging.info("Improved %s" % (gapName))
                fillSequence = allGapMetrics[gapName]["FillSequence"]
                fillQual = allGapMetrics[gapName]["FillQual"] 
                if allGapMetrics[gapName]['SpansGap']:
                    nGapType = "gap_filled"
                elif allGapMetrics[gapName].has_key("GapUnderestimated"):
                    nGapType = "gap_overfilled"
                else:
                    nGapType = "gap_reduced"
                fillShift = len(fillSequence) - gapSize
                
            """Here I know
            oStart, oEnd
            nStart, nEnd 
            if there are trims:
                fivePrimeTrim Amount & threePrimeTrim Amount
            nGapType - myFilling Type
            fillShift - how much what I put in is different from what was already there (Ns to non-Ns)
            """
            #write this contig's 5' trim
            if fivePrimeTrim != 0:
                shiftTable.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (originalNames[key], \
                             oStart, \
                             oStart + fivePrimeTrim, \
                             "na", \
                             "na", \
                             "trim"))
            #write contig
            shiftTable.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (originalNames[key], \
                             oStart, \
                             oEnd, \
                             nStart - fivePrimeTrim, \
                             nEnd - fivePrimeTrim + threePrimeTrim, \
                             "contig"))
            #write this contig's 3' trim
            if threePrimeTrim != 0:
                shiftTable.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (originalNames[key], \
                             oEnd + threePrimeTrim, \
                             oEnd, \
                             "na", \
                             "na", \
                             "trim"))
                threePrimeIndex = threePrimeTrim
            else:
                threePrimeIndex = len(groupedFasta[key][str(i)])

            #write the gap we're iterating
            shiftTable.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (originalNames[key], \
                             gapInfo[gapName].start, \
                             gapInfo[gapName].end, \
                             gapInfo[gapName].start + allShift - fivePrimeTrim + threePrimeTrim, \
                             gapInfo[gapName].end + allShift - fivePrimeTrim + threePrimeTrim + fillShift, \
                             nGapType))
            
            #output sequence
            curScafFasta.write(groupedFasta[key][str(i)][fivePrimeTrim:threePrimeIndex])
            curScafFasta.write(fillSequence)
            
            curScafQual.write(" ".join(map(str,groupedQual[key][str(i)][fivePrimeTrim:threePrimeIndex])))
            curScafQual.write(" " + " ".join(map(str,fillQual)) + " ")
            #output gap
            #update shifts
            oStart = gapInfo[gapName].end
            allShift += fillShift - fivePrimeTrim + threePrimeTrim

    
        #Finish Last Contig

        fivePrimeTrim = 0
        oEnd = oStart + len(groupedFasta[key][str(i+1)])
        nStart = oStart + allShift
        nEnd = oEnd + allShift
        #write trim
        if allGapMetrics.has_key(gapName) and allGapMetrics[gapName].has_key("RightTrim"):
            fivePrimeTrim = allGapMetrics[gapName]["RightTrim"]
            shiftTable.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (originalNames[key], \
                             oStart, \
                             oStart + fivePrimeTrim, \
                             "na", \
                             "na", \
                             "trim")) 
        #write contig
        shiftTable.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (originalNames[key], \
                         oStart, \
                         oEnd, \
                         nStart - fivePrimeTrim, \
                         nEnd - fivePrimeTrim, \
                         "contig"))
        curScafFasta.write(groupedFasta[key][str(i+1)][fivePrimeTrim:])
        outputFasta.write(wrap(curScafFasta.getvalue())+"\n")
        curScafQual.write(" ".join(map(str,groupedQual[key][str(i+1)][fivePrimeTrim:])))
        outputQual.write(curScafQual.getvalue()+"\n")
        logging.info("Finished Writing GroupedFasta")
    
    shiftTable.close()
    outputFasta.close()
    outputQual.close()
    logging.info("Finished Output")

def statsCollector(inputDir):
    folder = glob(os.path.join(inputDir,"ref*"))
    
    all = {}
    
    for f in folder:
        try:
            fh = open(os.path.join(f,"fillingMetrics.json"),'r')
        except IOError:
            all[f] = {"NoFillingMetrics":None}
            continue
        try:
            all[f] = json.load(fh)
            all[f]["FillQual"] = all[f]["FillQual"].split(' ')
        except ValueError:
            logging.error("WARNING! "+f+" didn't produce a valid JSON output in " + \
                          "fillingMetrics.json. Go check this file and if it is " + \
                          "not a plain text JSON file, try re-running the " + \
                          "assembly Process on this folder")
            exit(1)
        fh.close()
    
    noFillingMetrics = 0
    filled = 0
    overfilled = 0
    reduced = 0
    
    for i in all:
        if all[i].has_key("NoFillingMetrics"):
            noFillingMetrics += 1
            continue
        if all[i]["SpansGap"]:
            filled += 1
        elif all[i]["LeftContig"] or all[i]["RightContig"]:
            reduced += 1
            if all[i].has_key("GapUnderestimated"):
                overfilled += 1
    
    logging.info("Number of Gaps Addressed %d" % (len(all.keys())))
    logging.info("No Filling Metrics %d" % (noFillingMetrics))
    logging.info("Filled %d" % (filled))
    logging.info("Reduced %d" % (reduced))
    logging.info("Overfilled %d" % (overfilled))
    
    return all

def __initLog():
    """Logging"""
    logLevel = logging.DEBUG if DEBUGON else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))

def parseOpts():
    parser = OptionParser(USAGE)
    parser.add_option("--debug", action="store_true", 
                      help="Increases verbosity of logging")
    opts, args = parser.parse_args()
    
    if len(args) != 4:
        parser.error("Error! Incorrect number of arguments")
    DEBUGON = opts.debug
    inputDir = args[0]
    contigsFasta = args[1]
    contigsQual = args[2]
    gapInfoFile = args[3]
    
    return inputDir, contigsFasta, contigsQual, gapInfoFile

if __name__ == '__main__':
    inputDir, contigsFasta, contigsQual, gapInfoFile = parseOpts()
    __initLog()
    
    allFilling = statsCollector(inputDir)
    outputNewScaffold(allFilling, contigsFasta, contigsQual,gapInfoFile, inputDir)
    
