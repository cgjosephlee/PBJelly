#!/usr/bin/env python
import sys, os, logging, pickle, re, json
from glob import glob
from optparse import OptionParser
from collections import defaultdict
from FileHandlers import FastaFile, QualFile, GapInfoFile, wrap, LiftOverTable, LiftOverEntry
from Setup import refParser
from StringIO import StringIO

USAGE= ("Collection.py <inputDir> <reference.contigs.fasta> " \
        "<reference.contigs.qual> <gapInfo.bed>")

"""
TODO:

Need to Refactor all this code - 
especially since I've got LiftTable Objects.
"""
def outputNewScaffold(allFilling, contigsFasta, contigsQual, \
                      gapInfoFile, inputDir):
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
        contig = int(contig)
        groupedFasta[scafIndex][contig] = inputFasta[key]
        groupedQual[scafIndex][contig] = inputQual[key]
        originalNames[scafIndex] = origScaf+"|"+scafIndex

    groupedFasta = dict(groupedFasta)
    groupedQual = dict(groupedQual)
    
    #Build LiftTable
    liftTable = LiftOverTable()
    
    for key in groupedFasta:
        if len(groupedFasta[key]) == 1:
            liftEntry = LiftOverEntry(originalNames[key], 0, len(groupedFasta[key][0]),
                                       0, len(groupedFasta[key][0]),
                                       "contig")
            liftTable.addEntry(liftEntry)
            continue
        
        oStart = 0
        for i in xrange(len(groupedFasta[key])-1):
            gapName = key + "_" + str(i) + "_" + str(i+1)
            gap = gapInfo[gapName]
            #make contig
            oEnd = oStart + len(groupedFasta[key][i])
            nStart = oStart
            nEnd = oEnd
            contig = LiftOverEntry(originalNames[key], oStart, oEnd, nStart, nEnd, "contig")
            liftTable.addEntry(contig)
            
            #make gap
            liftGap = LiftOverEntry(originalNames[key], gap.start, gap.end,\
                                    gap.start, gap.end, "gap")
            liftTable.addEntry(liftGap)
            #Shift down.
            oStart = gap.end
        
        #make final contig
        oEnd = oStart + len(groupedFasta[key][i+1])
        nStart = oStart
        nEnd = oEnd
        contig = LiftOverEntry(originalNames[key], oStart, oEnd, nStart, nEnd, "contig")
        liftTable.addEntry(contig)
            
    #Consolidate sequence and update LiftTable
    for key in groupedFasta:
        if len(groupedFasta[key]) == 1:
            logging.info("Skipping Singleton %s" % (key))
            continue
        
        logging.debug("%d entries to parse for %s" % (len(groupedFasta[key])-1, key))
        for i in xrange(len(groupedFasta[key])-1):#pretty sure this is it
            gapName = key + "_" + str(i) + "_" + str(i+1)
            
            #curGap in the liftTable
            curGap = liftTable.getEntry(originalNames[key], gapInfo[gapName].start)

            if not allGapMetrics.has_key(gapName):
                logging.info("Didn't Address %s" % (gapName))
                curGap.gType = "gap_unaddressed"
                groupedFasta[key][i] += "N"*gapInfo[gapName].length
                groupedQual[key][i].extend([0] * gapInfo[gapName].length)
                #nothing to be done
                continue
            
            #filling metrics
            curMetrics = allGapMetrics[gapName]
                        
            #3' trim on previous contig
            if curMetrics.has_key("LeftTrim"):
                trimAmt = curMetrics["LeftTrim"]
                if trimAmt != 0:
                    tEnd = gapInfo[gapName].start
                    tStart = tEnd + trimAmt #Left Trims are negative numbers
                    trim = LiftOverEntry(originalNames[key], tStart, tEnd, 'na', 'na', "trim")
                    
                    #after this entry add a trim
                    liftTable.insertEntry(curGap, trim, after=False)
                    #take the bases off the contig
                    curContig = curGap.getPrev("contig")
                    curContig.nEnd += trimAmt#Shortening
                    groupedFasta[key][i] = groupedFasta[key][i][:trimAmt]
                    groupedQual[key][i] = groupedQual[key][i][:trimAmt]
                    
                    #update the liftTable
                    liftTable.updateScaffold(curContig, trimAmt)
                    logging.info("Trimmed %d from 3' end of %s" % (-trimAmt, key+"_"+str(i)))
                
            #Handle gap based on fill type
            if curMetrics.has_key("NoFillingMetrics"):
                logging.info("Didn't Pass Assembly %s" % (gapName))
                curGap.gType = "gap_failedAssembly"
                groupedFasta[key][i] += "N"*gapInfo[gapName].length
                groupedQual[key][i].extend([0] * gapInfo[gapName].length)

            #Same Safety Net From Assembly -- Shouldn't be necessary
            elif curMetrics["FillSequence"] == "complex":
                logging.info("Complex Fill on %s" % (gapName))
                curGap.gType = "gap_unimprovedC"
                nGapType = "gap_unimprovedC"
                groupedFasta[key][i] += "N"*gapInfo[gapName].length
                groupedQual[key][i].extend([0] * gapInfo[gapName].length)

            elif curMetrics['SpansGap']:
                logging.info("Closed %s" % (gapName))
                curGap.gType = "gap_closed"
                
                shift = curMetrics["FillLength"] - curMetrics["GapPredictedSize"]
                #Gap no longer exists in the assembly
                if curMetrics["FillLength"] != 0:#This happens around trims.
                    newSeq = LiftOverEntry(originalNames[key], 'na', 'na', curGap.nStart,
                                    curGap.nStart + curMetrics["FillLength"],
                                    "new_sequence")
                    liftTable.insertEntry(curGap, newSeq)
                    #update the scaffold for the difference
                    liftTable.updateScaffold(newSeq, shift)
                else:
                    liftTable.updateScaffold(curGap, shift)
                curGap.nStart = 'na'
                curGap.nEnd = 'na'
                                
                #add the sequence to the current contig
                groupedFasta[key][i] += curMetrics["FillSequence"]
                groupedQual[key][i].extend(curMetrics["FillQual"])

            elif curMetrics["FillSequence"] != "":
                fillSeq = curMetrics["FillSequence"]
                fillQual = curMetrics["FillQual"]
                
                if curMetrics.has_key("GapUnderestimated"):
                    logging.info("Overfilled %s" % (gapName))
                    curGap.gType = "gap_overfilled"
                    #shift all but what we leave for overfill
                    remainingSize = 25
                else:
                    logging.info("Reduced %s" % (gapName))
                    curGap.gType = "gap_reduced"
                    remainingSize = curMetrics["GapPredictedSize"] - \
                                    curMetrics["FillLength"]
                    if remainingSize == 0:
                        remainingSize = 25
                
                #Tell everything Downstream that there is no gap
                liftTable.updateScaffold(curGap, -curMetrics["GapPredictedSize"])
                #then we'll add the new stuff and push then where they should be
                
                #add the new 5' sequence
                fiveLen = fillSeq.index('N')
                fiveSeq = LiftOverEntry(originalNames[key], 'na', 'na', curGap.nStart, \
                                     curGap.nStart + fiveLen, "new_sequence")
                liftTable.insertEntry(curGap, fiveSeq, after=False)
                
                #update the gap's size
                curGap.nStart = curGap.prev.nEnd
                curGap.nEnd = curGap.nStart + remainingSize
                #nEnd = curGap.nStart + remainingSize
                #shift = nEnd - curGap.nEnd
                #curGap.nEnd = nEnd
                #liftTable.updateScaffold(curGap, shift)
                
                #add the new 3' sequence
                threeLen = len(fillSeq) - fillSeq.rindex('N') - 1
                threeSeq = LiftOverEntry(originalNames[key], 'na', 'na', curGap.nEnd, \
                                        curGap.nEnd + threeLen, "new_sequence")
                liftTable.insertEntry(curGap, threeSeq)
                
                #Then do the shifting
                liftTable.updateScaffold(threeSeq, fiveLen + remainingSize + threeLen)
                #liftTable.updateScaffold(curGap,remainingSize)
                #liftTable.updateScaffold(threeSeq, threeLen)
                #add the sequence to the current contig
                groupedFasta[key][i] += fillSeq
                groupedQual[key][i].extend(fillQual)
            else:
                logging.warning("EdgeCase!? " + gapName)
            
            #5' trim on the next contig
            if curMetrics.has_key("RightTrim"):
                trimAmt = curMetrics["RightTrim"]
                if trimAmt != 0:
                    tStart = gapInfo[gapName].end
                    tEnd = tStart + trimAmt
                    trim = LiftOverEntry(originalNames[key], tStart, tEnd, 'na', 'na', "trim")
                
                    #after this gap add a trim.
                    liftTable.insertEntry(curGap, trim)
                    #take the bases off the contig
                    curContig = curGap.getNext("contig")
                    curContig.nStart += trimAmt#Shortening
                    groupedFasta[key][i+1] = groupedFasta[key][i+1][trimAmt:]
                    groupedQual[key][i+1] = groupedQual[key][i+1][trimAmt:]
                
                    #update the liftTable - prev because I want to shift this contig, too.
                    liftTable.updateScaffold(curContig.prev, -trimAmt)
                    logging.info("Trimmed %d from 5' end of %s" % (trimAmt, key+"_"+str(i+1))) 
    
    #Finally output the new information
    logging.info("Writing Fasta & Qual")
    outputFasta = open(os.path.join(inputDir,"jellyOutput.fasta"),'w')
    outputQual  = open(os.path.join(inputDir,"jellyOutput.qual"),'w')
        
    for key in groupedFasta:
        outputFasta.write(">"+originalNames[key]+"\n")
        outputQual.write(">"+originalNames[key]+"\n")
        if len(groupedFasta[key]) == 1:#Fix indexing
            outputFasta.write(wrap(groupedFasta[key][0]) + "\n")
            outputQual.write(" ".join(map(str,groupedQual[key][0])) + "\n")
            continue
        
        for i in xrange(len(groupedFasta[key])):
            outputFasta.write(wrap(groupedFasta[key][i]))
            outputQual.write(" ".join(map(str,groupedQual[key][i])) + ' ')
        #Hackiest thing every. We need to clean up outputting qual information
        outputFasta.write('\n')
        outputQual.seek(outputQual.tell()-1)
        outputQual.write('\n')
    outputFasta.close()
    outputQual.close()

    logging.info("Writing LiftTable")
    outTable  = open(os.path.join(inputDir,"jellyLiftOverTable.txt"),'w')
    outTable.write("#scaffoldName\toStart\toEnd\tnStart\tnEnd\tfeatureType\n")
    outTable.write(str(liftTable))
    outTable.close()
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

def __initLog(debug):
    """Logging"""
    logLevel = logging.DEBUG if debug else logging.INFO
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
    inputDir = args[0]
    contigsFasta = args[1]
    contigsQual = args[2]
    gapInfoFile = args[3]
    
    return inputDir, contigsFasta, contigsQual, gapInfoFile, opts.debug

if __name__ == '__main__':
    inputDir, contigsFasta, contigsQual, gapInfoFile, debug = parseOpts()
    __initLog(debug)
    
    allFilling = statsCollector(inputDir)
    outputNewScaffold(allFilling, contigsFasta, contigsQual,gapInfoFile, inputDir)
    
