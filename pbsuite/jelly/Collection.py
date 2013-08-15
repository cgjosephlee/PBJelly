#!/usr/bin/env python

import sys, json
from glob import glob
from optparse import OptionParser

from pbsuite.utils.setupLogging import *
from pbsuite.utils.FileHandlers import *
from pbsuite.jelly.Jelly import JellyProtocol
from pbsuite.jelly.Support import SUPPORTFLAGS

import networkx as nx

USAGE= "Collection.py <protocol.xml>"

def makeFilMetName(a,b):
    """
    
    """
    a = a.replace('/','.')
    b = b.replace('/','.')
    j = [a,b]
    j.sort()
    return "_".join(j)

class FillingMetrics():
    
    def __init__(self, data, gapName):
        self.data = data
        self.gapName = gapName
        self.__parseData()

    def __parseData(self):
        logging.debug(self.gapName)
        data = self.data
        gapName = self.gapName
        g = gapName.split('_')
        if len(g) == 1:
            a = g[0]
            b = "-"
        else:
            a,b = g
        self.span = False
        self.singleExtend = False
        self.doubleExtend = False
        self.seed1Name = a
        self.seed1Strand = None
        self.seed2Name = b
        self.seed2Strand = None
        self.fillLength = 0
        self.seed1ExtendSeq = FastqEntry(None, "", "")
        self.seed2ExtendSeq = FastqEntry(None,"", "")
        self.sameStrand = self.seed1Name[-1] != self.seed2Name[-1]
        if data.has_key("predictedGapSize"):
            self.predictedGapSize = data["predictedGapSize"]
        else:
            self.predictedGapSize = None
        
        #Setting support types
        #Setting Strands
        if SUPPORTFLAGS.span in data["support"][0] and data["fillSeq"] is not None:
            self.span = True
            self.seed1Strand = '+' if data["spanSeedStrand1"] == '0' else '-'
            self.seed2Strand = '+' if data["spanSeedStrand2"] == '0' else '-'
        else:
            if data["extendF1Count"] != 0 and data["extendF2Count"] != 0:
                self.doubleExtend = True
                self.seed1Strand = '+' if data["extendF1SeedStrand"] == '0' else '-'
                self.seed2Strand = '+' if data["extendF2SeedStrand"] == '0' else '-'
                if data["extendSeq1"] != None:
                    self.seed1ExtendSeq = \
                        FastqEntry(b, data["extendSeq1"], "?"*len(data["extendSeq1"]))
                else:
                    self.seed1ExtenSeq = FastqEntry(b, "", "")
                if data["extendSeq2"] != None:
                    self.seed2ExtendSeq = \
                        FastqEntry(b, data["extendSeq2"], "?"*len(data["extendSeq2"]))
                else:
                    self.seed2ExtendSeq = FastqEntry(b, "", "")

            else:
                self.singleExtend = True
                if data["extendSeq1"] is not None:
                    self.seed1ExtendSeq = \
                        FastqEntry(b, data["extendSeq1"], "?"*len(data["extendSeq1"]))

                    self.seed1Strand = '+' if data["extendF1SeedStrand"] == '0' else '-'
                if data["extendSeq2"] is not None:
                    self.seed2ExtendSeq = \
                        FastqEntry(b, data["extendSeq2"], "?"*len(data["extendSeq2"]))

                    self.seed2Strand = '+' if data["extendF2SeedStrand"] == '0' else '-'
        self.fillLength  = data["fillBases"]
                            
    def getSequence(self):
        """
        Returns the full sequence this metric holds
        If it's a span, you'll get a full sequence.
        If it's a gap reduced, You'll get a single sequence
        with a gap ('N') between the pieces 

        Note:
            all gaps < 25bp (overfilled or otherwise) will be inflated
            to 25bp

        returns None if there is any issue with the metrics
        """
        logging.debug("Getting fill sequence for %s" % (self.gapName))
        if self.span:
            logging.debug('fill span')
            return FastqEntry(self.gapName, self.data["fillSeq"], "?"*len(self.data["fillSeq"]))
        
        if self.predictedGapSize is None:
            #We can't reduce a gap of unknown size
            #one would just get the extend seq
            logging.debug("No predicted gap size")
            return None
        
        gapLen = self.predictedGapSize - self.fillLength
        
        #No improvement at all
        if not (self.span or self.singleExtend or self.doubleExtend):
            logging.debug("Unimproved Gap - %s" % (self.gapName))
            return FastqEntry(self.gapName, 'N'*gapLen, '!'*gapLen)
            
        if gapLen < 25:
            gapLen = 25
        
        #Single end extension
        if self.singleExtend:
            if self.seed1Strand == '-':
                self.seed1ExtendSeq.reverseCompliment()

            logging.debug('single extend')
            return FastqEntry(self.gapName, 
                        self.seed1ExtendSeq.seq + \
                        ('N'*gapLen), \
                        self.seed1ExtendSeq.qual + \
                        ('!'*gapLen))
        
        #Stick them together!
        #Fill Sequence is on same strand as it should be
        if self.sameStrand and self.seed1Strand == self.seed2Strand:
            logging.debug('same strand success')
            if self.seed1Name.endswith('e5'):
                logging.debug('seed 1 first')
                return FastqEntry(self.gapName,
                            self.seed1ExtendSeq.seq + \
                            ('N'*gapLen) + \
                            self.seed2ExtendSeq.seq, \
                            self.seed1ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed2ExtendSeq.qual)
            elif self.seed2Name.endswith('e5'):
                logging.debug('seed 2 first')
                return FastqEntry(self.gapName,
                            self.seed2ExtendSeq.seq + \
                            ('N'*gapLen) + \
                            self.seed1ExtendSeq.seq, \
                            self.seed2ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed1ExtendSeq.qual)
                return self.seed2ExtendSeq + self.seed2ExtendSeq
            else:
                logging.error(("Huge Problem! This Should Never Happen!  "
                               "sameStrand strandsEqual"))
                exit(10)
        
        #one fill sequence needs to be flipped.
        elif self.sameStrand and self.seed1Strand != self.seed2Strand:
            logging.debug('same strand success - flip one')
            if self.seed1Name.endswith('e5'):
                logging.debug('seed 1 first, flip 2')
                self.seed2ExtendSeq.reverseCompliment()
                return FastqEntry(self.gapName,
                            self.seed1ExtendSeq.seq + \
                            ('N'*gapLen) + \
                            self.seed2ExtendSeq.seq, \
                            self.seed1ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed2ExtendSeq.qual)
            elif self.seed2Name.endswith('e5'):
                logging.debug('seed 2 first, flip 1')
                self.seed1ExtendSeq.reverseCompliment()
                return FastqEntry(self.gapName,
                            self.seed2ExtendSeq.seq + \
                            ('N'*gapLen) + \
                            self.seed1ExtendSeq.seq, \
                            self.seed2ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed1ExtendSeq.qual)
            else:
                logging.error(("Huge Problem! This Should Never Happen!  "
                               "sameStrand !strandsEqual"))
                exit(10)
        
        #One Sequence may be None and not sameStrand
        elif not self.sameStrand and (self.seed1Strand is None or self.seed2Strand is None):
            logging.debug("not same strand - one Null")
            if self.seed1Strand is None:
                logging.debug("seed1 None")
                #5' needs to be extended upstream -- 
                if self.seed2Name.endswith('e5'):
                    logging.debug("5' extend upstream")
                    return FastqEntry(self.gapName, \
                                ('N'*gapLen) + \
                                self.seed2ExtendSeq.seq, \
                                ('!'*gapLen) + \
                                self.seed2ExtendSeq.qual)
                else:
                    #3' needs to be extended downstream
                    logging.debug("3' extend upstream")
                    return FastqEntry(self.gapName, \
                                self.seed2ExtendSeq.seq + \
                                ('N'*gapLen), \
                                self.seed2ExtendSeq.qual + \
                                ('!'*gapLen)) 
            
            elif self.seed2Strand is None:
                logging.debug("seed2 None")
                if self.seed1Name.endswith('e5'):
                    logging.debug("5' extend upstream")
                    return FastqEntry(self.gapName, \
                                ('N'*gapLen) + \
                                self.seed1ExtendSeq.seq, \
                                ('!'*gapLen) + \
                                self.seed1ExtendSeq.qual)
                else:
                    logging.debug("3' extend upstream")
                    return FastqEntry(self.gapName, \
                                self.seed1ExtendSeq.seq + \
                                ('N'*gapLen), \
                                self.seed1ExtendSeq.qual + \
                                ('!'*gapLen)) 
        
        #one fill sequence may need to be filpped
        elif not self.sameStrand:
            logging.debug("not same strand")
            if self.seed1Strand == self.seed2Strand:
                logging.debug('strand equal, flip 2')
                self.seed2ExtendSeq.reverseCompliment()
                if self.seed2Strand == '-':#minus becomes plus
                    self.seed2Strand = '+'
                else:#Plus becomes minus
                    self.seed1Strand = '-'
            if self.seed1Strand == '+':
                logging.debug('strand nequal, 1 +')
                return FastqEntry(self.gapName,
                            self.seed1ExtendSeq.seq + \
                            ('N'*gapLen) + \
                            self.seed2ExtendSeq.seq, \
                            self.seed1ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed2ExtendSeq.qual)
            elif self.seed2Strand == '+':
                logging.debug('strand nequal, 2 +')
                return FastqEntry(self.gapName,
                            self.seed2ExtendSeq.seq + \
                            ('N'*gapLen) + \
                            self.seed1ExtendSeq.seq, \
                            self.seed2ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed1ExtendSeq.qual)
            elif self.seed1Strand == '+' and self.seed2Strand == '-':
                logging.error(("Huge Problem! I haven't done this, yet "
                               "Not sameStrand (+/-)"))
                exit(10)
            elif self.seed1Strand == '-' and self.seed2Strand == '+':
                logging.error(("Huge Problem! I haven't done this, yet "
                               "Not sameStrand (-/+)"))
                exit(10)
                return FastqEntry
            else:
                logging.error(("Huge Problem! This Should Never Happen!  "
                               "Not sameStrand"))
                exit(10)
            #I'm just going to take the directStrand sequence 
            # wherever I stich things together later will need to 
            #ensure this is correct
            
    def getExtendSequence(self, contigEnd):
        """
        Get the sequence that extends the specified contig end. Returns None if 
        this metric doesn't hold extending sequence for the contig end
        """
        if contigEnd == self.seed1Name:
            return self.seed1ExtendSeq
        if contigEnd == self.seed2Name:
            return self.seed2ExtendSeq
        return None
        pass
        
    def getSeedStrand(self, name):
        if name == self.seed1Name:
            return self.seed1Strand
        elif name == self.seed2Name:
            return self.seed2Strand
        
    def __str__(self):
        return json.dumps(self.data,indent=4)
    
class Collection():

    def __init__(self):
        self.parseOpts()
        setupLogging(self.debug)
    
    def parseOpts(self):
        parser = OptionParser(USAGE)
        parser.add_option("--debug", action="store_true", 
                        help="Increases verbosity of logging")
        opts, args = parser.parse_args()
        self.debug = opts.debug
        
        if len(args) != 1:
            parser.error("Error! Incorrect number of arguments")
        
        self.protocol = JellyProtocol(args[0])
    
    
    def loadReference(self):
        """
        Get the reference information 
        """
        myReference = FastaFile(self.protocol.reference)
        self.reference = {}
        if self.protocol.scaffoldQualName is not None:
            myQual = QualFile(self.protocol.scaffoldQualName)
            for key in myReference:
                self.reference[key.split('|')[-1]] = FastqEntry(key.split('|')[-1], myReference[key], myQual[key])
        else:
            for key in myReference:
                self.reference[key.split('|')[-1]] = FastqEntry(key.split('|')[-1], myReference[key], 'l'*len(myReference[key]))
        
        #Graph gml and what not.
        bmlFile = os.path.join(self.protocol.outDir, "assembly", "masterSupport.bml")
        if not os.path.exists(bmlFile):
            logging.error("Consolidated support graph file %s not found!" \
                            % (bmlFile))
            exit(20)
        
        self.gapInfo = GapInfoFile(self.protocol.gapTable)
        self.gapGraph = GapGraph(self.gapInfo)
        self.gapGraph.loadBML(bmlFile)
        self.inputGml = self.gapGraph.graph

    def metricsCollector(self):
        folder = glob(os.path.join(self.protocol.outDir,"assembly","ref*"))
        gapStats = open(os.path.join(self.protocol.outDir, "gap_fill_status.txt"), 'w')
        self.allMetrics = {}
        numGapsAddressed = len(folder)
        noFillingMetrics = 0
        nfilled = 0
        filled = 0
        overfilled = 0
        reduced = 0
        extended = 0
        
        for f in folder:
            gapName = f.split('/')[-1]
            try:
                fh = open(os.path.join(f,"fillingMetrics.json"),'r')
            except IOError:
                noFillingMetrics += 1
                gapStats.write("%s\tnofillmetrics\n" % gapName)
                continue
            
            try:
                myMetrics = FillingMetrics(json.load(fh), gapName)
                self.allMetrics[gapName] = myMetrics
            except ValueError:
                logging.error("WARNING! "+f+" didn't produce a valid JSON output in " + \
                            "fillingMetrics.json. Go check this file and if it is " + \
                            "not a plain text JSON file, try re-running the " + \
                            "assembly Process on this folder")
                exit(1)
            fh.close()
        
            if myMetrics.span:
                if myMetrics.data["spanSeedName"] == "tooShortNs":
                    gapStats.write("%s\tn_filled\n" % gapName)
                    nfilled += 1
                else:
                    gapStats.write("%s\tfilled\n" % gapName)
                    filled += 1
                    
            elif myMetrics.predictedGapSize is not None and myMetrics.predictedGapSize < myMetrics.fillLength:
                gapStats.write("%s\toverfilled\n" % gapName)
                overfilled += 1
            elif myMetrics.singleExtend:
                gapStats.write("%s\tsingleextend\n" % gapName)
                extended += 1
            elif myMetrics.doubleExtend:
                gapStats.write("%s\tdoubleextend\n" % gapName)
                reduced += 1        
        gapStats.close()
        logging.info("Number of Gaps Addressed %d" % (numGapsAddressed))
        logging.info("No Filling Metrics %d" % (noFillingMetrics))
        logging.info("Filled %d" % (filled))
        logging.info("NFilled %d" % (nfilled))
        logging.info("Single-End Reduced %d" % (extended))
        logging.info("Double-End Reduced %d" % (reduced))
        logging.info("Overfilled %d" % (overfilled))
    
    def cleanGraph(self):
        
        #Need fully connected graphs to get the diameter
        subG = nx.connected_component_subgraphs(self.inputGml)
        logging.info("Prefilter: %d subGraphs" % (len(subG)))
        
        self.subGraphs = []
        
        for myGraph in subG:
            #Break any edge that didn't make a correct assembly
            # or doesn't span (excluding scaffold/contig evidence)
            for a,b in myGraph.edges():
                name = makeFilMetName(a,b)
                if "Contig" in myGraph.edge[a][b]['evidence'] \
                    or "Scaffolding" in myGraph.edge[a][b]['evidence']:
                    pass
                elif name not in self.allMetrics.keys():
                    logging.debug("Breaking %s due to assembly failure" %  name)
                    myGraph.remove_edge(a, b)
                elif not self.allMetrics[name].span:
                    logging.debug("Breaking %s due to non-span" %  name)
                    myGraph.remove_edge(a, b)
    
            #Resolving "forked" nodes
            # Usually some repeat. Right now, it's all about the fill quality
            for node in myGraph.nodes_iter():
                if len(myGraph.edge[node]) > 2:
                    best = None
                    bestScore = 0
                    for edge in myGraph.edge[node]:
                        name = makeFilMetName(node,edge)
                        if name in self.allMetrics.keys():
                            data = self.allMetrics[name]
                            seq = data.getSequence()
                            if seq is None:
                                #I think I fixed this
                                logging.debug("About to Fail on %s" % (node))
                            
                            if len(seq.qual) == 0:
                                logging.info("NoFilling %s " % name)
                                myScore = 15
                            else:
                                myScore = sum([ord(y)-33 for y in seq.qual])/float(len(seq.qual))
                            
                            if myScore > bestScore:
                                bestScore = myScore
                                best = name
                            
                    logging.debug("Resolved fork to be %s" % best)
                    if best is None:
                        #Again, think I fixed
                        logging.debug("I don't know how this doesn't get set")
                        logging.debug("Node %s" % (node))
                        logging.debug(json.dumps(myGraph.edge[node], indent=4))
                        
                    for edge in list(myGraph.edge[node]):
                        name = makeFilMetName(node, edge)
                        if "Contig" in myGraph.edge[node][edge]['evidence'] \
                            or "Scaffolding" in myGraph.edge[node][edge]['evidence']:
                            pass
                        elif name != best:
                            logging.debug("Removed edge %s" % name)
                            myGraph.remove_edge(node,edge)
            
            #Getting the contig paths
            for i,s in enumerate(nx.connected_component_subgraphs(myGraph)):
                #print "prefilter diameter of testSub piece %d == %d" % (i, nx.diameter(s))
                #I may get an error here if my above cleaning work isn't good enough 
                #for every case
                try:
                    #I don't understand the error I'm getting here...
                    ends = nx.periphery(s)#Singletons...
                except AttributeError:
                    logging.debug("Weird error!")
                    logging.debug("Nodes " + json.dumps(myGraph.node))
                    logging.debug("edges " + json.dumps(myGraph.edge))
                    logging.debug("types " + str(type(s)) + " " + str(i) + " " + str(type(i)))
                    logging.debug("Trying again??")
                    ends = nx.periphery(s)
                    
                if len(ends) > 2 and len(ends) == s.number_of_nodes():
                    logging.warning("Circular graph detected. Breaking weakest edge")
                    worstScore = sys.maxint
                    worstEdge = None
                    for a,b in s.edges():
                        if "Contig" in myGraph.edge[a][b]['evidence'] \
                            or "Scaffolding" in myGraph.edge[a][b]['evidence']:
                            continue
                        name = makeFilMetName(a,b)
                        data = self.allMetrics[name]
                        if data.fillLength > 0:
                            myScore = sum([ord(y)-33 \
                                for y in data.getSequence().qual])/float(data.fillLength)
                        else:
                            myScore = 0
                        
                        logging.warning(name + " " + str(myScore))
                        if myScore < worstScore:
                            worstScore = myScore
                            worstEdge = (a,b)
                    logging.info("breaking at %s" % (str(worstEdge)))
                    s.remove_edge(*worstEdge)
                #if the above didn't if didn't fix periphery, we'll get
                #a value error and a problem parsing the graph...
                try:
                    a,b = nx.periphery(s)
                except ValueError:
                    logging.error("Graph doesn't have ends. Check it's repeats in collectionErr.gml")
                    logging.error(nx.periphery(s))
                    nx.write_gml(s, "collectionError.gml")
                    exit(1)
                #What I should default to here is just breaking all non captured links

                #nx.shortest_path(s,a,b)
                self.subGraphs.append(s)
    
        logging.info("PostFilter: %d subGraphs" % (len(self.subGraphs)))
    
    def grabContig(self, nodeA, nodeB):
        """
        grabs the contig from the reference that exists
        between nodes A and B
        """
        nodeA, nodeB = orderSeeds(nodeA, nodeB)
        logging.debug("Grabbing contig between nodes %s & %s" % (nodeA, nodeB))
        scafName = nodeA[:10]
        seq = self.reference[scafName]
        #let's get the start
        
        if nodeA.count('.') == 1:
            #find gap with /0 name
            gid = int(nodeA[nodeA.rindex('.')+1:-2])
            if nodeA.endswith('e3'):
                gapName = "%s_%d_%d" % (scafName, gid, gid+1)
            else:
                gapName = "%s_%d_%d" % (scafName, gid-1, gid)
            gap = self.gapInfo[gapName]
            start = gap.end
        else:#no / means it's got to be the beginning
            start = 0

        if nodeB.count('.') == 1:
            gid = int(nodeB[nodeB.rindex('.')+1:-2])
            try:
                gapName = "%s_%d_%d" % (scafName, gid-1, gid)
                gap = self.gapInfo[gapName]
                end = gap.start
            except KeyError:
                gapName = "%s_%d_%d" % (scafName, gid, gid+1)
                gap = self.gapInfo[gapName]
                end = gap.start
            end = gap.start
        else:# no/ means it's got to be the end
            end = None
            
        logging.debug("contig %s to %s" % (str(start), str(end)))
        return seq.subSeq(start,end)
            
    def outputContigs(self):
        """
        output all the contigs, use the span and stuff get
        """
        fout = open(os.path.join(self.protocol.outDir, "pbjelly.out.fasta"), 'w')
        qout = open(os.path.join(self.protocol.outDir, "pbjelly.out.qual" ), 'w')
        liftOverTable = {}#ContigName: [(piece, strand, size), ]
        for part,graph in enumerate(self.subGraphs):
            logging.debug("Contig %d" % part)
            liftTracker = []
            start, end = nx.periphery(graph)
            path = nx.shortest_path(graph, start, end)
            curFasta = []
            curQual = []
            name = makeFilMetName(start, "")
            #First guy's extender
            if name in self.allMetrics.keys():
                logging.debug("First Extender %s" % name)
                data = self.allMetrics[name]
                seq = data.getExtendSequence(name)
                strand = '+'
                if data.seed1Strand == '-':
                    seq.reverseCompliment()
                    strand = '-'
                liftTracker.append((name, strand, len(seq.seq)))
                curFasta.append(seq.seq)
                curQual.append(seq.qual)
            
            #Did we filp the previous sequence
            pFlip = 1 
            FirstFlip = None
            
            for i, nodeA in enumerate(path[:-1]):
                nodeB = path[i+1]
                logging.debug("Moving from %s to %s (p=%d)" % (nodeA, nodeB, pFlip))
                name = makeFilMetName(nodeA, nodeB)
                lookupName = isCapturedGap(name)
                logging.debug(lookupName)
                #Existing sequence -- put in A
                if "Contig" in graph.edge[nodeA][nodeB]['evidence']:
                    logging.debug("contig")
                    #need to output the contig seq
                    seq = self.grabContig(nodeA, nodeB)
                    strand = '+'
                    if pFlip == -1:
                        strand = '-'
                        seq.reverseCompliment()
                    curFasta.append(seq.seq)
                    curQual.append(seq.qual)
                    liftTracker.append((name, strand, len(seq.seq)))
                #We have to, at the very least, keep a gap in the sequence
                elif "Scaffolding" in graph.edge[nodeA][nodeB]['evidence'] and \
                                lookupName and name not in self.allMetrics.keys():
                    logging.debug("unimproved gap")
                    #keep mat orientation the same
                    a = nodeA[nodeA.rindex('.')+1:-2]
                    b = nodeB[nodeB.rindex('.')+1:-2]
                    j = [a,b]; j.sort(); a,b = j
                    gapName = "%s_%s_%s" % (nodeA[:10], a, b)
                    gapSize = gapInfo[gapName].length
                    curFasta.append("N"*gapSize)
                    curQual.append("!"*gapSize)
                    liftTracker.append((name, '?', gapSize))
                elif "Scaffolding" in graph.edge[nodeA][nodeB]['evidence'] and \
                                              name in self.allMetrics.keys():
                    logging.debug("improved gap")
                    data = self.allMetrics[name]
                    if data.span or data.doubleExtend:
                        seq = data.getSequence()
                    else:
                        seq = data.getExtendSeq(data.seed1Name)
                    
                    logging.debug('samestrand')
                    logging.debug(data.sameStrand)
                    if not data.sameStrand:
                        logging.error(("Gap %s has opposite strand "
                                         "fillers even though they're "
                                         "within scaffold gaps") % name)
                        exit(10)#never happens
                    strand = '+'
                    logging.debug(name)
                    #logging.debug(str(seq))
                    #if data.getSeedStrand(nodeA) == '-':
                    if (pFlip == -1 and data.getSeedStrand(nodeA) == '+') or \
                       (pFlip == 1 and data.getSeedStrand(nodeA) == '-'):
                        logging.debug("flippin")
                        strand = '+'
                        seq.reverseCompliment()
                    #logging.debug(str(seq))
                    liftTracker.append((name, strand, len(seq.seq)))
                    curFasta.append(seq.seq)
                    curQual.append(seq.qual)
                else:
                    logging.debug("new sequence")
                    #Else we have new sequence. 
                    #All of the previous filters removed non-span stuff
                    if lookupName:
                        data = self.allMetrics[lookupName]
                    else:
                        data = self.allMetrics[name]
                    #except KeyError:#apparently Scaffold isn't being placed in evidence correctly
                        #n = isCapturedGap(name)
                        #if n:
                            #size = self.gapInfo[n].length
                            #curFasta.append("N"*size)
                            #curQual.append("!"*size)
                        #else:
                            #raise KeyError("gap %s wasn't filtered or didn't produce sequence correctly" % name)
                        #continue
                    seq = data.getSequence()
                    
                    a = 1 if data.getSeedStrand(nodeA) == '+' else -1
                    b = 1 if data.getSeedStrand(nodeB) == '+' else -1

                    if FirstFlip is None:
                        #FirstFlip = a == -1
                        FirstFlip = nodeA.endswith('e5') 
                        #building up the 5' end
                    
                    if pFlip == a:
                        m = 1
                    else:
                        m = -1
                    
                    strand = '+'
                    if m == -1:
                        strand = '-'
                        seq.reverseCompliment()
                    liftTracker.append((name, strand, len(seq.seq)))
                    
                    curFasta.append(seq.seq)
                    curQual.append(seq.qual)
                    
                    pFlip = b * m
                    
            name = makeFilMetName(end, "")
            #Final guy's extender
            if name in self.allMetrics.keys():
                data = self.allMetrics[name]
                seq = data.getExtendSequence(name)
                strand = '+'
                if pFlip:
                    strand = '-'
                    seq.reverseCompliment()
                liftTracker.append((name, strand, len(seq.seq)))
                curFasta.append(seq.seq)
                curQual.append(seq.qual)
            
            #We may have been assembling - strand this whole time and we need
            # revcomp it
            if FirstFlip:
                tF = []
                tQ = []
                tL = []
                logging.debug("FirstFlipping %d" % (part))
                for i in curFasta:
                    tF.append(i.translate(revComp)[::-1])
                for i in curQual:
                    tQ.append(i[::-1])
                for i in liftTracker:
                    name, strand, size = i
                    if strand == '+':
                        strand = '-'
                    elif strand == '-':
                        strand = '+'
                    tL.append((name, strand, size))
                curFasta = tF
                curQual = tQ
                liftTracker = tL
                
            fout.write(">Contig%d\n%s\n" %  (part, "".join(curFasta)))
            qout.write(">Contig%d\n%s\n" % (part, "".join(curQual)))
            liftOverTable["Contig%d" % (part)] = liftTracker
        
        fout.close()
        qout.close()
        lout = open(os.path.join(self.protocol.outDir, 'liftOverTable.json'),'w')
        json.dump(liftOverTable, lout)
        lout.close()
    
    def run(self):
        logging.info("Grabbing Filling Metrics")
        self.metricsCollector()
        logging.info("Loading Reference")
        self.loadReference()
        logging.info("Removing Poorly Supported Edges")
        self.cleanGraph()
        logging.info("Outputting new contigs")
        self.outputContigs()
        logging.info("Finished!")
        #g = nx.Graph()
        #for i in self.subGraphs:
            #for n in i.nodes():
                #g.add_node(n)
            #for e in i.edges():
                #g.add_edge(*e)
        #nx.write_gml(g,"output.gml")

def orderSeeds(a, b):
    if a.count('.') == 1 and b.count('.') == 0:
        if a.endswith('e3'):
            return b, a
        else:
            return a, b
    if b.count('.') == 1 and a.count('.') == 0:
        if b.endswith('e5'):
            return b, a
        else:
            return a, b
    if a.count('.') == 0 and b.count('.') == 0:
        return a, b
    aint = int(a.split('.')[1][:-2])
    bint = int(b.split('.')[1][:-2])
    if aint < bint:
        return a, b
    return b,a

def isCapturedGap(name):
    """
    checks if assembly folder's name holds captured gap information.
    returns False if it doesn't. 
    returns the gap's name (in a gapinfofile) if it does
    """
    nodes = name.split('_')
    if len(nodes) != 2:
        return False

    a, b = nodes
    try:
        aref, aend = a.split('.')
        bref, bend = b.split('.')
    except ValueError:
        return False
    if aref != bref:
        return False
    if not aend.endswith("e3") or not bend.endswith("e5"):
        return False
    c = int(aend[:-2])
    d = int(bend[:-2])

    if not c + 1 == d:
        return False
    return "%s_%d_%d" % (aref, c, d)

if __name__ == '__main__':
    c = Collection()
    c.run()
    """
    Load the GapGraph.gml and the reference sequence
    
    for every successful ref_ref, add a link.
    
    resolve redundant links
        
    try to get the maximum diameter and whatever 
    
    traverse from end to end - marking where each piece is/was
        
    out
    """
