#!/usr/bin/env python
"""
document the protocol
"""

import os, sys, subprocess, logging, time
from optparse import OptionParser
from string import Template
from xml.etree import ElementTree
from glob import glob
import Stages
from CommandRunner import * 

USAGE = """USAGE: Jelly.py <stage> <protocol.xml> [-x \"--options for stage\"]
    
    Jelly is the driver for each stage of the 
    reference genome upgrade. 
    
    <stage> is one of
        %s
    <protocol.xml> contains the information about the
    data and parameters Jelly will run. See README.txt
    or the documentation for the Protocol's format.""" % ("\n\t".join(STAGES))

class JellyRunner():
    """
    Take a JellyProtocol and loads in the variables in a way that 
    JellyRunner can use it easily!
    """
    def __init__(self):
        """
        Given a protocol fn, load it up so we are ready to run. 
        """
        self.parseArgs()
        self.__initLog()
        self.parseProtocol()
        
    def __initLog(self):
        """Logging"""
        logLevel = logging.DEBUG if self.options.debug else logging.INFO
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
        logging.info("Running %s" % " ".join(sys.argv) )
    
    def parseArgs(self):
        """
        Uses OptionParser to parse out input
        Jelly.py <stage> <protocol>
        """
        parser = OptionParser(USAGE)
        parser.remove_option("-h")
        parser.add_option("-h", "--help", action="store_true", default=False)
        
        parser.add_option("--debug",action="store_true",default=False)
        parser.add_option("-x", dest="extras", type="string", default="", 
                help="-x \"<options>\" are options to pass into the stage you're running")
        
        self.options, args = parser.parse_args()

        if self.options.help == True:
            if len(args) == 1:
                if args[0] in STAGES:
                    print exe(Stages.PRINT_HELPS[args[0]])[1]
                    sys.exit(0)
                #Else, this will drop down to the next parser.error
            else:
                print parser.format_help()
                sys.exit(0)
        if len(args) != 2 or args[0] not in STAGES:
            parser.error("Invalid Arguments. Expected one of\n'%s'" % "', '".join(STAGES))
            sys.exit(1)
        self.executeStage = args[0]
        self.protocolName = args[1]
        
    def parseProtocol(self):
        """
        Parses a Jelly Protocol
        """
        try:
            file = ElementTree.parse(self.protocolName)
        except Exception:
            logging.error(("If you're actually sure the input Protocol is "\
                           "where you said it is, then the Protocol doesn't "\
                           "have valid XML Format. Use an XML validator."))
            sys.exit(1)
        
        root = file.getroot()
        
        refNode = root.find("reference")
        
        if refNode == None:
            logging.error("Protocol doesn't have <reference> element.")
            sys.exit(1)
        else:
            self.scaffoldName = refNode.text
            if not os.path.exists(self.scaffoldName):
                logging.error("Reference %s Does Not Exist!" % (self.scaffoldName))
                sys.exit(1)
            self.referenceNameBase = self.scaffoldName[:self.scaffoldName.rindex(".fasta")]
            self.scaffoldQualName = self.referenceNameBase+".qual"
            self.gapTable = self.referenceNameBase+".gapInfo.bed"
            
            self.contigName = self.referenceNameBase+".contigs.fasta"
            self.contigQualName = self.referenceNameBase+".contigs.qual"
            #Use the Scaffolding
            """Disabling for now
            if refNode.attrib["useContigs"] == "False":
                self.useContigs = False
                self.reference = self.scaffoldName
            else:
                self.useContigs = True
                self.reference = self.contigName
            """
            self.reference = self.scaffoldName
            self.useContigs = False
             
            self.referenceIndex = self.reference+".sa"
            
        inputNode = root.find("input")
        if inputNode.attrib.has_key("baseDir"):
            self.baseDir = inputNode.attrib["baseDir"]
        else:
            self.baseDir = ""
        
        self.inputs = []
        for input in inputNode:
            self.inputs.append(os.path.join(self.baseDir, input.text))
        if len(self.inputs) == 0:
            logging.error("Protocol doesn't specifiy any input inputs!")
            sys.exit(1)

        
        outputNode = root.find("outputDir")
        if outputNode == None:
            logging.warning("Output directory not specified. Using pwd. Hope you're cool with that...")
            self.outDir = os.getcwd()
        else:
            self.outDir = outputNode.text
        
        #Command Runner is going to sort all this out.
        #And create an object to do the calls.
        self.runCmd = CommandRunner(root.find("cluster"))
        
        blasrNode = root.find("blasr")
        if blasrNode == None:
            logging.warning("No blasr parameters!?")
            self.blasrParams = ""
        else:
            self.blasrParams = blasrNode.text
    
    
    def run(self):
        logging.info("Executing Stage: %s" % self.executeStage)
        
        if self.options.debug:
            Stages.DEBUG = "--debug"
        
        #Setup before a stage and run
        if self.executeStage == "setup":
            wDir = os.path.dirname(self.scaffoldName)
            myCommands = [Stages.setup(self.scaffoldName, \
                        self.scaffoldQualName, self.gapTable, \
                        self.options.extras)]
        
        elif self.executeStage == "mapping":
            wDir = os.path.join(self.outDir, "mapping")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
            
            myCommands = Stages.mapping(self.inputs, wDir, \
                        self.reference, self.referenceIndex, \
                        self.blasrParams, self.options.extras)
        
        elif self.executeStage == "support":
            wDir = os.path.join(self.outDir, "support")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
       
            myCommands = Stages.support(self.outDir, self.gapTable, \
                    wDir, self.options.extras) 
        
        elif self.executeStage == "extraction":
            wDir = os.path.join(self.outDir, "assembly")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
            
            myCommands = [Stages.extraction(self.scaffoldName, \
                        self.scaffoldQualName , self.gapTable, \
                        self.outDir, wDir, self.inputs, \
                        self.options.extras)]
        
        elif self.executeStage == "assembly":
            wDir = os.path.join(self.outDir, "assembly")
            
            myCommands = Stages.assembly(wDir, self.gapTable, \
                                         self.options.extras)
        
        elif self.executeStage == "output":
            wDir = os.path.join(self.outDir, "assembly")
            
            myCommands = [Stages.collection(wDir, \
                    self.contigName, \
                    self.contigQualName,\
                    self.gapTable, \
                    self.options.extras) ]
        
        logging.debug("CommandRunner Returned: " + 
            str(self.runCmd(myCommands, wDir, self.executeStage )) )
            
                  
        logging.info("Finished %s Stage: %s" % (self.runCmd.runType, self.executeStage))
        
if __name__ == '__main__':
    prog = JellyRunner()
    prog.run()
