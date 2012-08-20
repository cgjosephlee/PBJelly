
import glob, os, sys, logging
from string import Template
from CommandRunner import CommandSetup
from FileHandlers import GapInfoFile

DEBUG = ""#change to --debug for use
"""
This code is more about setting up commands to run other things, not actual computations.
"""

#Get Jelly SourceDir from the environment:
if os.environ.has_key("JELLYPATH"):
    SRCDIR = os.environ["JELLYPATH"]
else:
    sys.stderr.write("Error! JELLYPATH not found in environment variables. Did you source 'SetupPaths.sh'?\n")
    sys.exit(1)

PRINT_HELPS = {"setup": os.path.join(SRCDIR, "Setup.py --help"), \
               "mapping": "blasr -h", \
               "support": os.path.join(SRCDIR, "Support.py --help"), \
               "extraction": os.path.join(SRCDIR, "Extraction.py --help"), \
               "assembly": os.path.join(SRCDIR, "WrapAssembly.py --help"), \
               "output": os.path.join(SRCDIR, "Collection.py --help")}

def setup( scaffoldName, scaffoldQualName, gapInfoName , extras):
    """
    Generate all the information we need from the input scaffolding
    """
    command = Template(os.path.join(SRCDIR, \
                "Setup.py ${scaf} ${scafQual} -g ${gap} -i ${debug} ${extras}")).substitute( \
                    {"scaf":scaffoldName, \
                    "scafQual":scaffoldQualName, \
                    "gap":gapInfoName, \
                    "debug":DEBUG, \
                    "extras":extras})
    baseName = os.path.dirname(scaffoldName)
    ret = CommandSetup(command, "setup", os.path.join(baseName,"setup.out"), \
                        os.path.join(baseName,"setup.err"))
    return ret

def mapping(jobDirs, outDir, reference, referenceSa, parameters, extras):
    """
    Input:
        - a list of fasta inputs
        - an output directory
        - a reference - (should be contigs.. see scaffoldIntakeSetup
        - a pacbio indexed reference
    Task:
        - map each fasta to reference
    Output:
        - m4 alignments of the input sequence
    """
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    level = "DEBUG" if DEBUG != "" else "INFO"
    logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
    logging.info("Running blasr")
    
    mappingTemplate = Template("blasr ${fasta} ${ref} -m 4 -sa ${sa} -out ${outFile} ${parameters} ${extras}")
    ret = []
    for fasta in jobDirs:
        name = fasta[fasta.rindex('/')+1:]
        
        if not os.path.exists(fasta):
            logging.error("%s doesn't exist." % fasta)
            exit(1)
        
        outFile = os.path.join(outDir,name+".m4")
        if os.path.isfile(outFile):
            logging.warning("Output File %s already exists and will be overwritten." % (outFile))
        
        #Build Blasr Command 
        cmd = mappingTemplate.substitute( {"fasta":fasta, 
                           "ref":reference, 
                           "sa":referenceSa, 
                           "outFile":outFile, 
                           "parameters":parameters,
                           "extras":extras} )
        
        #Build CommandSetup to send to CommandRunner 
        jobname = name+".mapping"
        stdout = os.path.join(outDir, name+".out")
        stderr = os.path.join(outDir, name+".err")
        ret.append( CommandSetup(cmd, jobname, stdout, stderr) )
    
    return ret


def support(inputDir, gapTable, outputDir, extras):
    ret = []
    command = Template(os.path.join(SRCDIR, \
                    "SupportGaps.py ${inputm4} ${gapTable} ${outFile} ${debug} ${extras}"))
    mappingFiles = glob.glob(os.path.join(inputDir, "mapping/*.m4"))
    
    if len(mappingFiles) == 0:
        logging.warning("No mapping files found!")
        return ret
        
    
    for inputm4 in  mappingFiles:
        baseName = inputm4[inputm4.rindex('/')+1:inputm4.rindex(".m4")]
        outFile = os.path.join(outputDir, baseName+".gapCans")
        if os.path.isfile(outFile):
            logging.warning("Overwriting %s" % outFile)
        myCommand = command.substitute( {"inputm4": inputm4,\
                         "gapTable": gapTable,\
                         "outFile": outFile,\
                         "debug": DEBUG,\
                         "extras":extras} )
        
        ret.append( CommandSetup(myCommand,\
                     baseName+".support",\
                     os.path.join(outputDir,baseName+".out"),\
                                 os.path.join(outputDir,baseName+".err")) )

    return ret


def extraction(ref, qual, gapInfo, inputDir, outputDir, inputs, extras):
    #Refactor -- All the inputs is bad form.
    inputDir = os.path.join(inputDir, "support")
    if not os.path.exists(inputDir):
        logging.warning("%s was not found! Have you run 'support' yet?" % inputDir)
        sys.exit(1)

    command = Template(os.path.join(SRCDIR, \
                "Extraction.py ${ref} ${qual} ${gapInfo} ${inputDir} ${outputDir} ${inputs} ${debug} ${extras}"))
    myCommand = command.substitute({"ref": ref,
                    "qual": qual,
                    "gapInfo": gapInfo,
                    "inputDir": inputDir,\
                    "outputDir": outputDir,\
                    "inputs":" ".join(map(lambda x: x, inputs)),\
                    "debug":DEBUG,\
                    "extras":extras})
    
    return CommandSetup(myCommand, "extraction", \
                os.path.join(outputDir,"extraction.out"), \
                os.path.join(outputDir,"extraction.err"))

def assembly(inputDir, gapInfoFn, extras):
    gapInfo = GapInfoFile(gapInfoFn)
    command = Template(os.path.join(SRCDIR, \
                "WrapAssembly.py ${inputDir} ${gapInfo} ${debug} ${extras}"))
    ret = []
    allInputs = glob.glob(os.path.join(inputDir,"ref*"))
    if len(allInputs) == 0:
        logging.warning("No gaps to be assembled were found in %s! Have you run 'extraction' yet?" % inputDir)
        sys.exit(1)
    for input in allInputs:
        #Skip Contig Ends Assembly... Though we could extend?
        if input.split('/')[-1] not in gapInfo.keys():
            logging.info("Skipping "+input+" -- Potential Contig Extension")
            continue
        myCommand = command.substitute({"inputDir":input,\
                        "gapInfo":gapInfoFn,\
                        "debug":DEBUG,\
                        "extras":extras})
        
        ret.append(CommandSetup(myCommand, 
               os.path.join(input.split('/')[-1],"assembly"), \
               os.path.join(input,"assembly.out"), \
               os.path.join(input,"assembly.err")) )
    return ret

def collection(inputDir, fasta, qual, gapInfo, extras):
    command = Template(os.path.join(SRCDIR, \
                "Collection.py ${inputDir} ${fasta} ${qual} ${gapInfo} ${debug} ${extras}"))
    
    myCommand = command.substitute({"inputDir":inputDir,\
                    "fasta": fasta,\
                    "qual": qual,\
                    "gapInfo":gapInfo,\
                    "debug":DEBUG,\
                    "extras":extras})

    return CommandSetup(myCommand, \
            os.path.join(inputDir,"collectingOutput"), \
            os.path.join(inputDir,"output.out"), \
            os.path.join(inputDir,"output.err"))

