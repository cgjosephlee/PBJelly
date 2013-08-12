
import glob, os, sys, logging
from string import Template


from pbsuite.utils.CommandRunner import Command
from pbsuite.utils.FileHandlers import GapInfoFile

DEBUG = ""#change to --debug for use
"""
This code is more about setting up commands to run other things, not actual computations.
"""


PRINT_HELPS = {"setup": os.path.join("Setup.py --help"), \
               "mapping": "blasr -h", \
               "support": os.path.join("Support.py --help"), \
               "extraction": os.path.join("Extraction.py --help"), \
               "assembly": os.path.join("WrapAssembly.py --help"), \
               "output": os.path.join("Collection.py --help")}

def setup( scaffoldName, scaffoldQualName, gapInfoName , extras):
    """
    Generate all the information we need from the input scaffolding
    """
    #Error
    if scaffoldQualName is None:
        scaffoldQualName = ""
    
    #"Setup.py ${scaf} ${scafQual} -g ${gap} -i ${debug} ${extras}")).substitute( \
    command = Template("Setup.py ${scaf} -g ${gap} -i ${debug} ${extras}")\
                    .substitute( \
                    {"scaf":scaffoldName, \
                    "scafQual":scaffoldQualName, \
                    "gap":gapInfoName, \
                    "debug":DEBUG, \
                    "extras":extras})
    baseName = os.path.dirname(scaffoldName)
    ret = Command(command, "setup", os.path.join(baseName,"setup.out"), \
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
    #Sometimes, you can't use the sa index. It's a problem, I know
    mappingTemplate = Template("blasr ${fasta} ${ref} -m 4 -out ${outFile} ${parameters} ${extras}")
    #mappingTemplate = Template("blasr ${fasta} ${ref} -m 4 -sa ${sa} -out ${outFile} ${parameters} ${extras}")
    #Sam for misassembly ident
    #mappingTemplate = Template("blasr ${fasta} ${ref} -sam -sa ${sa} -out ${outFile} ${parameters} ${extras}")
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
        
        #Build Command to send to CommandRunner 
        jobname = name+".mapping"
        stdout = os.path.join(outDir, name+".out")
        stderr = os.path.join(outDir, name+".err")
        ret.append( Command(cmd, jobname, stdout, stderr) )
    
    return ret


def support(inputDir, gapTable, outputDir, extras):
    ret = []
    command = Template("Support.py ${inputm4} ${gapTable} ${outFile} ${debug} ${extras}")
    mappingFiles = glob.glob(os.path.join(inputDir, "mapping/*.m4"))
    
    if len(mappingFiles) == 0:
        logging.warning("No mapping files found!")
        return ret
    
    for inputm4 in  mappingFiles:
        baseName = inputm4[inputm4.rindex('/')+1:inputm4.rindex(".m4")]
        outFile = os.path.join(outputDir, baseName+".gml")
        if os.path.isfile(outFile):
            logging.warning("Overwriting %s" % outFile)
        myCommand = command.substitute( {"inputm4": inputm4,\
                         "gapTable": gapTable,\
                         "outFile": outFile,\
                         "debug": DEBUG,\
                         "extras":extras} )
        
        ret.append( Command(myCommand,\
                     baseName+".support",\
                     os.path.join(outputDir,baseName+".out"),\
                                 os.path.join(outputDir,baseName+".err")) )

    return ret

def extraction(outputDir, protocol, extras):
    command = Template("Extraction.py ${protocol} ${debug} ${extras}")
    myCommand = command.substitute({"protocol": protocol, \
                    "debug":DEBUG, \
                    "extras":extras})
    
    return Command(myCommand, "extraction", \
                os.path.join(outputDir,"extraction.out"), \
                os.path.join(outputDir,"extraction.err"))

def assembly(inputDir, gapInfoFn, extras):
        
    gapInfo = GapInfoFile(gapInfoFn)
    command = Template("Assembly.py ${inputDir} ${size} ${debug} ${extras}")
    ret = []
    allInputs = glob.glob(os.path.join(inputDir,"ref*"))
    if len(allInputs) == 0:
        logging.warning("No gaps to be assembled were found in %s! Have you run 'extraction' yet?" % inputDir)
        sys.exit(1)
        
    for inputDir in allInputs:
        #get The predicted size if exists
        mySize = ""
        gapName = inputDir.split('/')[-1]
        if gapName.count('.') > 0:
            g = gapName.split('_')
            if len(g) == 1:
                ref,cnam = g[0].split('.')
                cn = int(cnam[:-2])
                if cnam.endswith('e5'):
                    ca = cn-1
                    cb = cn
                elif cnam.endswith('e3'):
                    ca = cn
                    cb = cn+1
                else:
                    logging.error("Single End Extension is misFormatted for %s" % inputDir)
                    exit(1)
                size = gapInfo["%s_%d_%d" % (ref, ca, cb)].length
                mySize = "-p %d" % (size)
            elif len(g) == 2:
                ref,ca = g[0].split('.')
                ref,cb = g[1].split('.')
                #Hackey Shack. - Affect of sorting node names
                # to prevent redundancies during graph building
                ca = int(ca[:-2])
                cb = int(cb[:-2])
                j = [ca,cb]
                j.sort()
                ca, cb = j
                size = gapInfo["%s_%d_%d" % (ref, ca, cb)].length
                mySize = "-p %d" % (size)
            else:
                logging.error("Couldn't recreate gapName from refDir for %s" % inputDir)
                exit(1)
        myCommand = command.substitute({"inputDir":inputDir,\
                                        "size":mySize,\
                                        "debug":DEBUG,\
                                        "extras":extras})
        
        ret.append(Command(myCommand, 
               os.path.join(inputDir.split('/')[-1],"assembly"), \
               os.path.join(inputDir,"assembly.out"), \
               os.path.join(inputDir,"assembly.err")) )
    
    return ret


def collection(inputDir, protocol):
    command = Template("Collection.py ${protocol}")
    
    myCommand = command.substitute({"protocol": protocol.protocolName})
    
    return Command(myCommand, \
            os.path.join(inputDir,"collectingOutput"), \
            os.path.join(inputDir,"output.out"), \
            os.path.join(inputDir,"output.err"))

