from string import Template
import subprocess, signal, logging, os, stat, sys

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm

STAGES = ["setup", "mapping", "support", "extraction", "assembly", "output"]

def exe(cmd, timeout=1440):
    """
    Executes a command through the shell.
    timeout in minutes!
    """
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, \
                            stderr=subprocess.STDOUT, close_fds=True)
    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(int(timeout*60))  # 5 minutes
    try:
        stdoutVal, stderrVal =  proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d" % (timeout)))

        return 214,None,None
    
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal

class CommandSetup():
    #CommandSetup Replacement for namedtuple
    def __init__(self, cmd, jobname, stdout, stderr):
        self.cmd = cmd
        self.jobname = jobname
        self.stdout = stdout
        self.stderr = stderr
    
class CommandRunner():
    """
    Uses a command template to run stuff. This is helpful for cluster commands
    """
    def __init__(self, xmlNode=None):
        """
        if xmlNode:
            build how to run the command
        else:
            running the command will be easy
        """
        self.paramDict = {}
        self.buildTemplate(xmlNode)
    
    def buildTemplate(self, xmlNode):
        """
        Takes an xmlNode and buils the template..
        """
        if xmlNode == None:
            self.cmdTemplate = Template("${CMD} > ${STDOUT} 2> ${STDERR}")
            self.nJobs = 0
            for s in STAGES:
                self.paramDict[s] = {}
            self.runType = "Running"
        else:
            self.runType = "Submitting"
            command = xmlNode.find("command")
            if command == None:
                logging.error(("You're trying to use a cluster " \
                                  "template, but you didn't specify the " \
                                  "template. Please read the documentation." \
                                  "Exiting.\n"""))
                sys.exit(1)
            nJobs = xmlNode.find("nJobs")
            
            if nJobs == None or nJobs.text == '0':
                logging.warning(("Not Specifying number of jobs may " \
                                  "overload clusters."))
                self.nJobs = 0
            else:
                self.nJobs = int(nJobs.text)
            
            self.cmdTemplate = Template(command.text)
            """
            All of the below code was created to allow specifying
            resources on a per-stage basis.
            # We'll implement this in future versions of Jelly
            """
            myStages = xmlNode.findall("stage")
            #All -- this means everyone gets the same parameters
            if len(myStages) == 0:
                for s in STAGES:
                    self.paramDict[s] = {}
            #Whoops, didn't say all and don't have the stages specified
            elif len(myStages) != len(STAGES):
                logging.error(("You don't have parameters for "\
                                  "'all' or for each stage [ %s ] " \
                                  "specified in your cluster setup! " \
                                  "Please read the Documentation. " \
                                  "Exiting\n" % ", ".join(STAGES)))
                sys.exit(1) 
            #Okay, let's build our template 
            else:   
                for s in myStages:
                    if s.attrib["name"] not in STAGES:
                        logging.error(("%s is an unknown stage."\
                                          "Only use %s. (case sensitive)"\
                                          "Exiting\n" % \
                                          (s["name"], ", ".join(STAGES)) ))
                        sys.exit(1)
                    self.paramDict[s.attrib["name"]] = {}
                    for item in s:
                        self.paramDict[s.attrib["name"]][item.tag] = item.text  
        
        #if we've gotten this far, we can do a test substitue:
        for s in self.paramDict.keys():
            temp = dict(self.paramDict[s])
            
            temp.update({"CMD":"test", \
                         "STDOUT":"testo", \
                         "STDERR":"teste", \
                         "JOBNAME":"testn"})
            try:
                w = self.cmdTemplate.substitute(temp)
            except KeyError:
                logging.error(("Your Cluster Command " \
                                  "Doesn't Have enough arguments " \
                                  "to fill in the template\n" \
                                  "Please Fix this.\n"))
                sys.exit(1)
            
    def __call__(self, cmdSetups, wDir, stage):
        """
        Executes a list of command setup object
        if stdout/stderr are not specified, we write to /dev/null
        """
        outputRet =  []
        if self.nJobs == 0:
            for cs in cmdSetups:
                cmd = self.buildCommand(cs, stage)
                outputRet.append(exe(cmd))
        else:
            chunk = 0
            for commands in partition(cmdSetups, self.nJobs):
                outScript = open(os.path.join(wDir,stage+"_chunk%d.sh" % chunk),'w')
                outScript.write("#!/bin/bash\n\n")
                for j in commands:
                    outScript.write(j.cmd+"\n")
                outScript.close()
                #Add executeable 
                existing_permissions = stat.S_IMODE(os.stat(outScript.name).st_mode)
                if not os.access(outScript.name, os.X_OK):
                    new_permissions = existing_permissions | stat.S_IXUSR
                    os.chmod(outScript.name, new_permissions)
                
                cs = CommandSetup(outScript.name, \
                                stage+"_chunk%d" % chunk, \
                                os.path.join(wDir,stage+"_chunk%d.outlog" % chunk), \
                                os.path.join(wDir,stage+"_chunk%d.errlog" % chunk))
                cmd = self.buildCommand(cs, stage)
                outputRet.append(exe(cmd))
                chunk += 1
        
        return outputRet
        
    def buildCommand(self, cmdSetup, stage):
        #This is for cluster stages command setup
        runDict = dict(self.paramDict[stage]) 
        
        runDict.update({"CMD":cmdSetup.cmd, \
                        "STDOUT":cmdSetup.stdout, \
                        "STDERR":cmdSetup.stderr, \
                        "JOBNAME":cmdSetup.jobname})
        
        return self.cmdTemplate.substitute(runDict)

def partition(n,m):
    """
    Helper function. splits list n into m partitions
    """
    p = map(lambda x: list(), range(m))
    index = 0
    for item in n:
        p[index].append(item)
        if index < m-1:
            index += 1
        else:
            index = 0
    return filter(lambda x: len(x)>0, p)
