import re, sys, os, bisect
from collections import defaultdict
from StringIO import StringIO
from string import maketrans

"""
An object for loading Fasta entries from a file into a
dictionary {SequenceName:sequence}
"""
def wrap(string, width=100):
    return os.linesep.join( \
        [ string[i:i+width] for i in xrange(0,len(string),width) ] )

#Refactor this -- It's too slow 
def qwrap(lst, width = 40):
    newStr = ""
    
    for i in xrange(0,len(lst),width):
        newStr += " ".join(map(str,lst[i:i+width]))+" \n"
    
    return newStr.strip()

class FastaFile(dict):
    
    def __init__(self, fileName):
        super(dict)
        self.fileHandler = open(fileName,'r')
        self.__parse()
        self.fileHandler.close()
    
    def __parse(self):
        for line in self.fileHandler.readlines():
            if line.startswith('>'):
                curName = line.strip()[1:]
                self[curName] = StringIO()
                continue
            self[curName].write(line.strip())
    
        for key in self:
            self[key] = self[key].getvalue()
    
    def toString(self, key):
        return ">" + key + "\n" + wrap(self[key]).strip()

class QualFile(dict):
    
    def __init__(self, fileName, convert=True):
        super(dict)
        self.fileHandler = open(fileName,'r')
        self.convert = convert
        self.__parse()
        self.fileHandler.close()
    
    def __parse(self):
        splRE = re.compile("\s+")
        for line in self.fileHandler.readlines():
            if line.startswith('>'):
                curName = line.strip()[1:]
                if self.convert:
                    self[curName] = []
                else:
                    self[curName] = StringIO()
                continue
            
            if self.convert:
                self[curName].extend(map(int, splRE.split(line.strip())))
            else:
                self[curName].write(line.strip()+" ")
        
        if not self.convert:
            for key in self:
                self[key] = self[key].getvalue().strip()
    
    def __parse2(self):
        for line in self.fileHandler.readlines():
            if line.startswith('>'):
                curName = line.strip()[1:]
                continue
            self[curName].write(line.strip() + " ")
        
        for key in self:
            self[key] = map(int, re.split("\s+", self[key].getvalue().strip()))
    

    def toString(self, key):
        return ">" + key + "\n" + qwrap(self[key]).strip()
    
    def valToString(self,key):
        return qwrap(self[key]).strip()

#make an index class that fasta and qual can inherit
class IndexFasta():
    def __init__(self, file):
        self.file = open(file)
        self.buildindex()
    
    def __getitem__(self, key):
        try:
            pos = self.index[key]
        except KeyError:
            raise IndexError("no such item")
        f = self.file
        f.seek(pos)
        header = f.readline()
        data = []
        while True:
            line = f.readline().strip()
            if not line or line[0] == ">":
                break
            data.append(line.strip())
        return "".join(data)
    
    def buildindex(self):
        index = {}
        # scan the file
        f = self.file
        f.seek(0)
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                break
            if line[0] == ">":
                # save offset to header line
                header = line[1:].strip()
                index[header] = pos
        self.index = index

class IndexQual():
    def __init__(self, file, convert=True):
        self.file = open(file)
        self.convert = convert
        self.buildindex()
    
    def __getitem__(self, key):
        try:
            pos = self.index[key]
        except KeyError:
            raise IndexError("no such item")
        f = self.file
        f.seek(pos)
        header = f.readline()
        data = [] if self.convert else ""
        while True:
            line = f.readline()
            if not line or line[0] == ">":
                break
            if self.convert:
                data.extend(map(int, re.split("\s+", line.strip())))
            else:
                data += line
        return data
    
    def buildindex(self):
        index = {}
        # scan the file
        f = self.file
        f.seek(0)
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                break
            if line[0] == ">":
                # save offset to header line
                header = line[1:].strip()
                index[header] = pos
        self.index = index


"""
An object for handling gapInfo files, which are in bed format:

targetName \t start \t end \t gapName
"""

class GapInfoFile(dict):
    
    def __init__(self, fileName):
        super(dict)
        self.fileHandler = open(fileName,'r')
        #inorder lift of gaps
        for line in self.fileHandler.readlines():
            curGap = Gap(*line.strip().split('\t'))
            self[curGap.name] = curGap
        
        self.fileHandler.close()
    def getSortedGaps(self):
        """
        Returns dictionary of gaps partitioned by scaffold and sorted by
        location
        """
        ret = defaultdict(list)
        for key in self:
            #print str(self[key])
            #print ret[self[key].scaffold]
            bisect.insort(ret[self[key].scaffoldId], self[key])
            #print ret[self[key].scaffold]
            #raw_input()
            
        return dict(ret)
        
    def getScaffFromIndex(self, key):
        """
        Looks through file for scaffold name
        based on scaffold index
        """
        for x in self.keys():
            if x.startswith(key):
                return self[x].scaffold
        
        raise KeyError(key)

class Gap():
    
    def __init__(self, scaff, start, end, name):
        #scaff, start, end, name = line.strip().split(delim)
        ref, lcontig, rcontig = name.split('_')
        self.name = name
        self.scaffold = scaff
        self.scaffoldId = ref
        self.leftContig = self.scaffold+"/"+lcontig
        self.rightContig = self.scaffold+"/"+rcontig
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start
    def __str__(self):
        return "\t".join([self.scaffold, str(self.start), str(self.end), self.name])

    def __lt__(self, other):
        if type(other) == type(self):
            return self.start - other.start < 0
        elif type(other) == int:
            return self.start - other < 0
        else:
            raise AttributeError
    def __gt__(self, other):
        if type(other) == type(self):
            return self.start - other.start > 0
        elif type(other) == int:
            return self.start - other > 0
        else:
            raise AttributeError
    
"""
An object for handling m4 format mapped reads information efficently.

m4 format:

qname tname score pctsimilarity qstrand qstart qend qseqlength 
tstrand tstart tend tseqlength [mapqv ncells npaths]
"""

subreadGrab = re.compile(".*/(\d+)_(\d+)$")

class M4File(list):

    def __init__(self, file):
        super(list)
        if type(file) != type(StringIO()):
            file = open(file,'r')
        self.fileHandler = file
        self.__parse()
        self.fileHandler.close()
        
    def __parse(self):
        for line in self.fileHandler.readlines():
            try:
                self.append(M4Line(line.strip()))
            except TypeError, err:
                sys.stderr.write("BadM4Line! \n%s\n%s\n" % (line, str(err)) )
                sys.exit(1)
            

class M4Line():
    
    def __init__(self, line):
        data = re.split("\s+",line)
        
        self.qname           = data[0]
        self.tname           = data[1]
        self.score           = int(data[2])
        self.pctsimilarity   = float(data[3]) 
        self.qstrand         = data[4]
        self.qstart          = int(data[5])
        self.qend            = int(data[6])
        self.qseqlength      = int(data[7])
        self.tstrand         = data[8]
        self.tstart          = int(data[9])
        self.tend            = int(data[10])
        self.tseqlength      = int(data[11]) 
        """Old Blasr -- In SMRTAnalysis 1.3.0
        if self.qstrand == '1':
            self.qstart, self.qend = self.qseqlength - self.qend, \
                                     self.qseqlength - self.qstart 
        """
        if self.tstrand == '1':
            self.tstart, self.tend = self.tseqlength - self.tend, \
                                     self.tseqlength - self.tstart 
        #"""
        #Collect subread Information
        try:
            subStart, subEnd     = subreadGrab.match(self.qname).groups()
        except AttributeError:
            subStart, subEnd = self.qstart, self.qend
        
        self.qsubstart       = int(subStart)
        self.qsubend         = int(subEnd)
        self.qsublength      = self.qsubend - self.qsubstart
        self.queryPctAligned = (self.qend - self.qstart) \
                               / float(self.qsubend - self.qsubstart)
            
        
        self.flag = 0
        self.trim = False
    
    def __str__(self, convert=False):
        #if convert:#Put back into original space
        if self.tstrand == '1':
            self.tstart, self.tend = self.tseqlength - self.tend, \
                                     self.tseqlength - self.tstart 

        return " ".join(map(str, [self.qname, self.tname, self.score, \
                                  self.pctsimilarity, self.qstrand,   \
                                  self.qstart, self.qend, self.qseqlength, \
                                  self.tstrand, self.tstart, self.tend, \
                                  self.tseqlength]))
    
    def toBed(self):
        """
        Returns the M4Line to a bed file
        """
        if self.tstrand == '1':
            strand = "-"
            chromStart = str(self.tstart - (self.qseqlength - self.qend))
            chromEnd = str(self.tend + self.qstart)
        else:
            strand = "+"
            chromStart = str(self.tstart - self.qstart)
            chromEnd = str(self.tend + (self.qseqlength - self.qend))
        #chromStart = str(self.tstart - self.qstart)
        #chromEnd = str(self.tend + (self.qseqlength - self.qend))
        chrom = self.tname
        name = self.qname
        score = str(self.score)
        thickStart = str(self.tstart)
        thickEnd = str(self.tend)
        itemRgb = str(int(self.pctsimilarity))
        return "\t".join([chrom, chromStart, chromEnd, name, score, strand, \
                          thickStart, thickEnd, itemRgb])
"""
Parses M5 alignment object

"""

revComp = maketrans("ATCGNatcgn","TAGCNtagcn")

class M5File(list):
    
    def __init__(self, file):
        super(list)
        if type(file) != type(StringIO()):
            file = open(file,'r')
        self.fileHandler = file
        self.__parse()
        self.fileHandler.close()
    
    def __parse(self):
        """
        Returns a list of all the M5 lines
        """
        for line in self.fileHandler.readlines():
            try:
                self.append(M5Line(line.strip()))
            except TypeError, err:
                sys.stderr.write("BadM5Line! \n%s\n%s\n" % (line, str(err)) )
                sys.exit(1)
            
class M5Line():
    
    def __init__(self, line):
        data = re.split("\s+",line)
        
        self.qname  = data[0]
        self.qseqlength = int(data[1])
        self.qstart = int(data[2])
        self.qend   = int(data[3])
        self.qstrand    = data[4]
        self.tname  = data[5]
        self.tseqlength = int(data[6])
        self.tstart = int(data[7])
        self.tend   = int(data[8])
        self.tstrand    = data[9]
        self.score  = int(data[10])
        self.nMatch = int(data[11])
        self.nMismatch  = int(data[12])
        self.nInsert    = int(data[13])
        self.nDelete    = int(data[14])
        self.weirdUnknown    = data[15]
        self.querySeq   = data[16]
        self.compSeq    = data[17]
        self.targetSeq  = data[18]
        
        self.queryPctAligned = (self.qend - self.qstart)/float(self.qseqlength)
        self.pctsimilarity = self.nMatch / float(self.qend - self.qstart)
        #"""newBlasr
        if self.tstrand == '-':
            self.negStrand = True
            #translating to + strand.
            self.targetSeq = self.targetSeq.translate(revComp)[::-1]
            self.querySeq = self.querySeq.translate(revComp)[::-1]
            self.compSeq = self.compSeq[::-1]
            self.tstrand = '1'
        else:
            self.negStrand = False
            self.tstrand = '0'
        
        #M5 is now always in + strand orientation
        self.qstrand = '0' if self.qstrand == '+' else '1'        
        self.flag = 0
        self.trim = False
    
    def toBed(self):
        """
        Returns the M5Line to a bed file
        """
        if self.negStrand:
            strand = "-"
        else:
            strand = "+"

        chrom = self.tname
        chromStart = str(self.tstart - self.qstart)
        chromEnd = str(self.tend + (self.qseqlength - self.qend))
        name = self.qname
        score = str(self.score)
        thickStart = str(self.tstart)
        thickEnd = str(self.tend)
        itemRgb = str(int(self.pctsimilarity))
        return "\t".join([chrom, chromStart, chromEnd, name, score, strand, \
                          thickStart, thickEnd, itemRgb])

    
    def __str__(self):
        """
        Undo changes 
        """
        #"""newBlasr
        if self.negStrand:
            self.targetSeq = self.targetSeq.translate(revComp)[::-1]
            self.querySeq = self.querySeq.translate(revComp)[::-1]
            self.compSeq = self.compSeq[::-1]
            self.tstrand = '-'
        else:
            self.tstrand = '+'
        
        self.qstrand = '+' if self.qstrand == '0' else '-'        
        return " ".join(map(str, [self.qname, self.qseqlength, self.qstart, \
                                  self.qend, self.qstrand, self.tname, \
                                  self.tseqlength, self.tstart, self.tend, \
                                  self.tstrand, self.score, self.nMatch, \
                                  self.nMismatch, self.nInsert, self.nDelete, \
                                  self.weirdUnknown, \
                                  self.querySeq, self.compSeq, self.targetSeq]))

class GapCans(dict):    
    
    def __init__(self, default=list):
        super(dict)
        self.update({"LeftContig": default(),
                 "RightContig": default(),
                 "SpansGap": default()})
    
    def consolidate(self):
        """
        Helper method to try and see if a single read supports
        the left and right contig of a gap, put it into the SpansGap category
        Uses the dictionaries to speed up time.
        """
        tmpl = {}
        for i in self["LeftContig"]:
            tmpl[i.qname] = i
        
        for pos,i in enumerate(self["RightContig"]):
            try:
                both = tmpl[i.qname]#They're in both
                #Add to spans
                self["SpansGap"].append(i)
                #remove from left and right
                del(tmpl[i.qname])
                del(self["RightContig"][pos])
            except KeyError:
                #Not in both
                pass
        
        self["LeftContig"] = tmpl.values() 
    
    def extend(self, gapCans):
        """
        Appends gapCans to self
        """
        if gapCans == None:
            return
        self["LeftContig"].extend(gapCans["LeftContig"])
        self["RightContig"].extend(gapCans["RightContig"])
        self["SpansGap"].extend(gapCans["SpansGap"])
        
    def getAllReads(self):
        ret = []
        for key in self.keys():
            ret.extend(self[key])
        return ret
    
    def isEmpty(self):
        for key in self.keys():
            if len(self[key]) > 0:
                return False
        return True
    
    def __str__(self):
        ret = ""
        for key in self:
            ret += key +":["+",".join(self[key])+"],"
        ret = ret[:-1]
        return ret

def GapCansDecode(obj):
    if obj.trim:
        extra = "##%d#%d##" % obj.trim
    else:
        extra = ""
    
    return obj.qname + extra

class GapSeqs(GapCans):
    def __init__(self):
        super(GapCans, GapCans.__init__(self, default=foo))

class LiftOverTable():
    """
    TODO: 
        Make the entry take care of updating it's own coordinates
        instead of relying on the calling code to calculate shifts.
    """
    def __init__(self, fn=None):
        #Quick Look on a per entry basis
        self.hash = {}
        #For outputting all scaffolds later
        #And keeping contigs contained within 
        #scaffolds for updates
        self.scaffoldRoots = {}
        #For adding
        self.curRoot = None
        if fn != None:
            self.__parse(fn)

    def __parse(self, fn):
        fh = open(fn,'r')
        head =  fh.readline()
        for line in fh.readlines():
            entry = LiftOverEntry(*line.strip().split('\t'))
            self.addEntry(entry)
        fh.close()

    def addEntry(self, entry):
        if not self.scaffoldRoots.has_key(entry.scaffold) or self.scaffoldRoots[entry.scaffold] == None:
            self.scaffoldRoots[entry.scaffold] = entry
            self.curRoot = entry
        else:
            self.curRoot.next = entry
            entry.prev = self.curRoot
            self.curRoot = entry
        
        key = entry.scaffold+str(entry.oStart)
        self.hash[key] = entry
            
    def getEntry(self, scaffold, oStart):
        return self.hash[scaffold+str(oStart)]
    
    def removeEntry(self, entry):
        if not self.scaffoldRoots.has_key(entry.scaffold):
            raise KeyError("Scaffold %s Not Found" % entry.scaffold)
        
        if entry.prev != None:
            entry.prev.next = entry.next
        
        if entry.next != None:
            entry.next.prev = entry.prev
        #I need a better hashing function
        #key = entry.scaffold+str(entry.oStart)
        #del(self.hash[key])
        
        #if this is the only guy in the entire scaffolding
        #Get rid of him
        if self.scaffoldRoots[entry.scaffold].next == None:
            del(self.scaffoldRoots[entry.scaffold])
    
    def updateScaffold(self, entry, shift):
        """
        Takes the key(entry) and changes all
        downstream coordinates as applicable
        """
        #entry.nStart += startShift
        #The previous end changes. 
        #With the start's shift
        #entry.prev.nEnd += startShift
        
        while entry.next != None:
            entry = entry.next
            if type(entry.nStart) == str:
                continue
            entry.nStart += shift
            entry.nEnd += shift

    def insertEntry(self, existingEntry, newEntry, after=True):
        """
        Inserts a new LiftOverEntry into the Table
        #Note! Changes made to the scaffold before entry is inserted
        are not applied to this entry - so get everything inserted
        before you do this
        """
        key = newEntry.scaffold+str(newEntry.oStart)
        self.hash[key] = newEntry
        
        if after == True:
            newEntry.prev = existingEntry 
            newEntry.next = existingEntry.next 
            if newEntry.next != None:
                newEntry.next.prev = newEntry
            existingEntry.next = newEntry
        else:
            newEntry.next = existingEntry
            newEntry.prev = existingEntry.prev
            if newEntry.prev != None:
                newEntry.prev.next = newEntry
            existingEntry.prev = newEntry
            
        
    def __str__(self):
        """
        returns the table as a string
        """
        ret = ""
        for entry in self:
            ret += str(entry) + "\n"
        
        return ret

    def __iter__(self):
        for key in self.scaffoldRoots.keys():
            root = self.scaffoldRoots[key]
            while root.next != None:
                yield root
                root = root.next
            yield root
        
class LiftOverEntry():
    def __init__(self, scaffold, oStart, oEnd, nStart, nEnd, gType, \
                 prev = None, next = None):
        self.scaffold = scaffold
        
        if oStart == 'na':
            self.oStart = oStart
            self.oEnd = oEnd
        else:
            self.oStart = int(oStart)
            self.oEnd = int(oEnd)

        if nStart == 'na':
            self.nStart = nStart
            self.nEnd = nEnd
        else:
            self.nStart = int(nStart)
            self.nEnd = int(nEnd)
        
        self.next = next
        self.gType = gType
        self.prev = prev
        
    def getNext(self, gType):
        """
        return the next feature that is gType
        """
        if self.next == None:
            return None
        if self.next.gType == gType:
            return self.next
        else:
            return self.next.getNext(gType)

    def getPrev(self, gType):
        """
        return the next feature that is gType
        """
        if self.prev == None:
            return None
        if self.prev.gType == gType:
            return self.prev
        else:
            return self.prev.getPrev(gType)
            
    def __str__(self):
        return "\t".join([self.scaffold, str(self.oStart), str(self.oEnd), \
                          str(self.nStart), str(self.nEnd), self.gType])
