UsageInfo="""BiggsetParentMass.py

Required Options:
 -m [Directory] Directory (not file) where spectra files are stored


"""


import os
import sys
import getopt
import xml.sax.handler
import xml.sax

## auxiliary for the mzxml files
class XMLHandler(xml.sax.handler.ContentHandler):
    def __init__(self):
        self.inOffset = 0
        self.ParentMasses = []
        self.CID2PQDMapping = {} # key = CID ScanNum, value = PQD Scan
    def startElement(self, name, attributes):
        if name == "precursorMz":
            self.buffer = ""
            self.inOffset = 1
                
    def characters(self, data):
        if self.inOffset:
            self.buffer += data
    def endElement(self, name):
        if name == "precursorMz":
            self.inOffset = 0
            self.ParentMasses.append(float( self.buffer))

class AbacusClass():
    def __init__ (self):
        self.SpectraDir = None
        self.MaxPM = 0
        
        
    def Main(self):
        self.ParseSpectraDir()
        print "The biggest is %s"%self.MaxPM


    def ParseSpectraDir(self):
        for File in os.listdir(self.SpectraDir):
            (FileName, Extension) = os.path.splitext(File)
            Extension = Extension.lower()
            if not Extension == ".mzxml":
                continue
            print File
            FullPath = os.path.join(self.SpectraDir, File)
            self.ParseSpectrumFile(FullPath)

    def ParseSpectrumFile(self, FilePath):
        """Given a mzxml file, read through and get the spectrum pairs
        of neighboring scans with the same parent mass.  Put them in a
        dictionary of some kind.  Also make sure to get the file offset for such files.
        """
        Parser = xml.sax.make_parser()
        Handler = XMLHandler()
        Parser.setContentHandler(Handler)
        Parser.parse(FilePath)
        AllParentMasses = Handler.ParentMasses
        AllParentMasses.sort()
        print "First %f and last %f"%(AllParentMasses[0], AllParentMasses[-1])
        if AllParentMasses[-1] > self.MaxPM:
            self.MaxPM = AllParentMasses[-1]

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "m:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-m":
                # -m spectra file(s)
                self.SpectraDir = Value
        if not OptionsSeen.has_key("-m") :
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    DoStuff = AbacusClass()
    DoStuff.ParseCommandLine(sys.argv[1:])
    DoStuff.Main()        
