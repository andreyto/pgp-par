UsageInfo = """ PeptidomeCrossRefs.py
This program is designed to create a tab delimited
bunch of data for the cross references in a peptidome
meta data file.  It assumes nice consistent naming
Required Options
 -r [Directory] Directory of spectra (recursive parse)
 -w [FileName] OutFile
"""


import os
import getopt
import sys


class AbacusClass:
    def __init__(self):
        self.SpectraDir = None
        self.OutputFile = None
        self.AllSpectra = []
        self.ApprovedExtensions = [".mzXML"] # you can add more if you like
        
    def Main(self):
        self.CollectAllSpectra(self.SpectraDir)
        self.WriteItOut()
               
        
    def CollectAllSpectra(self, Directory):
        """This is a wrapper for directory listing
        """
        for Entity in os.listdir(Directory):
            FullPath = os.path.join(Directory, Entity)
            if os.path.isdir(FullPath):
                self.CollectAllSpectra(FullPath)
                continue
            #now I'm guessing it's a file.
            (Path, FileName) = os.path.split(FullPath)
            #print "Before %s"%FileName
            FileName = FileName.replace(".gz", "")
            #print "After %s"%FileName
            (Stub, Ext) = os.path.splitext(FileName)
            if Ext in self.ApprovedExtensions:
                self.AllSpectra.append(FileName)
 
    def WriteItOut(self):
        """ put it in the format """
        Handle = open(self.OutputFile, "wb")
        Algorithm = "Inspect, PepNovo"
        for File in self.AllSpectra:
            (Stub, Ext) = os.path.splitext(File)
            PSMFile = "%s.rescore.pepxml"%Stub # i put rescore in here because all of my stuff is rescored with Ari's pepnovo
            Fraction = Stub
            Type = Ext[1:] # to get rid of the '.' that come in the extension, like .mzXML
            Handle.write("%s\t%s\t%s\t%s\t%s\n"%(Fraction, File, Type, PSMFile, Algorithm))
        Handle.close()

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SpectraDir = Value
            elif Option == "-w":
                self.OutputFile = Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w"):
            print UsageInfo
            sys.exit(1)
            

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    BeanCounter = AbacusClass()
    BeanCounter.ParseCommandLine(sys.argv[1:])
    BeanCounter.Main()                
    
