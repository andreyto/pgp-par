UsageInfo = """ GelProteinAnalysis.py
This program is designed for processing the results of 2D gels, in that
it looks for the main protein found in each mzXML file, assuming that 
each mzXML file corresponds to a single protein.
Required Options
 -r [FileName] Filename of Inspect Results
 -d [Directory] Where to put individual results files
"""


import os
import getopt
import sys
import ResultsParser


class AbacusClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFile = None
        self.OutputDir = None
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProcessResultsFiles(self.InputFile, self.ParseFile)
               
        
    def ParseFile(self, FileName):
        """This is a wrapper for the Assess ProteinComponents
        """
        CurrentOutFile = None
        OutHandle = None
        Handle = open(FileName, "rb")
        LineCount = 0
        for Line in Handle.xreadlines():
            LineCount += 1
            if LineCount % 1000 == 0:
                print "Parsing Line %s"%LineCount
            if not Line.strip():
                continue
            Line = Line.strip()
            if Line[0] == "#":
                continue
            Bits = Line.split("\t")
            SpectrumFilePath = Bits[self.Columns.SpectrumFile]
            OutFile = self.GetFileName(SpectrumFilePath)
            print "%s to %s"%(SpectrumFilePath, OutFile)
            continue
            #check to see if we're in the same file
            if OutFile == CurrentOutFile:
                #still in the same place, write this line
                OutHandle.write("%s\n"%Line)
            else: 
                #close up and start over
                if CurrentOutFile: # as long as its not 'None'
                    OutHandle.close()
                #make a new full path
                CurrentOutFile = OutFile
                FullPath = os.path.join(self.OutputDir, OutFile)
                OutHandle = open (FullPath, "wb")
                #now write
                OutHandle.write("%s\n"%Line)
        Handle.close()

    def GetFileName(self, SpectrumFilePath):
        """Parameters: path to a spectrum file
        Return: the name of an appropriate output file 
        Description: Given a spectrum path /path/to/f1.mzxml, we
        return what a normal inspect output filename would be, f1.txt
        """
        (Path, FileName) = os.path.split(SpectrumFilePath)
        (Stub, Ext) = os.path.splitext(FileName)
        NewOut = "%s.txt"%Stub
        return NewOut

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFile = Value
            elif Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find output directory '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.OutputDir = Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-d"):
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
    
