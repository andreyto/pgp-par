UsageInfo="""ParseSignalP.py
THis script parses signalp output and returns a list of proteins
with their score, and predicted site of cleavage

Required Options:
 -d [Directory] directory with a bunch of .signalp.raw files 
 -w [FileName] Output of the program
"""

import os
import sys
import getopt
import SelectProteins

class WrapperClass:
    def __init__ (self):
        self.Directory = None
        self.OutFileName = "RenameYourOutput.txt" 
        
    def Main(self):
        OutHandle = open(self.OutFileName, "wb")
        FileNames = os.listdir(self.Directory)
        for FileName in FileNames:
            (Stub, Ext) = os.path.splitext(FileName)
            if not Ext == ".raw":
                continue
            FullPath = os.path.join(self.Directory, FileName)
            (Name, Score, PlusOneResidue) = self.ParseSingleFile(FullPath)
            if not Name:
                continue
            #print "%s\t%s\t%s"%(Name, Score, PlusOneResidue)
            OutHandle.write("%s\t%s\t%s\n"%(Name, Score, PlusOneResidue))

    def ParseSingleFile(self, FileName):
        """Looking for specific lines.  Very forced parsing
        """
        SignalPHandle = open (FileName, "rb")
        FastaLine = None
        DString = "       D"
        DLine = None
        for Line in SignalPHandle.xreadlines():
            if Line[0] == ">" and (not FastaLine): # not set yet
                FastaLine = Line
            #now looking for a line about the neural net prediction from D
            First8 = Line[:8]
            if First8 == DString:
                DLine = Line
                break
        DBits = DLine.split()
        # "D", "length of signal peptide e.g. 1-34", "probability", "cutoff", "YES/NO"
        #if DBits[-1] == "NO":
        #    return (None, None, None) #return junk
        SignalPeptideAminos = DBits[1]
        LastResidue = SignalPeptideAminos.split("-")[1]
        PlusOneResidue = int(LastResidue) + 1
        Name = FastaLine.split(" ")[0]
        Name = Name.replace(">", "")
        return (Name, DBits[2], PlusOneResidue)

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:w:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                # -r results file(s)
                if not os.path.isdir(Value):
                    print "** Error: '%s' is not a directory\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.Directory = Value
            if Option == "-w":
                # -r results file(s)
                self.OutFileName = Value
                
        if not OptionsSeen.has_key("-d")  :
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    DoStuff = WrapperClass()
    DoStuff.ParseCommandLine(sys.argv[1:])
    DoStuff.Main()        