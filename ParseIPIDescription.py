UsageInfo="""ParseIPIDescription.py
IPI database compiles proteins from all over, and then puts them into
a single place, sometimes the descriptions don't make it into the
Inspect output, so here we are going to parse that from the .dat file
which is in the format of swiss prot/trembl

ID   IPI123578906433789
...
DE   Some free text description.

Required Options:
 -r [FileName] .dat file for an IPI species
 -w [FileName] Output from program


"""

import os
import sys
import getopt

class AbacusClass():
    def __init__ (self):
        self.IPIPath = None
        self.OutputPath = None
    def Main(self):
        self.ParseDatabase()


    def ParseDatabase(self):
        """ID   IPI00516214.1         IPI;      PRT;   253 AA.
        it appears that there is a consistent number of spaces after ID.
        so Line[5:] should return the rest, which can then be easily split
        DE   PATHOGENESIS-RELATED THAUMATIN FAMILY PROTEIN.
        """
        Accession = None
        Description = ""
        Handle = open(self.IPIPath, "rb")
        OutHandle = open(self.OutputPath, "wb") 
        for Line in Handle.xreadlines():
            if Line[:2] == "ID":
                #new ID found, store old and reset description
                if Accession:
                    OutHandle.write("%s\t%s\n"%(Accession, Description))
                    Description = ""                    
                RestOfLine = Line.strip()[5:]
                Bits = RestOfLine.split(" ")
                #print Bits[0] # should be the accession
                Accession = Bits[0]
                if not Accession[:3] == "IPI":
                    print "SNAFU parsing Line:",Line
            if Line[:2] == "DE":
                RestOfLine = Line.strip()[5:]
                Description += RestOfLine
        Handle.close()
        OutHandle.close()

                
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "w:r:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-w":
                # -r results file(s)
                self.OutputPath = Value
            if Option == "-r":
                self.IPIPath = Value
        if not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-r"):
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
