"""ParseClustalForIdentity.py
This script is a parser for the ClustalW program.  It
looks for sequence identity, marked as * in the .aln file, and checks to
see if a given alignment has X% identity.  it reports alignments
that fail the curoff
"""

UsageInfo="""QuantifySequenceIdentity.py
Find the sequence identity between all proteins in a database.
Keep track of pairs whose % identity exceeds a cutoff

Required Options:
 -c [Directory] Directory of clustal results
 -d [TrieFile] Reference database (S. cerevisiae)
 -w [FileName] Output from program

Additional Options:
 -m [Float] minimum percent identity, defaults to 0.9
"""

import os
import sys
import getopt
import SelectProteins

class WrapperClass:
    def __init__ (self):
        self.OutputPath = None
        self.ClustalResultsDir = None
        self.DBPath = []
        self.Cutoff = 0.5
        self.Histogram = {} # key => count
    def Main(self):
        ## load the reference database
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        ## parse clustal output
        self.RecurseOnDirectory(self.ClustalResultsDir, self.ParseAlignmentResults)
        ## print off results
        for Bin in range(101):
            Count = self.Histogram.get(Bin, 0)
            print "%s\t%s"%(Bin, Count)

    def RecurseOnDirectory(self, FilePath, Callback):
        if os.path.isdir(FilePath):
            FileNames = os.listdir(FilePath)
            for FileNameIndex in range(len(FileNames)):
                FileName = FileNames[FileNameIndex]
                (Stub, Extension) = os.path.splitext(FileName)
                if Extension.lower() not in (".aln"):
                    continue
                SubFilePath = os.path.join(FilePath, FileName)
                apply(Callback, (SubFilePath,))
        else:
            apply(Callback, (FilePath,))

    def GetReferenceLength(self, LocusToMatch):
        """Given a locus name, like YAL012C, we find the length of the protein
        """
        #lower case everything for comparison
        LocusToMatch.lower()
        for (ID, Fasta) in self.ProteinPicker.ProteinNames.items():
            Bits = Fasta.split(" ")
            Locus = Bits[0].lower()
            if Locus == LocusToMatch:
                return len(self.ProteinPicker.ProteinSequences[ID])
        print "No match for %s"%LocusToMatch
        return 0            


    def ParseAlignmentResults(self, FilePath):
        """basically counts the '*' in a file, returns
        Gets the length of the reference protein (the name of the file)
        MAJOR Assumption: you have no '*' in the protein name or sequence.  Doing such would be DUMB, and make the method fail
        """
        (PathStub, FileName) = os.path.split(FilePath)
        Locus = FileName.replace(".aln", "")
        ReferenceLength = self.GetReferenceLength(Locus)
        if not ReferenceLength:
            return
        Handle = open(FilePath, "rb")
        Alignment = Handle.read()
        Handle.close()
        Count = Alignment.count("*")
        FractionIdentical = Count / float(ReferenceLength)
        Bin = int(FractionIdentical * 100)
        #if Bin < 5:
        #    print Locus, FractionIdentical
        if not self.Histogram.has_key(Bin):
            self.Histogram[Bin] =0
        self.Histogram[Bin] += 1


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:w:c:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-c":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ClustalResultsDir = Value
            if Option == "-d":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath.append(Value)
            if Option == "-w":
                # -r results file(s)
                self.OutputPath = Value
            if Option == "-m":
                self.Cutoff = float (Value)
        if not OptionsSeen.has_key("-d")  or not OptionsSeen.has_key("-w"):
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