"""QuantifySequenceIdentity.py
This script is a wrapper and parser for the ClustalW program.  It takes
two proteins from a database, makes a temp.fasta and has clustal align it.  Then
it looks for sequence identity, marked as * in the .aln file, and checks to
see if a given alignment has X% identity.  it keeps a hash of all such alignments
the exceed the cutoff
"""

UsageInfo="""QuantifySequenceIdentity.py
Find the sequence identity between all proteins in a database.
Keep track of pairs whose % identity exceeds a cutoff

Required Options:
 -d [TrieFile] Database of protein sequences, use multiple times if needed
 -w [FileName] Output from program

Additional Options:
 -c [Float] Cutoff for percent identity, defaults to 0.9
"""

import os
import sys
import getopt
import SelectProteins

class WrapperClass:
    def __init__ (self):
        self.DBPath = []
        self.DB = None
        self.OutputPath = None
        self.Cutoff = 0.9
        self.TempFileName = "temp.delme.fasta"
        self.ProteinAssociations = {} # protein1 ->[homolog1, homolog2,..]
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        self.OutHandle = open(self.OutputPath, "wb")
        self.PerfectHandle = open("%s.perfect.txt"%self.OutputPath, "wb")
        self.CompareAllProteins()
        self.OutputResults()

    def OutputResults(self):
        for Name in self.ProteinAssociations.keys():
            self.OutHandle.write("%s is at least %f identical to %d proteins \n"%(Name, self.Cutoff, len(self.ProteinAssociations[Name])))
            self.OutHandle.write("\t%s\n"%self.ProteinAssociations[Name])
    def CompareAllProteins(self):
        """double for loop, get two proteins, make a temp.fasta
        call clustal and parse results
        """
        NumProteins = len(self.ProteinPicker.ProteinSequences)
        for Index in range(NumProteins):
            Sequence1 = self.ProteinPicker.ProteinSequences[Index]
            for Jndex in range(Index+1, NumProteins):
                Sequence2 = self.ProteinPicker.ProteinSequences[Jndex]
                self.MakeTempFastaFile(Sequence1, Sequence2)
                Command = "clustalw.exe %s"%self.TempFileName
                os.system(Command)
                IdenticalResidues = float(self.ParseAlignmentResults()) # cast to float for the division 
                FractionSeq1 = IdenticalResidues / len(Sequence1)
                FractionSeq2 = IdenticalResidues / len(Sequence2)
                Name1 = self.ProteinPicker.ProteinNames[Index]
                Name2 = self.ProteinPicker.ProteinNames[Jndex]
                if FractionSeq1 >= self.Cutoff:
                    #make a listing in the results
                    if FractionSeq1 == 1.0:
                        self.PerfectHandle.write("%f identity %s to %s\n"%(FractionSeq1, Name1, Name2))
                    if not self.ProteinAssociations.has_key(Name1):
                        self.ProteinAssociations[Name1] = []
                    self.ProteinAssociations[Name1].append(Name2)
                if FractionSeq2 >= self.Cutoff:
                    #make a listing in the results
                    if FractionSeq2 == 1.0:
                        self.PerfectHandle.write("%f identity %s to %s\n"%(FractionSeq2, Name2, Name1))
                    if not self.ProteinAssociations.has_key(Name2):
                        self.ProteinAssociations[Name2] = []
                    self.ProteinAssociations[Name2].append(Name1)
                self.CleanUp()
            break

    def CleanUp(self):
        """removes the files from the alignment"""
        Fasta = self.TempFileName
        Align = self.TempFileName.replace(".fasta", ".aln")
        Aux = self.TempFileName.replace(".fasta", ".dnd")
        os.remove(Fasta)
        os.remove(Align)
        os.remove(Aux)


    def MakeTempFastaFile(self, Sequence1, Sequence2):
        """Makes a temp.fasta file from the two sequences given"""
        Handle = open(self.TempFileName, "wb")
        Handle.write(">Sequence1\n%s\n>Sequence2\n%s\n"%(Sequence1, Sequence2))
        Handle.close()

    def ParseAlignmentResults(self):
        """basically counts the '*' in a file, returns
        MAJOR Assumption: you have no '*' in the protein name or sequence.  Doing such would be DUMB, and make the method fail
        """
        AlignmentResult = self.TempFileName.replace(".fasta", ".aln")
        Handle = open(AlignmentResult, "rb")
        Alignment = Handle.read()
        Handle.close()
        Count = Alignment.count("*")
        return Count


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:w:c:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
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
            if Option == "-c":
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