UsageInfo = """CharacterizeAlternativeSplicing.py
Looks at the sequence level differences between alternative isoforms from the
same locus for TAIR.  Then print out a histogram of the differences.

Required Options:
 -d [TrieFile] Database
 -c [Clustal] Path to the clustal Executable
 -r [FileName] Results to compare against the whole

"""

import os
import getopt
import sys
import SelectProteins

class AnalysisClass:
    def __init__(self):
        self.DatabasePaths = [] #allow multiple
        self.Histogram = {}
        self.ObservedHistogram = {}
        self.AlternativelySplicedLoci = {} # Integer accessions => dictionary of Seqeunce => IDs
        self.ClustalPath = None
        self.ClustalOutputPath = None
        self.OutputFilePath = "RenameYourOutput.txt"
        self.ResultsPath = None
        self.ObservedAltSpliceLoci = {} # key = locus, value = dummy
        self.KuldgeString = ""

    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DatabasePaths)
        self.FindLoci()
        self.ParseObservedLoci()
        #self.DebugChecker()
        ## make a path for all this output.
        print self.ClustalPath
        (Path, Exe) = os.path.split(self.ClustalPath)
        self.ClustalOutputPath = os.path.join(Path, "OutputDump")
        for Locus in self.AlternativelySplicedLoci.keys():
            #for Locus in ["AT1G65440", "AT1G78870", "AT2G16430", "AT3G59430", "AT4G21540", "AT5G60940" ]:
            self.WrapClustalW(Locus)
        self.PrintHistogram()
        print self.KuldgeString


    def ParseObservedLoci(self):
        """Given an input file, where the first thing is the name of a gene which
        is observed in two protein isoforms, make a dictionary of these things
        """
        if not self.ResultsPath:
            return
        Handle = open(self.ResultsPath, "rb")
        for Line in Handle.xreadlines():
            Bits = Line.strip().split(" ")
            if len(Bits) < 3:
                continue
            if Bits[1] == "must" and Bits[2] == "be":
                #a good line
                self.ObservedAltSpliceLoci[Bits[0]] = 1
        

    def PrintHistogram(self):
        Keys = self.Histogram.keys()
        Keys.sort()
        MaxKey = Keys[-1]
        Handle = open(self.OutputFilePath, "wb")
        Handle.write("Amino Acid Differences\tNumber of Isoform Pairs\tObservedPairs\n")
        for Key in range(MaxKey + 1):
            Handle.write( "%s\t%s\t%s\n"%(Key,self.Histogram.get(Key,0),self.ObservedHistogram.get(Key,0)))
        Handle.close()
        


    def WrapClustalW(self, Locus):
        """Get the dictionary of sequences corresponding to the locus, and then
        make a quicky file and clustal it.  Then we have to parse the output
        and figure out how much of the sequence it shares.
        """
        ## 1. create the fasta file
        Counter = 0
        Dictionary = self.AlternativelySplicedLoci[Locus]
        Length = len(Dictionary)
        Keys = Dictionary.keys()
        ##double for loop here to make sure that we only have pairwise compares
        for Index in range(Length):
            for Jndex in range(Index+1, Length):
                Counter += 1
                Seq1 = Keys[Index]
                ID1 = Dictionary[Seq1]
                Seq2 = Keys[Jndex]
                ID2 = Dictionary[Seq2]
                FilePath = os.path.join(self.ClustalOutputPath, "%s.%s.txt"%(Locus,Counter))
                Handle = open(FilePath, "wb")
                Handle.write(">%s\n%s\n"%(ID1, Seq1))
                Handle.write(">%s\n%s\n"%(ID2, Seq2))
                Handle.close()
                ## 2. Call clustal
                Command = "%s %s"%(self.ClustalPath, FilePath)
                print Command
                os.system(Command)
                ## 3. Parse the output
                AlignmentFile = os.path.join(self.ClustalOutputPath, "%s.%s.aln"%(Locus,Counter))
                NumAligned = self.ParseAlignment(AlignmentFile)
                ## 4. add to the histogram.  here the counting is actually difficult
                ## if we had    TTTTT and TTTT, the number of stars would be 4, and
                ## we could say that the real difference between the two is 1 amino acid
                ## so we add to the histogram the number: max(len(seq1), len(seq2)) - numaligned
                Max = len(Seq1)
                if len(Seq2) > Max:
                    Max = len(Seq2)
                AminoAcidDifferences = Max - NumAligned
                if AminoAcidDifferences > 200:
                    print "$$$$$$$$$$$ %s, %s"%(AminoAcidDifferences, Locus)
                if not self.Histogram.has_key(AminoAcidDifferences):
                    self.Histogram[AminoAcidDifferences] = 0
                    self.ObservedHistogram[AminoAcidDifferences] = 0
                self.Histogram[AminoAcidDifferences] += 1
                if self.ObservedAltSpliceLoci.has_key(Locus):
                    self.KuldgeString += "%s\t%s\n"%(Locus, AminoAcidDifferences)
                    self.ObservedHistogram[AminoAcidDifferences] += 1
                    

    def ParseAlignment(self, FilePath):
        """open a clustal Alignment file and parse it, trying to get out
        the number of identical amino acids between two sequences.  basically counting the * of the alignment
        Note.  this is imperfect, an underestimate of the actual differences.  just because of parsing hassle
        see this sequence
            AT1G01020.2      IGVLSANAAFIISFAIATKGLLNEVSRRREIMLGIFISSYFKIFLLAMLVCCS-------
            AT1G01020.1      IGVLSANAAFIISFAIATKGLLNEVSRRREIMLGIFISSYFKIFLLAMLVWEFPMSVIFF
                             **************************************************          

            AT1G01020.2      -----FTS---------------------------HLIPNIEVPNFLSIP----------
            AT1G01020.1      VDILLLTSNSMALKVMTESTMTRCIAVCLIAHLIRFLVGQIFEPTIFLIQIGSLLQYMSY
                                  :**                           .*: :*  *.:: *           
        all that star stuff on the bottom row is total crap, but that's to hard to parse out that
        it's not a real alignment, so I'm sticking with counting stars.        
        """
        Handle = open(FilePath, "rb")
        NumStars = 0
        for Line in Handle.xreadlines():
            NumStars += Line.count("*")
            #print Line
            #print Line.count("*")
        Handle.close()
        return NumStars
                         

    def DebugChecker(self):
        "do whatever I want"
        for Accession in self.AlternativelySplicedLoci.keys():
            Dictionary = self.AlternativelySplicedLoci[Accession]
            print "I have multiple protein isoforms at %s"%Accession
            for (Sequence, ID) in Dictionary.items():
                print self.ProteinPicker.ProteinNames[ID]
                print Sequence
            print "\n\n"

    def FindLoci(self):
        """look through the database and find loci with multiple
        protein isoforms, and then save it out to an array.
        """
        SimpleHash = {} # key = AT accession (integer part only)  value = ID
        for (ID, Fasta) in self.ProteinPicker.ProteinNames.items():
            FullAccession = Fasta.split()[0]
            Dot = FullAccession.find(".")
            IntegerAccession = FullAccession[:Dot]
            #print IntegerAccession
            if not SimpleHash.has_key(IntegerAccession):
                SimpleHash[IntegerAccession] = []
            SimpleHash[IntegerAccession].append(ID)
        ## 2. now go and find the distinct protein isoforms
        for (IntegerAccession, List) in SimpleHash.items():
            # we have a list of IDs, let's get their proteins, and see if it's in there
            SequenceHash = {} # key = sequence, value = ID
            for ID in List:
                Sequence = self.ProteinPicker.ProteinSequences[ID]
                if not SequenceHash.has_key(Sequence):
                    SequenceHash[Sequence] = ID
            if len(SequenceHash) > 1:
                ## here we have two things.  Hooray
                self.AlternativelySplicedLoci[IntegerAccession] = SequenceHash
                # i hope the above assignment doesn't over write stuff on multiple iterations.
                # we'll dump stuff as a check.
            

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:c:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ResultsPath = Value
            if Option == "-d":
                self.DatabasePaths.append(Value)
            if Option == "-c":
                self.ClustalPath = os.path.abspath(Value)
            if Option == "-w":
                self.OutputFilePath = Value
        if not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-c"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Hunt = AnalysisClass()
    Hunt.ParseCommandLine(sys.argv[1:])
    Hunt.Main()            