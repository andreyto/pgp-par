UsageInfo = """EvaluateSignalP.py
Compare the predictions made by signalp to proteomic data.
See what proteomics confirms and what it rejects.

Required Options
 -r [FileName] File with the first obseved peptide for each protein
 -s [FileName] Signal peptide predictions from signalp
 -c [FileName] File with the proteomically called signal cleavage sites
 -w [FileName]  Output file
 
 
"""

import sys
import os
import getopt
import traceback
import ResultsParser
import BasicStats
from Utils import *
Initialize()


class FinderClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.CleavageResults = None
        self.PredictedSignalPFile = None
        self.OutputPath = "RenameYourOutput.txt"
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        CleavageHash = self.ParseProteomicsResults(self.CleavageResults)
        ProteomicsHash = self.ParseProteomicsResults(self.ReferenceResults)
        SignalPHash = self.ParseSignalP(self.PredictedSignalPFile)
        print "I have %s cleavage"%len(CleavageHash)
        print "I have %s proteomics"%len(ProteomicsHash)
        print "I have %s signalp"%len(SignalPHash)
        self.GraphItOut(CleavageHash, ProteomicsHash, SignalPHash)
        
    def GraphItOut(self, CleavageHash, ProteomicsHash, SignalPHash):
        """Parameters: three hashes gi->number
        """
        Accept = {} #score =>count
        Reject = {} #score => count
        AllSignalP = {} #score => count
        for GI in SignalPHash.keys():
            (SignalPStart, SignalPProb) = SignalPHash[GI]
            ScoreBin = round(SignalPProb, 2) #gets me two decimals of precision
            if not AllSignalP.has_key(ScoreBin):
                #set all of them up
                AllSignalP[ScoreBin] = 0
                Accept[ScoreBin] = 0
                Reject[ScoreBin] = 0
            AllSignalP[ScoreBin] += 1
            #do we have this in our cleavage set (proteomically confirmed signal peptide cleavage)
            if CleavageHash.has_key(GI):
                ProteomicsStart = CleavageHash[GI]
                print "%s verified by proteomics signal peptide cleavage." %GI
                print  "\tProteomics Start %s, SignalP start %s (score %s)"%(ProteomicsStart, SignalPStart, SignalPProb)
                Accept[ScoreBin] += 1
            elif ProteomicsHash.has_key(GI):
                #we don't have a clear signal peptide cleavage, but we did find something
                ProteomicsStart = ProteomicsHash[GI]
                
                if (SignalPStart - ProteomicsStart) > 3:
                    if SignalPProb > 0.4:
                        print "%s verified protein existence." %GI
                        print "\t##We disagree"
                        print "\tProteomics Start %s, SignalP start %s (score %s)"%(ProteomicsStart, SignalPStart, SignalPProb)
                    Reject[ScoreBin] += 1
                
        #end of the for
        HistogramList = []
        HistogramList.append(AllSignalP)
        HistogramList.append(Accept)
        HistogramList.append(Reject)
        Handle = open(self.OutputPath, "wb")
        BasicStats.PrettyPrintMultiHistogram(HistogramList, Handle)
        
    def ParseProteomicsResults(self, Path):
        """Parameters: Path to a file
        Return: Dictionary with Accession->(firstobservedresidue)
        Description: These files are from the MapNTerminalPeptides.py and
        massaged gently by me.  
        """
        Handle = open(Path, "rb")
        Dictionary  = {}
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            Bits = Line.split("\t")
            ProteinName = Bits[0]
            GI = self.JustGI(ProteinName)
            LenOfPrefix = int (Bits[1])
            if not Dictionary.has_key (GI):
                Dictionary[GI] = LenOfPrefix
            else:
                CurrLen = Dictionary[GI] 
                if LenOfPrefix < CurrLen:
                   Dictionary[GI] = LenOfPrefix
            
        return Dictionary 
         
    def ParseSignalP(self, Path):
        """This function parses the signal p prediction file, which is a list of gi numbers.
        I just have to map those the the proteins in our database, and then return the IDs.
        gi|22125314|ref|NP_668737.1|    0.554    30
        gi|22127555|ref|NP_670978.1|    0.831    21
        gi|22127683|ref|NP_671106.1|    0.792    20
        gi|22125472|ref|NP_668895.1|    0.866    28
        ## warning this number is 1 based, not zero based like we often use.  that's the reason for
        subtracting zero
        """
        Dictionary = {}
        Handle = open(Path, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            Line = Line.strip()
            (NameStub, Probability, ProteinStart) = Line.split("\t") #and some other crap that should also match.
            ProteinStart = int (ProteinStart)
            Probability = float (Probability)
            GI = self.JustGI(NameStub)
            Dictionary[GI] = (ProteinStart, Probability)
        return Dictionary


    def JustGI(self, Accessions):
        """Parameters: string with some | separated stuff
        Return: the gi number (2nd thing)
        Description: SignalP does not give me the full name, so I'm taking this
        """
        Bits = Accessions.split("|")
        return Bits[1]
 
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:c:w:s:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceResults = Value
            if Option == "-c":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.CleavageResults = Value
            if Option == "-s":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.PredictedSignalPFile = Value
                    
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-c") or not OptionsSeen.has_key("-s"):
            print UsageInfo
            sys.exit(1)




if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gumshoe = FinderClass()
    Gumshoe.ParseCommandLine(sys.argv[1:])
    Gumshoe.Main()
    