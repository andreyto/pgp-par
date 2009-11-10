#!/usr/bin/env python

UsageInfo = """DiffInspect.py
Find the differences between two Inspect runs that should be the 
same.  We assume the same databases and spectral input, and are
merely trying to register the raw Inspect output. 

Required Parameters
 -a [FileName] File or directory of Inspect results.
 -b [FileName] File or directory of Inspect results.
 -w [FileName] Output of the test
 
Additional Parameters
 -1 Only check the top hit per spectrum. Default to check all answers.

"""

import sys
import os
import getopt
import traceback
import ResultsParser
from Utils import *
Initialize()


class SpectrumResult(ResultsParser.ResultsParser):
    """This tiny object is meant to hold inspect results
    and facilitate their comparison.
    """
    def __init__(self):
        self.PSMs = {} #key = rank of PSM, value = bits of line
        ResultsParser.ResultsParser.__init__(self)
    def AddPSM(self, Rank, Bits):
        if self.PSMs.has_key(Rank):
            print "ERROR: trying to add a PSM with a rank already in use."
            return
        self.PSMs[Rank] = Bits
    def GetAminos(self):
        """Parameters: none
        Return: a list of Aminos for all these PSMs
        Description: This gives you quick access to some of the result.
        if it so happens that this is all you care about
        """
        List = []
        for (Rank, Bits) in self.PSMs.items():
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            List.append( Peptide.Aminos)
        return List
    
    def GetAminosToRank(self):
        """Parameters: none
        Return: dictionary Aminos->Rank
        Description: make this dictionary for you if you find it useful
        """
        Dict = {}
        for (Rank, Bits) in self.PSMs.items():
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            Dict[Peptide.Aminos] = Rank
        return Dict



class FinderClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ResultsAPath = None
        self.ResultsBPath = None 
        self.ResultsADictionary = {} #key = (file, spectrumNum) value = SpectrumResult
        self.ResultsBDictionary = {}
        self.CurrentResults = None
        self.OutputPath = "RenameYourOutput.txt"
        self.TopPSMOnly = 0 #whether to parse only the best psm
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.CurrentResults = self.ResultsADictionary #totally cheesy flag, because
        # I can't pass anything else into the parse method
        self.ProcessResultsFiles(self.ResultsAPath, self.ParseInspectCallback)
        self.CurrentResults = self.ResultsBDictionary
        self.ProcessResultsFiles(self.ResultsBPath, self.ParseInspectCallback)
        print "I got %s from A and %s from B"%(len(self.ResultsADictionary), len(self.ResultsBDictionary))
        self.CompareResultSets(self.ResultsADictionary, self.ResultsBDictionary)

    def CompareResultSets(self, DictionaryA, DictionaryB):
        """Parameters: Two dictionaries with SpectrumResults objects as values
        Return: Nothing
        Description: go through all the results and figure out if they differ
        """ 
        DebugFile = "A623A003.mzXML"
        DebugSpectrum = "786"
        Handle = open(self.OutputPath, "wb")
        for Key in DictionaryA.keys():
            ResultA = DictionaryA[Key]
            (File, Spectrum) = Key
            if not DictionaryB.has_key(Key):
                Handle.write("Results set B has no values for %s, %s\n"%(File, Spectrum))
                continue # no more can be done on this spectrum
            ResultB = DictionaryB[Key]
            #now how deep do we want to go.
            AminosA = ResultA.GetAminos()
            AminosB = ResultB.GetAminos()
            LosersOfA = self.ComparePeptideLists(AminosA, AminosB) #those in A without a match in B
            NumDisagreements = len(LosersOfA)
            if NumDisagreements: #shortcut because I expect zero
                Handle.write( "I got %s disagreements for %s, %s\n"%(NumDisagreements, File, Spectrum))
                LosersOfB = self.ComparePeptideLists(AminosB, AminosA) #those in B without a match in A
                self.DeepComparison(ResultA, ResultB, LosersOfA, LosersOfB, Handle)
                
        Handle.close()
  
    def DeepComparison(self, ResultA, ResultB, LosersOfA, LosersOfB, Handle):
        """Parameters: Two SpectrumResults Objects, two lists of peptides absent from the comparative set
            and an open file handle
        Return: None
        Description: we got two sets of PSMs that are not the same. Here
        we want to take a better look at the differences
        """
        #Bits that I care about: MQScore Vector Components cols 7-12, not in the ResultsParser.Columns
        #also columns 5 and 6 which are score and length
        ColumnsOfConcern = range(5, 13)
        #first what may seem duplicitous, just to get a handle on what is different
        AminosToRankA = ResultA.GetAminosToRank()
        AminosToRankB = ResultB.GetAminosToRank()
        #so - LosersOfA are those peptides in A without a match in B, let's print out everything
        ##relevant on those peptides.
        StringA = self.ListToString(LosersOfA)
        StringB = self.ListToString(LosersOfB)
        Handle.write("The following peptides do not overlap between the results sets\n")
        Handle.write("From A: %s\n"%StringA)
        Handle.write("From B: %s\n"%StringB)
        #now let's check out the rest of the hits to make sure that everything is same-same
        #in this loop we are looking at aminos which are the same between the two results sets
        # having already printed out our loosers.
        for (RankA, BitsA) in ResultA.PSMs.items():
            AnnotationA = BitsA[self.Columns.Annotation]
            PeptideA = GetPeptideFromModdedName(AnnotationA)
            AminosA = PeptideA.Aminos
            if AminosA in LosersOfA:
                continue #we've done what we need with these guys
            # a totally chicken hack.  I just can't deal with this now.
            if not AminosToRankB.has_key(AminosA):
                continue
            RankB = AminosToRankB[AminosA]
            BitsB = ResultB.PSMs[RankB]
            #now let's check the bits that we care about
            Okay = 1
            for Col in ColumnsOfConcern:
                BitA = BitsA[Col]
                BitB = BitsB[Col]
                if not BitA == BitB:
                    Okay = 0
                    break
            #now print if there's an issue
            if not Okay:
                StringBitsA = self.ListToString(BitsA[5:13])
                StringBitsB = self.ListToString(BitsB[5:13])
                Handle.write("Discrepencies for the annotation of %s\n"%AminosA)
                Handle.write("%s\n"%StringBitsA)
                Handle.write("%s\n"%StringBitsB)
               
                    

    def AmbiguousAminoAcidReplacement(self, String):
        NewString = String.replace("I", "L")
        NewString = NewString.replace("Q", "K")
        return NewString
            
            
    def ListToString(self, List):
        """Paramters: a List with simple values
        Return: a string 
        Description: turns it into a string for printing
        """
        String = ""
        for Item in List:
            String += "%s, "%Item
        String = String[:-2] #chomp off the 2 last characters
        return String

        
    def DictionaryToString(self, Dictionary):
        """Paramters: a dictionary with simple key value pairs
        Return: a string 
        Description: turns it into a string for printing
        """
        Keys = Dictionary.keys()
        Keys.sort()
        String = ""
        for Key in Keys:
            Value = Dictionary[Key]
            String += "\t%s -> %s\n"%(Key, Value)
        String = String[:-1] #chomp off the last character
        return String
    
    def ComparePeptideLists(self, ListA, ListB):
        """Parameters: Two lists of amino acid strings, no order required
        Return: List of disagreeing string
        Description: simply compare the two lists.  I expect exactly
        the same thing.
        """
        Disagreements = []
        for Item in ListA:
            if not Item in ListB:
                #now we don't find it.  WE should try some L->I and K->Q conversions
                FoundWithReplacements = self.FindWithReplacements(Item, ListB)
                if not FoundWithReplacements:
                    Disagreements.append(Item)
                
        return Disagreements


    def FindWithReplacements(self, Item, List):
        """Parameters: A string, a list of strings
        Return: 0/1 found or not
        Description: WE can't find the string in the list, but are
        curious about if any isomers of it may exist.  You see Inspect
        can't tell the difference between I and L, or Q and K.  So let's do
        some replacement, and see how things work.
        """
        NewList = []
        NewItem  = self.AmbiguousAminoAcidReplacement(Item)
        for ListMember in List:
            NewString = self.AmbiguousAminoAcidReplacement(ListMember)
            NewList.append(NewString)
        #now check to see if it's there
        if NewItem in NewList:
            return 1
        return 0
       
    def ParseInspectCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out Inspect Results to get peptide annotations, 
        Putting them in a hash for safe keeping
        """
        Handle = open(FilePath, "rb")
        Rank = 0 #we use zero based counting, like any good programmer
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            if Line[0] == "#":
                continue # comment line
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                SpectrumPath = Bits[self.Columns.SpectrumFile]
                ScanNumber = Bits[self.Columns.ScanNumber]
            except:
                traceback.print_exc()
                continue # SNAFU
            #now we see where to put it
            (Path, File) = os.path.split(SpectrumPath)
            Key = (File, ScanNumber)
            if not self.CurrentResults.has_key(Key):
                #this should be the first PSM of a spectrum.
                #let's make it and set the clocks back to zero.
                self.CurrentResults[Key] = SpectrumResult()
                Rank = 0
            #now let's put this in
            self.CurrentResults[Key].AddPSM(Rank, Bits)
            Rank += 1
        Handle.close()

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "a:b:w:1")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-a":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ResultsAPath = Value
            elif Option == "-b":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ResultsBPath = Value
            elif Option == "-1":
                self.TopPSMOnly = 1
            elif Option == "-w":
                self.OutputPath = Value
            else:
                print "Option %s not recognized.  Fatal Error."%Option
                print UsageInfo
                sys.exit(1)
        if not OptionsSeen.has_key("-a") or not OptionsSeen.has_key("-b"):
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
