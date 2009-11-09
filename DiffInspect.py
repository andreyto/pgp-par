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
    def GetPeptides(self):
        """Parameters: none
        Return: a list of peptides for all these PSMs
        Description: This gives you quick access to some of the result.
        if it so happens that this is all you care about
        """
        List = []
        for (Rank, Bits) in self.PSMs.items():
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            List.append( Peptide.Aminos)
        return List



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
            PeptidesA = ResultA.GetPeptides()
            PeptidesB = ResultB.GetPeptides()
            Disagreements = self.ComparePeptideLists(PeptidesA, PeptidesB)
            if Disagreements: #shortcut because I expect zero
                Handle.write( "I got %s disagreements for %s, %s\n"%(Disagreements, File, Spectrum))
                #now that we have disagreements, we should look more closely
                PeptidesA.sort()# I sort here because I don't want to waste the time otherwise
                PeptidesB.sort()
                StringA = self.ListToString(PeptidesA)
                StringB = self.ListToString(PeptidesB)
                Handle.write("%s\n"%StringA)
                Handle.write("%s\n"%StringB)
        Handle.close()
        

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
        Return: Int - the number of disagreements
        Description: simply compare the two lists.  I expect exactly
        the same thing.
        """
        Disagreements = 0
        for Item in ListA:
            if not Item in ListB:
                #now we don't find it.  WE should try some L->I and K->Q conversions
                FoundWithReplacements = self.FindWithReplacements(Item, ListB)
                if not FoundWithReplacements:
                    Disagreements += 1
                
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
        NewItem = Item.replace("I", "L")
        NewItem = NewItem.replace("Q", "K")
        for ListMember in List:
            NewString = ListMember.replace("I", "L")
            NewString = NewString.replace("Q", "K")
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
