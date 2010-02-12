UsageInfo = """BacterialProteinInference.py
Given a set of annotation for some spectra, and a protein database
we report the proteins present.  We do very little in the way of
inference.  We do not resolve non-uniquely mapping peptides. Thus
this is best used for bacteria, where the problem is minimal.

Required Options
 -r [FileName] File or directory of results.
 -d [Trie file] Proteome Database used for search
 -w [FileName]  Output file
 
Additional Options
 -v   Flag to print a more verbose trace of the program progression
     through various methods and whatnot
 -p [float] Pvalue cutoff (Default 0.05)
 
"""


import getopt
import traceback
import InspectResults
import SelectProteins
from Utils import *
Initialize()


class FinderClass():
    def __init__(self):
        self.ReferenceResults = None
        self.OutputPath = "RenameYourOutput.txt"
        self.DatabasePaths = [] #possibly multiple
        self.AllPeptides = {} # AminoSequence -> sPeptide Object
        self.PValueLimit = 0.05 #pvalues are 0.00 (good) to 1.0 (bad) 
        self.Verbose = 0
        self.NameToID = {}
        self.SpectrumCount = 0
        self.MinPeptides = 2
        
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DatabasePaths)
        self.ParseInspect(self.ReferenceResults)
        print "I found %s peptides from %s spectra"%(len(self.AllPeptides), self.SpectrumCount)
        Proteins = self.PutPeptidesOnProteins()
        print "I found %s proteins (with 1 or more peptides)"%len(Proteins)
        self.PrettyPrintProteins(Proteins)
        
    def PrettyPrintProteins(self, Proteins):
        """Parameters: Dictionary of proteins (name->list of peptide objects)
        Return:None
        Description: WE print out first the proteins with unique peptides, and then those than
        only have shared peptides
        """
        ThoseWithUnique = []#simple list of protein names, that we use to key the dict later
        ThoseWithoutUnique = []
        ToPrintCount = 0
        
        for ProteinName in Proteins.keys():
            ListOfPeptides = Proteins[ProteinName]
            if len(ListOfPeptides) < self.MinPeptides:
                continue
            ToPrintCount += 1
            HasUnique = 0
            for Peptide in ListOfPeptides:
                if Peptide.Unique:
                    HasUnique = 1
                    #break
            if HasUnique:
                ThoseWithUnique.append(ProteinName)
            else:
                ThoseWithoutUnique.append(ProteinName)
        print "%s proteins pass the minimum peptide cutoff (%s with a unique peptide)"%(ToPrintCount, len(ThoseWithUnique))
        #now we go and print stuff off
        Handle = open(self.OutputPath, "wb")
        Header = "ProteinName\tNum Peptides\tUnique Count\tUniquePeptides\tSharedPeptides"
        Handle.write("%s\n"%Header)
        for ProteinName in ThoseWithUnique:
            UniqueString = ""
            NonUniqueString = ""
            Count = 0
            UniqueCount = 0
            for Peptide in Proteins[ProteinName]:
                Count += 1
                if Peptide.Unique:
                    UniqueString += "%s, "%Peptide.Aminos
                    UniqueCount += 1
                else:
                    NonUniqueString += "%s, "%Peptide.Aminos
            #now dump crap
            Line = "%s\t%s\t%s\t%s\t%s\n"%(ProteinName, Count, UniqueCount, UniqueString, NonUniqueString)
            Handle.write(Line)
        #now those unfortunate 
        for ProteinName in ThoseWithoutUnique:
            NonUniqueString = ""
            Count = 0
            for Peptide in Proteins[ProteinName]:
                Count += 1
                NonUniqueString += "%s,"%Peptide.Aminos
            #now dump crap
            Line = "%s\t%s\t%s\t%s\t%s\n"%(ProteinName, Count, 0, "", NonUniqueString)
            Handle.write(Line)
        Handle.close()
       
        
        

    def PutPeptidesOnProteins(self):
        """Get the location for the start of the peptide, and then just keep track of the smallest one
        per protein
        """
        PeptideObjectInProtein = {} #protein->list of objects
        for Aminos in self.AllPeptides.keys():
            DBLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
            #now some true hackery.  WE make a variable in the peptide object
            #which is not defined in the object itself, but we make it here and
            #use it here.  So it's not generally useful.  Work on that Sam.
            if len(DBLocations) > 1:
                self.AllPeptides[Aminos].Unique = 0
            else:
                self.AllPeptides[Aminos].Unique = 1
            for (ProteinID, PeptideStartAA) in DBLocations:
                ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                self.NameToID[ProteinName] = ProteinID
                if not PeptideObjectInProtein.has_key(ProteinName):
                    PeptideObjectInProtein[ProteinName] = []
                PeptideObjectInProtein[ProteinName].append(self.AllPeptides[Aminos])
        return PeptideObjectInProtein
        
    def ParseInspect(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out Inspect Results to get peptide annotations, 
        Putting them in a hash for safe keeping
        """
        inspectParser = InspectResults.Parser( FilePath )
        for result in inspectParser:
            try:
                Annotation = result.Annotation
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
                PValue = result.PValue
            except:
                traceback.print_exc()
                continue # SNAFU
            #here's some hacking that needs to be fixed.  I currently want to filter stuff by lfdr, but
            #that may not always be present
            if result.LFDR == None: #meaning that I don't have the column in question
                if PValue > self.PValueLimit: #substitute pvalue for LFDR
                    continue
            else:
                LFDR = result.LFDR
                PValue = LFDR
                if LFDR > self.PValueLimit:
                    continue
            if not self.AllPeptides.has_key(Aminos):
                self.AllPeptides[Aminos] = Peptide
            else:
                if PValue <  self.AllPeptides[Aminos]:
                    self.AllPeptides[Aminos] = Peptide
            self.SpectrumCount += 1



    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:vp:")
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
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DatabasePaths.append( Value)
            if Option == "-w":
                self.OutputPath = Value
            if Option == "-p":
                #undocumented pvalue filter
                self.PValueLimit = float(Value)
            if Option == "-v":
                self.Verbose = 1
                    
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d"):
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
    