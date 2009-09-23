UsageInfo = """ GelProteinAnalysis.py
This program is designed for processing the results of 2D gels, in that
it looks for the main protein found in each mzXML file, assuming that 
each mzXML file corresponds to a single protein.
Required Options
 -r [FileName] Filename or directory of annotations
 -d [TrieFile] Database for the Search
 -w [FileName] Output filename
 
 -R     Flag to remove sticky proteins from the report (default no).
"""


import os
import getopt
import sys
import ResultsParser
import SelectProteins
import ProteinStatistics
from Utils import *
Initialize()


class QuickyProteinObject:
    def __init__(self, Name, Sequence):
        self.Name = Name
        self.Sequence = Sequence
        self.Peptides = {} #key = (aminos, filename) value = spectrumcount
    def SpectraForFile(self, TargetFileName):
        """For a given filename, we count the spectra found
        """
        Count = 0
        for (Aminos, File) in self.Peptides.keys():
            if TargetFileName == File:
                Count += self.Peptides[(Aminos, File)]
        return Count

    def PrintMe(self, TargetFileName=None):
        """can weed out for a filename, if offered
        """
        print "I am protein %s"%self.Name
        for (Aminos, File) in self.Peptides.keys():
            if TargetFileName == File:
                print "%s\t%s"%(Aminos, self.Peptides[(Aminos, File)])

    def RemoveSpotFromMe(self, Spot):
        DeleteThese = []
        for Key in self.Peptides.keys():
            (Aminos, File) = Key
            if Spot == File:
                DeleteThese.append(Key)
        #now remove
        for Key in DeleteThese:
            del self.Peptides[Key]
        
    


class AbacusClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFile = None
        self.OutputFile = None
        self.Header = None
        self.UniquePeptides = {} # Peptide => (MQScore, Bits)
        self.UniqueSites = {} #site in DB => (PLS, Bits)
        self.DBPath = []  # array for possible multiple files
        self.TrypsinCount = 0
        self.KeratinCount = 0
        self.ProteinObjects = {} # key = proteinID (like returned from find peptide location), value = object
        self.Spots = {} #key = filename (spot) value = [proteinIDs found, 2, 3, ..]
        self.MostProteinsInSpot = 0
        self.RemoveSticky = 0
        self.ProteinPairs ={} #key = (ID1, ID2) value = [spot1, spot2, spot3, ...]
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        self.ProcessResultsFiles(self.InputFile, self.ParseFile)
        
        self.FlagStickyProteins()
        self.PrintDominantProteins()
        #self.FindPairs()
        #self.PrintProteinPairs()
        #self.FindMultiSpotProteins()
        
    def FindMultiSpotProteins(self):
        """For this method, we are trying to find proteins which appear in multiple spots on a gel
        THus we leverage the naming conventions for gels.  First pass is just to print them out so that I can 
        see them, and then perhaps do somehting fancier like looking for separate parts of the protein in 
        separate spots.
        """
        for (ProteinID, ProteinObject) in self.ProteinObjects.items():
            #get all my spots, not a good container for this yet, but we'll make due
            MyGels = {}
            print "\n"
            for Key in ProteinObject.Peptides.keys():
                (Aminos, File) = Key
                GelName = self.GetGelName(File)
                SpotName = self.GetSpotName(File)
                if not MyGels.has_key(GelName):
                    MyGels[GelName] = []
                SpotsInGel = MyGels[GelName]
                if not SpotName in SpotsInGel:
                    MyGels[GelName].append(SpotName)
            #now go through and find out if a gel has multiple files
            for Gel in MyGels.keys():
                if len(MyGels[Gel]) > 1:
                    self.DetermineMultiSpotOverlap(ProteinObject, Gel)
                    
                    
    def DetermineMultiSpotOverlap(self, ProteinObject, Gel):
        """Now I have a protein, and the name of a gel that has multiple spots.
        I want to see if there is good overlap in the peptides from these spots
        """
        PeptidesBySpot = {} # spot=>list of aminos
        for Key in ProteinObject.Peptides.keys():
            (Aminos, File) = Key
            IteratorGelName = self.GetGelName(File)
            if not IteratorGelName == Gel:
                continue # wrong gel
            ## file is like /export/Projects/YPestis/6framesearch/test/A412A031.rescore.txt
            Spot = self.GetSpotName(File)
            if not PeptidesBySpot.has_key(Spot):
                PeptidesBySpot[Spot] = []
            PeptidesBySpot[Spot].append(Aminos)
        ##now go through and get the bounds of these aminos on the protein
        CoverageBySpot = {}
        for (Spot, PeptideList) in PeptidesBySpot.items():
            for Peptide in PeptideList:
                Start = ProteinObject.Sequence.find(Peptide)
                Stop = Start + len(Peptide)
                if not CoverageBySpot.has_key(Spot):
                    CoverageBySpot[Spot] = (Start, Stop)
                else:
                    (CurrStart, CurrStop) = CoverageBySpot[Spot]
                    if Start < CurrStart:
                        CurrStart = Start
                    if Stop > CurrStop:
                        CurrStop = Stop
                    #now reset
                    CoverageBySpot[Spot] = (CurrStart, CurrStop)
        #now perhaps print it all out.  then look by eye?
        #print ProteinObject.Name
        #print CoverageBySpot
        #do I have spots that have zero overlap?
        AllSpots = CoverageBySpot.keys()
        for Index in range (len(AllSpots)):
            for Jndex in range( Index, len(AllSpots)):
                Spot_I = AllSpots[Index]
                Spot_J = AllSpots[Jndex]
                (Start_I, Stop_I) = CoverageBySpot[Spot_I]
                (Start_J, Stop_J) = CoverageBySpot[Spot_J]
                #mutually exclusive ?
                Exclusive = 1
                RangeJ = range (Start_J, Stop_J + 1)
                RangeI = range (Start_I, Stop_I + 1)
                for IterI in RangeI:
                    if IterI in RangeJ:
                        Exclusive = 0
                        break
                if Exclusive:
                    print ProteinObject.Name
                    print CoverageBySpot
        
    def FlagStickyProteins(self):
        """Parameters: none
        Return: None
        Description: for 2D gels, sometimes proteins are *sticky* in the column.  
        For example.  Protein A is very abundant, and present in spot 1.  Spot 1 is
        picked and sent through the MS/MS and identified.  Then spot 2 is picked.  
        Spot 2 is not located near spot 1 on the gel, but protein A shows up in
        the MS/MS analysis of spot 2.  This is most likely caused by some of protein A
        being left in the column and coming out during the second round.
        The purpose of this method is to flag, and optionally remove, these sticky
        protein identifications.
        """
        Spots = self.Spots.keys()
        Spots.sort() # the BASIC assumption here is that the ordering of filenames
        ##             follows the ordering of spots picked from the gel.
        PreviousSpotProteins = []
        SticksFound = [] # list of tuple (spot, proteinID) 
        for Spot in Spots:
            CurrentFullList = self.Spots[Spot]
            CurrentPassList = [] # has at least 2, not contaminants
            for CurrSpotProtein in CurrentFullList:
                ## first check to see if it's a contaminant.  I don't want to deal with those
                CurrName = self.ProteinPicker.ProteinNames[CurrSpotProtein]
                if self.CheckContaminants(CurrName):
                    continue
                # I also want to check that the number of spectra for this ID is better than one-hit wonder
                Spectra = self.ProteinObjects[CurrSpotProtein].SpectraForFile(Spot)
                if Spectra < 2:
                    continue
                CurrentPassList.append(CurrSpotProtein)
                if CurrSpotProtein in PreviousSpotProteins:
                    SticksFound.append((Spot, CurrSpotProtein))
                    #print "Sticky: %s in %s"%(CurrName, Spot)
            PreviousSpotProteins = CurrentPassList
        ##now here we can remove these sticky ids, if we are directed to
        if not self.RemoveSticky:
            return
        for Tuple in SticksFound:
            #get rid of both references
            (Spot, ProteinID) = Tuple
            self.Spots[Spot].remove(ProteinID)
            #now the proteinObject
            self.ProteinObjects[ProteinID].RemoveSpotFromMe(Spot)
        
    def ParseFile(self, FileName):
        """This is a wrapper for the Assess ProteinComponents
        """
        PeptidesSeen = {}
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
            SpectrumFile = Bits[self.Columns.SpectrumFile]
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            if not PeptidesSeen.has_key(Peptide.Aminos):
                PeptidesSeen[Peptide.Aminos] =0
            PeptidesSeen[Peptide.Aminos] +=1
        Handle.close()
        # We've finished our parsing, now put things in the proper structures
        self.FillObjects(FileName, PeptidesSeen)


    def FillObjects(self, FileName, PeptidesDict):
        """Given the results corresponding to a single gelspot(single mzXML file), 
        we try and find the most common protein.  We will also keep track of minor components. 
        I think that I will also keep track of hits to common contaminants separately.  perhaps
        reporting those as a QC for the experiment.
        """
        LineCount = 0
        for (Aminos, Spectra) in PeptidesDict.items():
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            #for all of the locations, put data into the quicky protein object
            for (ProteinID, StartResidue) in Locations:
                # StartResidue is ZERO based, most people think 1 based
                # just remember that for printing out and stuff
                if not self.ProteinObjects.has_key(ProteinID):
                    ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                    ProteinSequence = self.ProteinPicker.ProteinSequences[ProteinID]
                    self.ProteinObjects[ProteinID] = QuickyProteinObject(ProteinName, ProteinSequence)
                # now put in the data points
                Key = (Aminos, FileName)
                if not self.ProteinObjects[ProteinID].Peptides.has_key(Key):
                    self.ProteinObjects[ProteinID].Peptides[Key] = 0 # initialize
                self.ProteinObjects[ProteinID].Peptides[Key] += Spectra # count the spectra
                # now keep track of the other pointer, the spot to protein
                if not self.Spots.has_key(FileName):
                    self.Spots[FileName] = [] #empty list
                ProteinsInSpot = self.Spots[FileName]
                if not ProteinID in ProteinsInSpot:
                    self.Spots[FileName].append(ProteinID)

    def PrintDominantProteins(self):
        """After we've assembled our data, let's go and print off the dominant proteins
        for each spot
        """
        Handle = open(self.OutputFile, "wb")
        PrintHeader = "FileName\tTotalSpectra\tContaminating spectra\tDominantProtein\tSpectra\tNextProtein\tSpectraCount\n"
        Handle.write(PrintHeader)
        for (FileName, ProteinList) in self.Spots.items():
            #print "I'm Spot %s"%FileName
            PrintLine = "%s\t"%FileName
            #now we get the number of spectra for each protein in the spot
            #note we only are looking for real proteins, not contaminants
            SpectrumCount = {} # count
            AllSpectra = 0
            ContaminatingSpectra = 0
            for ProteinID in ProteinList: #remmeber these are IDs, not names
                ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                Object = self.ProteinObjects[ProteinID]
                #Object.PrintMe(FileName)
                SpectraForID = Object.SpectraForFile(FileName)
                AllSpectra += SpectraForID #just a global counter  
                if self.CheckContaminants(ProteinName):
                    ContaminatingSpectra += SpectraForID
                    continue
                if SpectraForID < 2: # we don't care about these proteins
                    continue
                SpectrumCount[SpectraForID] = ProteinName # the ones we actually care about
            
            PrintLine += "%s\t%s\t"%(AllSpectra, ContaminatingSpectra)
            #Now find the top 3
            CountList = SpectrumCount.keys()
            CountList.sort()
            CountList.reverse()
            WinningScores = []
            if len(CountList) >= 3:
                WinningScores.extend(CountList[:3])
            elif len(CountList) in [1,2]:
                WinningScores.extend(CountList) # just get the all
            ##I'm a huge idiot, I can't figure out how to sort this stupid thing by spectrum count
            ## because it may not be unique, I can't use it as the key, so I have to sort by values
            ## which I can't figure out how to do.  so who cares.  I'll just do it differently
            ## it's a 3N loop, which is dumb, but not a big deal         
            for Score in WinningScores:
                for (Count, ProteinName) in SpectrumCount.items():
                    if Count == Score: #print it out
                        PrintLine += "%s\t%s\t"%(ProteinName, Count)
            Handle.write("%s\n"%PrintLine)
        Handle.close()
                    
    def FindPairs(self):
        """Parameters: None
        Return: None
        Description: Now that we've identified protein components of spots, we want to find pairs of proteins
        that are often in the same spot.  We plot these out.  In the future we are also going to check
        to see if they are the same molecular weight, and thus expected to be in the same spot.
        NOTE: this is most effective when trailers have been removed (see FlagStickyProteins)
        """
        for Spot in self.Spots.keys():
            CurrentFullList = self.Spots[Spot]
            CurrentPassList = [] # has at least 2 spectra, not contaminants
            for CurrSpotProtein in CurrentFullList:
                ## first check to see if it's a contaminant.  I don't want to deal with those
                CurrName = self.ProteinPicker.ProteinNames[CurrSpotProtein]
                if self.CheckContaminants(CurrName):
                    continue
                # I also want to check that the number of spectra for this ID is better than one-hit wonder
                Spectra = self.ProteinObjects[CurrSpotProtein].SpectraForFile(Spot)
                if Spectra < 2:
                    continue
                CurrentPassList.append(CurrSpotProtein)
            #now do the double loop and match up everyone
            if len(CurrentPassList) < 2:
                continue #to next spot
            # I think that I only want to have the major component matched up with the 
            # assorted list of minor components.  then i can have some amount of information
            #encoded in the key (Major, minor)
            ## the major one shoul dbe the first one in the list (added first wit the append above)
            MajorProteinID = CurrentPassList.pop(0) # pop element 0
            for ID in CurrentPassList:
                Tuple = (MajorProteinID, ID)
                if not self.ProteinPairs.has_key(Tuple):
                    self.ProteinPairs[Tuple] = []
                self.ProteinPairs[Tuple].append(Spot)
                
                
    def PrintProteinPairs(self, MinGelOccurances = 3):
        """Parameters: Minimum number of replicate gels that this pair shows up on (as judged by the filename)
        Return: None
        Description: print off protein pairs
        """
        Handle = open ("ProteinPairsInGel.txt", "wb")
        HeaderString = "Major Protein\tMW\tMinor Protein\tMW\t Difference MW\n"
        Handle.write(HeaderString)
        for Pair in self.ProteinPairs.keys():
            (MajorID, MinorID) = Pair
            MajorSequence = self.ProteinPicker.ProteinSequences[MajorID]
            MinorSequence = self.ProteinPicker.ProteinSequences[MinorID]
            MajorName = self.ProteinPicker.ProteinNames[MajorID]
            MinorName = self.ProteinPicker.ProteinNames[MinorID]
            GelCounts = {} # filename, count
            ListOfOccurances = self.ProteinPairs[Pair]
            for Occurance in ListOfOccurances:
                GelName = self.GetGelName(Occurance)
                if not GelCounts.has_key(GelName):
                    GelCounts[GelName] = 0
                GelCounts[GelName] += 1
            # now we check for min occurances
            if not len(GelCounts) >= MinGelOccurances:
                continue #to next pair of proteins
            # now we compute the mw difference (and soon the pI difference
            MajorMW = ProteinStatistics.GetMW(MajorSequence)
            MinorMW = ProteinStatistics.GetMW(MinorSequence)
            DiffMW = MajorMW-MinorMW
            String = "%s\t%s\t%s\t%s\t%s"%(MajorName, MajorMW, MinorName, MinorMW, DiffMW)
            # now put on all the gels and crap
            for (Gel, Count) in GelCounts.items():
                String += "\t%s(%s)"%(Gel, Count)
            Handle.write("%s\n"%String)
        Handle.close()
            
    def GetGelName(self, String):
        """Parameters: string of the filename
        Return: Gel Name
        Description: a function that only works if you name gels like me
        /export/Projects/YPestis/6framesearch/test/A412A031.rescore.txt
        """
        FileName = os.path.split(String)[-1]
        #now we get the spot from the multiple extensions
        Spot = FileName.split(".")[0]  # get rid of all the multiple extensions
        # now just the gel name (lose all the numbers at the end)
        ## total hack for today.
        GelName = Spot[:-3]
        #print "%s, %s"%(FileName, GelName)
        return GelName

    def GetSpotName(self, String):
        """Parameters: string of the filename
        Return: Spot Name
        Description: a function that only works if you name gels like me
        /export/Projects/YPestis/6framesearch/test/A412A031.rescore.txt
        """
        FileName = os.path.split(String)[-1]
        #now we get the spot from the multiple extensions
        Spot = FileName.split(".")[0]  # get rid of all the multiple extensions
        return Spot

    def CheckContaminants(self, ProteinName):
        """For a protein name, check to see if its a common contaminant
        You can add things to the list here for customization
        """
        CommonContaminants = ["trypsin", "Keratin"]
        IsContaminant = 0
        for Contaminant in CommonContaminants:
            if ProteinName.find(Contaminant) > -1:
                IsContaminant = 1
                break
        return IsContaminant

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:R")
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
                self.DBPath.append(Value)
            elif Option == "-w":
                self.OutputFile = Value
            elif Option == "-R":
                self.RemoveSticky = 1
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
    
