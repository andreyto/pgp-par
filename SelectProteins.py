"""
Helper class for PValue.py and FDR.py:
Given an f-score cutoff, select a parsimonious collection of proteins which
account for most / all of the annotations.
"""
import os
import sys
import traceback
import struct
import InspectResults
from Utils import *
Initialize()

class ProteinSelector():
    def __init__(self):
        self.PeptideDict = {} # aminos -> location list
        self.ProteinPeptideCounts = {}
        self.ProteinSpectrumCounts = {}
        #self.SharedProteinPeptides = {}
        #self.SharedProteinSpectra = {}
        self.PeptideSpectrumCounts = {}
        self.ProteinPeptides = {} # Protein -> List of aminos
        self.ProteinNames = {}
        self.ProteinSequences = {}
        self.MQScoreWeight = 0.3
        self.DeltaScoreWeight = 1.5
        self.MinimumPeptideLength = 7
        self.BestScoresByPeptide = {}
        self.PValueCutoff = None
        self.MaxFileLines = None
        # if RetainRepresentativeCount is set, then we remember the
        # best n spectra for a particular annotation in the dictionary
        # self.BestRepresentatives
        self.RetainRepresentativeCount = None
        self.BestRepresentatives = {}
        self.AnnotationSpectrumCounts = {}
        self.FScoreCutoff2 = None
        self.FScoreCutoff3 = None

    def FindPeptideLocations(self, Aminos):
        PrevPos = -1
        LocationList = []
        while (1):
            Pos = self.DB.find(Aminos, PrevPos + 1)
            if Pos == -1:
                break
            # Which protein does Pos lie in?
            LowIndex = 0
            HighIndex = len(self.ProteinPos) - 1
            # Pos >= ProteinPos[LowIndex] and Pos < ProteinPos[HighIndex]
            # Special case - last protein:
            if Pos > self.ProteinPos[HighIndex]:
                ProteinID = HighIndex
                ResidueNumber = Pos - self.ProteinPos[HighIndex]
            else:
                while (1):
                    if LowIndex+1==HighIndex:
                        ProteinID = LowIndex
                        ResidueNumber = Pos - self.ProteinPos[LowIndex]
                        break
                    MidIndex = (LowIndex + HighIndex) / 2
                    if Pos >= self.ProteinPos[MidIndex]:
                        LowIndex = MidIndex
                    else:
                        HighIndex = MidIndex
            LocationList.append((ProteinID, ResidueNumber))
            PrevPos = Pos
        return LocationList
    def OldFindPeptideLocations(self, Aminos):
        LocationList = []
        #print "Find locations for %s..."%Aminos
        for (ID, Sequence) in self.ProteinSequences.items():
            Pos = Sequence.find(Aminos)
            if Pos != -1:
                LocationList.append((ID, Pos))
                #print "Found at pos %s in %s"%(Pos, ID)
        if len(LocationList) == 0:
            print "*** WARNING: Peptide '%s' not found in the database."%Aminos
        return LocationList
    def LoadDB(self, DBPath):
        DBFile = open(DBPath, "rb")
        self.DB = DBFile.read()
        DBFile.close()
        IndexPath = os.path.splitext(DBPath)[0] + ".index"
        IndexFile = open(IndexPath, "rb")
        BlockSize = struct.calcsize("<qi80s")
        ID = 0
        PrevID = None
        self.ProteinPos = []
        while (1):
            Block = IndexFile.read(BlockSize)
            if not Block:
                break
            Info = struct.unpack("<qi80s", Block)
            Name = Info[2]
            NullPos = Name.find("\0")
            if NullPos !=- 1:
                Name = Name[:NullPos]
            self.ProteinNames[ID]= Name
            StartPos = Info[1]
            self.ProteinPos.append(StartPos)
            if PrevID != None:
                self.ProteinSequences[PrevID] = self.DB[self.ProteinPos[PrevID]:StartPos - 1]
            PrevID = ID
            ID += 1
        self.ProteinSequences[PrevID] = self.DB[self.ProteinPos[PrevID]:]
    def LoadMultipleDB(self, DBPathList):
        """" Given a list of DB pathnames, load all the corresponding DB """
        ID = 0
        self.DB = ""
        self.ProteinPos = []
        for DBPath in DBPathList:
            print "loading %s"%DBPath
            DBFile = open(DBPath, "rb")
            OldDB = self.DB
            self.DB += DBFile.read()    # concatenate all DBs sequentially
            DBFile.close()
            IndexPath = os.path.splitext(DBPath)[0] + ".index"
            IndexFile = open(IndexPath, "rb")
            BlockSize = struct.calcsize("<qi80s")
            PrevID = None
            while (1):
                Block = IndexFile.read(BlockSize)
                if not Block:
                    break
                Info = struct.unpack("<qi80s", Block)
                Name = Info[2]
                NullPos = Name.find("\0")
                if NullPos !=- 1:
                    Name = Name[:NullPos]
                self.ProteinNames[ID]= Name
                StartPos = Info[1] + len(OldDB) # adjust StartPos for adding a new DB                
                self.ProteinPos.append(StartPos)
                if PrevID != None:
                    self.ProteinSequences[PrevID] = self.DB[self.ProteinPos[PrevID]:StartPos - 1]
                PrevID = ID
                ID += 1
            self.ProteinSequences[PrevID] = self.DB[self.ProteinPos[PrevID]:]
    def OldLoadDB(self, DBPath):
        """
        Load the database, popluating self.ProteinSequences
        """
        print "LoadDB(%s)"%DBPath
        IndexPath = os.path.splitext(DBPath)[0] + ".index"
        IndexFile = open(IndexPath, "rb")
        DBFile = open(DBPath, "rb")
        BlockSize = struct.calcsize("<qi80s")
        PrevName = None
        PrevID = None
        PrevStartPos = None
        ID = 0
        while (1):
            Block = IndexFile.read(BlockSize)
            if not Block:
                break
            Info = struct.unpack("<qi80s", Block)
            Name = Info[2]
            NullPos = Name.find("\0")
            if NullPos !=- 1:
                Name = Name[:NullPos]
            StartPos = Info[1]
            self.ProteinNames[ID] = Name
            if PrevName != None:
                DBFile.seek(PrevStartPos)
                Sequence = DBFile.read(StartPos - PrevStartPos)
                Sequence = Sequence.replace("*", "")
                self.ProteinSequences[PrevID] = Sequence
            PrevName = Name
            PrevID = ID
            PrevStartPos = StartPos
            ID += 1
        if PrevName != None:
            DBFile.seek(StartPos)
            Sequence = DBFile.read()
            self.ProteinSequences[PrevID] = Sequence
            #self.ProteinNames[PrevID] = Name

        DBFile.close()
        IndexFile.close()
    def ChooseProteins(self):
        """
        Iteratively select proteins which account for all the peptides.
        """
        self.SelectedProteins = {} # Protein -> (Peptides, Spectra)
        self.PeptideProteins = {} # Peptide -> final selection of protein
        print "\n\n\n"
        print "CHOOSE PROTEINS:"
        for (Peptide, SpectrumCount) in self.PeptideSpectrumCounts.items():
            for (ProteinID, Pos) in self.PeptideDict[Peptide]:
                self.ProteinSpectrumCounts[ProteinID] = self.ProteinSpectrumCounts.get(ProteinID, 0) + SpectrumCount
        while (1):
            BestCandidate = None
            BestScore = None
            for Protein in self.ProteinPeptideCounts.keys():
                if self.SelectedProteins.has_key(Protein):
                    continue
                PeptideCount = self.ProteinPeptideCounts[Protein]
                SpectrumCount = self.ProteinSpectrumCounts.get(Protein, 0)
                Score = (PeptideCount, SpectrumCount)
                #print Protein, Score
                if Score > BestScore:
                    BestScore = Score
                    BestCandidate = Protein
            if not BestScore:
                break
            (PeptideCount, SpectrumCount) = BestScore
            if PeptideCount == 0:
                break
            #%%%
            print "Accept protein %s (%s)\n  Gets %s peptides, %s spectra"%(BestCandidate, self.ProteinNames[BestCandidate], PeptideCount, SpectrumCount)
            self.SelectedProteins[BestCandidate] = BestScore
            # Lay claim to all the (not-yet-claimed) peptides:
            for Peptide in self.ProteinPeptides[BestCandidate]:
                if not self.PeptideProteins.has_key(Peptide):
                    #print "  Grab %s spectra from peptide %s"%(self.PeptideSpectrumCounts[Peptide], Peptide)
                    self.PeptideProteins[Peptide] = BestCandidate
                    # Other proteins (if not already accepted) lose a peptide, and some spectra:
                    for (OtherProtein, Pos) in self.PeptideDict[Peptide]:
                        if self.SelectedProteins.has_key(OtherProtein):
                            continue
                        self.ProteinPeptideCounts[OtherProtein] -= 1
                        self.ProteinSpectrumCounts[OtherProtein] = self.ProteinSpectrumCounts.get(OtherProtein, 0) - self.PeptideSpectrumCounts[Peptide]
        # Sanity check - the selected proteins have peptides, the unselected proteins have 0
        for Protein in self.ProteinPeptideCounts.keys():
            ProteinName = self.ProteinNames[Protein]
            PeptideCount = self.ProteinPeptideCounts[Protein]
            SpectrumCount = self.ProteinSpectrumCounts.get(Protein, 0)
            if self.SelectedProteins.has_key(Protein) and PeptideCount <= 0:
                print "** Warning: Selected protein %s (%s) has %s peptides!"%(Protein, ProteinName, PeptideCount)
            if not self.SelectedProteins.has_key(Protein) and PeptideCount != 0:
                print "** Warning: Unelected protein %s (%s) has %s peptides!"%(Protein, ProteinName, PeptideCount)
    def ParseAnnotations(self, FileName):
        """
        Parse annotations, remembering all protein locations for each peptide.
        """
        print "Parse %s..."%FileName
        File = InspectResults.Parser(FileName)
        OldSpectrum = None
        Stub = os.path.split(FileName)[1]
        LineNumber = 0
        for row in File:
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "%s %s..."%(Stub, LineNumber)
                if self.MaxFileLines != None and LineNumber >= self.MaxFileLines:
                    return # Quick-parse, for debugging only!
            try:
                Spectrum = (row.SpectrumFile, row.ScanNumber)
            except:
                continue # header line
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            try:
                MQScore = row.MQScore
                DeltaScore = row.DeltaScore
                Charge = row.Charge
            except:
                traceback.print_exc()
                print row
                continue
            # Apply a threshold: EITHER f-score cutoff (default) OR p-value cutoff
            if self.PValueCutoff != None:
                try:
                    PValue = row.PValue
                except:
                    traceback.print_exc()
                    print row
                    continue
                PeptideScore = (-PValue, MQScore)
                if PValue > self.PValueCutoff:
                    continue
            else:
                if Charge < 3:
                    WeightedScore = self.MQScoreWeight * MQScore + self.DeltaScoreWeight * (DeltaScore / self.MeanDeltaScore2)
                    if WeightedScore < self.FScoreCutoff2:
                        continue
                else:
                    WeightedScore = self.MQScoreWeight * MQScore + self.DeltaScoreWeight * (DeltaScore / self.MeanDeltaScore3)
                    if WeightedScore < self.FScoreCutoff3:
                        continue
                PeptideScore = WeightedScore
                
            try:
                Peptide = GetPeptideFromModdedName(row.Annotation)
            except:
                continue
            if len(Peptide.Aminos) < self.MinimumPeptideLength:
                continue
            # Remember this peptide:
            if not self.PeptideDict.get(Peptide.Aminos):
                #print "FindPeptideLocations(%s)"%Peptide.Aminos
                # It's a new peptide!  Figure out where it falls in the database:
                LocationList = self.FindPeptideLocations(Peptide.Aminos)
                for (Protein, Pos) in LocationList:
                    if not self.ProteinPeptides.has_key(Protein):
                        self.ProteinPeptides[Protein] = []
                    self.ProteinPeptides[Protein].append(Peptide.Aminos)
                self.PeptideDict[Peptide.Aminos] = LocationList
                for (ProteinNumber, Dummy) in LocationList:
                    self.ProteinPeptideCounts[ProteinNumber] = self.ProteinPeptideCounts.get(ProteinNumber, 0) + 1
            else:
                # We've seen this peptide before:
                LocationList = self.PeptideDict[Peptide.Aminos]
            OldScore = self.BestScoresByPeptide.get(Peptide.Aminos, -9999)
            self.BestScoresByPeptide[Peptide.Aminos] = max(PeptideScore, OldScore)
            self.PeptideSpectrumCounts[Peptide.Aminos] = self.PeptideSpectrumCounts.get(Peptide.Aminos, 0) + 1
            ##############################################################
            # Populate self.BestRepresentative, if requested:
            if self.RetainRepresentativeCount:
                Peptide.MQScore = MQScore
                Peptide.PValue = PValue
                Peptide.SpectrumFilePath = row.SpectrumFile
                Peptide.ScanNumber = row.ScanNumber
                Peptide.SpectrumFilePos = row.FileOffset
                Key = Peptide.GetFullModdedName()
                RepresentativeList = self.BestRepresentatives.get(Key, [])
                Tuple = (PeptideScore, Peptide)
                RepresentativeList.append(Tuple)
                RepresentativeList.sort()
                self.BestRepresentatives[Key] = RepresentativeList[-self.RetainRepresentativeCount:]
                self.AnnotationSpectrumCounts[Key] = self.AnnotationSpectrumCounts.get(Key, 0) + 1
                

if __name__ == "__main__":
    # Test
    Bob = ProteinSelector()
    Bob.LoadDB("database\DictyCommon.Aug28.FS2.trie")
    print Bob.FindPeptideLocations("GTVESEMAEQDSLLNKLNK")
    print Bob.FindPeptideLocations("TSEGDFTLLLGQIVDNQIGDLNKSG")
    print Bob.FindPeptideLocations("YAVFAPGLADVVIEVVAK")
