"""
Analyze peptide endpoints from a large peptide MS/MS search experiment.
Most peptides are generated from a larger protein from two cuts by trypsin, or by
one tryptic cut and the N- or C-terminus of the full protein.  Some peptide endpoints
are produced by cleavage of an N-terminal methionine, or a signal peptide.  Some
peptide matches are incorrect.  Some peptide endpoints result from non-specific chemical
cleavage (or fragmentation during electrospray); for instance, breakage before proline.
And some peptide endpoints are produced by in vivo cleavage by some other protease; these
cases are interesting.
"""
RUN_FAST = 0
RUN_FAST_LIMIT = 1000
from Utils import *
Initialize()
import os
import struct
import sys
import traceback
import cPickle
#import ShewORF

PValueCutoff = 0.01
#FScoreCutoff = 3.12
FScoreCutoff = 5.27
CHROMOSOME_SIZE = 4969804
PLASMID_SIZE = 161614

MAXIMUM_GROUP_GAP = 200

MAXIMUM_SIGNAL_LENGTH = 60

SpectrumCountCutoff = 2
CleavageRateCutoff = 0.1

class ShewGene:
    def __init__(self):
        self.Strand = 1
        self.Start = 0
        self.End = 0
        self.Length = 0
        self.AALength = 0
        self.PeptideCount = 0

class CleavageClass:
    def __init__(self):
        self.Strand = 1 # 1 forward, -1 reverse, 2 plasmid, -2 reverse-plasmid
        self.DBPos = 0
        self.ChromosomePos = 0
        # Dictionary for finding a "representative" peptide for the cleavage.
        # self.CPeptides[Aminos] = count of occurrences of that peptide;
        # self.CPeptides[None] = count of all C-terminal coverage.
        self.CPeptides = {}
        self.NPeptides = {}
        self.CleavageRate = 0
        self.CleavageCount = 0
        self.GeneName = ""
        self.AAPosition = "" # offset from the N-terminus
        self.Category = None
        
class CleavageCategorizer:
    def __init__(self):
        self.ShewGenes = {}
        self.ParseDatabase()
        self.DBLength = len(self.DB)
        self.CleavageCounts = [0]*(self.DBLength + 1)
        self.Coverage = [0]*(self.DBLength + 1) # uncleaved coverage of residues
        #self.UncleavedCoverage = [0]*(self.DBLength + 1)
        self.Cleavages = {} # Keys of the form (Strand, Position)
        self.GenesSeenAlready = {} # keys of the form (start, end)
    def LoadSEEDIDs(self):
        """
        Load a file which maps between TIGR gene IDs (SO#####) and SEED IDs.
        Sample line:
        fig|211586.1.peg.1    SO0001
        """
        File = open("Shew\\so_ids.txt", "rb")
        self.SEEDIDs = {}
        for FileLine in File.xreadlines():
            Bits = FileLine.strip().split("\t")
            if len(Bits) > 1:
                self.SEEDIDs[Bits[1]] = Bits[0]
        File.close()
        # Assign Seed IDs to our genes:
        for Gene in self.ShewGenes.values():
            Gene.SEEDID = self.SEEDIDs.get(Gene.Name, "")
    def GetSEEDURL(self, GeneName):
        SEEDID = self.SEEDIDs.get(GeneName, "")
        return "http://theseed.uchicago.edu/FIG/protein.cgi?prot=%s"%SEEDID
    def LoadShewanellaSequences(self):
        File = open("Database\\Shew.fasta", "rb")
        self.GeneSequences = {}
        CurrentSequence = ""
        for FileLine in File.xreadlines():
            if FileLine[0] == ">":
                # finish previous:
                if CurrentSequence:
                    self.GeneSequences[CurrentName] = CurrentSequence
                CurrentSequence = ""
                CurrentName = FileLine[1:].split()[0]
            else:
                CurrentSequence += FileLine.strip()
        if CurrentSequence != "":
            self.GeneSequences[CurrentName] = CurrentSequence
        File.close()
    def LoadShewanellaGenes(self, InputFileName):
        """
        Parse a file of shewanella genes.  We have one from TIGR and one from GeneMark.
        """
        # Sample line:
        # SO0001\t774\t334\t441\t146
        File = open(InputFileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if len(Bits) < 4:
                continue
            try:
                Gene = ShewGene()
                Gene.Name = Bits[0]
                Start = int(Bits[1])
                End = int(Bits[2])
            except:
                continue
            if self.GenesSeenAlready.has_key((Start, End)):
                continue
            self.GenesSeenAlready[(Start, End)] = 1
            if Start < End:
                Gene.Strand = 1
                Gene.Start = Start
                Gene.End = End
            else:
                Gene.Strand = -1
                Gene.Start = End 
                Gene.End = Start
            Gene.End += 1 # The end is NOT inclusive
            Gene.Length = Gene.End - Gene.Start
            Gene.Coverage = [0] * Gene.Length # flags for coverage by peptides
            Gene.AALength = Gene.Length / 3
            # Make note of genes from the PLASMID:
            if Gene.Name[:3] == "SOA" or Gene.Name[0] == "P":
                if Gene.Strand == 1:
                    Gene.Strand = 2
                else:
                    Gene.Strand = -2
            # Get the DATABASE BYTE POSITION for the gene!
            Start = Gene.Start
            End = Gene.End
            if Gene.Strand == 1:
                Gene.ReadingFrame = Start % 3
                if Gene.ReadingFrame == 0:
                    RecordNumber = 2
                elif Gene.ReadingFrame == 1:
                    RecordNumber = 0
                elif Gene.ReadingFrame == 2:
                    RecordNumber = 1
                Gene.DBPos = self.DBRecordPositions[RecordNumber] + (Start / 3)
                if Gene.ReadingFrame == 0:
                    Gene.DBPos -= 1
            elif Gene.Strand == -1:
                Gene.ReadingFrame = (CHROMOSOME_SIZE - End) % 3
                if Gene.ReadingFrame == 0:
                    RecordNumber = 3
                elif Gene.ReadingFrame == 1:
                    RecordNumber = 4
                elif Gene.ReadingFrame == 2:
                    RecordNumber = 5
                Gene.DBPos = self.DBRecordPositions[RecordNumber] + (CHROMOSOME_SIZE - End) / 3
            elif Gene.Strand == 2:
                Gene.ReadingFrame = Start % 3
                if Gene.ReadingFrame == 0:
                    RecordNumber = 8
                elif Gene.ReadingFrame == 1:
                    RecordNumber = 6
                elif Gene.ReadingFrame == 2:
                    RecordNumber = 7
                Gene.DBPos = self.DBRecordPositions[RecordNumber] + (Start / 3)
                if Gene.ReadingFrame == 0:
                    Gene.DBPos -= 1                                    
            elif Gene.Strand == -2:
                Gene.ReadingFrame = (PLASMID_SIZE - End) % 3
                if Gene.ReadingFrame == 0:
                    RecordNumber = 9
                elif Gene.ReadingFrame == 1:
                    RecordNumber = 10
                elif Gene.ReadingFrame == 2:
                    RecordNumber = 11
                Gene.DBPos = self.DBRecordPositions[RecordNumber] + (PLASMID_SIZE - End) / 3
##            print "Gene %s: strand %s frame %s %s-%s DBPos %s"%(Gene.Name, Gene.Strand, Gene.ReadingFrame, Gene.Start, Gene.End, Gene.DBPos)
##            print "  %s:%s"%(self.DB[Gene.DBPos - 10:Gene.DBPos], self.DB[Gene.DBPos:Gene.DBPos + 60])
            #GeneDict[Gene.Name] = Gene
            self.ShewGenes[Gene.Name] = Gene
##            Sequence = self.GeneSequences.get(Gene.Name, "")
##            Gene.TruePos = self.DB.find(Sequence[1:]) - 1
##            Str = "%s\t%s\t%s\t%s\t%s\t%s\t"%(Gene.Name, Gene.Strand, Gene.Start, Gene.End, Gene.DBPos, Gene.TruePos)
##            TruePositionsFile.write(Str + "\n")
    def ParseDatabase(self):
        """
        Read the amino acid database.  Populate self.DB, and self.DBRecordPositions
        """
        #DBFile = open("Database\\shew.trie", "rb")
        DBFile = open("Database\\ShewanellaWG.trie", "rb")
        self.DB = DBFile.read()
        DBFile.close()
        self.DBRecordPositions = [0]
        Pos = 0
        while (1):
            StarPos = self.DB.find("*", Pos + 1)
            if StarPos == -1:
                break
            self.DBRecordPositions.append(StarPos + 1)
            Pos = StarPos
        print "DBRecordPositions:", self.DBRecordPositions
    def GetChromosomePos(self, DBBytePos):
        """
        Given a byte-position in the sequence database, return a tuple of the
        form (Strand, Position).
        """
        RecordNumber = None 
        for RecordIndex in range(len(self.DBRecordPositions)):
            if DBBytePos < self.DBRecordPositions[RecordIndex]:
                RecordNumber = RecordIndex - 1
                RecordByteOffset = DBBytePos - self.DBRecordPositions[RecordIndex - 1]
                break
        if RecordNumber == None:
            # The match is in the LAST record 
            RecordNumber = len(self.DBRecordPositions) - 1 #default is Last Record
            RecordByteOffset = DBBytePos - self.DBRecordPositions[-1]
            return (0, RecordByteOffset, 0)
        if RecordNumber in (0, 1, 2):
            if RecordNumber == 0:
                Frame = 1
            elif RecordNumber == 1:
                Frame = 2
            else:
                Frame = 0
            ChromPos = Frame + (RecordByteOffset * 3)
            if Frame == 0:
                ChromPos += 3
            return (+1, ChromPos, Frame)
        elif RecordNumber in (3, 4, 5):
            Frame = RecordNumber - 3
            return (-1, CHROMOSOME_SIZE - Frame - (RecordByteOffset * 3), Frame)
        elif RecordNumber in (6, 7, 8):
            if RecordNumber == 6:
                Frame = 1
            elif RecordNumber == 7:
                Frame = 2
            else:
                Frame = 0
            ChromPos = Frame + (RecordByteOffset * 3)
            if Frame == 0:
                ChromPos += 3                                    
            return (+2, ChromPos, Frame)
        elif RecordNumber in (9, 10, 11):
            Frame = RecordNumber - 9
            return (-2, PLASMID_SIZE - Frame - (RecordByteOffset * 3), Frame) 
        else:
            Name = "Contaminant"
            return (0, RecordByteOffset, 0)
    def GetNiceChromosomeName(self, DBBytePos):
        RecordNumber = None 
        for RecordIndex in range(len(self.DBRecordPositions)):
            if DBBytePos < self.DBRecordPositions[RecordIndex]:
                RecordNumber = RecordIndex - 1
                RecordByteOffset = DBBytePos - self.DBRecordPositions[RecordIndex - 1]
                break
        if RecordNumber == None:
            # The match is in the LAST record 
            RecordNumber = len(self.DBRecordPositions) - 1 #default is Last Record
            RecordByteOffset = DBBytePos - self.DBRecordPositions[-1]
        print "DBBytePos %s -> record %s byte offset %s"%(DBBytePos, RecordNumber, RecordByteOffset)
        # Record number is 0-2 for forward strand, 3-5 for reverse,
        # 6-11 for plasmid, 12-13 for trypsin, then keratin stuff
        if RecordNumber in (0, 1, 2):
            Frame = RecordNumber
            ChromosomePos = "c+:%d"%(Frame + (RecordByteOffset * 3))
        elif RecordNumber in (3, 4, 5):
            Frame = RecordNumber - 3
            ChromosomePos = "c-:%d"%(CHROMOSOME_SIZE - Frame - (RecordByteOffset * 3))
        elif RecordNumber in (6, 7, 8):
            Frame = RecordNumber - 6
            ChromosomePos = "p+:%d"%(Frame  + (RecordByteOffset * 3))
        elif RecordNumber in (9, 10, 11):
            Frame = RecordNumber - 9
            ChromosomePos = "p-:%d"%(PLASMID_SIZE - Frame - (RecordByteOffset * 3))
        else:
            Name = "Contaminant"
            ChromosomePos = "%s.%d.%d"%(Name, RecordNumber, RecordByteOffset)
        return ChromosomePos
    def ParseTidyCleavages(self):
        self.PeptideDict = {}
        ORFFile = open("ShewORFPeptides.txt","rb")
        for FileLine in ORFFile.xreadlines():
            Bits = FileLine.split("\t")
            try:
                Peptide = Bits[0]
                Count = int(Bits[1])
            except:
                continue
            DBPos = self.PeptideDict.get(Peptide, None)
            if DBPos == None:
                DBPos = self.DB.find(Peptide)
                self.PeptideDict[Peptide] = DBPos
            if DBPos != -1:
                #print Peptide.Aminos
                self.CleavageCounts[DBPos] += Count
                CPos = DBPos + len(Peptide)
                self.CleavageCounts[CPos] += Count
                for X in range(DBPos, CPos):
                    self.Coverage[X] += Count
                    #print "Coverage %s: %s"%(X, self.Coverage[X])
                for X in range(DBPos + 1, CPos):
                    self.UncleavedCoverage[X] += Count
            else:
                print "Warning: Peptide not in database!", Peptide
        ORFFile.close()
        self.SaveCleavageEvents()
                
    def xParseCleavageEvents(self):
        self.PeptideDict = {}
        ResultsDirectory = r"F:\ftproot\Shewanella\ShewPosFixed"
        FileCount = 0
        for FileName in os.listdir(ResultsDirectory):
            print "Parse %s..."%FileName
            FilePath = os.path.join(ResultsDirectory, FileName)
            try:
                File = open(FilePath, "rb")
            except:
                traceback.print_exc()
                continue
            OldSpectrum = None
            for FileLine in File.xreadlines():
                Bits = FileLine.split("\t")
                try:
                    PValue = float(Bits[10])
                    Spectrum = (Bits[0], Bits[1])
                    FScore = float(Bits[16])
                except:
                    continue
                if Spectrum == OldSpectrum:
                    continue
                if FScore < FScoreCutoff:
                    continue
                #if PValue > PValueCutoff:
                #    continue
                OldSpectrum = Spectrum
                Peptide = GetPeptideFromModdedName(Bits[2][2:-2])
                DBPos = self.PeptideDict.get(Peptide.Aminos, None)
                if DBPos == None:
                    DBPos = self.DB.find(Peptide.Aminos)
                    self.PeptideDict[Peptide.Aminos] = DBPos
                if DBPos != -1:
                    #print Peptide.Aminos
                    self.CleavageCounts[DBPos] += 1
                    CPos = DBPos + len(Peptide.Aminos)
                    self.CleavageCounts[CPos] += 1
                    for X in range(DBPos, CPos):
                        self.Coverage[X] += 1
                        #print "Coverage %s: %s"%(X, self.Coverage[X])
                    for X in range(DBPos + 1, CPos):
                        self.UncleavedCoverage[X] += 1
                else:
                    print "Warning: Peptide not in database!", Peptide.Aminos
            File.close()
            FileCount += 1
            if FileCount > 1 and RUN_FAST:
                break # finish the parse phase quickly, for debugging
        self.SaveCleavageEvents()
    def SaveCleavageEvents(self):
        # Save all the cleavage events for further analysis:
        print "Saving cleavage events..."
        OutputFile = open("ShewCleavage.dat", "wb")
        cPickle.dump(self.Coverage, OutputFile)
        cPickle.dump(self.UncleavedCoverage, OutputFile)
        cPickle.dump(self.CleavageCounts, OutputFile)
        OutputFile.close()
    def LoadCleavageEvents(self):
        InputFile = open("ShewCleavage.dat", "rb")
        self.Coverage = cPickle.load(InputFile)
        self.UncleavedCoverage = cPickle.load(InputFile)
        self.CleavageCounts = cPickle.load(InputFile)
        InputFile.close()
    def CategorizeCleavages(self):
        """
        Iterate over all cleavage events observed in the database.  Assign each one
        to a category:
        - Tryptic terminus
        - C-terminus
        - proteolysis withut decay
        - start site / signal peptide cleavage / proteolysis without witnessed N-terminal flanking sequence
        - decay of another endpoint
        """
        CleavageReportFile = open("CleavageReport.txt", "wb")
        for Pos in range(self.DBLength):
            CleaveCount = self.CleavageCounts[Pos]
            NCoverage = self.Coverage[Pos - 1]
            CCoverage = self.Coverage[Pos]
            CleaveRate = self.CleavageCounts[Pos] / max(1.0, float(self.CleavageCounts[Pos] + self.UncleavedCoverage[Pos]))
            NextCleaveCount = self.CleavageCounts[Pos + 1]
            NextCleaveRate = self.CleavageCounts[Pos + 1] / max(1.0, float(self.CleavageCounts[Pos + 1] + self.UncleavedCoverage[Pos + 1]))
            if Pos == 0:
                PrevCleaveCount = 0
                PrevCleaveRate = 0
                PrevAA = "X"
            else:
                PrevCleaveCount = self.CleavageCounts[Pos - 1]
                PrevCleaveRate = self.CleavageCounts[Pos - 1] / max(1.0, float(self.CleavageCounts[Pos - 1] + self.UncleavedCoverage[Pos - 1]))
                PrevAA = self.DB[Pos - 1]
            if Pos < len(self.DB) - 1:
                NextAA = self.DB[Pos + 1]
            else:
                NextAA = "X"
            AA = self.DB[Pos]
            if CleaveCount < 5:
                continue
            PrefixAminos = self.DB[max(0, Pos - 10):Pos]
            SuffixAminos = self.DB[Pos:Pos + 10]
            Aminos = "%s:%s"%(PrefixAminos, SuffixAminos)
            if PrevAA in ("R", "K") and AA != "P":
                Category = "Tryptic"
            elif AA in ("X", "*"):
                Category = "StopCodon"
            elif (AA in ("U","O") and CCoverage < 5):
                Category = "StopCodon"
            elif (NextCleaveCount > CleaveCount or PrevCleaveCount > CleaveCount):
                Category = "Decay"
            elif CCoverage < NCoverage * 0.1:
                Category = "Suffix cut"
            elif NCoverage < CCoverage * 0.1:
                Category = "Prefix cut"
            else:
                Category = "Internal cut"
            PosName = self.GetNiceChromosomeName(Pos)
            Str = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(Pos, PosName, Aminos, CleaveCount, NCoverage, CCoverage, CleaveRate, Category)
            CleavageReportFile.write(Str + "\n")
    def ParseTIGRCoverage(self):
        print "Measure coverage of residues by a search of TIGR proteins..."
        if os.path.exists("TIGRCoverage.dat"):
            PickleFile = open("TIGRCoverage.dat", "rb")
            self.TIGRCoverage = cPickle.load(PickleFile)
            self.TIGRAllCoverage = cPickle.load(PickleFile)
            PickleFile.close()
            return        
        self.TIGRCoverage = [0]*(self.DBLength + 1) # *uncleaved* coverage of residues
        self.TIGRAllCoverage = [0]*(self.DBLength + 1)
        File = open("ShewTIGRPeptides001.txt", "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            Bits = FileLine.split()
            try:
                PeptideHitCount = int(Bits[1])
            except:
                continue
            LineNumber += 1
            Peptide = Bits[0]
            # Verify that the peptide is found exactly once in the database:
            DBPos = self.DB.find(Peptide)
            if DBPos == -1:
                # Try all but the first residue:
                DBPos = self.DB.find(Peptide[1:])
                if DBPos == -1:
                    print "Warning: Unknown peptide encountered '%s'"%Peptide
                    continue
                DBPos -= 1
            ########################################
            # Note uncleaved coverage of database positions:
            for X in range(DBPos + 1, DBPos + len(Peptide)):
                self.TIGRCoverage[X] += PeptideHitCount
            for X in range(DBPos, DBPos + len(Peptide)):
                self.TIGRAllCoverage[X] += PeptideHitCount
        File.close()
        PickleFile = open("TIGRCoverage.dat", "wb")
        cPickle.dump(self.TIGRCoverage, PickleFile)
        cPickle.dump(self.TIGRAllCoverage, PickleFile)
        PickleFile.close()
    def ParseCleavageEvents(self):
        """
        Parse peptides from a file.  File lines have the format "Peptide\tCount".
        First determine whether the peptide has multiple occurrences in the genome
        and if so, discared them.  Then find/create a cleavage record for the peptide's endpoints.
        """
        self.Coverage = [0]*(self.DBLength + 1) # *uncleaved* coverage of residues
        self.AllCoverage = [0]*(self.DBLength + 1)
        self.Cleavages = {}
        #File = open("ShewORFPeptides.txt", "rb")
        File = open("Shew\\PValueLong.1.8.txt", "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            Bits = FileLine.split()
            try:
                PeptideHitCount = int(Bits[1])
            except:
                continue
            LineNumber += 1
            if RUN_FAST and LineNumber > RUN_FAST_LIMIT:
                break
            if LineNumber % 1000 == 0:
                print "Line %s..."%LineNumber
            Peptide = Bits[0]
            # Verify that the peptide is found exactly once in the database:
            DBPos = self.DB.find(Peptide)
            if DBPos == -1:
                print "Warning: Unknown peptide encountered '%s'"%Peptide
                continue
            SecondPos = self.DB.find(Peptide, DBPos + 1)
            if SecondPos != -1:
                print "Skip peptide with multiple database positions: %s"%Peptide
                continue
            ########################################
            # Note uncleaved coverage of database positions:
            for X in range(DBPos + 1, DBPos + len(Peptide)):
                self.Coverage[X] += PeptideHitCount
            for X in range(DBPos, DBPos + len(Peptide)):
                self.AllCoverage[X] += PeptideHitCount
            ########################################
            # Get the strand and position of the N-terminus, record the cleavage:
            (StrandN, PosN, FrameN) = self.GetChromosomePos(DBPos)
            if StrandN:
                Key = (StrandN, PosN)
                Cleavage = self.Cleavages.get(Key, None)
                if not Cleavage:
                    Cleavage = CleavageClass()
                    Cleavage.Strand = StrandN
                    Cleavage.ChromosomePos = PosN
                    Cleavage.DBPos = DBPos
                    self.Cleavages[Key] = Cleavage
                Cleavage.CPeptides[Peptide] = Cleavage.CPeptides.get(Peptide, 0) + PeptideHitCount
                Cleavage.CPeptides[None] = Cleavage.CPeptides.get(None, 0) + PeptideHitCount
            ########################################
            # Get the strand and position of the C-terminus, record the cleavage:
            (StrandC, PosC, FrameC) = self.GetChromosomePos(DBPos + len(Peptide))
            if StrandC:
                Key = (StrandC, PosC)
                Cleavage = self.Cleavages.get(Key, None)
                if not Cleavage:
                    Cleavage = CleavageClass()
                    Cleavage.Strand = StrandC
                    Cleavage.ChromosomePos = PosC
                    Cleavage.DBPos = DBPos + len(Peptide)
                    self.Cleavages[Key] = Cleavage
                Cleavage.NPeptides[Peptide] = Cleavage.NPeptides.get(Peptide, 0) + PeptideHitCount
                Cleavage.NPeptides[None] = Cleavage.NPeptides.get(None, 0) + PeptideHitCount
                
    def FindGenesForEndpoints(self):
        """
        For each cleavage that we've observed:
        - Compute the cleavage rate
        - Find the containing gene, if any
        """
        for (Key, Cleavage) in self.Cleavages.items():
            Cleavage.CleavageCount = Cleavage.CPeptides.get(None, 0) + Cleavage.NPeptides.get(None, 0)
            UncleavedCount = self.Coverage[Cleavage.DBPos]
            Cleavage.CleavageRate = Cleavage.CleavageCount / float(Cleavage.CleavageCount + UncleavedCount)
            # If the cleavage is rare, or very rare relative to spectrum coverage, skip reporting it:
            if Cleavage.CleavageRate < CleavageRateCutoff or Cleavage.CleavageCount < SpectrumCountCutoff:
                continue
            for Gene in self.ShewGenes.values():
                GeneDBEnd = Gene.DBPos + Gene.Length / 3
                if Cleavage.DBPos >= Gene.DBPos - 10 and Cleavage.DBPos <= GeneDBEnd + 10:
                    Cleavage.GeneName = Gene.Name
                    Cleavage.AAPosition = Cleavage.DBPos - Gene.DBPos
                    #if Cleavage.Strand > 0:
                    #    Cleavage.GenePosition = Cleavage.ChromosomePos - Gene.Start
                    #else:
                    #    Cleavage.GenePosition = Gene.End - 1 - Cleavage.ChromosomePos
    def CategorizeEndpoints(self):
        for (Key, Cleavage) in self.Cleavages.items():
            # If the cleavage is rare, or very rare relative to spectrum coverage, skip reporting it:
            if Cleavage.CleavageRate < CleavageRateCutoff or Cleavage.CleavageCount < SpectrumCountCutoff:
                continue
            PrefixAA = self.DB[Cleavage.DBPos - 1]
            SuffixAA = self.DB[Cleavage.DBPos]
            NextAA = self.DB[Cleavage.DBPos + 1]
            CoverageN = self.AllCoverage[Cleavage.DBPos - 1]
            CoverageC = self.AllCoverage[Cleavage.DBPos]
            #CoverageN = Cleavage.NPeptides.get(None, 0)
            #CoverageC = Cleavage.CPeptides.get(None, 0)
            if PrefixAA in ("K", "R") and SuffixAA != "P":
                Cleavage.Category = "Trypsin"
            elif SuffixAA in ("X", "*", "O", "U"):
                Cleavage.Category = "CTerminus"
            elif Cleavage.AAPosition == 0:
                Cleavage.Category = "NTerminus"
            elif Cleavage.AAPosition == 1:
                Cleavage.Category = "NTMCleavage"
            elif Cleavage.AAPosition > 0 and Cleavage.AAPosition < MAXIMUM_SIGNAL_LENGTH and CoverageN < CoverageC * 0.1:
                Cleavage.Category = "SignalCleavage"
            else:
                Cleavage.Category = "Other"
    def ReportEndpoints(self):
        OutputFile = open("EndpointReport.txt", "wb")
        Keys = self.Cleavages.keys()
        Keys.sort()
        Header = "DBPos\tCount\tCleavageRate\tCategory\tStrand\tPosition\tGeneName\tGenePos\tCoverageN\tCoverageC\tTIGRCoverN\tTIGRCoverC\tFlankingAA\tNTerminal\tCTerminal\tGeneSpectra\tGeneCoverFrac\tPredisi\tPScore\tSignalP\tSPScore\t"
        OutputFile.write(Header + "\n")
        for Key in Keys:
            Cleavage = self.Cleavages[Key]
            # If the cleavage is rare, or very rare relative to spectrum coverage, skip reporting it:
            if Cleavage.CleavageRate < CleavageRateCutoff or Cleavage.CleavageCount < SpectrumCountCutoff:
                continue
            Str = "%s\t%s\t%s\t%s\t"%(Cleavage.DBPos, Cleavage.CleavageCount, Cleavage.CleavageRate, Cleavage.Category)
            Str += "%s\t%s\t%s\t%s\t"%(Cleavage.Strand, Cleavage.ChromosomePos, Cleavage.GeneName, Cleavage.AAPosition)
            Str += "%s\t%s\t"%(self.AllCoverage[Cleavage.DBPos - 1], self.AllCoverage[Cleavage.DBPos])
            Str += "%s\t%s\t"%(self.TIGRAllCoverage[Cleavage.DBPos - 1], self.TIGRAllCoverage[Cleavage.DBPos])
            # Get the most common n-terminal and c-terminal peptides for this cleavage:
            BestNTerminal = ""
            BestNCount = 0
            BestCTerminal = ""
            BestCCount = 0
            for (Peptide, Count) in Cleavage.NPeptides.items():
                if Peptide!=None and Count > BestNCount:
                    BestNCount = Count
                    BestNTerminal = Peptide
            for (Peptide, Count) in Cleavage.CPeptides.items():
                if Peptide!=None and Count > BestCCount:
                    BestCCount = Count
                    BestCTerminal = Peptide
            FlankingAminos = "%s:%s"%(self.DB[max(0, Cleavage.DBPos - 10):Cleavage.DBPos],
                                      self.DB[Cleavage.DBPos:Cleavage.DBPos + 10])
            Str += "%s\t%s\t%s\t"%(FlankingAminos, BestNTerminal, BestCTerminal)
            GeneCoverage = self.GeneCoverage.get(Cleavage.GeneName, (0,0))
            if Cleavage.GeneName:
                Str += "%s\t%s\t"%(GeneCoverage[0], GeneCoverage[1])
                for SignalDict in (self.SignalPredictionsA, self.SignalPredictionsB):
                    SignalPrediction = SignalDict.get(Cleavage.GeneName, None)
                    if SignalPrediction == None:
                        Str += "\t\t"
                    else:
                        Str += "%s\t%s\t"%(SignalPrediction[0],SignalPrediction[2])
            print Str
            OutputFile.write(Str + "\n")
        OutputFile.close()
    def LoadShewanellaGeneCoverage(self):
        #File = open("ShewProteinCoverage.txt", "rb")
        File = open("ShewCoverage.TIGR.txt", "rb")
        self.GeneCoverage = {}
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                PepCount = int(Bits[1])
                SpecCount = int(Bits[2])
                CoverageRate = float(Bits[3])
            except:
                continue # header
            GeneName = Bits[0]
            self.GeneCoverage[GeneName] = (SpecCount, CoverageRate)
        File.close()
    def LoadSignalPredictions(self):
        self.SignalPredictionsA = {} # Gene -> Locus, Flag, Score
        self.SignalPredictionsB = {}
        File = open("PredisiResults.txt", "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.split("\t")
            if len(Bits) < 3:
                continue
            GeneName = Bits[0].split()[0]
            if not self.ShewGenes.has_key(GeneName):
                print "* Warning: Unknown gene '%s' in line %s of PredisiResults"%(GeneName, LineNumber)
                continue
            Locus = int(Bits[2])
            Score = float(Bits[1])
            if Bits[3] == "Y":
                Flag = 1
            else:
                Flag = 0
            self.SignalPredictionsA[GeneName] = (Locus, Flag, Score)
        File.close()
        ######################################################
        File = open("SignalP_NN_CleavagePrediction.txt", "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.split("\t")
            if len(Bits)<3:
                continue
            GeneName = Bits[0].split()[0]
            if not self.ShewGenes.has_key(GeneName):
                print "* Warning: Unknown gene '%s' in line %s of PredisiResults"%(GeneName, LineNumber)
                continue
            Locus = int(Bits[2]) - 1
            Score = float(Bits[1])
            if Bits[3] == "Y":
                Flag = 1
            else:
                Flag = 0
            self.SignalPredictionsB[GeneName] = (Locus, Flag, Score)
        File.close()
    def RefuteSignalPeptides(self):
        GeneNames = self.ShewGenes.keys()
        GeneNames.sort()
        OutputFile = open("SignalPredictionsConfirmDeny.txt", "wb")
        OutputFile.write("Gene\tPredisiSite\tPredisiScore\tNCoverage\tCoverage\tTIGR-N\tTIGR-C\tCleavageRate\tSignalPSite\tSignalPScore\tNCoverage\tCoverage\tTIGR-N\tTIGR-C\tCleavageRate\t\n")
        for GeneName in GeneNames:
            Gene = self.ShewGenes[GeneName]
            Str = "%s\t"%GeneName
            for SignalDict in [self.SignalPredictionsA, self.SignalPredictionsB]:
                SignalPrediction = SignalDict.get(GeneName, None)
                if SignalPrediction == None:
                    Str += "\t\t\t\t"
                else:
                    AAPos = SignalPrediction[0]
                    DBPos = Gene.DBPos + AAPos
                    NCoverage = self.AllCoverage[DBPos - 1]
                    CCoverage = self.AllCoverage[DBPos]
                    TNCoverage = self.TIGRAllCoverage[DBPos - 1]
                    TCCoverage = self.TIGRAllCoverage[DBPos]
                    if CCoverage > 0:
                        CleavageRate = (self.AllCoverage[DBPos] - self.Coverage[DBPos]) / float(self.AllCoverage[DBPos])
                    else:
                        CleavageRate = "None(%s/%s/%s/%s)"%(self.AllCoverage[DBPos], self.AllCoverage[DBPos+1], self.Coverage[DBPos], self.Coverage[DBPos+1])
                    Str += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(AAPos, SignalPrediction[2], NCoverage, CCoverage, TNCoverage, TCCoverage, CleavageRate)
            OutputFile.write(Str + "\n")
    def ReadNovelPeptides(self, AcceptablePeptideFileName):
        "Helper function for CategorizeNovelPeptides: Parse all the peptides for 5% false discovery rate."
        File = open(AcceptablePeptideFileName, "rb")
        self.FeaturesByStrand = {1:[], -1:[], 2:[], -2:[], 0:[]}
        Strands = [1, -1, 2, -2]
        #############################################
        # Read all the peptides that were identified.  self.FeaturesByStrand[Strand#] is a list of
        # tuples of the form (StrandPosition, DBPosition, Feature), where Feature is either
        # a Peptide instance or a ShewGene instance.
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 100 == 0:
                print "Line %s..."%LineNumber
            if RUN_FAST and LineNumber > RUN_FAST_LIMIT:
                break
            Bits = FileLine.split("\t")
            try:
                SpecCount = int(Bits[1])
                ORFSize = int(Bits[2])
            except:
                continue # comment line
            Peptide = GetPeptideFromModdedName(Bits[0])
            Peptide.ORFSize = ORFSize
            DBPos = self.DB.find(Peptide.Aminos)
            if DBPos == -1:
                print "** Warning: Peptide '%s' not found in whole-genome database."%Peptide.Aminos
                continue
            DBPos2 = self.DB.find(Peptide.Aminos, DBPos + 1)
            if DBPos2 != -1:
                print "Skipping peptide '%s', found in two or more genome positions."%Peptide.Aminos
                continue
            (Strand, Pos, Frame) = self.GetChromosomePos(DBPos)
            Peptide.ReadingFrame = Frame
            Peptide.HitCount = SpecCount
            if Strand > 0:
                Peptide.Start = Pos
                Peptide.End = Pos + len(Peptide.Aminos) * 3
            else:
                Peptide.End = Pos
                Peptide.Start = Pos - (len(Peptide.Aminos) * 3)
            self.FeaturesByStrand[Strand].append((Pos, DBPos, Peptide))
        File.close()
        print "All peptides have been parsed!"        
    def DetermineNovelPeptideCoverage(self):
        """
        Helper function for CategorizeNovelPeptides.
        Determine the level of coverage of each gene (by peptides), and the level of coverage
        of each peptide (by genes):
        """
        VerboseFlag = 0
        Strands = [1, -1, 2, -2]
        for Strand in Strands:
            print "Check overlap on strand %s..."%Strand
            Features = self.FeaturesByStrand[Strand]
            for FeatureIndex in range(len(Features)):
                (StartPos, DBPos, Item) = Features[FeatureIndex]
                if isinstance(Item, ShewGene):
                    continue
                Item.Coverage = "UNCOVERED"
                Item.CoveredFlag = 0
                for NearIndex in range(FeatureIndex - 800, FeatureIndex + 800):
                    if NearIndex < 0:
                        continue
                    if NearIndex >= len(Features):
                        continue
                    if NearIndex == FeatureIndex:
                        continue
                    NearItem = Features[NearIndex][-1]
                    if not isinstance(NearItem, ShewGene):
                        continue
                    # Determine overlap:
                    if VerboseFlag:
                        print "Check for overlap between peptide %s %s-%s and gene %s %s-%s"%(Item.Aminos, Item.Start, Item.End, NearItem.Name, NearItem.Start, NearItem.End)
                    if Strand > 0:
                        GenePosStart = Item.Start - NearItem.Start
                        #GenePosEnd = GenePosStart + (len(Item.Aminos) * 3) #Item.EndPos - NearItem.Start
                    else:
                        GenePosStart = NearItem.End - Item.End  #NearItem.Start - Item.StartPos
                        #GenePosEnd = GenePosStart + (len(Item.Aminos) * 3) #NearItem.Start - Item.EndPos
                    GenePosEnd = GenePosStart + (len(Item.Aminos) * 3) #Item.EndPos - NearItem.Start
                    if VerboseFlag:
                        print "-> Gene pos %s..%s"%(GenePosStart, GenePosEnd)
                    if GenePosStart < 0 and GenePosEnd < 0:
                        continue
                    if GenePosStart >= NearItem.Length and GenePosEnd >= NearItem.Length:
                        continue
                    if VerboseFlag:
                        print "Peptide %s Gene %s GenePosStart %s GenePosEnd %s"%(Item.Aminos, NearItem.Name, GenePosStart, GenePosEnd)
                    # The peptide Item overlaps with the gene NearItem!
                    # Flag coverage in the gene:
                    for Pos in range(GenePosStart, GenePosEnd):
                        if Pos >=0 and Pos < NearItem.Length:
                            NearItem.Coverage[Pos] = 1
                    # Flag coverage of the peptide.  Note that we MIGHT have two genes which both
                    # overlap it, in which case we should take the best coverage we can get.
                    if not Item.CoveredFlag:
                        if GenePosStart < 0 or GenePosEnd < 0:
                            Item.Coverage = "OverlapStart %s"%(NearItem.Name)
                        elif GenePosStart >= NearItem.Length or GenePosEnd >= NearItem.Length:
                            Item.Coverage = "OverlapEnd %s"%(NearItem.Name)
                            #print "Peptide %s-%s OVERLAPS gene %s-%s"%(Item.StartPos, Item.EndPos, NearItem.Start, NearItem.End)
                        elif NearItem.ReadingFrame != Item.ReadingFrame:
                            Item.Coverage = "WrongFrameCoverage %s"%(NearItem.Name)
                        else:
                            Item.Coverage = "%s"%NearItem.Name
                            NearItem.PeptideCount += 1
                            # Hack: If gene name starts with "S" it's a TIGR gene and we consider this to be coverage
                            if NearItem.Name[0] == "S":
                                Item.CoveredFlag = 1
                            #print "Peptide %s-%s covered by gene %s-%s"%(Item.StartPos, Item.EndPos, NearItem.Start, NearItem.End)
                if VerboseFlag:
                    print "Peptide %s coverage: %s"%(Item.Aminos, Item.Coverage)
    def GroupNovelPeptides(self):
        """
        Helper function for CategorizeNovelPeptides.  Organizes un-covered peptides into groups.
        """
        self.AllGroups = []
        Strands = [1, -1, 2, -2]
        class PeptideGroup:
            def __init__(self):
                self.End = None
                self.PeptideCount = 0
                self.SpectrumCount = 0
                self.Strand = 0
                self.Entries = []
                self.BestORFSize = 0
            def AddFeature(self, Feature):
                self.Entries.append(Feature)
                Item = Feature[-1]
                if isinstance(Item, PeptideClass):
                    if not Item.CoveredFlag:
                        self.End = max(Item.End, self.End)
                        self.PeptideCount += 1
                        self.SpectrumCount += Item.HitCount
                        self.BestORFSize = max(self.BestORFSize, Item.ORFSize)
        for Strand in Strands:
            CurrentGroup = None
            Features = self.FeaturesByStrand[Strand]
            for FeatureIndex in range(len(Features)):
                FeatureTuple = Features[FeatureIndex]
                Feature = FeatureTuple[-1]
                Grouped = 0
                # If we are in a group, then extend it, unless we gapped:
                if CurrentGroup != None:
                    GapSize = Feature.Start - CurrentGroup.End
                    if GapSize > MAXIMUM_GROUP_GAP:
                        # Time to end the group.  Extend it up to the first gene:
                        for NearIndex in range(FeatureIndex, len(Features)):
                            NearFeature = Features[PrevIndex]
                            CurrentGroup.AddFeature(NearFeature)
                            if isinstance(NearFeature[-1], ShewGene):
                                break
                        #CurrentGroup.AddFeature(FeatureTuple)
                        #CurrentGroup.append(FeatureTuple) # flanking entry
                        CurrentGroup = None
                    else:
                        CurrentGroup.AddFeature(FeatureTuple)
                        Grouped = 1
                # Check to see whether this is an uncovered peptide, which should start a new group:
                if Grouped:
                    continue
                if isinstance(Feature, ShewGene):
                    continue
                if Feature.CoveredFlag:
                    continue
                # Aha!  This is an uncovered peptide!
                CurrentGroup = PeptideGroup()
                CurrentGroup.Strand = Strand
                # Add preceding entries, up to the first gene:
                for PrevIndex in range(FeatureIndex - 1, -1, -1):
                    PrevFeature = Features[PrevIndex]
                    CurrentGroup.AddFeature(PrevFeature)
                    if isinstance(PrevFeature[-1], ShewGene):
                        break
                CurrentGroup.Entries.reverse()
                CurrentGroup.AddFeature(FeatureTuple)
                self.AllGroups.append(CurrentGroup)
    def CategorizeNovelPeptides(self, AcceptablePeptideFileName):
        Strands = [1, -1, 2, -2]
        self.ReadNovelPeptides(AcceptablePeptideFileName)
        OutputFile = open("Shew\\Features.txt", "wb")
        # Write column headings:
        OutputFile.write("Gene\tStrand\tDBPos\tStart\tEnd\tFrame\tName\tSEEDID\tCoverage\tPeptides\t\n")
        OutputFile.write("Peptide\tStrand\tDBPos\tStart\tEnd\tFrame\tSequence\tSpectrumCount\tCoveringGene\t\n")
        # Add genes to the feature lists, and sort the feature lists:
        for Gene in self.ShewGenes.values():
            self.FeaturesByStrand[Gene.Strand].append((Gene.Start, Gene.DBPos, Gene))
        for List in self.FeaturesByStrand.values():
            List.sort()
        self.DetermineNovelPeptideCoverage()
        ############################################################################################
        # Set the coverage-fraction of all genes:
        for Strand in Strands:
            for FeatureTuple in self.FeaturesByStrand[Strand]:
                (StartPos, DBPos, Feature) = FeatureTuple
                if isinstance(Feature, ShewGene):
                    CoverCount = 0
                    for CoverageFlag in Feature.Coverage:
                        if CoverageFlag:
                            CoverCount += 1
                    Feature.CoverageRate = CoverCount / float(max(Feature.Length, 1))
        ##############################################
        # Get novel-peptide groups:
        self.GroupNovelPeptides()
        GroupFile = open("Shew\\UncoveredBlocks.txt", "wb")
        # Write column headings:
        GroupFile.write("Strand\tStart\tEnd\tBasePairs\tPeptideCount\tSpectrumCount\tORFSize\t\n")
        GroupFile.write("\tGene\tStrand\tDBPos\tStart\tEnd\tFrame\tName\tSEEDID\tCoverage\tPeptides\t\n")
        GroupFile.write("\tPeptide\tStrand\tDBPos\tStart\tEnd\tFrame\tSequence\tSpectrumCount\tCoveringGene\tORFSize\t\n")
        # Report the groups in order by number of peptides:
        SortedGroups = []
        for Group in self.AllGroups:
            SortedGroups.append((Group.PeptideCount, Group))
        SortedGroups.sort()
        SortedGroups.reverse()
        for (PeptideCount, Group) in SortedGroups:
            UncoveredStart = 999999999999999
            UncoveredEnd = 0
            Strand = Group.Strand
            for Tuple in Group.Entries:
                Feature = Tuple[-1]
                if isinstance(Feature, ShewGene):
                    continue
                if Feature.CoveredFlag:
                    continue
                UncoveredStart = min(UncoveredStart, Feature.Start)
                UncoveredEnd = max(UncoveredEnd, Feature.End)
            # First line: General group info:
            Str = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(Strand, UncoveredStart, UncoveredEnd,
                UncoveredEnd - UncoveredStart, Group.PeptideCount, Group.SpectrumCount,Group.BestORFSize)
            GroupFile.write(Str + "\n")
            # Now write one line for each feature in the group:
            for (StartPos, DBPos, Item) in Group.Entries:
                if isinstance(Item, ShewGene):
                    Str = "\tGene\t%s\t%s\t%s\t%s\t%s\t"%(Strand, DBPos, getattr(Item, "Start", ""), getattr(Item, "End", ""), Item.ReadingFrame)
                    Str += "%s\t%s\t%s\t%s\t"%(Item.Name, Item.SEEDID, Item.CoverageRate, Item.PeptideCount)
                else:
                    # It's a peptide instance:
                    Str = "\tPeptide\t%s\t%s\t%s\t%s\t%s\t"%(Strand, DBPos, Item.Start, Item.End, Item.ReadingFrame)
                    Str += "%s\t%s\t%s\t%s\t"%(Item.Aminos, Item.HitCount, Item.Coverage, Item.ORFSize)
                print Str
                GroupFile.write(Str + "\n")
        GroupFile.close()
        ############################################################################################    
        # Print all the peptides and known-genes from each strand:
        for Strand in Strands:
            Features = self.FeaturesByStrand[Strand]
            for (StartPos, DBPos, Item) in Features:
                if isinstance(Item, ShewGene):
                    Str = "Gene\t%s\t%s\t%s\t%s\t%s\t"%(Strand, DBPos, getattr(Item, "Start", ""), getattr(Item, "End", ""), Item.ReadingFrame)
                    Str += "%s\t%s\t%s\t%s\t"%(Item.Name, Item.SEEDID, Item.CoverageRate, Item.PeptideCount)
                else:
                    # It's a peptide instance:
                    Str = "Peptide\t%s\t%s\t%s\t%s\t%s\t"%(Strand, DBPos, Item.Start, Item.End, Item.ReadingFrame)
                    Str += "%s\t%s\t%s\t"%(Item.Aminos, Item.HitCount, Item.Coverage)
                print Str
                OutputFile.write(Str + "\n")
        OutputFile.close()

def CleavageReportMain():
    EndPointer = CleavageCategorizer()
    #EndPointer.ParseDatabase()
    EndPointer.LoadShewanellaGenes("ShewGenePositions.txt")
    EndPointer.LoadShewanellaGenes("GeneMarkGenePositions.txt")
    EndPointer.LoadShewanellaGeneCoverage()
    EndPointer.LoadSignalPredictions()
    EndPointer.ParseTIGRCoverage()
    print "Parse peptides..."
    EndPointer.ParseCleavageEvents()
    print "Find cleavage rates and containing genes:"
    EndPointer.FindGenesForEndpoints()
    print "Categorize endpoints:"
    EndPointer.CategorizeEndpoints()
    print "Report endpoints:"
    EndPointer.ReportEndpoints()
    EndPointer.RefuteSignalPeptides()
    #Angband.ParseCleavageEvents() # do this ONCE
    #Angband.ParseTidyCleavages()
    #Angband.LoadCleavageEvents()
    #Angband.CategorizeCleavages()

def FindNovelPeptidesMain():
    PeptideFinder = CleavageCategorizer()
    PeptideFinder.LoadShewanellaSequences()
    PeptideFinder.LoadShewanellaGenes("ShewGenePositions.txt")
    PeptideFinder.LoadShewanellaGenes("GeneMarkGenePositions.txt")
    PeptideFinder.LoadSEEDIDs()
    #PeptideFinder.ParseDatabase() # called in constructor
    #PeptideFinder.ParseTIGRCoverage()
    #PeptideFinder.CategorizeNovelPeptides("Shew\\NewPValue.05.txt")
    PeptideFinder.CategorizeNovelPeptides("Shew\\PValueLong.1.8.txt")
    #PeptideFinder.CategorizeNovelPeptides("Test.txt")
    #PeptideFinder.CategorizeNovelPeptides("Test.txt")
    
    
if __name__ =="__main__":
    import psyco
    psyco.full()
    if len(sys.argv) > 1:
        Command = sys.argv[1].lower()
    else:
        Command = "cleavage"
    if Command == "novel":
        FindNovelPeptidesMain()
    elif Command == "cleavage":
        CleavageReportMain()
    else:
        print "** UNKNOWN command %s"%Command
