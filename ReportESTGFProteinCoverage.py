"""
This script reads the output of AlignProteomeToGenome, as well as the database
of EST intervals and GF intervals.  It reports how much of the protein
is covered by ESTs and by GFs - how many residues, how many exons are touched
or covered entirely, how many edges are covered.  
"""
import os
import sys
import string
import traceback
import struct
import AlignProteomeToGenome
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

# Fonts are broken on Linux.  (Tried pdf, pcf, and pil formats...but no luck)
# So, we'll content ourselves with a hideous default font if we must:
try:
    TheFont = ImageFont.truetype("Times.ttf", 12)
except:
    TheFont = ImageFont.load_default()

class Colors:
    "Bag of constants specifying a color scheme"
    Background = (255, 255, 255)
    Exon = (155, 155, 155)
    SNP = (55, 55, 55)
    Mismatch = (255, 100, 100)
    ESTCorrect = (200, 0, 0)
    ESTWrong = (100, 255, 100)
    LinkCorrect = (200, 0, 0)
    LinkWrong = (125, 125, 255)
    GF = (0, 0, 100)
    Label = (0, 0, 0)
    Border = (200, 200, 200)
    SpanningEdge = (0, 155, 155)
    

class ExonClass:
    "On exon, as determined by aligning protein sequence against genomic sequence."
    def __init__(self):
        self.Start = None
        self.End = None
        self.ProteinStart = None
        self.ProteinEnd = None
        self.EdgeResidue = None
        self.Mismatches = []
        self.SNPs = [] 
    def __cmp__(self, Other):
        if not Other:
            return 1
        if self.Start < Other.Start:
            return -1
        elif self.Start > Other.Start:
            return 1
        if self.End < Other.End:
            return -1
        if self.End > Other.End:
            return 1
        return 0
    def __str__(self):
        return "<x %s-%s>"%(self.Start, self.End)
    
class GeneMapper:
    """
    Main class for GeneMapper: Reads in ESTs and SNPs and gene predictions, and reads
    a neighborhood along the genome, maps a protein to the genome, and then reports
    how much of the protein is covered by ESTs and GF and SNPs.
    """
    def __init__(self, ChromosomeNumber):
        self.ChromosomeNumber = ChromosomeNumber
        self.OutputFile = None
        self.GenerateImagesFlag = 0 # OFF by default!
    def ProcessBatch(self, Aligner, ChromosomeNumber, OutputFileName, SpecialTarget = None):
        self.OutputFile = open(OutputFileName, "w")
        self.GFLinkFile = open("TrueGFLinks.txt", "wb")
        HeaderLine = "#ProteinID\tProteinName\tChromosome\tStrand\tLength\tCoverage\tExonCount\tExons\tSNPs\tMismatches\tCoveragePercent\t\t"
        #HeaderLine = "#ProteinIndex\tProteinName\tChromosome\tStrand\tSeed\tLength\tCoverage\tExonCount\tExons\tMismatches\tSNPs\tCoveragePercent\t\t"
        HeaderLine += "ESTCoverage\tESTExonsCovered\tESTExonsTouched\tESTJunctions\tESTGaps\tESTFirst\tESTLast\t"
        HeaderLine += "GFCoverage\tGFExonsCovered\tGFExonsTouched\tGFJunctions\tGFGaps\tGFFirst\tGFLast\t"
        HeaderLine += "EitherCoverage\tEitherExonsCovered\tEitherExonsTouched\tEitherJunctions\tEitherGaps\tEitherFirst\tEitherLast\t"
        self.OutputFile.write(HeaderLine+"\n")
        Keys = Aligner.Proteins.keys()
        Keys.sort()
        if SpecialTarget:
            print "SPECIAL TARGET IS:", SpecialTarget
        for Key in Keys:
            Protein = Aligner.Proteins[Key]
            if SpecialTarget != None and Protein.ProteinID != SpecialTarget:
                continue
            # Skip proteins with no alignment:
            print "Protein %s has %s seeds"%(Key, len(Protein.Seeds))
            if not Protein.Seeds:
                continue
            # Skip proteins aligned to another chromosome:
            Seed = Protein.Seeds[0]
            print "Seed is on chromosome %s"%Seed.Chromosome
            if Seed.Chromosome != ChromosomeNumber:
                continue
            self.TrueStrand = Seed.Strand
            if self.TrueStrand == 0:
                self.TrueStrand = -1
            ProteinLength = Protein.Length
            ProteinCoverage = Seed.ResultCoverage
            ExonString = AlignProteomeToGenome.UnparseExonString(Seed.ResultExons)
            self.ParseExons(ExonString)
            Mismatches = Seed.ResultMismatches
            if Mismatches == None:
                Mismatches = []
            SNPs = Seed.ResultSNPs
            if SNPs == None:
                SNPs = []
            ProteinName = Protein.Name
            ProteinID = Key
            if SpecialTarget != None and ProteinID != SpecialTarget:
                continue
            # %%%
            if ProteinCoverage < ProteinLength:
                continue # Skip all records not PERFECTLY COVERED.
            print Mismatches
            self.MeasureCoverage(ProteinLength, self.TrueStrand, Mismatches, 1)
            self.MeasureCoverage(ProteinLength, self.TrueStrand, Mismatches, 0)
            self.MeasureCoverage(ProteinLength, self.TrueStrand, Mismatches, 2)
            self.WriteGFLinkages(ProteinLength, self.TrueStrand, Mismatches) #%%%
            # Construct a list of fields, for output:
            #HeaderLine = "#ProteinID\tProteinName\tChromosome\tStrand\tLength\tCoverage\tExonCount\t
            #Exons\tSNPs\tMismatches\tCoveragePercent\t\t"
            Bits = [Key, Protein.Name, Seed.Chromosome, Seed.Strand, Protein.Length, Seed.ResultCoverage,
                    len(Seed.ResultExons), ExonString, SNPs, Mismatches, Seed.ResultCoveragePercent]
            Bits.append("")
            Bits.append("%s"%self.ESTCoverage)
            Bits.append("%s"%self.ESTExonsCovered)
            Bits.append("%s"%self.ESTExonsTouched)
            Bits.append("%s"%self.ESTEdgesCovered)
            Bits.append("%s"%self.ESTEdgesGapped)
            Bits.append("%s"%self.ESTFirstFlag)
            Bits.append("%s"%self.ESTLastFlag)
            Bits.append("%s"%self.GFCoverage)
            Bits.append("%s"%self.GFExonsCovered)
            Bits.append("%s"%self.GFExonsTouched)
            Bits.append("%s"%self.GFEdgesCovered)
            Bits.append("%s"%self.GFEdgesGapped)
            Bits.append("%s"%self.GFFirstFlag)
            Bits.append("%s"%self.GFLastFlag)
            Bits.append("%s"%self.EitherCoverage)
            Bits.append("%s"%self.EitherExonsCovered)
            Bits.append("%s"%self.EitherExonsTouched)
            Bits.append("%s"%self.EitherEdgesCovered)
            Bits.append("%s"%self.EitherEdgesGapped)
            Bits.append("%s"%(self.ESTFirstFlag or self.GFFirstFlag))
            Bits.append("%s"%(self.ESTLastFlag or self.GFLastFlag))
            Bits = map(str, Bits) # some entries might be ints, string.join doesn't like that
            # Write to output file:
            Str = string.join(Bits, "\t")
            self.OutputFile.write(Str)
            self.OutputFile.write("\n")
            if not self.GenerateImagesFlag:
                print Str
                continue
            print "\n"*10
            print Str
            self.AllMismatches = Mismatches
            # Note mismatches and snps in the "owning" exons:
            for Exon in self.ExonList:
                for Tuple in Mismatches:
                    (ProteinPos, Dummy) = Tuple
                    if ProteinPos >= Exon.Start and ProteinPos < Exon.End:
                        Exon.Mismatches.append(Tuple)
                for Tuple in SNPs:
                    (ProteinPos, Dummy) = Tuple
                    if ProteinPos >= Exon.Start and ProteinPos < Exon.End:
                        Exon.SNPs.append(Tuple)
            # Count total exon coverage:
            self.TotalExonCoverage = 0
            for Exon in self.ExonList:
                self.TotalExonCoverage += (Exon.ProteinEnd - Exon.ProteinStart)
                if Exon.EdgeResidue:
                    self.TotalExonCoverage += 1
            self.TotalExonCoverage -= len(self.AllMismatches)
            # Load the protein sequence:
            self.LoadProteinSequence(ProteinID)
            File = open("GeneMap\\%s.txt"%ProteinID, "w")
            File.write(self.ProteinSequence)
            File.close()
            self.GenerateImage(self.ExonList, ProteinID)
    def LoadProteinSequence(self, ProteinID):
        IndexFilePath = "Database\\IPIv314.index"
        DBFilePath = "Database\\IPIv314.trie"
        IndexFile = open(IndexFilePath, "rb")
        DBFile = open(DBFilePath, "rb")
        IndexFile.seek(92 * ProteinID + 8)
        DBFilePos = struct.unpack("<i", IndexFile.read(4))[0]
        DBFile.seek(DBFilePos)
        self.ProteinSequence = ""
        while (1):
            Stuff = DBFile.read(1024)
            if not Stuff:
                break
            StarPos = Stuff.find("*")
            if StarPos == -1:
                self.ProteinSequence += Stuff
                continue
            self.ProteinSequence += Stuff[:StarPos]
            break
        IndexFile.close()
        DBFile.close()
    def ParseExons(self, ExonString):
        "Parse self.ExonList and self.JunctionList from a string."
        # 114963200-114963317-24-30-31,114964231-114964391-32-55,
        self.ExonList = []
        self.JunctionList = []
        if ExonString[0] == '"':
            ExonString = ExonString[1:-1]
        for ExonBit in ExonString.split(","):
            DashBits = ExonBit.split("-")
            Exon = ExonClass()
            self.ExonList.append(Exon)
            Exon.Start = int(DashBits[0])
            Exon.End = int(DashBits[1])
            Exon.ProteinStart = int(DashBits[2])
            Exon.ProteinEnd = int(DashBits[3])
            if len(DashBits) > 4:
                Exon.EdgeResidue = int(DashBits[4])
        for Index in range(len(self.ExonList) - 1):
            self.JunctionList.append((self.ExonList[Index].End, self.ExonList[Index + 1].Start))
    def WriteGFLinkages(self, ProteinLength, Strand, Mismatches):
        if Strand == 1:
            IntervalDict = self.GFIntervalDictForward
        else:
            IntervalDict = self.GFIntervalDictReverse
        for ExonIndex in range(len(self.ExonList) - 1):
            Exon = self.ExonList[ExonIndex]
            NextExon = self.ExonList[ExonIndex + 1]
            ThisInterval = IntervalDict.get((Exon.Start, Exon.End), None)
            if ThisInterval == None:
                continue
            NextInterval = IntervalDict.get((NextExon.Start, NextExon.End), None)
            if NextInterval == None:
                continue
            IntronLength = NextExon.Start - Exon.End
            Str = "%s\t%s\t%s\t%s\t"%(ThisInterval + NextInterval, 
                                      min(ThisInterval, NextInterval),
                                      max(ThisInterval, NextInterval),
                                      IntronLength)
            self.GFLinkFile.write(Str + "\n")
    def MeasureCoverage(self, ProteinLength, Strand, Mismatches, ESTFlag):
        if ESTFlag == 2:
            DataSourceCount = 2
            if Strand == 1:
                IntervalDicts = [self.ESTIntervalDictForward, self.GFIntervalDictForward]
                IntervalLists = [self.ESTIntervalsForward, self.GFIntervalsForward]
                JunctionDicts = [self.ESTJunctionDictForward, self.GFJunctionDictForward]
                JunctionLists = [self.ESTJunctionsForward, self.GFJunctionsForward]
            else:
                IntervalDicts = [self.ESTIntervalDictReverse, self.GFIntervalDictReverse]
                IntervalLists = [self.ESTIntervalsReverse, self.GFIntervalsReverse]
                JunctionDicts = [self.ESTJunctionDictReverse, self.GFJunctionDictReverse]
                JunctionLists = [self.ESTJunctionsReverse, self.GFJunctionsReverse]
        elif ESTFlag == 1:
            DataSourceCount = 1
            self.ESTFirstFlag = 0
            self.ESTLastFlag = 0
            if Strand == 1:
                IntervalDicts = [self.ESTIntervalDictForward]
                IntervalLists = [self.ESTIntervalsForward]
                JunctionDicts = [self.ESTJunctionDictForward]
                JunctionLists = [self.ESTJunctionsForward]
            else:
                IntervalDicts = [self.ESTIntervalDictReverse]
                IntervalLists = [self.ESTIntervalsReverse]
                JunctionDicts = [self.ESTJunctionDictReverse]
                JunctionLists = [self.ESTJunctionsReverse]      
        else:
            DataSourceCount = 1
            self.GFFirstFlag = 0
            self.GFLastFlag = 0
            if Strand == 1:
                IntervalDicts = [self.GFIntervalDictForward]
                IntervalLists = [self.GFIntervalsForward]
                JunctionDicts = [self.GFJunctionDictForward]
                JunctionLists = [self.GFJunctionsForward]
            else:
                IntervalDicts = [self.GFIntervalDictReverse]
                IntervalLists = [self.GFIntervalsReverse]
                JunctionDicts = [self.GFJunctionDictReverse]
                JunctionLists = [self.GFJunctionsReverse]
        ExonsTouched = 0
        ExonsCovered = 0
        EdgesCovered = 0
        EdgesGapped = 0 # splice junctions s.t. the flanking exons are both hit but the junction's missed
        CoverageFlags = [0] * ProteinLength
        #print ProteinLength, CoverageFlags
        ExonHits = []
        EdgeHits = []
        for ExonIndex in range(len(self.ExonList)):
            CoverFlag = 0
            TouchFlag = 0
            BestOverlap = 0
            JunctionFlag = 0
            BestOverlapInterval = None
            for DataSourceType in range(DataSourceCount):
                IntervalDict = IntervalDicts[DataSourceType]
                IntervalList = IntervalLists[DataSourceType]
                JunctionDict = JunctionDicts[DataSourceType]
                JunctionList = JunctionLists[DataSourceType]
                Exon = self.ExonList[ExonIndex]
                if IntervalDict.has_key((Exon.Start, Exon.End)):
                    print "Interval (%s-%s) found: Exon %s ok!"%(Exon.Start, Exon.End, ExonIndex)
                    CoverFlag = 1
                    TouchFlag = 1
                    ##for X in range(Exon.ProteinStart, Exon.ProteinEnd):
                    ##    CoverageFlags[X] = 1
                    ##ExonsCovered += 1
                    ##ExonsTouched += 1
                else:
                    # The exon isn't exactly present.  find the best overlap:
                    for (Start, End) in IntervalList:
                        if Start >= Exon.End:
                            break
                        if End < Exon.Start:
                            continue
                        Overlap = min(Exon.End, End) - max(Exon.Start, Start)
                        if Overlap > BestOverlap:
                            BestOverlap = Overlap
                            BestOverlapInterval = (Start, End)
                if ExonIndex == len(self.ExonList) - 1:
                    continue
                NextExon = self.ExonList[ExonIndex + 1]
                # Check for the junction:
                Start = Exon.End
                End = NextExon.Start
                if JunctionDict.has_key((Start, End)):
                    print "[%s] Link from %s to %s found - junction %s ok!"%(ESTFlag, Start, End, ExonIndex)
                    #EdgesCovered += 1
                    JunctionFlag = 1
                else:
                    BestDistance = 9999
                    for (JStart, JEnd) in JunctionList:
                        if (JStart > Start + 20):
                            break
                        if (JStart < Start - 20):
                            continue
                        Dist = abs(Start - JStart) + abs(End - JEnd)
                        if (Dist < BestDistance):
                            BestDistance = Dist
                            BestEdge = (JStart, JEnd)
                    if BestDistance < 9999:
                        Str = "Junction %s-%s not found (closest was %s-%s) "%(Start, End, BestEdge[0], BestEdge[1])
                    else:
                        Str = "Junction %s-%s not found (nothing near) "%(Start, End)
                    print "**", Str
            if JunctionFlag:
                EdgeHits.append(1)
                EdgesCovered += 1
                if (Exon.EdgeResidue):
                    print "Junction residue %s ok!"%(Exon.EdgeResidue)
                    CoverageFlags[Exon.EdgeResidue] = 1
            else:
                EdgeHits.append(0)
            if CoverFlag:
                ExonsCovered += 1
                ExonsTouched += 1
                for X in range(Exon.ProteinStart, Exon.ProteinEnd):
                    CoverageFlags[X] = 1
                ExonHits.append(1)
                ####################################################################
                # Track whether we saw the first / last exon in the gene:
                if ExonIndex == 0:
                    if Strand == 1:
                        if ESTFlag == 1:
                            self.ESTFirstFlag = 1
                        elif ESTFlag == 0:
                            self.GFFirstFlag = 1
                    else:
                        if ESTFlag == 1:
                            self.ESTLastFlag = 1
                        elif ESTFlag == 0:
                            self.GFLastFlag = 1
                if ExonIndex == len(self.ExonList) - 1:
                    if Strand == 1:
                        if ESTFlag == 1:
                            self.ESTLastFlag = 1
                        elif ESTFlag == 0:
                            self.GFLastFlag = 1
                    else:
                        if ESTFlag == 1:
                            self.ESTFirstFlag = 1
                        elif ESTFlag == 0:
                            self.GFFirstFlag = 1
            elif BestOverlap == 0:
                print "[%s] Exon %s-%s missing! "%(ESTFlag, Exon.Start, Exon.End)
                ExonHits.append(0)
            elif BestOverlap < (Exon.End - Exon.Start):
                ExonHits.append(0)
                print "[%s] Exon %s-%s (%s-%s) truncated (best overlap %s-%s) "%(ESTFlag, Exon.Start, Exon.End,
                    Exon.ProteinStart, Exon.ProteinEnd, BestOverlapInterval[0], BestOverlapInterval[1])
                ExonsTouched += 1
                if Strand == 1:
                    FirstIndex = Exon.ProteinStart + (max(BestOverlapInterval[0], Exon.Start) - Exon.Start) / 3
                    LastIndex = Exon.ProteinStart + (min(BestOverlapInterval[1], Exon.End) - Exon.Start) / 3
                else:
                    FirstIndex = Exon.ProteinStart + (Exon.End - min(Exon.End, BestOverlapInterval[1])) / 3
                    LastIndex = Exon.ProteinStart + (Exon.End - max(Exon.Start, BestOverlapInterval[0])) / 3
                FirstIndex = max(FirstIndex, Exon.ProteinStart)
                LastIndex = min(LastIndex, Exon.ProteinEnd)
                CoverRange = (FirstIndex, LastIndex)                    
                print "Overlap cover-range is:",CoverRange
                for X in range(CoverRange[0], CoverRange[1]):
                    if X>=0 and X<len(CoverageFlags):
                        CoverageFlags[X] = 1
            else:
                ExonHits.append(0)
                print "[%s] Exon #%d (%s-%s) contained in interval %s-%s"%(ESTFlag, ExonIndex, Exon.Start, Exon.End,
                    BestOverlapInterval[0], BestOverlapInterval[1])
                ExonsCovered += 1
                ExonsTouched += 1 # full counts as touched, too
                for X in range(Exon.ProteinStart, Exon.ProteinEnd):
                    if X>=0 and X<len(CoverageFlags):
                        CoverageFlags[X] = 1
                ####################################################################
                # (copypasta) Track whether we saw the first / last exon in the gene:
                if ExonIndex == 0:
                    if Strand == 1:
                        if ESTFlag == 1:
                            self.ESTFirstFlag = 1
                        elif ESTFlag == 0:
                            self.GFFirstFlag = 1
                    else:
                        if ESTFlag == 1:
                            self.ESTLastFlag = 1
                        elif ESTFlag == 0:
                            self.GFLastFlag = 1
                if ExonIndex == len(self.ExonList) - 1:
                    if Strand == 1:
                        if ESTFlag == 1:
                            self.ESTLastFlag = 1
                        elif ESTFlag == 0:
                            self.GFLastFlag = 1
                    else:
                        if ESTFlag == 1:
                            self.ESTFirstFlag = 1
                        elif ESTFlag == 0:
                            self.GFFirstFlag = 1
                            
        for Index in range(len(EdgeHits) - 1):
            if ExonHits[Index] and ExonHits[Index+1] and not EdgeHits[Index]:
                EdgesGapped += 1
        # End of exon loop.  Un-cover the mismatch residues:
        for (ProteinPos, Dummy) in Mismatches:
            try:
                CoverageFlags[ProteinPos] = 0
            except:
                continue # out of range - ignore
        ResiduesCovered = 0
        for Flag in CoverageFlags:
            ResiduesCovered += Flag
        # Save results:
        if ESTFlag == 2:
            self.EitherCoverage = ResiduesCovered
            self.EitherExonsCovered = ExonsCovered
            self.EitherExonsTouched = ExonsTouched
            self.EitherEdgesCovered = EdgesCovered
            self.EitherEdgesGapped = EdgesGapped
        elif ESTFlag == 1:
            self.ESTCoverage = ResiduesCovered
            self.ESTExonsCovered = ExonsCovered
            self.ESTExonsTouched = ExonsTouched
            self.ESTEdgesCovered = EdgesCovered
            self.ESTEdgesGapped = EdgesGapped
        else:
            self.GFCoverage = ResiduesCovered
            self.GFExonsCovered = ExonsCovered
            self.GFExonsTouched = ExonsTouched
            self.GFEdgesCovered = EdgesCovered
            self.GFEdgesGapped = EdgesGapped
    def ReadESTs(self):
        "Read EST intervals for the chromosome, so we can know whether there's EST evidence for exons."
        self.ESTIntervalDictForward = {}
        self.ESTIntervalsForward = []
        self.ESTJunctionDictForward = {}
        self.ESTJunctionsForward = []
        self.ESTIntervalDictReverse = {}
        self.ESTIntervalsReverse = []
        self.ESTJunctionDictReverse = {}
        self.ESTJunctionsReverse = []
        # Read ESTs for both strands, ignoring orientation for now
        for (Strand, ReverseChar) in [(1, "+"), (-1, "-")]:
            #FileName = "EST\\%s%s.filtered"%(self.ChromosomeNumber, ReverseChar)
            FileName = "ESTREF\\%s%s.filtered"%(self.ChromosomeNumber, ReverseChar) #%%% TEMP: use ESTREF
            print "Read ESTs from %s..."%FileName
            File = open(FileName, "rb")
            while (1):
                Data = File.read(12)
                if not Data:
                    break
                (Start, End, Count) = struct.unpack("<iii", Data)
                Interval = (Start, End)
                if Strand == 1:
                    self.ESTIntervalDictForward[Interval] = 1
                    self.ESTIntervalsForward.append(Interval)
                else:
                    self.ESTIntervalDictReverse[Interval] = 1
                    self.ESTIntervalsReverse.append(Interval)
                JunctionCount = struct.unpack("<i", File.read(4))[0]
                for X in range(JunctionCount):
                    Data = File.read(12)
                    (JunctionStart, Dummy, Dummy2) = struct.unpack("<iif", Data)
                    Key = (JunctionStart, Start)
                    if Strand == 1:
                        if not self.ESTJunctionDictForward.has_key(Key):
                            self.ESTJunctionDictForward[Key] = 1 # don't store count, we weight it equally with GFJunctionDict
                            self.ESTJunctionsForward.append(Key)
                    else:
                        if not self.ESTJunctionDictReverse.has_key(Key):
                            self.ESTJunctionDictReverse[Key] = 1 
                            self.ESTJunctionsReverse.append(Key)
                            
            File.close()
        print "Sort EST intervals..."
        self.ESTIntervalsForward.sort()
        self.ESTIntervalsReverse.sort()
        self.ESTJunctionsForward.sort()
        self.ESTJunctionsReverse.sort()
    def ReadGeneFinding(self):
        """
        Read gene-finder output, to get a lot of putative intervals.
        """
        self.GFIntervalDictForward = {}
        self.GFIntervalsForward = []
        self.GFJunctionDictForward = {}
        self.GFJunctionsForward = []
        self.GFIntervalDictReverse = {}
        self.GFIntervalsReverse = []
        self.GFJunctionDictReverse = {}
        self.GFJunctionsReverse = []
        # Read intervals for both strands, ignoring orientation for now
        for (Strand, ReverseChar) in [(1, "+"), (-1, "-")]:
            # %%%% TEMP, to use new geneid linking procedure
            #FileName = "GeneFindDB\\%s%s.dat"%(self.ChromosomeNumber, ReverseChar)
            FileName = "NGeneFindDB\\%s%s.dat"%(self.ChromosomeNumber, ReverseChar)
            print "Read GF intervals from %s..."%FileName
            File = open(FileName, "rb")
            while (1):
                Data = File.read(16)
                if not Data:
                    break
                (Start, End, Flags, Score) = struct.unpack("<iiii", Data)
                Interval = (Start, End)
                if Strand == 1:
                    self.GFIntervalDictForward[Interval] = Score
                    self.GFIntervalsForward.append(Interval)
                else:
                    self.GFIntervalDictReverse[Interval] = Score
                    self.GFIntervalsReverse.append(Interval)
                JunctionCount = struct.unpack("<i", File.read(4))[0]
                for X in range(JunctionCount):
                    Data = File.read(8)
                    (JunctionStart, Dummy) = struct.unpack("<if", Data)
                    Key = (JunctionStart, Start)
                    if Strand == 1:
                        if not self.GFJunctionDictForward.has_key(Key):
                            self.GFJunctionDictForward[Key] = 1
                            self.GFJunctionsForward.append(Key)
                    else:
                        if not self.GFJunctionDictReverse.has_key(Key):
                            self.GFJunctionDictReverse[Key] = 1
                            self.GFJunctionsReverse.append(Key)
                            
        print "Sort GF intervals..."
        self.GFIntervalsForward.sort()
        self.GFIntervalsReverse.sort()
        self.GFJunctionsForward.sort()
        self.GFJunctionsReverse.sort()
    def GenerateImage(self, ExonList, ProteinID):
        """
        Let's generate a .png image, diagramming our exons (with mismatches)
        together with the ESTS and gene-predictions nearby.  Our image will
        have one column for each nucleotide, and we trim introns.  Even so,
        it could be a rather wide image.
        """
        # First, let's decide how to group the exons, so that we know how
        # wide the image will be:
        ExonGroups = []
        NextExon = 0
        CurrentGroupEnd = None
        while (1):
            if NextExon >= len(ExonList):
                break
            Exon = ExonList[NextExon]
            if CurrentGroupEnd != None and (Exon.Start - CurrentGroupEnd) < 150:
                ExonGroups[-1].append(Exon)
            else:
                ExonGroups.append([Exon])
            CurrentGroupEnd = Exon.End
            NextExon += 1
        Width = 0
        for ExonGroup in ExonGroups:
            print "Draw exon group with %s members."%len(ExonGroup)
            for Exon in ExonGroup:
                print "  %s-%s"%(Exon.Start, Exon.End)
            
            GroupWidth = ExonGroup[-1].End - ExonGroup[0].Start + 100 # padding on both sides
            Width += GroupWidth + 1
        Height = 250
        self.GenomeY = 125
        self.Image = Image.new("RGB", (Width, Height), Colors.Background)
        print "->-> Created image %s x %s"%(Width, Height)
        self.Draw = ImageDraw.Draw(self.Image)
        # Caption
        if self.TrueStrand == 1:
            StrandName = "forward"
        else:
            StrandName = "reverse"
        Str = "Protein %s maps to %s exons on chromosome %s (%s strand)"%(ProteinID, len(ExonList), self.ChromosomeNumber, StrandName)
        self.Draw.text((5, 2), Str, Colors.Label)
        Str = "Exons %s EST %s GF %s"%(self.TotalExonCoverage, self.ESTCoverage, self.GFCoverage)
        if len(self.AllMismatches):
            Str += " %s mismatches"%(len(self.AllMismatches))
        self.Draw.text((5, 10), Str, Colors.Label)
        self.GroupLeft = 0
        for GroupIndex in range(len(ExonGroups)):
            ExonGroup = ExonGroups[GroupIndex]
            print "Draw exon group with %s members."%len(ExonGroup)
            for Exon in ExonGroup:
                print "  %s-%s"%(Exon.Start, Exon.End)
            self.GroupGenomeStart = ExonGroup[0].Start - 50
            self.GroupGenomeEnd = ExonGroup[-1].End + 50
            self.GroupWidth = self.GroupGenomeEnd - self.GroupGenomeStart
            Right = self.GroupLeft + self.GroupWidth
            ########################################################
            self.GroupRight = self.GroupLeft + self.GroupWidth
            self.Draw.line((self.GroupLeft + 2, self.GenomeY, self.GroupRight - 2, self.GenomeY), Colors.Label)
            # Sigil: ----//----
            self.Draw.line((self.GroupLeft, self.GenomeY + 2, self.GroupLeft + 4, self.GenomeY - 2), Colors.Label)
            self.Draw.line((self.GroupRight - 4, self.GenomeY + 2, self.GroupRight, self.GenomeY - 2), Colors.Label)
            # Draw genome, and label the genome positions
            LabelGenomePos = ((self.GroupGenomeStart + 51) / 100) * 100
            LabelY = 20 + (GroupIndex % 2)*10
            while LabelGenomePos < self.GroupGenomeEnd:
                print LabelGenomePos, self.GroupGenomeStart, self.GroupGenomeEnd
                X = self.GetX(LabelGenomePos)
                # Don't draw the label if we're close to the edge of the exon group.
                # (We don't want the label to spill over *too* far into the adjoining group)
                if (X > self.GroupLeft + 25 and X < self.GroupRight - 25):
                    self.Draw.line((X, self.GenomeY, X, LabelY + 10), Colors.Border)
                    Str = str(LabelGenomePos)
                    self.Draw.text((X - 3*len(Str), LabelY), Str, Colors.Label)
                LabelGenomePos += 100
            ########################################################
            # Draw the exons, and mismatches:
            for Exon in ExonGroup:
                LeftX = self.GetX(Exon.Start)
                RightX = self.GetX(Exon.End)
                Exon.LeftX = LeftX
                Exon.RightX = RightX
                for Y in range(self.GenomeY - 4, self.GenomeY + 5):
                    self.Draw.line((LeftX, Y, RightX, Y), Colors.Exon)
                for Mismatch in Exon.Mismatches:
                    X = self.GetX(Mismatch[1])
                    self.Draw.line((X, self.GenomeY - 4, X, self.GenomeY + 5), Colors.Mismatch)
                    self.Draw.line((X-1, self.GenomeY - 4, X-1, self.GenomeY + 5), Colors.Mismatch)
                    self.Draw.line((X+1, self.GenomeY - 4, X+1, self.GenomeY + 5), Colors.Mismatch)
                for SNP in Exon.SNPs:
                    X = self.GetX(SNP[1])
                    self.Draw.line((X, self.GenomeY - 4, X, self.GenomeY + 5), Colors.SNP)
                    self.Draw.line((X-1, self.GenomeY - 4, X-1, self.GenomeY + 5), Colors.SNP)
                    self.Draw.line((X+1, self.GenomeY - 4, X+1, self.GenomeY + 5), Colors.SNP)
                    
                Str = "%s-%s"%(Exon.ProteinStart, Exon.ProteinEnd)
                self.Draw.text((LeftX, self.GenomeY + 8), Str, Colors.Label)
                Sequence = self.ProteinSequence[Exon.ProteinStart:Exon.ProteinEnd]
                MaxPrintLen = max(8, (RightX - LeftX) / 7)
                if len(Sequence) < MaxPrintLen:
                    Str = Sequence
                else:
                    Half = (MaxPrintLen / 2) - 1
                    Str = "%s..%s"%(Sequence[:Half], Sequence[-Half:])
                self.Draw.text((LeftX, self.GenomeY + 16), Str, Colors.Label)
            # Draw ESTs and GFs:
            for ExonIndex in range(len(ExonList)):
                Exon = ExonList[ExonIndex]
                if Exon not in ExonGroup:
                    continue
                if ExonIndex < len(ExonList) - 1:
                    NextExon = ExonList[ExonIndex + 1]
                    if NextExon in ExonGroup:
                        NextExonInGroup = 1
                    else:
                        NextExonInGroup = 0
                else:
                    NextExon = None
                    NextExonInGroup = 0
                self.DrawIntervalsAndEdges(Exon, NextExon, NextExonInGroup, 0)
                self.DrawIntervalsAndEdges(Exon, NextExon, NextExonInGroup, 1)
            # Draw border with the next exon-group:
            if (GroupIndex < len(ExonGroups)):
                self.Draw.line((Right, 0, Right, Height), Colors.Border)
            self.GroupLeft += self.GroupWidth + 1
        # Finally, draw any edges that skip exons:
        self.DrawExonSkippingEdges(ExonList, 0)
        self.DrawExonSkippingEdges(ExonList, 1)
        OutputFileName = "GeneMap\\%s.png"%ProteinID
        self.Image.save(OutputFileName, "png")
    def GetX(self, GenomePosition):
        X = (GenomePosition - self.GroupGenomeStart) # could scale here
        X = max(0, min(X, self.GroupWidth))
        X += self.GroupLeft
        return X
    def DrawExonSkippingEdges(self, ExonList, GFFlag):
        if GFFlag:
            if self.TrueStrand == 1:
                IntervalList = self.GFIntervalsForward
                JunctionList = self.GFJunctionsForward
                JunctionDict = self.GFJunctionDictForward
            else:
                IntervalList = self.GFIntervalsReverse
                JunctionList = self.GFJunctionsReverse
                JunctionDict = self.GFJunctionDictReverse
            StartY = self.GenomeY + 75
            YDir = 1
        else:
            if self.TrueStrand == 1:
                IntervalList = self.ESTIntervalsForward
                JunctionList = self.ESTJunctionsForward
                JunctionDict = self.ESTJunctionDictForward
            else:
                IntervalList = self.ESTIntervalsReverse
                JunctionList = self.ESTJunctionsReverse
                JunctionDict = self.ESTJunctionDictReverse
            IntervalStartY = self.GenomeY - 8
            StartY = self.GenomeY - 47
            YDir = -1
        Y = StartY
        for ExonIndexA in range(len(ExonList) - 2):
            for ExonIndexB in range(ExonIndexA + 2, len(ExonList)):
                Junction = (ExonList[ExonIndexA].End, ExonList[ExonIndexB].Start)
                if JunctionDict.has_key(Junction):
                    # It's an exon-spanning edge!
                    print "Exon-spanning edge!  Exon %s to exon %s"%(ExonIndexA, ExonIndexB)
                    X1 = ExonList[ExonIndexA].RightX
                    X2 = ExonList[ExonIndexB].LeftX
                    while (1):
                        Pixel = self.Image.getpixel((X1, Y))
                        if Pixel[0]!=255:
                            Y += YDir*3
                            if (Y > 248 or Y < 1):
                                Y = StartY
                                break
                            continue
                        self.Draw.line((X1, Y, X2, Y), Colors.SpanningEdge)
                        self.Draw.line((X1, Y+YDir, X2, Y+YDir), Colors.SpanningEdge)
                        break
    def DrawIntervalsAndEdges(self, Exon, NextExon, NextExonInGroup, GFFlag):
        if GFFlag:
            if self.TrueStrand == 1:
                IntervalList = self.GFIntervalsForward
                JunctionList = self.GFJunctionsForward
            else:
                IntervalList = self.GFIntervalsReverse
                JunctionList = self.GFJunctionsReverse
            IntervalStartY = self.GenomeY + 28
            EdgeStartY = self.GenomeY + 30
            YDir = 1
        else:
            if self.TrueStrand == 1:
                IntervalList = self.ESTIntervalsForward
                JunctionList = self.ESTJunctionsForward
            else:
                IntervalList = self.ESTIntervalsReverse
                JunctionList = self.ESTJunctionsReverse
            IntervalStartY = self.GenomeY - 8
            EdgeStartY = self.GenomeY - 10
            YDir = -1
        BestGuys = []
        FarLeft = Exon.Start - 10
        for (Start, End) in IntervalList:
            if End < Exon.Start:
                continue
            if Start > Exon.End:
                break
            Distance = abs(Exon.Start - Start) + abs(Exon.End - End)
            BestGuys.append((Distance, Start, End))
        BestGuys.sort()
        Y = IntervalStartY
        for (Score, Start, End) in BestGuys[:10]:
            if Score == 0:
                Color = Colors.ESTCorrect
            else:
                Color = Colors.ESTWrong
            X1 = self.GetX(Start)
            X2 = self.GetX(End)
            self.Draw.line((X1, Y, X2, Y), Color)
            self.Draw.line((X1, Y+YDir, X2, Y+YDir), Color)
            if Score > 0:
                # Color the tips of the interval if endpoints are correct.
                if (Start == Exon.Start):
                    self.Draw.line((X1, Y, X1+3, Y), Colors.ESTCorrect)
                    self.Draw.line((X1, Y+YDir, X1+3, Y+YDir), Colors.ESTCorrect)
                if (End == Exon.End):
                    self.Draw.line((X2 - 3, Y, X2, Y), Colors.ESTCorrect)
                    self.Draw.line((X2 - 3, Y+YDir, X2, Y+YDir), Colors.ESTCorrect)
            Y += 4*YDir
        if not NextExon:
            return            
        #######################
        FarLeft = Exon.End - 25
        FarRight = Exon.End + 25
        EndFarRight = NextExon.End - 25
        BestGuys = []
        for (Start, End) in JunctionList:
            if Start < FarLeft:
                continue
            if Start > FarRight:
                break
            if End > EndFarRight:
                continue
            Distance = abs(Exon.End - Start) + abs(NextExon.Start - End)
            BestGuys.append((Distance, Start, End))
        BestGuys.sort()
        Y = EdgeStartY
        for (Distance, Start, End) in BestGuys[:10]:
            if Distance == 0:
                Color = Colors.LinkCorrect
            else:
                Color = Colors.LinkWrong
            X1 = self.GetX(Start)
            X2 = self.GetX(End)
            if not NextExonInGroup:
                # Special case: This edge extends forward to the next exon-group!
                if X2 >= self.GroupRight:
                    X2 = self.GroupRight + 50
                    X2 += (End - NextExon.Start)
                    X2 = max(self.GroupRight, X2)
            print "Exon edge from %s to %s (true %s to %s),\n  that means x %s to %s"%(Start, End, Exon.End, NextExon.Start, X1, X2)
            print X1, Y, X2, Y, Color
            self.Draw.line((X1, Y, X2, Y), Color)
            self.Draw.line((X1, Y+YDir, X2, Y+YDir), Color)
            Y += 4*YDir
    def ReportCoverageByExon(self, Aligner, ChromosomeNumber, OutputFileName):
        """
        Alternate pathway to ProcessBatch:
        - Iterate over ALL the proteins in the input file (with correct chromosome)
        - Keep a list of ALL exons...but ignore extremely short exons
        - Report coverage at the exon level
        """
        OutputFile = open(OutputFileName, "wb")
        Exons = {} # interval -> occurrence count
        ReverseExons = {} # interval -> occurrence count
        Introns = {}
        ReverseIntrons = {}
        Keys = Aligner.Proteins.keys()
        Keys.sort()
        MinimumExonLength = 50 # exons must be at least this many nucleotides in length.
        for Key in Keys:
            Protein = Aligner.Proteins[Key]
            # Skip proteins with no alignment:
            print "Protein %s has %s seeds"%(Key, len(Protein.Seeds))
            if not Protein.Seeds:
                continue
            # Skip proteins aligned to another chromosome:
            Seed = Protein.Seeds[0]
            print "Seed is on chromosome %s"%Seed.Chromosome
            if Seed.Chromosome != ChromosomeNumber:
                continue
            # Skip all proteins that aren't PERFECTLY covered!
            if Seed.ResultCoveragePercent < 1:
                continue
            # Remember all the exons:
            for ExonIndex in range(len(Seed.ResultExons)):
                ExonTuple = Seed.ResultExons[ExonIndex]
                Tuple = (ExonTuple[0], ExonTuple[1])
                ExonLength = ExonTuple[1] - ExonTuple[0]
                if Seed.Strand == 1:
                    Exons[Tuple] = Exons.get(Tuple, 0) + 1
                else:
                    ReverseExons[Tuple] = ReverseExons.get(Tuple, 0) + 1
                if ExonIndex < len(Seed.ResultExons) - 1:
                    NextExonTuple = Seed.ResultExons[ExonIndex + 1]
                    NextExonLength = NextExonTuple[1] - NextExonTuple[0]
                    if ExonLength > MinimumExonLength and NextExonLength > MinimumExonLength:
                        Tuple = (ExonTuple[1], NextExonTuple[0])
                        if Seed.Strand == 1:
                            Introns[Tuple] = Introns.get(Tuple, 0) + 1
                        else:
                            ReverseIntrons[Tuple] = ReverseIntrons.get(Tuple, 0) + 1
        ################################################
        # Iterate over exons, measure coverage at exon and at nt level:
        AllBaseCount = 0
        AllExonCount = 0
        TooShortExonCount = 0
        AllExonMultiCount = 0
        ESTBaseCoverCount = 0
        ESTExonCoverCount = 0
        EitherBaseCoverCount = 0
        EitherExonCoverCount = 0
        GFBaseCoverCount = 0
        GFExonCoverCount = 0
        for (Strand, ExonDict) in [(1, Exons), (-1, ReverseExons)]:
            if Strand == 1:
                ESTIntervalDict = self.ESTIntervalDictForward
                GFIntervalDict = self.GFIntervalDictForward
                ESTIntervalList = self.ESTIntervalsForward
                GFIntervalList = self.GFIntervalsForward
            else:
                ESTIntervalDict = self.ESTIntervalDictReverse
                GFIntervalDict = self.GFIntervalDictReverse
                ESTIntervalList = self.ESTIntervalsReverse
                GFIntervalList = self.GFIntervalsReverse
            for (Interval, Count) in ExonDict.items():
                (Start, End) = Interval
                ExonLength = (End - Start)
                if ExonLength < MinimumExonLength:
                    TooShortExonCount += 1
                    continue
                AllExonMultiCount += Count
                AllExonCount += 1
                AllBaseCount += ExonLength
                ESTCoverage = [0] * ExonLength
                GFCoverage = [0] * ExonLength
                EitherCoverage = [0] * ExonLength
                ExonCoveredFlag = 0
                ####################################################################
                # EST coverage:
                if ESTIntervalDict.has_key(Interval):
                    ESTExonCoverCount += 1
                    ESTCoverage = [1] * ExonLength
                    ExonCoveredFlag = 1
                else:
                    for (IntervalStart, IntervalEnd) in ESTIntervalList:
                        if (IntervalEnd < Start):
                            continue
                        if (IntervalStart > End):
                            break
                        OverlapStart = max(IntervalStart, Start)
                        OverlapEnd = min(IntervalEnd, End)
                        for X in range(OverlapStart, OverlapEnd):
                            ESTCoverage[X - Start] = 1
                for CoverFlag in ESTCoverage:
                    ESTBaseCoverCount += CoverFlag
                ####################################################################
                # GF coverage (copypasta):
                if GFIntervalDict.has_key(Interval):
                    GFExonCoverCount += 1
                    GFCoverage = [1] * ExonLength
                    ExonCoveredFlag = 1
                else:
                    for (IntervalStart, IntervalEnd) in GFIntervalList:
                        if (IntervalEnd < Start):
                            continue
                        if (IntervalStart > End):
                            break
                        OverlapStart = max(IntervalStart, Start)
                        OverlapEnd = min(IntervalEnd, End)
                        for X in range(OverlapStart, OverlapEnd):
                            GFCoverage[X - Start] = 1
                for CoverFlag in GFCoverage:
                    GFBaseCoverCount += CoverFlag
                if (ExonCoveredFlag):
                    EitherExonCoverCount += 1
                for Index in range(len(ESTCoverage)):
                    EitherBaseCoverCount += (ESTCoverage[Index] or GFCoverage[Index])
        OutputFile.write("Exons: %s (%s occurrences, %s too short)\n"%(AllExonCount, AllExonMultiCount, TooShortExonCount))
        OutputFile.write("Bases: %s\n"%(AllBaseCount))
        OutputFile.write("EST Coverage: %s exons, %s bases\n"%(ESTExonCoverCount, ESTBaseCoverCount))
        OutputFile.write("GF Coverage: %s exons, %s bases\n"%(GFExonCoverCount, GFBaseCoverCount))
        ########################################################################
        # Now, measure intron coverage.
        IntronCount = 0
        IntronMultiCount = 0
        ESTIntronCoverage = 0
        GFIntronCoverage = 0
        EitherIntronCoverage = 0
        for (Strand, IntronDict) in [(1, Introns), (-1, ReverseIntrons)]:
            if Strand == 1:
                ESTJunctionDict = self.ESTJunctionDictForward
                ESTJunctionList = self.ESTJunctionsForward
                GFJunctionDict = self.GFJunctionDictForward
                GFJunctionList = self.GFJunctionsForward
            else:
                ESTJunctionDict = self.ESTJunctionDictReverse
                ESTJunctionList = self.ESTJunctionsReverse
                GFJunctionDict = self.GFJunctionDictReverse
                GFJunctionList = self.GFJunctionsReverse
            for (Junction, Count) in IntronDict.items():
                IntronCount += 1
                IntronMultiCount += Count
                (Start, End) = Junction
                CoverFlag = 0
                if ESTJunctionDict.has_key(Junction):
                    ESTIntronCoverage += 1
                    CoverFlag = 1
                if GFJunctionDict.has_key(Junction):
                    GFIntronCoverage += 1
                    CoverFlag = 1
                if CoverFlag:
                    EitherIntronCoverage += 1
        OutputFile.write("Introns: %s (%s occurrences)\n"%(IntronCount, IntronMultiCount))
        OutputFile.write("EST intron coverage: %s\n"%ESTIntronCoverage)
        OutputFile.write("GF intron coverage: %s\n"%GFIntronCoverage)
        Str = "Final\t%s\t"%ChromosomeNumber
        Str += "%s\t%s\t%s\t%s\t"%(AllExonCount, AllExonMultiCount, TooShortExonCount, AllBaseCount)
        Str += "%s\t%s\t"%(IntronCount, IntronMultiCount)
        Str += "%s\t%s\t%s\t"%(ESTBaseCoverCount, ESTExonCoverCount, ESTIntronCoverage)
        Str += "%s\t%s\t%s\t"%(GFBaseCoverCount, GFExonCoverCount, GFIntronCoverage)
        Str += "%s\t%s\t%s\t"%(EitherBaseCoverCount, EitherExonCoverCount, EitherIntronCoverage)
        print Str
        OutputFile.write(Str + "\n")
class SpliceFindColumnNumbers:
    ChromosomeNumber = 20
    Strand = 21
    ProteinLength = 22
    Exons = 25
    Mismatches = 26
    SNPs = 27
    ProteinID = 0
    ProteinName = 1

class GeneMappingColumnNumbers:
    ChromosomeNumber = 2
    Strand = 3
    ProteinLength = 5
    Coverage = 6
    Exons = 8
    Mismatches = 9
    SNPs = 10
    ProteinID = 0
    ProteinName = 1    

if __name__ == "__main__":
    import psyco
    psyco.full()
    if len(sys.argv) > 1:
        ChromosomeNumber = int(sys.argv[1])
        SpecialTarget = None
        Mapper = GeneMapper(ChromosomeNumber)
        OutputFileName = "ESTGFCoverage.%s.txt"%ChromosomeNumber
        if len(sys.argv) > 2:
            SpecialTarget = int(sys.argv[2])
            OutputFileName = "ESTGFCoverage.Special.txt"
            Mapper.GenerateImagesFlag = 1
        #OutputFileName = "%s.Coverage"%Stub
        print "Read known gene mappings..."
        Aligner = AlignProteomeToGenome.GenomeAligner()
        Aligner.Load("GeneMappings.Best.txt")
        Mapper.ReadGeneFinding()
        Mapper.ReadESTs()
        Mapper.ProcessBatch(Aligner, ChromosomeNumber, OutputFileName, SpecialTarget)
        #Mapper.ReportCoverageByExon(Aligner, ChromosomeNumber, OutputFileName)
        sys.exit(1)
