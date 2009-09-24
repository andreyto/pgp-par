import os
import sys
import math
import time
import traceback
import struct

try:
    from PIL import Image
    from PIL import ImageDraw
    from PIL import ImageFont
except:
    print "No PIL found."
    pass

USE_QUICK_EXONS = 1 # default
# Add current directory to the path, for use on the grid:
sys.path.append(".")
try:
    import QuickExons
except:
    print "* Warning * QuickExons not available."
    USE_QUICK_EXONS = 0
VERBOSE = 0

##GeneMapper:
##This script takes, as input, a protein sequence and a "seed" chromosome and position.
##As output, it provides a list of genomic intervals which encode the protein - or
##at least, as much as the protein as possible.  Intervals are selected with attention
##to splice boundaries - often there are multiple exon boundaries with the same
##translation, but typically only one will have a decent splice boundary.
##ALL INTERVALS are half-open, including the lower but not the upper index.

def PrintObjects():
    "For debugging of excess memory usage"
    import gc
    TypeCounts = {}
    for Object in gc.get_objects():
        TypeCounts[type(Object)] = TypeCounts.get(type(Object), 0) + 1
    SortList = []
    for (Type, Count) in TypeCounts.items():
        SortList.append((Count, Type))
    SortList.sort()
    SortList.reverse()
    print "\nGC report:"
    for (Count, Type) in SortList:
        print Count, Type

# Fonts are broken on Linux.  (Tried pdf, pcf, and pil formats...but no luck)
# So, we'll content ourselves with a hideous default font if we must:
try:
    TheFont = ImageFont.truetype("Times.ttf", 12)
except:
    try:
        TheFont = ImageFont.load_default()
    except:
        TheFont = None

GeneticCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
               "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
               "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
               "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
               "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
               "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
               "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
               "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
               "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
               }
RCDict = {"A":"T", "G":"C", "T":"A", "C":"G",
          "a":"T", "g":"C", "t":"A", "c":"G"}

#A    |    x
#C A G|G U x A G U . . . . . . . . C A G|G 
#3 4 5 6 7 8 9 1011                15161718
SpliceBoundaryProfile = [
      {"A":0.3501  ,"C":0.3420  ,"G":0.1873  ,"T":0.1207},
      {"A":0.6085  ,"C":0.1094  ,"G":0.1220  ,"T":0.1601},
      {"A":0.0956  ,"C":0.0381  ,"G":0.7933  ,"T":0.0730},
      {"A":0.0028  ,"C":0.0010  ,"G":0.9949  ,"T":0.0014},
      {"A":0.0028  ,"C":0.0124  ,"G":0.0012  ,"T":0.9836},
      {"A":0.5512  ,"C":0.0304  ,"G":0.3846  ,"T":0.0338},
      {"A":0.7011  ,"C":0.0764  ,"G":0.1142  ,"T":0.1084},
      {"A":0.0787  ,"C":0.0545  ,"G":0.8026  ,"T":0.0642},
      {"A":0.1570  ,"C":0.1569  ,"G":0.1926  ,"T":0.4935},
      {"A":0.0560  ,"C":0.6400  ,"G":0.0094  ,"T":0.2947},
      {"A":0.9937  ,"C":0.0011  ,"G":0.0013  ,"T":0.0039},
      {"A":0.0011  ,"C":0.0014  ,"G":0.9943  ,"T":0.0032},
      {"A":0.2453  ,"C":0.1438  ,"G":0.4948  ,"T":0.1160},
    ]

def GetSpliceSignalScore(SpliceSignal):
    Score = 0
    for Pos in range(len(SpliceSignal)):
        Frequency = SpliceBoundaryProfile[Pos].get(SpliceSignal[Pos], 0.25)
        Score += math.log(Frequency)
    return Score

def ReverseComplement(DNA):
    "Returns the reverse complement of a DNA sequence."
    Str = ""
    for Index in range(len(DNA) - 1, -1, -1):
        Str += RCDict.get(DNA[Index], DNA[Index])
    return Str
def Translate(DNA):
    "Returns the peptide translation of a sequence."
    Peptide = ""
    for Index in range(0, len(DNA) - 2, 3):
        Codon = DNA[Index:Index+3].upper()
        AA = GeneticCode.get(Codon, "")
        Peptide += AA
        # Include a character for stop codons, and then end it:
        #if AA == "X":
        #    break
    return Peptide

class ExonClass:
    "An exon, as determined by aligning protein sequence against genomic sequence."
    def __init__(self):
        self.Start = None
        self.End = None
        self.Offset = None
        self.ProteinStart = None
        self.ProteinEnd = None
        self.Pos = None
        self.Mismatches = [] # list of codon-centers where mismatches were seen
        self.SNPs = [] # list of tuples of the form (GenomePosition, ProteinPos) for SNPs
        self.MismatchVerbose = []
        # If the edge leading forward from this exon produces an amino acid, then
        # EdgeResidue is the index of that AA within the protein.
        self.EdgeResidue = None
        self.ReadStart = None
        self.ReadEnd = None
        self.PrevExon = None
    def __cmp__(self, Other):
        "Exons are sorted by ascending genomic position"
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
    def DropLostFeatures(self):
        """
        If we've shifted the boundaries of the exon, then we can forget about any
        mismatches or SNPs which lie outside the exon:
        """
        SNPIndex = 0
        while SNPIndex < len(self.SNPs):
            (ProteinPos, GenomePos) = self.SNPs[SNPIndex]
            if ProteinPos < self.ProteinStart or ProteinPos >= self.ProteinEnd:
                del self.SNPs[SNPIndex]
            else:
                SNPIndex += 1
        MMIndex = 0
        while MMIndex < len(self.Mismatches):
            (ProteinPos, GenomePos) = self.Mismatches[MMIndex]
            if ProteinPos < self.ProteinStart or ProteinPos >= self.ProteinEnd:
                del self.Mismatches[MMIndex]
            else:
                MMIndex += 1
        
class GeneMapper:
    """
    Main application class: Reads in a genomic neighborhood, maps a protein
    to it, then reports the intervals and coverage.
    """
    Span = 4000000
    HalfSpan = 2000000
    UseLongWords = {33202:1, 36298: 1, 17453: 1, 11639: 1} # long repetitive runs
    def __init__(self, ChromosomeNumber):
        self.ChromosomeNumber = ChromosomeNumber
        self.OutputFile = None
        self.SpecialTarget = None
        self.WordSize = 6
    def ReadSNPs(self):
        """
        Parse SNPs from binary file.  If we find mismatches, we'll check whether the amino acid
        can be 'fixed' by picking a good SNP.
        """
        self.SNPDict = {}
        File = open(os.path.join("SNP", "%s.snp"%self.ChromosomeNumber), "rb")
        SNPCount = 0
        while (1):
            Data = File.read(4)
            if not Data:
                break
            SNPCount += 1
            try:
                Pos = struct.unpack("<i", Data)[0]
                Type = ord(File.read(1))
                # type 0 -> 2 forms, type 1 -> 3 forms, type 2 -> all 4 forms
                Forms = File.read(Type + 2)
                self.SNPDict[Pos] = Forms
            except:
                print SNPCount, len(Data), Data, Pos, Forms
                raise
        File.close()
    def ReadGenome(self, Seed):
        "Read in the chunk of the genome including position Seed"
        self.Seed = Seed
        self.Start = max(0, Seed - self.HalfSpan)
        self.Start -= (self.Start % 3) # start should be multiple of three, for easier codons
        self.GenomeFilePath = "e:\\chromosome\\Chr%s.trie"%self.ChromosomeNumber
        if not os.path.exists(self.GenomeFilePath):
            self.GenomeFilePath = "/scratch/stanner/blast/Chromosome/Chr%s.trie"%self.ChromosomeNumber
        File = open(self.GenomeFilePath,"rb")
        if VERBOSE:
            print "Seed %s Seek to %s"%(Seed, self.Start)
        File.seek(self.Start)
        self.GenomeSequence = File.read(self.Span).upper()
        self.Length = len(self.GenomeSequence)
        if VERBOSE:
            print "Genome sequence: '%s' ... '%s'"%(self.GenomeSequence[:20], self.GenomeSequence[-20:]) #%%%
            print "Forward0:", Translate(self.GenomeSequence[:90])
            print "Forward1:", Translate(self.GenomeSequence[1:91])
            print "Forward2:", Translate(self.GenomeSequence[2:92])
        File.close()
        self.End = self.Start + len(self.GenomeSequence) # non-inclusive
    def GenerateExonList(self):
        """
        Input: self.WordPositions, a list mapping from a protein residue number to a
        list of tuples of the form (ChromosomePosition, Strand).
        Output: A list of genome-disjoint exons, all from the same strand, such that the number of
        WordPositions contained in the exons in maximal.
        Algorithm: Build up (greedily) possible exons based on WordPositions.  For each exon,
        compute PrefixScore (the maximum number of amino acids covered by any exon chain
        strictly preceding the exon) and Score (PrefixScore plus the exon's length).
        Then take the exon with maximal score, and read back along its chain.  
        """
        ForwardExons = []
        ReverseExons = []
        for ProteinPos in range(len(self.WordPositions)):
            WordList = self.WordPositions[ProteinPos]
            #print "%s: %s %s words (%s exons so far)"%(ProteinPos, self.ProteinSequence[ProteinPos:ProteinPos + self.WordSize], len(WordList), len(ForwardExons))
            for (ChromosomePosition, Strand) in WordList:
                if Strand == 1:
                    # Try to enlarge an existing exon:
                    FoundFlag = 0
                    for OldExon in ForwardExons:
                        # If we have almost-full overlap with the old exon, then assimilate into it:
                        if OldExon.End == ChromosomePosition + (self.WordSize - 1) * 3 and OldExon.ProteinEnd == ProteinPos + (self.WordSize - 1):
                            OldExon.End += 3
                            OldExon.ProteinEnd += 1
                            FoundFlag = 1
                            #break
                    if not FoundFlag:
                        # Build and add a new exon:
                        Exon = ExonClass()
                        Exon.Start = ChromosomePosition
                        Exon.End = ChromosomePosition + (self.WordSize * 3)
                        Exon.ProteinStart = ProteinPos
                        Exon.ProteinEnd = ProteinPos + self.WordSize
                        Exon.ReadStart = Exon.Start
                        Exon.ReadEnd = Exon.End
                        ForwardExons.append(Exon)
                else:
                    # Try to enlarge an existing exon:
                    FoundFlag = 0
                    for OldExon in ReverseExons:
                        if OldExon.Start == ChromosomePosition - ((self.WordSize - 1) * 3) + 1 and OldExon.ProteinEnd == ProteinPos + (self.WordSize - 1):
                            OldExon.Start -= 3
                            OldExon.ProteinEnd += 1
                            FoundFlag = 1
                            break
                    if not FoundFlag:
                        # Build and add a new exon:
                        Exon = ExonClass()
                        Exon.Start = ChromosomePosition - (self.WordSize * 3) + 1
                        Exon.End = ChromosomePosition + 1
                        Exon.ProteinStart = ProteinPos
                        Exon.ProteinEnd = ProteinPos + self.WordSize
                        Exon.ReadStart = Exon.Start
                        Exon.ReadEnd = Exon.End                        
                        ReverseExons.append(Exon)
        ####################################################################################
        # Find best chain on the forward strand.
        if len(ForwardExons) > 1000:
            # PRUNE unreasonably short exons:
            for Exon in ForwardExons[:]:
                if Exon.ProteinEnd - Exon.ProteinStart < 10:
                    ForwardExons.remove(Exon)        
        ForwardExons.sort()
        if VERBOSE:
            for Exon in ForwardExons:
                print "%s %s-%s (%s-%s) %s"%(Exon.ProteinEnd - Exon.ProteinStart, Exon.ProteinStart, Exon.ProteinEnd, Exon.Start, Exon.End,
                                          self.ProteinSequence[Exon.ProteinStart:Exon.ProteinEnd])
        
        # Keep track of the best chain endpoint seen:
        BestChainEnd = None
        for X in range(len(ForwardExons)):
            # Compute PrefixScore, PrefixExon, and Score for each exon.
            Exon = ForwardExons[X]
            Exon.PrefixScore = 0
            for PrevExonIndex in range(X):
                PrevExon = ForwardExons[PrevExonIndex]
                # We can link back to an exon which begins earlier in the protein *and* earlier in
                # the genome.  We MAY overlap somewhat with the preceding exon, in protein positions
                # (if exons stretched over the splice boundary) and possibly in genomic positions (if
                # the reference sequence contains an insertion).  
                if PrevExon.Start < Exon.Start and PrevExon.ProteinStart < Exon.ProteinStart:
                    End = min(Exon.ProteinStart, PrevExon.ProteinEnd)
                    if (PrevExon.End > Exon.Start):
                        OverlapSize = (PrevExon.End - Exon.Start + 2) / 3
                        End = min(End, PrevExon.ProteinEnd - OverlapSize)
                    Score = PrevExon.PrefixScore + (End - PrevExon.ProteinStart)
                    if Score >= Exon.PrefixScore:
                        Exon.PrefixScore = Score
                        Exon.PrevExon = PrevExon
            Exon.Score = Exon.PrefixScore + (Exon.ProteinEnd - Exon.ProteinStart)
            if (BestChainEnd == None or Exon.Score > BestChainEnd.Score):
                BestChainEnd = Exon
        BestForwardChainScore = 0
        BestForwardChain = []
        if BestChainEnd:
            # Put the best-chain exons into a list:
            BestForwardChainScore = BestChainEnd.Score
            Exon = BestChainEnd
            while Exon:
                BestForwardChain.append(Exon)
                PrevExon = Exon.PrevExon
                # Eliminate any genomic overlap now:
                if (PrevExon and PrevExon.End > Exon.Start):
                    OverlapSize = (PrevExon.End - Exon.Start + 2) / 3
                    PrevExon.End = Exon.Start
                    PrevExon.ProteinEnd -= OverlapSize
                    if VERBOSE:
                        print "REPAIRED genomic overlap, dropped %s amino acids"%OverlapSize
                Exon = PrevExon
        ####################################################################################
        # Find best chain on the reverse strand.
        if len(ReverseExons) > 1000:
            # PRUNE unreasonably short exons:
            for Exon in ReverseExons[:]:
                if Exon.ProteinEnd - Exon.ProteinStart < 10:
                    ReverseExons.remove(Exon)        
        ReverseExons.sort()
        ReverseExons.reverse() # oriented forward along protein, backward along chromosome
        if VERBOSE:
            for Exon in ReverseExons:
                print "%s %s-%s (%s-%s) %s"%(Exon.ProteinEnd - Exon.ProteinStart, Exon.ProteinStart, Exon.ProteinEnd, Exon.Start, Exon.End,
                                          self.ProteinSequence[Exon.ProteinStart:Exon.ProteinEnd])
        BestChainEnd = None
        for X in range(len(ReverseExons)):
            # Compute PrefixScore, PrefixExon, and Score for each exon.
            # Keep track of the best end seen to the chain:            
            Exon = ReverseExons[X]
            Exon.PrefixScore = 0
            for PrevExonIndex in range(X):
                PrevExon = ReverseExons[PrevExonIndex]
                if PrevExon.End > Exon.End and PrevExon.ProteinStart < Exon.ProteinStart:
                    End = min(PrevExon.ProteinEnd, Exon.ProteinStart)
                    if (PrevExon.Start < Exon.End):
                        OverlapSize = (Exon.End - PrevExon.Start + 2) / 3
                        End = min(End, PrevExon.ProteinEnd - OverlapSize)
                    Score = PrevExon.PrefixScore + (End - PrevExon.ProteinStart)
                    if Score >= Exon.PrefixScore:
                        Exon.PrefixScore = Score
                        Exon.PrevExon = PrevExon
            Exon.Score = Exon.PrefixScore + (Exon.ProteinEnd - Exon.ProteinStart)
            #print "Exon %s-%s pref %s score %s"%(Exon.ProteinStart, Exon.ProteinEnd, Exon.PrefixScore, Exon.Score)
            if (BestChainEnd == None or Exon.Score > BestChainEnd.Score):
                BestChainEnd = Exon
        BestReverseChainScore = 0
        BestReverseChain = []
        if BestChainEnd:
            # Put the best-chain exons into a list:
            BestReverseChainScore = BestChainEnd.Score
            Exon = BestChainEnd
            while Exon:
                BestReverseChain.append(Exon)
                PrevExon = Exon.PrevExon
                if VERBOSE:
                    if PrevExon:
                        print "Exon %s-%s (%s-%s) prev %s-%s (%s-%s)"%(Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd,
                            PrevExon.Start, PrevExon.End, PrevExon.ProteinStart, PrevExon.ProteinEnd)
                # Eliminate any genomic overlap now:
                if (PrevExon and PrevExon.Start < Exon.End):
                    OverlapSize = (Exon.End - PrevExon.Start + 2) / 3
                    PrevExon.Start = Exon.End
                    PrevExon.ProteinEnd -= OverlapSize
                    if VERBOSE:
                        print "REPAIRED genomic overlap, dropped %s amino acids"%OverlapSize
                        print "PrevExon corrected to %s-%s (%s-%s)"%(PrevExon.Start, PrevExon.End, PrevExon.ProteinStart, PrevExon.ProteinEnd)
                Exon = PrevExon
        ############################################################################
        # Return the best forward chain OR best reverse chain, whichever scored highest.
        if BestForwardChainScore >= BestReverseChainScore:
            BestForwardChain.sort()
            ExonList = BestForwardChain
            self.TrueStrand = 1
        else:
            BestReverseChain.sort()
            ExonList = BestReverseChain
            self.TrueStrand = -1
        # Remove pruned-away exons:
        for Exon in ExonList[:]:
            if Exon.ProteinEnd == Exon.ProteinStart:
                ExonList.remove(Exon)
        self.ExonCoverageFlags = [0] * self.ProteinLen
        for Exon in ExonList:
            for X in range(Exon.ProteinStart, Exon.ProteinEnd):
                self.ExonCoverageFlags[X] = 1
        return ExonList
    def DebugPrintExonList(self, ExonList):
        print "Exons (%d):"%len(ExonList)
        for Index in range(len(ExonList)):
            Exon = ExonList[Index]
            AACount = Exon.ProteinEnd - Exon.ProteinStart
            if Exon.EdgeResidue:
                AACount += 1
            print " Exon#%d: P%s..%s %s..%s (len %s) AA%s %s"%(Index, Exon.ProteinStart, Exon.ProteinEnd, Exon.Start, Exon.End, Exon.End - Exon.Start, AACount, Exon.EdgeResidue)
            print "    %s"%(self.ProteinSequence[Exon.ProteinStart:Exon.ProteinEnd])
            if self.TrueStrand == 1:
                Prefix = self.GenomeSequence[Exon.ReadStart - 3 - self.Start:Exon.ReadStart - self.Start]
                Suffix = self.GenomeSequence[Exon.ReadEnd - self.Start:Exon.ReadEnd + 3 - self.Start]
            else:
                Prefix = self.GenomeSequence[Exon.ReadEnd + 1 - self.Start:Exon.ReadEnd + 4 - self.Start]
                Suffix = self.GenomeSequence[Exon.ReadStart - 3 - self.Start:Exon.ReadStart - self.Start]
            print "    Prefix %s suffix %s"%(Prefix, Suffix)
            if len(Exon.Mismatches):
                print "    %s mismatches: %s"%(len(Exon.Mismatches), Exon.Mismatches)
            if Index < len(ExonList) - 1:
                print "     Junction %s->%s"%(Exon.End, ExonList[Index + 1].Start)    
    def GetWordPositions(self):
        # WordDict[Word] is a list of protein positions where the seed Word occurs.
        WordDict = {}        
        self.WordPositions = []
        for X in range(self.ProteinLen):
            self.WordPositions.append([])
        for Pos in range(self.ProteinLen - self.WordSize + 1):
            Word = self.ProteinSequence[Pos:Pos + self.WordSize]
            if not WordDict.has_key(Word):
                WordDict[Word] = []
            WordDict[Word].append(Pos)
        ##################################################################
        # Find, and make note of, all the words coded on the forward strand:
        for Offset in range(3):
            print "Reading frame %s..."%Offset
            CodedProtein = ""
            GenomePos = Offset + self.Start
            SequencePos = Offset
            # Step forward by one codon at a time, remembering the most recent 8 amino acids:
            while (SequencePos < self.Length - 2):
                Codon = self.GenomeSequence[SequencePos:SequencePos + 3]
                AA = GeneticCode.get(Codon, "X")
                CodedProtein = (CodedProtein + AA)[-self.WordSize:]
                PosList = WordDict.get(CodedProtein, [])
                for Pos in PosList:
                    CodingStart = GenomePos - (self.WordSize - 1)*3
                    self.WordPositions[Pos].append((CodingStart, 1))
                    if VERBOSE:
                        print "Word %s '%s' encoded by %s-%s: %s"%(Pos, self.ProteinSequence[Pos:Pos + self.WordSize],
                            CodingStart, GenomePos + 3,
                            self.GenomeSequence[CodingStart - self.Start:SequencePos + 3])
                SequencePos += 3
                GenomePos += 3
        ##################################################################
        # Find, and make note of, all the words coded on the reverse strand:    
        for Offset in range(3):
            print "Reverse strand, reading frame %s..."%Offset
            CodedProtein = ""
            GenomePos = self.End - 1 - Offset
            SequencePos = self.Length - 1 - Offset
            while (SequencePos > 1): # last codon uses indices 2,1,0
                Codon = self.GenomeSequence[SequencePos - 2:SequencePos + 1]
                Codon = ReverseComplement(Codon)
                AA = GeneticCode.get(Codon, "X")
                CodedProtein = (CodedProtein + AA)[-self.WordSize:]
                PosList = WordDict.get(CodedProtein, [])
                for Pos in PosList:
                    CodingStart = GenomePos + (self.WordSize - 1)*3
                    self.WordPositions[Pos].append((CodingStart, -1))
                    if VERBOSE:
                        print "Word %s '%s' encoded by %s-%s: %s"%(Pos, self.ProteinSequence[Pos:Pos+self.WordSize],
                            GenomePos - 2, CodingStart + 1,
                            ReverseComplement(self.GenomeSequence[SequencePos - 2:CodingStart + 1 - self.Start]))
                SequencePos -= 3
                GenomePos -= 3

    def GenerateExonSeeds(self):
        """
        Populate the list WordPositions, where WordPositions[X] is a list of tuples (Location, Strand)
        giving genomic locations where the word beginning at protein position X is encoded.  
        """
        self.WordSize = 6
        if USE_QUICK_EXONS:
            # New way:
            if self.UseLongWords.has_key(self.ProteinID):
                self.WordPositions = QuickExons.GetWordPositions(self.ProteinSequence, self.GenomeFilePath, self.Seed, 10)
                self.WordSize = 10
            else:
                self.WordPositions = QuickExons.GetWordPositions(self.ProteinSequence, self.GenomeFilePath, self.Seed, 6)
                self.WordSize = 6
        else:
            # Old way:
            self.GetWordPositions()
            
        ##################################################################
        # Note the total word-coverage of the protein.  It should be 100%
        # other than mismatches and codons split across splice junctions
        CoverageFlags = [0]*self.ProteinLen
        ReverseCoverageFlags = [0]*self.ProteinLen
        for WordStart in range(self.ProteinLen):
            for (GenomePos, Strand) in self.WordPositions[WordStart]:
                for Pos in range(WordStart, WordStart + self.WordSize):
                    if Strand == 1:
                        CoverageFlags[Pos] = 1
                    else:
                        ReverseCoverageFlags[Pos] = 1
        TotalCoverage = 0
        TotalReverseCoverage = 0
        for X in range(len(self.ProteinSequence)):
            if CoverageFlags[X]:
                TotalCoverage += 1
            if ReverseCoverageFlags[X]:
                TotalReverseCoverage += 1
        WordCoverage = max(TotalCoverage, TotalReverseCoverage)
        if TotalCoverage >= TotalReverseCoverage:
            WordCoverage = TotalCoverage
        else:
            WordCoverage = TotalReverseCoverage
        CoveragePercent = WordCoverage / float(self.ProteinLen)
        print "Word-coverage: %s of %s residues (%.1f%%)"%(WordCoverage, self.ProteinLen, 100*CoveragePercent)
    def OutputExonList(self, ExonList):
        """
        Produce output, to self.OutputFile and/or stdout.  Provide:
        ProteinID, Chromosome, Strand, self.ProteinLen, Coverage, ExonList
        """
        Str = "%s\t%s\t%s\t%s\t%s\t"%(self.ProteinID, self.ProteinName, self.ChromosomeNumber, self.TrueStrand, self.Seed)
        Coverage = 0
        for Exon in ExonList:
            Coverage += (Exon.ProteinEnd - Exon.ProteinStart)
            if Exon.EdgeResidue != None:
                Coverage += 1
        Str += "%s\t%s\t%s\t"%(self.ProteinLen, Coverage, len(ExonList))
        ExonStr = ""
        MismatchList = []
        SNPList = []
        for Exon in ExonList:
            if Exon.EdgeResidue != None:
                ExonStr += "%s-%s-%s-%s-%s,"%(Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd, Exon.EdgeResidue)
            else:
                ExonStr += "%s-%s-%s-%s,"%(Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd)
            MismatchList.extend(Exon.Mismatches)
            SNPList.extend(Exon.SNPs)
        Str += "%s\t"%(ExonStr[:-1]) # remove trailing comma
        Str += "%s\t%s\t"%(MismatchList, SNPList)
        print Str
        sys.stdout.flush()
        if self.OutputFile:
            self.OutputFile.write(Str)
            self.OutputFile.write("\n")
            self.OutputFile.flush()
    def RefineExonEdges(self, ExonList):
        """
        This is where we 'polish' the edges between exons, so that junctions are splice
        boundaries wherever possible, and so that exons don't cover protein sequence
        reduntantly
        """
        for ExonIndex in range(len(ExonList) - 1):
            Exon = ExonList[ExonIndex]
            NextExon = ExonList[ExonIndex + 1]
            if self.TrueStrand == 1:
                # FORWARD strand:
                AAOverlap = Exon.ProteinEnd - NextExon.ProteinStart
                if AAOverlap < -3:
                    # These two exons have a GAP between them in the protein sequence which
                    # can't be bridged by an edge.
                    continue
                if AAOverlap <= 0:
                    # There's a gap of one or two amino acids.  So, extend the exons step closer together,
                    # even if that introduces mismatches, so that we can reach a splice boundary.
                    ExtendAA = max(1, -AAOverlap)
                    for X in range(1, ExtendAA + 1):
                        # Extend Exon forward:
                        ProteinPos = Exon.ProteinEnd + X - 1
                        AA = self.ProteinSequence[ProteinPos]
                        CodonStart = Exon.End + (X - 1) * 3
                        (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                        if AAStr and AA == AAStr[0]:
                            pass
                        elif AA in AAStr:
                            Exon.SNPs.append((ProteinPos, CodonStart))
                        else:
                            Exon.Mismatches.append((ProteinPos, CodonStart))
                        # Extend NextExon backward:
                        ProteinPos = NextExon.ProteinStart - X
                        AA = self.ProteinSequence[ProteinPos]
                        CodonStart = NextExon.Start - X*3
                        (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                        if AAStr and AA == AAStr[0]:
                            pass
                        elif AA in AAStr:
                            NextExon.SNPs.append((ProteinPos, CodonStart))
                        else:
                            NextExon.Mismatches.append((ProteinPos, CodonStart))
                    Exon.ProteinEnd += ExtendAA
                    Exon.End += ExtendAA*3
                    NextExon.ProteinStart -= ExtendAA
                    NextExon.Start -= ExtendAA*3
                    AAOverlap = Exon.ProteinEnd - NextExon.ProteinStart
                MaxLeftShift = AAOverlap*3 + 2                
                if VERBOSE:
                    print
                    print "AAOverlap %s MaxLeftShift %s"%(AAOverlap, MaxLeftShift)
                    print "Exon A (%d): %d-%d (%d-%d)"%(ExonIndex, Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd)
                    print "Exon B (%d): %d-%d (%d-%d)"%(ExonIndex + 1, NextExon.Start, NextExon.End, NextExon.ProteinStart, NextExon.ProteinEnd)
                    print "Exon A mismatches:", Exon.Mismatches
                    print "Exon B mismatches:", NextExon.Mismatches
                BestShift = None
                BestShiftScore = (-9999, -9999, -9999)
                for LeftShift in range(-2, MaxLeftShift):
                    # Consider setting Exon.End = Exon.End - LeftShift,
                    # and NextExon.Start = NextExon.Start - LeftShift + 3*AAOverlap
                    LeftEdge = Exon.End - LeftShift
                    RightEdge = NextExon.Start - LeftShift + 3*AAOverlap
                    # First: Does it look like a decent splice signal?
                    SpliceSignal = self.GenomeSequence[LeftEdge - 3 - self.Start:LeftEdge + 6 - self.Start]
                    SpliceSignal += self.GenomeSequence[RightEdge - 3 - self.Start:RightEdge + 1 - self.Start]
                    SpliceSignalScore = GetSpliceSignalScore(SpliceSignal)
                    # Next: How many amino acids do we get with this exon pair?
                    # Count AAs in the left exon:
                    if LeftShift == -3:
                        ExonProteinEnd = Exon.ProteinEnd + 1
                    else:
                        ExonProteinEnd = Exon.ProteinEnd - (LeftShift + 2) / 3
                    if (ExonProteinEnd < 0):
                        continue # invalid shift!
                    AACountLeft = (ExonProteinEnd - Exon.ProteinStart)
                    for (ProteinPos, CodonPos) in Exon.Mismatches:
                        if ProteinPos < ExonProteinEnd:
                            AACountLeft -= 1
                    # Count AAs in the right exon:
                    NextExonProteinStart = NextExon.ProteinStart + AAOverlap - (LeftShift / 3)
                    AACountRight = (NextExon.ProteinEnd - NextExonProteinStart)
                    for (ProteinPos, CodonPos) in NextExon.Mismatches:
                        if ProteinPos >= NextExonProteinStart:
                            AACountRight -= 1
                    # Count one AA from the edge:
                    if (LeftShift % 3 == 0):
                        LinkCodon = ""
                    elif (LeftShift % 3 == 1):
                        # Get length-2 suffix from left exon, length-1 prefix from right exon
                        LinkCodon = self.GenomeSequence[LeftEdge - 2 - self.Start] + \
                            self.GenomeSequence[LeftEdge - 1 - self.Start] + \
                            self.GenomeSequence[RightEdge - self.Start]
                    else:
                        # Get length-1 suffix from left exon, length-2 prefix from right exon
                        LinkCodon = self.GenomeSequence[LeftEdge - 1 - self.Start] + \
                            self.GenomeSequence[RightEdge - self.Start] + \
                            self.GenomeSequence[RightEdge + 1 - self.Start]
                    if GeneticCode.get(LinkCodon, "X") == self.ProteinSequence[ExonProteinEnd]:
                        AACountLink = 1
                        LinkResidue = ExonProteinEnd
                    else:
                        AACountLink = 0
                        LinkResidue = None
                    AACount = AACountLeft + AACountLink + AACountRight
                    if VERBOSE:
                        print "LeftShift %d: Splice '%s' score %.3f AACount %d"%(LeftShift, SpliceSignal, SpliceSignalScore, AACount)
                        print "AACount: %d(%d) + %d + %d(%d)"%(AACountLeft, ExonProteinEnd - Exon.ProteinStart, AACountLink,
                            AACountRight, NextExon.ProteinEnd - NextExonProteinStart)
                        print "EPE %s seq %s cod %s %s"%(ExonProteinEnd, self.ProteinSequence[ExonProteinEnd], LinkCodon, GeneticCode.get(LinkCodon, "X"))
                    if SpliceSignalScore >= -13:
                        SpliceSignalFlag = 1
                    else:
                        SpliceSignalFlag = 0
                    ShiftScore = (SpliceSignalFlag, AACount, SpliceSignalScore)
                    if ShiftScore > BestShiftScore:
                        BestShiftScore = ShiftScore
                        BestShift = LeftShift
                        BestLinkResidue = LinkResidue
                #####################################
                # Apply the shift:
                Exon.End -= BestShift
                NextExon.Start = NextExon.Start + 3*AAOverlap - BestShift
                if BestShift == -3:
                    Exon.ProteinEnd = Exon.ProteinEnd + 1
                else:
                    Exon.ProteinEnd = Exon.ProteinEnd - (BestShift + 2) / 3
                Exon.EdgeResidue = BestLinkResidue
                Exon.SpliceBoundaryScore = BestShiftScore[-1] # splice score is the last member in the tuple
                NextExon.ProteinStart = NextExon.ProteinStart + AAOverlap - (BestShift / 3)
                Exon.DropLostFeatures()
                NextExon.DropLostFeatures()
                if VERBOSE:
                    print "After shift %s: %s-%s (%s-%s) and %s-%s (%s-%s), edgeres %s"%(BestShift, Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd,
                        NextExon.Start, NextExon.End, NextExon.ProteinStart, NextExon.ProteinEnd, Exon.EdgeResidue)
            else:
                # REVERSE strand: COPY-PASTA from the forward strand.
                AAOverlap = NextExon.ProteinEnd - Exon.ProteinStart 
                if AAOverlap < -3:
                    # These two exons have a GAP between them in the protein sequence which
                    # can't be bridged by an edge.
                    continue
                if AAOverlap <= 0:
                    # There's a gap of one or two amino acids.  So, extend the exons step closer together,
                    # even if that introduces mismatches, so that we can reach a splice boundary.
                    ExtendAA = max(1, -AAOverlap)
                    for X in range(1, ExtendAA + 1):
                        # Extend Exon forward:
                        ProteinPos = Exon.ProteinStart - X
                        AA = self.ProteinSequence[ProteinPos]
                        CodonStart = Exon.End - 1 + (X * 3)
                        (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                        if AAStr and AA == AAStr[0]:
                            pass
                        elif AA in AAStr:
                            Exon.SNPs.append((ProteinPos, CodonStart))
                        else:
                            Exon.Mismatches.append((ProteinPos, CodonStart))
                        # Extend NextExon backward:
                        ProteinPos = NextExon.ProteinEnd + X - 1
                        AA = self.ProteinSequence[ProteinPos]
                        CodonStart = NextExon.Start + 2 - (X * 3)
                        (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                        if AAStr and AA == AAStr[0]:
                            pass
                        elif AA in AAStr:
                            NextExon.SNPs.append((ProteinPos, CodonStart))
                        else:
                            NextExon.Mismatches.append((ProteinPos, CodonStart))
                    Exon.ProteinStart -= ExtendAA
                    Exon.End += ExtendAA*3
                    NextExon.ProteinEnd += ExtendAA
                    NextExon.Start -= ExtendAA*3
                    AAOverlap = NextExon.ProteinEnd - Exon.ProteinStart 
                MaxLeftShift = AAOverlap*3 + 2
                if VERBOSE:
                    print
                    print "AAOverlap %s MaxLeftShift %s"%(AAOverlap, MaxLeftShift)
                    print "Exon A (%d): %d-%d (%d-%d)"%(ExonIndex, Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd)
                    print "Exon B (%d): %d-%d (%d-%d)"%(ExonIndex + 1, NextExon.Start, NextExon.End, NextExon.ProteinStart, NextExon.ProteinEnd)
                    print "Exon A mismatches:", Exon.Mismatches
                    print "Exon B mismatches:", NextExon.Mismatches
                BestShift = None
                BestShiftScore = (-9999, -9999, -9999)
                for LeftShift in range(-2, MaxLeftShift):
                    # Consider setting Exon.End = Exon.End - LeftShift,
                    # and NextExon.Start = NextExon.Start - LeftShift + 3*AAOverlap
                    LeftEdge = Exon.End - LeftShift
                    RightEdge = NextExon.Start - LeftShift + 3*AAOverlap
                    # First: Does it look like a decent splice signal?
                    SpliceSignal = self.GenomeSequence[LeftEdge - 1 - self.Start:LeftEdge + 3 - self.Start]
                    SpliceSignal += self.GenomeSequence[RightEdge - 6 - self.Start:RightEdge + 3 - self.Start]
                    SpliceSignal = ReverseComplement(SpliceSignal)
                    SpliceSignalScore = GetSpliceSignalScore(SpliceSignal)
                    # Next: How many amino acids do we get with this exon pair?
                    # Count AAs in the left exon:
                    AACount = 0
                    if LeftShift == -3:
                        ExonProteinStart = Exon.ProteinStart - 1
                    else:
                        ExonProteinStart = Exon.ProteinStart + (LeftShift + 2) / 3
                    if (ExonProteinStart >= self.ProteinLen):
                        continue # invalid shift!
                    AACount += (Exon.ProteinEnd - Exon.ProteinStart)
                    for (ProteinPos, CodonPos) in Exon.Mismatches:
                        if ProteinPos >= ExonProteinStart:
                            AACount -= 1
                    # Count AAs in the right exon:
                    NextExonProteinEnd = NextExon.ProteinEnd - AAOverlap + (LeftShift / 3)
                    AACount += (NextExonProteinEnd - NextExon.ProteinStart)
                    for (ProteinPos, CodonPos) in NextExon.Mismatches:
                        if ProteinPos < NextExon.ProteinEnd:
                            AACount -= 1
                    # Count one AA from the edge:
                    if (LeftShift % 3 == 0):
                        LinkCodon = ""
                    elif (LeftShift % 3 == 1):
                        # Get length-2 suffix from left exon, length-1 prefix from right exon
                        LinkCodon = self.GenomeSequence[LeftEdge - 2 - self.Start] + \
                            self.GenomeSequence[LeftEdge - 1 - self.Start] + \
                            self.GenomeSequence[RightEdge - self.Start]
                    else:
                        # Get length-1 suffix from left exon, length-2 prefix from right exon
                        LinkCodon = self.GenomeSequence[LeftEdge - 1 - self.Start] + \
                            self.GenomeSequence[RightEdge - self.Start] + \
                            self.GenomeSequence[RightEdge + 1 - self.Start]
                    LinkCodon = ReverseComplement(LinkCodon)
                    try:
                        self.ProteinSequence[ExonProteinStart - 1]
                    except:
                        print "!???"
                        print ExonProteinStart
                        print LinkCodon
                        print LeftShift
                        print "ProtLen:", self.ProteinSequence
                        raise
                    if GeneticCode.get(LinkCodon, "X") == self.ProteinSequence[ExonProteinStart - 1]:
                        AACount += 1
                        LinkResidue = ExonProteinStart - 1
                    else:
                        LinkResidue = None
                    if VERBOSE:
                        print "LeftShift %d: Splice '%s' score %.3f AACount %d"%(LeftShift, SpliceSignal,
                            SpliceSignalScore, AACount)
                    if SpliceSignalScore >= -13:
                        SpliceSignalFlag = 1
                    else:
                        SpliceSignalFlag = 0
                    ShiftScore = (SpliceSignalFlag, AACount, SpliceSignalScore)
                    if ShiftScore > BestShiftScore:
                        BestShiftScore = ShiftScore
                        BestShift = LeftShift
                        BestLinkResidue = LinkResidue
                #####################################
                # Apply the shift:
                Exon.End -= BestShift
                NextExon.Start = NextExon.Start + 3*AAOverlap - BestShift 
                if BestShift == -3:
                    Exon.ProteinStart = Exon.ProteinStart - 1
                else:
                    Exon.ProteinStart = Exon.ProteinStart + (BestShift + 2) / 3
                Exon.EdgeResidue = BestLinkResidue
                Exon.SpliceBoundaryScore = BestShiftScore[-1] # splice score is the last member in the tuple
                NextExon.ProteinEnd = NextExon.ProteinEnd - AAOverlap + (BestShift / 3)
                Exon.DropLostFeatures()
                NextExon.DropLostFeatures()
                
    def ExtendExons(self, ExonList):
        """
        We've generated an initial exon list.  Now we'll extend exons.  We're allowed
        to cross n mismatches if we match n bases.  We stop extending if hit too many
        mismatches, or we come to the edge of the protein
        """
        for ExonIndex in range(len(ExonList)):
            Exon = ExonList[ExonIndex]
            if ExonIndex:
                PrevExon = ExonList[ExonIndex - 1]
            else:
                PrevExon = None
            if ExonIndex < len(ExonList) - 1:
                NextExon = ExonList[ExonIndex + 1]
            else:
                NextExon = None
            if self.TrueStrand == 1:
                ##################################
                # Extend back along the genome
                BestExtension = 0
                BestMatchCount = 0
                MaxExtension = Exon.ProteinStart
                # Don't overlap the genomic interval of the previous exon:
                if PrevExon:
                    MaxExtension = min(MaxExtension, (Exon.Start - PrevExon.End) / 3)
                    if VERBOSE:
                        print "Match extension: %s, %s, %s"%(Exon.Start, PrevExon.End, MaxExtension) #%%%
                MatchCount = 0
                SNPList = []
                MismatchList = []
                MismatchCount = 0                
                for Extension in range(1, MaxExtension + 1):
                    CodonStart = Exon.Start - Extension * 3
                    (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                    ProteinPos = Exon.ProteinStart - Extension
                    AA = self.ProteinSequence[ProteinPos]
                    if AAStr and AA == AAStr[0]:
                        MatchCount += 1
                    elif AA in AAStr:
                        MatchCount += 1
                        SNPList.append((ProteinPos, CodonStart))
                    else:
                        MismatchCount += 1
                        MismatchList.append((ProteinPos, CodonStart))
                    if (MatchCount >= MismatchCount and MatchCount > BestMatchCount):
                        BestExtension = Extension
                        BestMatchCount = MatchCount
                    if MatchCount + (MaxExtension - Extension) < MismatchCount:
                        break # branch cut: This extension can't possibly recover.
                # Move the exon's start and protein-start:
                Exon.Start -= BestExtension * 3
                Exon.ProteinStart -= BestExtension
                # Note any SNPs and mismatches that we required:
                for (ProteinPos, Pos) in SNPList:
                    if ProteinPos >= Exon.ProteinStart:
                        Exon.SNPs.append((ProteinPos, Pos))
                for (ProteinPos, Pos) in MismatchList:
                    if ProteinPos >= Exon.ProteinStart:
                        Exon.Mismatches.append((ProteinPos, Pos))
                ##################################
                # Extend FORWARD along the genome
                BestExtension = 0
                BestMatchCount = 0
                MaxExtension = self.ProteinLen - Exon.ProteinEnd
                # Don't overlap the genomic interval of the previous exon:
                if NextExon:
                    MaxExtension = min(MaxExtension, (NextExon.Start - Exon.End) / 3)
                    if VERBOSE:
                        print "Match extension: %s, %s, %s"%(NextExon.Start, Exon.End, MaxExtension) #%%%                
                MatchCount = 0
                SNPList = []
                MismatchList = []
                MismatchCount = 0                
                for Extension in range(1, MaxExtension + 1):
                    CodonStart = Exon.End + (Extension - 1) * 3
                    (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                    ProteinPos = Exon.ProteinEnd + Extension - 1
                    AA = self.ProteinSequence[ProteinPos]
                    if AAStr and AA == AAStr[0]:
                        MatchCount += 1
                    elif AA in AAStr:
                        MatchCount += 1
                        SNPList.append((ProteinPos, CodonStart))
                    else:
                        MismatchCount += 1
                        MismatchList.append((ProteinPos, CodonStart))
                    if (MatchCount >= MismatchCount and MatchCount > BestMatchCount):
                        BestExtension = Extension
                        BestMatchCount = MatchCount
                    if MatchCount + (MaxExtension - Extension) < MismatchCount:
                        break # branch cut: This extension can't possibly recover.
                        
                # Move the exon's end and protein-end:
                Exon.End += BestExtension * 3
                Exon.ProteinEnd += BestExtension
                # Note any SNPs and mismatches that we required:
                for (ProteinPos, Pos) in SNPList:
                    if ProteinPos < Exon.ProteinEnd:
                        Exon.SNPs.append((ProteinPos, Pos))
                for (ProteinPos, Pos) in MismatchList:
                    if ProteinPos < Exon.ProteinEnd:
                        Exon.Mismatches.append((ProteinPos, Pos))
            else:
                # Reverse strand (copy-pasta from forward strand)
                # Extend forward along the genome
                BestExtension = 0
                BestMatchCount = 0
                MaxExtension = Exon.ProteinStart
                # Don't overlap the genomic interval of the previous exon:
                if NextExon:
                    MaxExtension = min(MaxExtension, (NextExon.Start - Exon.End) / 3)
                    if VERBOSE:
                        print "Match extension: %s, %s, %s"%(NextExon.Start, Exon.End, MaxExtension) #%%%
                MatchCount = 0
                SNPList = []
                MismatchList = []
                MismatchCount = 0                
                for Extension in range(1, MaxExtension + 1):
                    CodonStart = Exon.End - 1 + (Extension * 3)
                    (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                    ProteinPos = Exon.ProteinStart - Extension
                    AA = self.ProteinSequence[ProteinPos]
                    #print "Extending exon %s (%s-%s) (%s-%s) forward step %s:"%(ExonIndex, Exon.Start, Exon.End, Exon.ProteinStart,Exon.ProteinEnd, Extension)
                    #print "Pos %s Codon %s AAstr%s desired %s"%(ProteinPos, Codon, AAStr, AA)
                    if AAStr and AA == AAStr[0]:
                        MatchCount += 1
                        #print "Pos%s match!"%ProteinPos
                    elif AA in AAStr:
                        MatchCount += 1
                        SNPList.append((ProteinPos, CodonStart))
                        #print "Pos%s SNP!"%ProteinPos
                    else:
                        MismatchCount += 1
                        MismatchList.append((ProteinPos, CodonStart))
                        #print "Pos%s mismatch!"%ProteinPos
                    if (MatchCount >= MismatchCount and MatchCount > BestMatchCount):
                        BestExtension = Extension
                        BestMatchCount = MatchCount
                    if MatchCount + (MaxExtension - Extension) < MismatchCount:
                        break # branch cut: This extension can't possibly recover.
                # Move the exon's end and protein-start:
                Exon.End += BestExtension * 3
                Exon.ProteinStart -= BestExtension
                if VERBOSE:
                    print "Extend %s AA"%BestExtension
                # Note any SNPs and mismatches that we required:
                for (ProteinPos, Pos) in SNPList:
                    if ProteinPos >= Exon.ProteinStart:
                        Exon.SNPs.append((ProteinPos, Pos))
                for (ProteinPos, Pos) in MismatchList:
                    if ProteinPos >= Exon.ProteinStart:
                        Exon.Mismatches.append((ProteinPos, Pos))
                ##################################
                # Extend BACKWARD along the genome
                BestExtension = 0
                BestMatchCount = 0
                MaxExtension = self.ProteinLen - Exon.ProteinEnd
                # Don't overlap the genomic interval of the previous exon:
                if PrevExon:
                    MaxExtension = min(MaxExtension, (Exon.Start - PrevExon.End) / 3)
                    if VERBOSE:
                        print "Match extension: %s, %s, %s"%(Exon.Start, PrevExon.End, MaxExtension) #%%%
                MatchCount = 0
                SNPList = []
                MismatchList = []
                MismatchCount = 0                
                for Extension in range(1, MaxExtension + 1):
                    CodonStart = Exon.Start + 2 - (Extension * 3)
                    (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                    ProteinPos = Exon.ProteinEnd + Extension - 1
                    AA = self.ProteinSequence[ProteinPos]
                    #print "Extending exon %s (%s-%s) (%s-%s) back step %s:"%(ExonIndex, Exon.Start, Exon.End, Exon.ProteinStart,Exon.ProteinEnd, Extension)
                    #print "Pos %s Codon %s AAstr%s desired %s"%(ProteinPos, Codon, AAStr, AA)
                    if AAStr and AA == AAStr[0]:
                        MatchCount += 1
                        #print "Pos%s match"%ProteinPos
                    elif AA in AAStr:
                        MatchCount += 1
                        SNPList.append((ProteinPos, CodonStart))
                        #print "Pos%s SNP"%ProteinPos
                    else:
                        MismatchCount += 1
                        MismatchList.append((ProteinPos, CodonStart))
                        #print "Pos%s mismatch"%ProteinPos
                    if (MatchCount >= MismatchCount and MatchCount > BestMatchCount):
                        BestExtension = Extension
                        BestMatchCount = MatchCount
                    if MatchCount + (MaxExtension - Extension) < MismatchCount:
                        break # branch cut: This extension can't possibly recover.
                # Move the exon's start and protein-end:
                Exon.Start -= BestExtension * 3
                Exon.ProteinEnd += BestExtension
                # Note any SNPs and mismatches that we required:
                for (ProteinPos, Pos) in SNPList:
                    if ProteinPos < Exon.ProteinEnd:
                        Exon.SNPs.append((ProteinPos, Pos))
                for (ProteinPos, Pos) in MismatchList:
                    if ProteinPos < Exon.ProteinEnd:
                        Exon.Mismatches.append((ProteinPos, Pos))
        # A cleanup step: It's very unlikely, but we may have completely covered
        # an exon's protein sequence using a neighboring exon.  If that's the case,
        # then we can dump the spanned exon:
        ExonIndex = 0
        while ExonIndex < len(ExonList):
            Exon = ExonList[ExonIndex]
            if ExonIndex > 0:
                PrevExon = ExonList[ExonIndex - 1]
                if self.TrueStrand == 1:
                    if PrevExon.ProteinEnd >= Exon.ProteinEnd:
                        del ExonList[ExonIndex]
                        continue
                else:
                    if PrevExon.ProteinStart <= Exon.ProteinStart:
                        print "*** Exon assimilated!  Prev %s-%s, exon %s %s-%s"%(PrevExon.ProteinStart, PrevExon.ProteinEnd,
                            ExonIndex, Exon.ProteinStart, Exon.ProteinEnd)
                        del ExonList[ExonIndex]
                        continue
            if ExonIndex < len(ExonList) - 1:
                NextExon = ExonList[ExonIndex + 1]
                if self.TrueStrand == 1:
                    if NextExon.ProteinStart <= Exon.ProteinStart:
                        del ExonList[ExonIndex]
                        ExonIndex = max(ExonIndex - 1, 0)
                        continue
                else:
                    if NextExon.ProteinEnd >= Exon.ProteinEnd:
                        print "*** Exon assimilated!  Next %s-%s, exon %s %s-%s"%(NextExon.ProteinStart, NextExon.ProteinEnd,
                            ExonIndex, Exon.ProteinStart, Exon.ProteinEnd)
                        del ExonList[ExonIndex]
                        ExonIndex = max(ExonIndex - 1, 0)
                        continue
            ExonIndex += 1
        return ExonList
    def GetPossibleTranslations(self, CodonStart):
        """
        Return a tuple (StandardCodon, PossibleAA), where PossibleAA is
        usually a string of length 1.  If there are SNP(s), then PossibleAA
        is a string of all possible encoded amino acids, STARTING with the
        standard one.
        """        
        if self.TrueStrand == 1:
            StandardCodon = self.GenomeSequence[CodonStart - self.Start:CodonStart - self.Start + 3]
            BasesA = self.SNPDict.get(CodonStart, StandardCodon[0])
            BasesB = self.SNPDict.get(CodonStart + 1, StandardCodon[1])
            BasesC = self.SNPDict.get(CodonStart + 2, StandardCodon[2])
            AAStr = Translate(StandardCodon)
            #print "GetPossibleTranslations:", BasesA, BasesB, BasesC
            for BaseA in BasesA:
                for BaseB in BasesB:
                    for BaseC in BasesC:
                        Codon = "%s%s%s"%(BaseA, BaseB, BaseC)
                        AA = Translate(Codon)
                        AAStr += AA
        else:
            StandardCodon = self.GenomeSequence[CodonStart - self.Start - 2:CodonStart - self.Start + 1]
            StandardCodon = ReverseComplement(StandardCodon)
            BasesA = self.SNPDict.get(CodonStart, StandardCodon[0])
            BasesA = RCDict.get(BasesA, BasesA)
            BasesB = self.SNPDict.get(CodonStart - 1, StandardCodon[1])
            BasesB = RCDict.get(BasesB, BasesB)
            BasesC = self.SNPDict.get(CodonStart - 2, StandardCodon[2])
            BasesC = RCDict.get(BasesC, BasesC)
            AAStr = Translate(StandardCodon)
            #print "GetPossibleTranslations:", BasesA, BasesB, BasesC
            for BaseA in BasesA:
                for BaseB in BasesB:
                    for BaseC in BasesC:
                        Codon = "%s%s%s"%(BaseA, BaseB, BaseC)
                        AA = Translate(Codon)
                        AAStr += AA
        return (StandardCodon, AAStr)
            
    def MergeCloseExons(self,ExonList):
        """
        If exons are close (3-27 base pairs) with compatible reading frame, then probably
        they just have mismatch(es) separating them.  Glue them together into one long
        exon, and remember the mismatch.
        """
        ExonIndex = 0
        while (ExonIndex < len(ExonList) - 1):
            Exon = ExonList[ExonIndex]
            NextExon = ExonList[ExonIndex + 1]
            Dist = NextExon.Start - Exon.End
            #print "** Exon %s and %s: end %s start %s gap %s"%(ExonIndex, ExonIndex + 1, Exon.End, NextExon.Start, Dist)            
            if Dist in (0, 3, 6, 9, 12, 15, 18, 21, 24, 27):
                AADist = Dist / 3
                if self.TrueStrand == 1:
                    # Glue together two exons on the forward strand:
                    for Index in range(AADist):
                        CodonStart = Exon.End + Index * 3 
                        (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                        print "Extend exon %s-%s step %s:"%(Exon.Start, Exon.End, Index)
                        print "CodonStart %s codon %s translations %s"%(Exon.End + Index*3, Codon, AAStr)
                        AAStrPos = AAStr.find(self.ProteinSequence[Exon.ProteinEnd + Index])
                        if AAStrPos == 0:
                            # Not a mismatch or a SNP, just a match that's not contained in a word
                            # because of flanking mismatches.
                            self.ExonCoverageFlags[Exon.ProteinEnd + Index] = 1
                        elif AAStrPos != -1:
                            self.ExonCoverageFlags[Exon.ProteinEnd + Index] = 1
                            Exon.SNPs.append((Exon.ProteinEnd + Index, CodonStart))
                        else:
                            Exon.Mismatches.append((Exon.ProteinEnd + Index, CodonStart))
                            #Exon.MismatchVerbose.append("%s (not %s) at protein pos %s"%(AAStr, self.ProteinSequence[Exon.ProteinEnd + Index], Exon.ProteinEnd + Index))
                    Exon.ProteinEnd = NextExon.ProteinEnd
                else:
                    # Reverse strand:
                    for Index in range(AADist):
                        CodonStart = Exon.End + 2 + Index * 3 # codon is Exon.End, Exon.End + 1, Exon.End + 2
                        (Codon, AAStr) = self.GetPossibleTranslations(CodonStart)
                        print "Extend exon %s-%s step %s: CodonStart %s codon %s translations %s"%(Exon.Start, Exon.End, Index, Exon.End + Index*3, Codon, AAStr)
                        AAStrPos = AAStr.find(self.ProteinSequence[Exon.ProteinStart - Index - 1])
                        if AAStrPos == 0:
                            # Not a mismatch or a SNP, just a match that's not in any word
                            # because of flanking mismaches
                            self.ExonCoverageFlags[Exon.ProteinStart - Index - 1] = 1
                        elif AAStrPos != -1:
                            self.ExonCoverageFlags[Exon.ProteinStart - Index - 1] = 1
                            Exon.SNPs.append((CodonStart, Exon.ProteinStart - Index - 1))
                        else:
                            #Codon = self.GenomeSequence[CodonStart:CodonStart + 3]
                            Exon.Mismatches.append((Exon.ProteinStart - Index - 1, CodonStart))
                            #Exon.MismatchVerbose.append("%s (not %s) at protein pos %s"%(AAStr, self.ProteinSequence[Exon.ProteinStart - Index], Exon.ProteinStart - Index))
                    Exon.ProteinStart = NextExon.ProteinStart
                Exon.End = NextExon.End
                del ExonList[ExonIndex + 1]
            elif Dist < 27:
                # The exons are separated by a short gap, but the gap size is not a multiple of three.
                # This may correspond to an indel of size 1 or so.  We'll call everything in the
                # seal a 'mismatch' for now because giving credit for SNPs here is a nuisance.
                print "MERGE across a STRANGE GAP dist %s"%Dist
                print "%s-%s (%s-%s) and %s-%s (%s-%s)"%(Exon.Start, Exon.End,
                    Exon.ProteinStart,Exon.ProteinEnd, NextExon.Start, NextExon.End,
                    NextExon.ProteinStart, NextExon.ProteinEnd)
                if self.TrueStrand == 1:
                    AADist = NextExon.ProteinStart - Exon.ProteinEnd
                    if abs(AADist*3 - Dist) < 6:
                        GenomePos = Exon.End
                        for AAPos in range(Exon.ProteinEnd, NextExon.ProteinStart):
                            Exon.Mismatches.append((AAPos, GenomePos))
                            GenomePos += 3
                        Exon.ProteinEnd = NextExon.ProteinEnd
                        Exon.End = NextExon.End                            
                        del ExonList[ExonIndex + 1]
                    else:
                        ExonIndex += 1
                else:
                    AADist = Exon.ProteinStart - NextExon.ProteinEnd
                    if abs(AADist*3 - Dist) < 6:
                        GenomePos = NextExon.Start - 3
                        for AAPos in range(NextExon.ProteinEnd, Exon.ProteinStart):
                            Exon.Mismatches.append((AAPos, GenomePos))
                            GenomePos -= 3
                        Exon.ProteinStart = NextExon.ProteinStart
                        Exon.End = NextExon.End                            
                        del ExonList[ExonIndex + 1]                        
                    else:
                        ExonIndex += 1
            else:
                ExonIndex += 1
        return ExonList
    def ProcessBatch(self, InputFileName, OutputFileName):
        """
        Call ReadGenome() and then AlignProtein() for a list of proteins, read from an input file.
        The input file should be tab-delimited.  Each line must have the following fields:
        PeptideID, PeptideName, ChromosomeLocation, self.ProteinSequence
        """
        self.OutputFile = open(OutputFileName, "w")
        HeaderLine = "#ProteinIndex\tProteinName\tChromosome\tStrand\tSeed\tLength\tCoverage\tExonCount\tExons\tMismatches\tSNPs\t\n"
        self.OutputFile.write(HeaderLine)
        InputFile = open(InputFileName, "rb")
        LineNumber = 0
        for FileLine in InputFile.xreadlines():
            LineNumber += 1
##            if LineNumber > 50:
##                break #%%% TEMP
            if FileLine[0] == "#":
                continue # comment!
            Bits = FileLine.split("\t")
            self.ProteinID = int(Bits[0])
            if self.SpecialTarget != None and int(Bits[0]) != self.SpecialTarget:
                continue
            self.ProteinName = Bits[1]
            ChromosomePosition = int(Bits[2])
            # Ensure that an X in the source sequence doesn't match ANYTHING:
            self.ProteinSequence = Bits[3].strip().replace("X", "$")
            print self.ProteinID, self.ProteinName
            try:
                self.ReadGenome(ChromosomePosition)
                self.AlignProtein(self.ProteinID, self.ProteinSequence)
            except:
                traceback.print_exc()
            print "\n"*10
        InputFile.close()
    def AlignProtein(self, ProteinID, ProteinSequence, Seed = None):
        """
        MAIN METHOD: Align a protein against our genome neighborhood.
        """
        self.ProteinSequence = ProteinSequence
        self.ProteinLen = len(self.ProteinSequence)
        self.ProteinID = ProteinID
        print "Handle protein %s (len %s)..."%(ProteinID, self.ProteinLen)
        # Prepare a list of exon seeds.  Exon seeds are genome positions encoding
        # short substrings of the protein, called "words".
        self.GenerateExonSeeds()
        # Find a chain of exons that captures as much of the protein as possible.
        # Some seeds are spurious matches, so we use d.p. to get the optimal chain.
        if VERBOSE:
            print "Generate exon list:"
        ExonList = self.GenerateExonList()
        if VERBOSE:        
            print "-=-"*20
            print "INITIAL exon list:"
            self.DebugPrintExonList(ExonList)
        ExonList = self.MergeCloseExons(ExonList)
        if VERBOSE:
            print "-=-"*20
            print "Exon list after MERGE:"
            self.DebugPrintExonList(ExonList)
        ExonList = self.ExtendExons(ExonList)
        if VERBOSE:
            print "-=-"*20
            print "Exon list after EXTEND:"
            self.DebugPrintExonList(ExonList)
        #ExonList = self.MergeCloseExons(ExonList)
        self.RefineExonEdges(ExonList)
        #ExonList = self.PruneShortExons(ExonList) # let's not, at least for now.
        if VERBOSE:
            print "-=-"*20
            print "Exon list after REFINE-EDGES:"
            self.DebugPrintExonList(ExonList)
        self.OutputExonList(ExonList)
        ####################################################
        # If we were passed a seed object, then store the results in the object:
        if Seed != None:
            Coverage = 0
            for Exon in ExonList:
                Coverage += (Exon.ProteinEnd - Exon.ProteinStart)
                if Exon.EdgeResidue != None:
                    Coverage += 1
            print "Update seed: Coverage is %s"%Coverage
            Seed.ResultCoverage = Coverage
            Seed.ResultCoveragePercent = Coverage / float(max(1, self.ProteinLen))
            Seed.ResultExons = []
            Seed.ResultSNPs = []
            Seed.ResultMismatches = []
            MismatchList = []
            SNPList = []
            for Exon in ExonList:
                if Exon.EdgeResidue != None:
                    Seed.ResultExons.append((Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd, Exon.EdgeResidue))
                    #ExonStr += "%s-%s-%s-%s-%s,"%(Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd, Exon.EdgeResidue)
                else:
                    #ExonStr += "%s-%s-%s-%s,"%(Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd)
                    Seed.ResultExons.append((Exon.Start, Exon.End, Exon.ProteinStart, Exon.ProteinEnd))
                Seed.ResultMismatches.extend(Exon.Mismatches)
                Seed.ResultSNPs.extend(Exon.SNPs)
            
    def PruneShortExons(self, ExonList):
        # Prune any extremely short exons, since a 6-mer exon is probably a spurious hit:
        for Exon in ExonList[:]:
            Length = Exon.ProteinEnd - Exon.ProteinStart
            if Length < 8:
                ExonList.remove(Exon)
        return ExonList
    
def TestMain():
    Mapper = GeneMapper(23)
    Mapper.ReadSNPs()
    Protein = "MSIVRLSVHAKWIMGKVTGTKMQKTAKVRVIRLVLDPHLLKYYNKQKTYFAHNALQQCTIGDIVLLKALP"
    Protein += "VPRTKHVKHELAEIVFKVGKLVDPVTGKPCAGTTYLESPLSSETTQGVDGASRPSRGPAPCRAGPGARHL"
    Protein += "RPWPESPRPEPRGLPGPGRGSMATWRRDGRLTGGQRLLCAGLAGTLSLSLTAPLELATVLAQVGVVRGHAR"
    Protein += "GPWATGHRVWRAEGLRALWKGNAVACLRLFPCSAVQLAAYRKFVVLFTDDLGHISQWSSIMAGSLAGMVSTI"
    Protein += "VTYPTDLIKTRLIMQNILEPSYRGLLHAFSTIYQQEGFLALYRGVSLTVVGALPFSAGSLLVYMNLEKIWNG"
    Protein += "PRDQFSLPQNFANVCLAAAVTQTLSFPFETVKRKMQSCCLFAYDGCVPLEAVMGRLAQWVTV"
    Mapper.ReadGenome(118315113)
    Mapper.ProteinName = "IPI00002036.3"
    Mapper.AlignProtein("IPI00002036.3", Protein)
def GetProteinSequence(ProteinID):
    IndexFile = open("Database\\IPIv315.index", "rb")
    DBFile = open("Database\\IPIv315.trie", "rb")                
    IndexFile.seek(92 * ProteinID + 8)
    DBFilePos = struct.unpack("<i", IndexFile.read(4))[0]
    DBFile.seek(DBFilePos)
    Sequence = ""
    while (1):
        Stuff = DBFile.read(1024)
        if not Stuff:
            break
        StarPos = Stuff.find("*")
        if StarPos == -1:
            Sequence += Stuff
            continue
        Sequence += Stuff[:StarPos]
        break
    IndexFile.close()
    DBFile.close()
    return Sequence

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        pass
    VERBOSE = 1
    ChromosomeNumber = int(sys.argv[1])
    Position = int(sys.argv[2])
    ProteinID = int(sys.argv[3])
    ProteinSequence = GetProteinSequence(ProteinID)
    ProteinSequence = ProteinSequence.replace("X", "$")
    Mapper = GeneMapper(ChromosomeNumber)
    Mapper.ReadSNPs()
    Mapper.ProteinID = ProteinID
    Mapper.ProteinName = ""
    Mapper.ReadGenome(Position)
    Mapper.AlignProtein(ProteinID, ProteinSequence)
    #TestMain()
    #import profile
    #profile.run("Main()")
    #sys.exit(1)
##    ################################################ %%%%%%%%%%%%%%%%
##    ChromosomeNumber = int(sys.argv[1])
##    Mapper = GeneMapper(ChromosomeNumber)
##    Mapper.ReadSNPs()
##    GeneMapperOutputFileName = "GeneMapperOutput.%s.txt"%ChromosomeNumber
##    InputFileName = "GeneMapperInput.%s.txt"%ChromosomeNumber
##    if len(sys.argv) > 2:
##        Mapper.SpecialTarget = int(sys.argv[2])
##        print "Special target:",Mapper.SpecialTarget
##        #InputFileName = "GeneMapperInput.Test.txt" #%%%
##        GeneMapperOutputFileName = "GeneMapperOutput.Test.txt"
##        VERBOSE = 1
##    Mapper.ProcessBatch(InputFileName, GeneMapperOutputFileName)
##    sys.exit(1)
##    ################################################ %%%%%%%%%%%%%%%%
##    if len(sys.argv) > 1:
##        InputFile = sys.argv[1]
##        ChromosomeNumber = int(sys.argv[2])
##        OutputFile = "%s.out"%InputFile
##        Mapper = GeneMapper(ChromosomeNumber)
##        Mapper.ReadGeneFinding()
##        Mapper.ReadESTs()
##        Mapper.ReadSNPs()
##        Mapper.ProcessBatch(InputFile, OutputFile)
##        sys.exit(1)
##    # Some test scaffolding for testing individual proteins
##    Mapper = GeneMapper(23)
##    Mapper.ReadSNPs()
##    Protein = "MSIVRLSVHAKWIMGKVTGTKMQKTAKVRVIRLVLDPHLLKYYNKQKTYFAHNALQQCTIGDIVLLKALP"
##    Protein += "VPRTKHVKHELAEIVFKVGKLVDPVTGKPCAGTTYLESPLSSETTQGVDGASRPSRGPAPCRAGPGARHL"
##    Protein += "RPWPESPRPEPRGLPGPGRGSMATWRRDGRLTGGQRLLCAGLAGTLSLSLTAPLELATVLAQVGVVRGHAR"
##    Protein += "GPWATGHRVWRAEGLRALWKGNAVACLRLFPCSAVQLAAYRKFVVLFTDDLGHISQWSSIMAGSLAGMVSTI"
##    Protein += "VTYPTDLIKTRLIMQNILEPSYRGLLHAFSTIYQQEGFLALYRGVSLTVVGALPFSAGSLLVYMNLEKIWNG"
##    Protein += "PRDQFSLPQNFANVCLAAAVTQTLSFPFETVKRKMQSCCLFAYDGCVPLEAVMGRLAQWVTV"
##    Mapper.ReadGenome(118315113)
##    Mapper.AlignProtein("IPI00002036.3", Protein)
