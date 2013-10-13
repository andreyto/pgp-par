"""
Parse the output of GeneID to get binary files of intervals.
"""
import os
import sys
import struct
import math
from Utils import *

#Species = "arabidopsis"
Species = "human"

# Add AT MOST this many splice junctions for a particular interval.
# (We don't want to have ridiculous numbers of edges, since that's not realistic,
# and it makes it very hard to disconnect break the exon graph into genes)
MaxLinks = 10

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

def GetSpliceAcceptorScore(Sequence):
    Score = 0
    for X in range(4):
        Freq = SpliceBoundaryProfile[9 + X].get(Sequence[X], 0.25)
        Score += math.log(Freq)
    return Score

def GetSpliceDonorScore(Sequence):
    Score = 0
    for X in range(9):
        Freq = SpliceBoundaryProfile[X].get(Sequence[X], 0.25)
        Score += math.log(Freq)
    return Score


class GeneFlags:
    """
    Flags for intervals.  We annotate the reading frame of an exon by where codons
    start, modulo 3.  That way, when we intersect and merge intervals in splicedb.c,
    we don't need to re-calculate the reading frame!
    """
    Frame0 = 1
    Frame1 = 2
    Frame2 = 4
    First = 8
    Last = 16
    AllReadingFrames = 7

class GeneIDInterval:
    def __init__(self):
        self.Start = None
        self.End = None
        self.Flags = 0
        self.Score = -99
        self.JunctionsBack = [] # entries of the form (EndPos, Score)
        # Prefix (and suffix) flags keep track whether we have a prefix (or suffix)
        # of the specified length.
        self.Prefixes = [0, 0, 0]
        self.Suffixes = [0, 0, 0]
    def AddLinkBack(self, Pos):
        for (EndPos, Score) in self.JunctionsBack:
            if Pos == EndPos:
                return
        self.JunctionsBack.append((Pos, 0))
    
class GeneIDParser:
    def __init__(self):
        # throw out all exon predictions below a cutoff:
        if Species == "human":
            self.ScoreCutoff = -1.0
        else:
            self.ScoreCutoff = 3.5 
        self.IntervalDict = {} # (Start, End) -> interval object
        self.LinkDict = {}
    def Parse(self, GeneIDPath, ReverseChar):
        LineNumber = 0
        self.ReverseChar = ReverseChar
        InFile = open(GeneIDPath, "rb")
        TotalBasePairs = 0
        for FileLine in InFile.xreadlines():
            LineNumber += 1
            if (LineNumber % 10000 == 0):
                print LineNumber, len(self.IntervalDict.keys())
                #break #%%%
            Bits = FileLine.split("\t")
            try:
                Start = int(Bits[3]) - 1
                End = int(Bits[4])
                Score = float(Bits[5])
                ReadingFrame = int(Bits[7])
            except:
                continue
            # change the meaning of reading frame from "first codon starts after prefix length n"
            # to "codons start at genomic position n modulo 3"
            PrefixLength = ReadingFrame
            SuffixLength = ((End - Start) - PrefixLength) % 3
            if ReverseFlag:
                # Important: The first codon of this putative exon starts at (End-1) if it's
                # in reading frame 0, at (End - 1 - ReadingFrame) in general.  
                ReadingFrame = (End - 1 - ReadingFrame) % 3
            else:
                ReadingFrame = (ReadingFrame + Start) % 3
            if ReadingFrame == 0:
                ReadingFrameFlag = GeneFlags.Frame0
            elif ReadingFrame == 1:
                ReadingFrameFlag = GeneFlags.Frame1
            elif ReadingFrame == 2:
                ReadingFrameFlag = GeneFlags.Frame2
            if Bits[6] != ReverseChar:
                continue # wrong strand!
            if Score < self.ScoreCutoff:
                continue # bleah, not very convincing
            if Bits[2] == "Internal":
                LinkageFlag = 0
            elif Bits[2] == "First":
                LinkageFlag = GeneFlags.First
            elif Bits[2] == "Terminal":
                LinkageFlag = GeneFlags.Last
            elif Bits[2] == "Single":
                LinkageFlag = GeneFlags.First & GeneFlags.Last
            else:
                print "Skip:", Bits[2]
                continue # Not an exon feature!
            Score = int(Score * 100) # convert to an int
            Key = (Start, End)
            OldInterval = self.IntervalDict.get(Key, None)
            # It's POSSIBLE (but rare) to see the same interval twice, with
            # different reading frames.  If so, flag both frames and take the
            # best score seen:
            if OldInterval:
                OldInterval.Score = max(OldInterval.Score, Score)
                OldInterval.Flags = OldInterval.Flags | ReadingFrameFlag
                OldInterval.Prefixes[PrefixLength] = 1
                OldInterval.Suffixes[SuffixLength] = 1
                continue
            Interval = GeneIDInterval()
            Interval.Start = Start # correct residue-numbering
            Interval.End = End # correct residue-numbering
            Interval.Flags = ReadingFrameFlag + LinkageFlag
            Interval.Prefixes[PrefixLength] = 1
            Interval.Suffixes[SuffixLength] = 1
            Interval.Score = Score
            self.IntervalDict[Key] = Interval
            TotalBasePairs += (End - Start)
        InFile.close()
        print "Accepted %s bases, %s intervals out of %s lines."%(TotalBasePairs, len(self.IntervalDict.keys()), LineNumber)
    def ScoreSpliceEndpoints(self, GenomeFilePath): #UNUSED
        """
        Score the splice donor and acceptor sites.  REDUNDANT, because GeneID already scored them.
        """
        self.SpliceDonors = {}
        self.SpliceAcceptors = {}
        GenomeFile = open(GenomeFilePath, "rb")
        # Process sites in order along the file, since that may be faster:
        Keys = self.IntervalDict.keys()
        Keys.sort()
        KeyCount = len(Keys)
        if self.ReverseChar == "+":
            KeyIndex = 0
            for Key in Keys:
                KeyIndex += 1
                if KeyIndex % 10000 == 0:
                    print "%s/%s..."%(KeyIndex, KeyCount)
                Interval = self.IntervalDict[Key]
                if not Interval.Flags & GeneFlags.First:
                    AcceptorScore = self.SpliceAcceptors.get(Interval.Start, None)
                    if AcceptorScore == None:
                        # Read the relevant sequence:
                        GenomeFile.seek(Interval.Start -  3)
                        Sequence = GenomeFile.read(4)
                        #print "Get splice acceptor score for %s-%s '%s'"%(Interval.Start, Interval.End, Sequence)
                        AcceptorScore = GetSpliceAcceptorScore(Sequence.upper())
                        self.SpliceAcceptors[Interval.Start] = AcceptorScore
                if not Interval.Flags & GeneFlags.Last:
                    DonorScore = self.SpliceDonors.get(Interval.End, None)
                    if DonorScore == None:
                        GenomeFile.seek(Interval.End - 3)
                        Sequence = GenomeFile.read(9)
                        #print "Get splice donor score for %s-%s '%s'"%(Interval.Start, Interval.End, Sequence)
                        DonorScore = GetSpliceDonorScore(Sequence.upper())
                        self.SpliceDonors[Interval.End] = DonorScore
##                print "Interval %s-%s: Acceptor score %.3f donor score %.3f"%(Interval.Start, Interval.End,
##                    self.SpliceAcceptors.get(Interval.Start, 0), self.SpliceDonors.get(Interval.End, 0))
        else:
            # Copy-pasta for the reverse strand
            for Key in Keys:
                Interval = self.IntervalDict[Key]
                if not Interval.Flags & GeneFlags.First:
                    AcceptorScore = self.SpliceAcceptors.get(Interval.End, None)
                    if AcceptorScore == None:
                        # Read the relevant sequence:
                        GenomeFile.seek(Interval.End)
                        Sequence = GenomeFile.read(4)
                        Sequence = ReverseComplement(Sequence)
                        AcceptorScore = GetSpliceAcceptorScore(Sequence.upper())
                        self.SpliceAcceptors[Interval.Start] = AcceptorScore
                if not Interval.Flags & GeneFlags.Last:
                    DonorScore = self.SpliceDonors.get(Interval.Start, None)
                    if DonorScore == None:
                        GenomeFile.seek(Interval.Start - 6)
                        Sequence = GenomeFile.read(9)
                        Sequence = ReverseComplement(Sequence)
                        DonorScore = GetSpliceDonorScore(Sequence.upper())
                        self.SpliceDonors[Interval.End] = DonorScore
        GenomeFile.close()
    def LinkExons(self, ReverseFlag):
        """
        We must decide which of our putative intervals are linked to each other.
        Assume for this discussion that we're on the forward strand:
        Interval A can link to B if:
        - A isn't flagged as terminal, and B isn't flagged as initial
        - B starts after A ends, but not TOO far after
        - B's reading frames include a legal reading frame after passing through A (e.g. if A has
          starting frame 0 and length 61, then B should be in reading frame 1)
        - We could examine the splice-signal strength as well, but splice-signal was already
          factored into the original exon scores.  
        We limit the degree of nodes in the graph: Each genomic coordinate can have at
        most ten edges linking in, at most ten edges linking out.  
        """
        
        # Record all the linkages, for evaluating our scoring model:
        AllLinkagesFile = open("GFLinkages.txt", "wb")
        
        #VerbosePoints = [47397674,47412711,47421273,47421575,47422292,47422540,47422899,47423106]
        VerbosePoints = [15919951,15944937,15948479,15981282,
                         15983102,15988080,15990852,15991318,
                         15992784,15994150,16008067,16009446,
                         16009807,16010677,16011225
                         ] #[171686859]
        VerbosePoints = []
        self.GenomicBackwardEdges = {}
        self.GenomicForwardEdges = {}
        MinIntronLength = 25
        MaxIntronLength = 20000
        Keys = self.IntervalDict.keys()
        Keys.sort()
        LinkCount = 0
        VerboseFlag = 0
        LenA = len(Keys)
        for IndexA in range(LenA):
            IntervalA = self.IntervalDict[Keys[IndexA]]
            if IntervalA.End in VerbosePoints:
                VerboseFlag = 1
                print
                print "Interval (%s-%s)"%(IntervalA.Start, IntervalA.End)
            else:
                VerboseFlag = 0
            if (IntervalA.Flags & GeneFlags.Last) and not ReverseFlag:
                continue
            if (IntervalA.Flags & GeneFlags.First) and ReverseFlag:
                continue
            if IndexA % 10000 == 0:
                print "%d/%d..."%(IndexA, LenA)
            # DesirableLinks is a list of entries of the form (score, intervalB), one
            # entry for each compatible interval we're considering linking forward to.
            DesirableLinks = {} # NextExonStart -> (Score, Interval)
            for IndexB in range(IndexA + 1, len(Keys)):
                IntervalB = self.IntervalDict[Keys[IndexB]]
                IntronLength = IntervalB.Start - IntervalA.End
                if IntronLength < MinIntronLength:
                    continue # keep looking for a mate                
                if VerboseFlag:
                    print "Consider linking to %s-%s intron %s"%(IntervalB.Start, IntervalB.End, IntronLength)
                    print "Suffixes %s prefixes %s"%(IntervalA.Suffixes, IntervalB.Prefixes)
                if IntronLength > MaxIntronLength:
                    break # stop iteration B, introns will get LONGER from now on                
                if (IntervalB.Flags & GeneFlags.First) and not ReverseFlag:
                    continue
                if (IntervalB.Flags & GeneFlags.Last) and ReverseFlag:
                    continue
                LegalLink = 0
                Complements = [0, 2, 1]
                if ReverseFlag:
                    for X in range(3):
                        if IntervalA.Prefixes[X] and IntervalB.Suffixes[Complements[X]]:
                            LegalLink = 1
                            break
                else:
                    for X in range(3):
                        if IntervalA.Suffixes[X] and IntervalB.Prefixes[Complements[X]]:
                            LegalLink = 1
                            break
                if not LegalLink:
                    continue
                # Given a choice, who should we link to?  We prefer
                # links to high-scoring exons, but also want a reasonable intron length:
                
##                Score = IntervalB.Score
##                if IntronLength < 70:
##                    IntronScore = 365
##                elif IntronLength < 150:
##                    pass # 70-150 is the 'sweet spot'
##                elif IntronLength < 500:
##                    IntronScore = 96
##                elif IntronLength < 1000:
##                    IntronScore = 141
##                elif IntronLength < 5000:
##                    IntronScore = 241
##                else:
##                    IntronScore = 324
##                #Score += IntronScore
##                # NEW APPROACH:
##                # Try the NEAREST TEN exons instead of the "best" ten.
##                IntronScore = -IntronLength
##                Score = IntronScore
                Score = 0
                SummedExonScore = IntervalA.Score + IntervalB.Score
                if SummedExonScore <= 700:
                    Score += -8e-6*SummedExonScore*SummedExonScore + 0.0058 * SummedExonScore - 0.5825
                elif SummedExonScore <= 1300:
                    Score += -4e-6*SummedExonScore*SummedExonScore + 0.0079 * SummedExonScore - 3.4403
                else:
                    Score += 0.0
                IntronScore = 9e-9*IntronLength*IntronLength - 0.0004*IntronLength + 1.3198
                Score += IntronScore
                if VerboseFlag:
                    print "-> SummedExonScore %s IntronScore %s Score %s"%(SummedExonScore, IntronScore, Score)
                if not DesirableLinks.has_key(IntervalB.Start) or DesirableLinks[IntervalB.Start][0] < Score:
                    DesirableLinks[IntervalB.Start] = (Score, IntronScore, IntervalB)
                #IntervalB.AddLinkBack(IntervalA.End)
                #self.LinkDict[(IntervalA.End, IntervalB.Start)] = 1
                #LinkCount += 1
            ##########################################################
            # After iterating over IntervalB, add up to MaxLinks links.
            if VerboseFlag:
                print "DesirableLinks forward from (%s-%s):"%(IntervalA.Start, IntervalA.End)
                for Tuple in DesirableLinks.values():
                    print "  ",Tuple[0], Tuple[1], Tuple[2].Start, Tuple[2].End
            PositionsSeen = {}
            for (Score, IntronScore, IntervalB) in DesirableLinks.values():
                BackScore = IntronScore #%%% IntervalA.Score + IntronScore
                if VerboseFlag:
                    print "Consider %s, %s-%s"%(Score, IntervalB.Start, IntervalB.End)
                    print "Forward edges:", self.GenomicForwardEdges.get(IntervalA.End, [])
                    print "Backward edges:", self.GenomicBackwardEdges.get(IntervalB.Start, [])
                if not PositionsSeen.has_key(IntervalB.Start):
                    Str = "%s\t%s\t%s\t%s\t%s\t"%(Score, IntervalA.Score + IntervalB.Score,
                                              min(IntervalA.Score, IntervalB.Score),
                                              max(IntervalA.Score, IntervalB.Score),
                                              IntervalB.Start - IntervalA.End)
                    AllLinkagesFile.write(Str + "\n")
                    PositionsSeen[IntervalB.Start] = 1
                # If the edge from IntervalA.End to IntervalB.Start is good enough
                # as a forward-edge and as a backward-edge, we'll add it.
                if not self.GenomicForwardEdges.has_key(IntervalA.End):
                    self.GenomicForwardEdges[IntervalA.End] = []
                ForwardList = self.GenomicForwardEdges[IntervalA.End]
                if len(ForwardList) >= 10 and Score <= ForwardList[0][0]:
                    # No link: There are already 10 better links extending forward
                    # from IntervalA.End
                    continue
                # Don't add this link if it's already there:
                FoundAlready = 0
                for X in range(len(ForwardList)):
                    (OldScore, OldPos) = ForwardList[X]
                    if OldPos == IntervalB.Start:
                        FoundAlready = 1
                        ForwardList[X] = (max(OldScore, Score), OldPos)
                        ForwardList.sort()
                        break
                if FoundAlready:
                    continue
                if not self.GenomicBackwardEdges.has_key(IntervalB.Start):
                    self.GenomicBackwardEdges[IntervalB.Start] = []
                BackwardList = self.GenomicBackwardEdges[IntervalB.Start]
                if len(BackwardList) >= 10 and BackScore <= BackwardList[0][0]:
                    # No link: There are already 10 better links extending backward
                    # from IntervalB.Start
                    continue
                # Don't add this link if it's already there:
                FoundAlready = 0
                for X in range(len(BackwardList)):
                    (OldScore, OldPos) = BackwardList[X]
                    if OldPos == IntervalA.End:
                        FoundAlready = 1
                        BackwardList[X] = (max(OldScore, BackScore), OldPos)
                        BackwardList.sort()
                        break
                if FoundAlready:
                    continue
                # It's a keeper!  Displace older links, if necessary.
                ForwardList.append((Score, IntervalB.Start))
                ForwardList.sort()
                if len(ForwardList) > 10:
                    (CrummyScore, CrummyPos) = ForwardList[0]
                    CorrespondingList = self.GenomicBackwardEdges.get(CrummyPos, [])
                    for Tuple in CorrespondingList:
                        if Tuple[1] == IntervalA.End:
                            CorrespondingList.remove(Tuple)
                            break
                self.GenomicForwardEdges[IntervalA.End] = ForwardList[-10:]
                BackwardList.append((BackScore, IntervalA.End))
                BackwardList.sort()
                if len(BackwardList) > 10:
                    (CrummyScore, CrummyPos) = BackwardList[0]
                    CorrespondingList = self.GenomicForwardEdges.get(CrummyPos, [])
                    for Tuple in CorrespondingList:
                        if Tuple[1] == IntervalB.Start:
                            CorrespondingList.remove(Tuple)
                            break
                self.GenomicBackwardEdges[IntervalB.Start] = BackwardList[-10:]
            ################################################
            if IntervalA.End in VerbosePoints:
                print "Chosen edges:"
                ForwardList = self.GenomicForwardEdges.get(IntervalA.End, [])
                for Tuple in ForwardList:
                    print "  ",Tuple[0], Tuple[1]
        # Now use self.GenomicBackwardEdges to set the real links:
        for Interval in self.IntervalDict.values():
            LinkList = self.GenomicBackwardEdges.get(Interval.Start, [])
            for (Score, BackPoint) in LinkList:
                #print "%s-%s links back to %s (score %d)"%(Interval.Start, Interval.End, BackPoint, Score)
                Interval.AddLinkBack(BackPoint)
                self.LinkDict[(BackPoint, Interval.Start)] = 1
                LinkCount += 1
        print "Created a total of %s distinct links (%s joins done) (%s intervals)"%(len(self.LinkDict.keys()), LinkCount, len(Keys))
    def WriteIntervals(self, OutputPath):
        print "Writing intervals out to %s..."%OutputPath
        Keys = self.IntervalDict.keys()
        Keys.sort()
        File = open(OutputPath, "wb")
        IntervalCount = 0
        JunctionCount = 0
        for Key in Keys:
            Interval = self.IntervalDict[Key]
            # Only pass along reading frame flags:
            Flags = Interval.Flags & GeneFlags.AllReadingFrames
            Str = struct.pack("<iiii", Interval.Start, Interval.End, Flags, Interval.Score)
            File.write(Str)
            IntervalCount += 1
            JunctionsBack = self.GenomicBackwardEdges.get(Interval.Start, [])
            File.write(struct.pack("<i", len(JunctionsBack)))
            for (Score, EndPos) in JunctionsBack:
                File.write(struct.pack("<if", EndPos, Score))
                JunctionCount += 1
            #for (EndPos, Score) in Interval.JunctionsBack:
            #    File.write(struct.pack("<if", EndPos, Score))
        File.close()
        print "Intervals read:", IntervalCount
        print "Junctions read:", JunctionCount
        GeneParseSummary = open("GeneParseSummary.txt", "a")
        GeneParseSummary.write("%s\t%s\t%s\t%s\t\n"%(ChromosomeNumber, ReverseFlag, IntervalCount, JunctionCount))
if __name__ == "__main__":
    import psyco
    psyco.full()
    ChromosomeNumber = int(sys.argv[1])
    ReverseFlag = int(sys.argv[2])
    if ReverseFlag:
        ReverseChar = "-"
    else:
        ReverseChar = "+"
    Parser = GeneIDParser()
    if Species == "human":
        GeneIDPath = os.path.join(r"E:\Chromosome\GeneIDOutput", "Exons%d.txt"%ChromosomeNumber)
        #OutputPath = os.path.join("GeneFindDB\\%s%s.dat"%(ChromosomeNumber, ReverseChar))
        #%%%%%%%%%%%%%%%%%%%%
        # Temp, for new intron linkage scheme:
        OutputPath = os.path.join("NGeneFindDB\\%s%s.dat"%(ChromosomeNumber, ReverseChar))
        GenomeDir = "e:\\chromosome"
    else:
        GeneIDPath = os.path.join(r"E:\Chromosome\arabidopsis\GeneIDOutput", "Exons%d.txt"%ChromosomeNumber)
        OutputPath = os.path.join("ATGeneFindDB\\%s%s.dat"%(ChromosomeNumber, ReverseChar))
        GenomeDir = "e:\\chromosome\\arabidopsis"
        
    print "Parse genomic intervals from %s..."%GeneIDPath
    Parser.Parse(GeneIDPath, ReverseChar)
    GenomeFileName = "%s.trie"%ReverseChromosomeMap[ChromosomeNumber]
    GenomeFilePath = os.path.join(GenomeDir, GenomeFileName)
    #print "Score all splice boundaries..."
    #Parser.ScoreSpliceEndpoints(GenomeFilePath)
    print "Link compatible exons with reasonable intron lengths..."
    Parser.LinkExons(ReverseFlag)
    Parser.WriteIntervals(OutputPath)
    