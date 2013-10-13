"""
Given a collection of alignments of ESTs and mRNAs against the human genome, construct
a database of putative exons.
We handle each chromosome separately, and we handle forward sequences separately from reverse-complement
sequences.  For each chromosome/orientation combo, we proceed as follows:
- Read in list of all distinct genomic intervals.  Save the start point, end point, GI-ID# of the first
  sequence witnessing the interval, interval occurrence count, and a list of intervals it links to (plus an occurrence-count
  for these linkages).  If the same
  interval is observed multiple times, keep incrementing its occurrence-count.  Optionally, compute a SCORE.
- Possibly, perform a filter, dropping all intervals with poor score / low occurrence-count

The remainder of the algorithm 'walks along the chromosome', by writing out groups of exons.
- Select the first (in chromosome position) interval not yet flagged as done.
- We want to group exons by gene, since a smaller and simpler exon graph is easier to work with (and we forbid
  any ludicrously long-distance exon bridging).  So, we perform a LIMITED amount of interval assimilation:
  - Take the "next" interval, as sorted by starting point, and put it in the active pool.  (Or, if 'unsatisfied' intervals remain from
    a previous step, throw them into the active pool)
  - The interval's starting point + 100000 serves as an extension limit, k
  Iterative:
   - Add all intervals which (a) intersect with an interval in the pool and (b) end before termination point k
   - If an interval I in the pool links to I2, and I2 ends before k, add I2 to the pool.
   - If an interval I in the pool links to I2, and I2 ends beyond k, add I to a list of 'unsatisfied' intervals not satisfied in this round
 Terminate iteration after the first cycle when nothing's added to the active pool.
- We now have a collection of intervals, which we process to remove overlap.  (Overlap is bad, it adds redundancy to the
  search, which makes running-time longer and exon occurrence count inaccurate).  Overlap removal is iterative.  For each interval:
    Find the first interval which overlaps with this one.  (If none, continue)
    If they share an endpoint, merge their prefix-lists and create two sub-intervals.  If they properly overlap (or one is properly
    contained), construct two child intervals for one.  The links between ADJACENT intervals are special (occurrence-count -1), and
    can be followed 'for free'.  

Now we have a collection of DISJOINT intervals.  We build (up to) three exons for each interval.  However, our troubles aren't over,
because intervals of length 1 or 2 aren't meaningful as exons!  The tricky point is that if your exon has a suffix of length 1, you can't
link to an interval of length 1.  Rather, you must step THROUGH that interval to each of its children, accumulating one base from each:
#1----
      \
#2---->A----->G
      / \
#3----   ---->GA

Exon #1 has suffix length 0.  It links to interval [A] in reading frame 0, with no AA
Exon #2 has suffix length 1.  It links to interval [G] in reading frame 1, with an AA.  And it links to interval [GA] in rf1, with an AA.
  These special 'composite links' receive the worst of the link-weights of the sub-links.
Exon #3 has suffix length 2.  It links to interval [A] in reading frame 0, with an AA.


- All these exons and their links are now written out as one "gene".  Any intervals not on the 'unsatisfied' list are flagged as done.
  Continue the iteration for the next unfinished interval.
"""
import sys
import os
import struct

import traceback

try:
    from PIL import Image
    from PIL import ImageDraw
    from PIL import ImageFont
except:
    print "<< no imaging libraries imported ... plotting code inactive >>"

VerboseFlag = 0

class Colors:
    "Color scheme for coverage plots"
    White = (255,255,255)
    Green = (0,255,0)
    Blue = (0,0,255)
    PaleBlue = (10,10,80)
    Red = (255,0,0)
    Grey = (155,155,155)
    ResidueNumber = (155, 155, 155)
    #Grey = (0,0,0)
    Background = (255, 255, 255)
    # Color schemes for different occurrence counts:
    Occurrences = [(255,0,0), (155, 155, 155), (145, 145, 145),
                   (135, 135, 135), (125, 125, 125), (115, 115, 115),
                   (105, 105, 105), (95, 95, 95), (85, 85, 85),
                   (75, 75, 75), (65, 65, 65), (55, 55, 55),
                   (45, 45, 45), (35, 35, 35), (25, 25, 25),
                   (15, 15, 15), (5, 5, 5), (0,0,0)]


# Map of chromosome names to numbers.  X is called 23, and Y is called 24.
ChromosomeMap = {"chr1":1, "chr2":2, "chr3":3, "chr4":4,
                 "chr5":5, "chr6":6, "chr7":7, "chr8":8,
                 "chr9":9, "chr10":10, "chr11":11, "chr12":12,
                 "chr13":13, "chr14":14, "chr15":15, "chr16":16,
                 "chr17":17, "chr18":18, "chr19":19, "chr20":20,
                 "chr21":21, "chr22":22, "chrX":23, "chrY":24,
                 }


SequenceLengthHistogram = {}
EdgeCountHistogram = {}
GIIDDict = None

def PruneExon(Sequence):
    """
    If an exon contains more than one stop codon, then drop everything between the first and last...but,
    don't drop an ORF of size >100
    """
    OldX = Sequence.find("X")
    while (1):
        NextX = Sequence.find("X", OldX + 1)
        if NextX == -1:
            break
        if NextX > OldX + 100:
            # A long open reading frame - keep it.  (Assuming uniform, independent distribution of
            # nucleotides, the odds of an ORF of length 100 is 0.82%)
            OldX = NextX
        else:
            Sequence = Sequence[:OldX] + Sequence[NextX:]
    return Sequence

class ExonClass:
    "Class representing an interval in a reading frame."
    def __init__(self, Start, End):
        self.Start = Start
        self.End = End
        self.Prefix = ""
        self.Suffix = ""
        self.Sequence = ""
        self.Shift = None
        self.ForwardLinks = [] # Edges are tuples of the form (AA, OtherExon, Flags)
        self.ForwardLinkPower = []
        self.ForwardLinkAA = []
        self.BackwardLinks = [] # Edges are tuples of the form (AA, OtherExon, Flags)
        self.BackwardLinkPower = []
        self.BackwardLinkAA = []
        self.Index = None
    def __cmp__(self, Other):
        if self.Start < Other.Start:
            return -1
        if self.Start > Other.Start:
            return 1
        if self.Shift < Other.Shift:
            return -1
        if self.Shift > Other.Shift:
            return 1
        return 1

# ChromosomePaths is a list, indexed by chrom-number, of the sequence-file paths.
ChromosomeDir = r"C:\source\Bafna\Splice\chromFa"
ChromosomePaths = [None]
for X in range(1, 23):
    ChromosomePaths.append("%s\\Chr%d.trie"%(ChromosomeDir, X))
ChromosomePaths.append("%s\\ChrX.trie"%ChromosomeDir)
ChromosomePaths.append("%s\\ChrY.trie"%ChromosomeDir)

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
RCDict = {"A":"T", "G":"C", "T":"A", "C":"G"}

def ReverseComplement(DNA):
    Str = ""
    for Index in range(len(DNA) - 1, -1, -1):
        Str += RCDict.get(DNA[Index], DNA[Index])
    return Str

def Translate(DNA):
    Peptide = ""
    for Index in range(0, len(DNA) - 2, 3):
        Codon = DNA[Index:Index+3]
        AA = GeneticCode.get(Codon, "")
        Peptide += AA
        # Include a character for stop codons, and then end it:
        #if AA == "X":
        #    break
    return Peptide

# Use the "no links out" and "no links in" flags to indicate that an interval was
# not matched in its entirety, but is a sub-interval of some larger match.  That
# means there's a boundary with no evidence it's a splice entry or (exit) site.
# Start with two intervals which overlap:
#           |--------------------| A
#                       |------------------| B

#           |----------||--------||--------| A B C
# A gets flag NoLinksOut (since we've never seen an intron start at the end of A), and
# C gets NoLinksIn (since we've never seen an intron end at the start of C). 
class IntervalClass:
    "A genomic interval, covered by one or more ESTs."
    FlagProcessed = 1
    FlagNoLinksOut = 2
    FlagNoLinksIn = 3
    def __init__(self, Start, End):
        self.Start = Start
        self.End = End
        self.Flags = 0
        self.Score = 0
        self.Occurrences = 0
        # The entries in ForwardLinks are other intervals.  The corresponding entries in ForwardLinksCounts
        # are positive integers (for the number of EST alignments supporting that splice junction), or -1 if the
        # interval is adjacent.
        self.ForwardLinks = []
        self.ForwardLinkAA = []
        self.ForwardLinkCounts = []
        self.BackwardLinks = []
        self.BackwardLinkAA = []
        self.BackwardLinkCounts = []
        # Most intervals have three exons, for three reading frames.  Exception: Length-1 intervals have reading frames 0 and 1 only.
        self.Exons = []
        self.IDList = []
    def __cmp__(self, Other):
        "For calling sort() on an interval list.  Sort by start position first, then end position."
        if self.Start < Other.Start:
            return -1
        if self.Start > Other.Start:
            return 1
        if self.End < Other.End:
            return -1
        if self.End > Other.End:
            return 1
        if id(self) < id(Other):
            return -1
        elif id(self) > id(Other):
            return 1
        return 0 

class ExonCollector:
    "Master object for constructing a splice-tolerant database"
    def __init__(self, Chromosome, ReverseFlag):
        self.Intervals = [] # sorted list
        self.IntervalDict = {} # start, end -> interval
        self.Chromosome = Chromosome
        self.ReverseFlag = ReverseFlag
        self.RecordNumber = 0
        print "Reverse flag:", self.ReverseFlag
    def GatherExons(self, FileName):
        """
        Parse intervals from an alignment file in sim4 format.  Accumulate the list self.Intervals and the
        dictionary IntervalDict.  Discard any intervals that don't come from self.Chromosome in the right orientation (self.ReverseFlag)
        """
        global GIIDDict
        File = open(FileName, "rb")
        State = 0 # skipping current record
        RecordLineNumber = 0
        LineNumber = 0
        PreviousInterval = None
        while (1):
            FilePos = File.tell()
            FileLine = File.readline()
            if not FileLine:
                break
            FileLine = FileLine.strip()
            LineNumber += 1
            if LineNumber%1000 == 0:
                print "Line %s..."%LineNumber
            if FileLine == "sim4begin":
                RecordLineNumber = 0
                PreviousInterval = None
                RecordStartFilePos = FilePos
                State = 1 # innocent, until proven guilty :)
                continue
            if FileLine == "sim4end":
                continue
            if State == 0:
                continue
            RecordLineNumber += 1
            #print State, RecordLineNumber, LineNumber, FileLine[:60].strip()
            if RecordLineNumber == 1: # Line 1: Chromosome position and orientation
                # Line of the form "2[339-0-0] 0[227615646-227622915] <279-0-96-forward-unknown>"
                Bits = FileLine.split("-")
                if len(Bits)<7:
                    # a bogus-looking line.
                    State = 0
                    continue
                # Reject, if the direction is wrong:
                Direction = Bits[6]
                if self.ReverseFlag!=None:
                    if Bits[6] == "forward":
                        if self.ReverseFlag:
                            State = 0
                    elif Bits[6] == "complement":
                        if not self.ReverseFlag:
                            State = 0
                Bracket = FileLine.find("[")
                Bracket = FileLine.find("[", Bracket + 1)
                EndBracket = FileLine.find("]", Bracket)
                RangeBits = FileLine[Bracket+1:EndBracket].split("-")
                OverallStart = int(RangeBits[0])
                OverallEnd = int(RangeBits[1])
            elif RecordLineNumber == 2: #edef, the GI number
                GINumber = int(FileLine.split("|")[1])
                #print "GINumber %s corresponds to file position %s"%(GINumber, RecordStartFilePos)
                if GIIDDict!=None:
                    if GIIDDict.has_key(GINumber):
                        OldEntry = GIIDDict[GINumber]
                        if type(OldEntry)==type([]):
                            GIIDDict[GINumber].append(RecordStartFilePos)
                        else:
                            GIIDDict[GINumber] = [OldEntry, RecordStartFilePos]
                    else:
                        GIIDDict[GINumber] = RecordStartFilePos
                #print "  GI number parsed:", GINumber
            elif RecordLineNumber == 3: #ddef, the chromosome number
                Chrom = FileLine.split()[0]
                Chrom = Chrom[6:]
                Number = ChromosomeMap.get(Chrom, 0)
                if self.Chromosome!=None and Number != self.Chromosome:
                    State = 0
                #print "  Chromosome %s"%Chrom
            else: # an exon line
                # Line format: 377-622 (117076-117321) <246-0-100> ->
                Paren = FileLine.find("(")
                Paren2 = FileLine.find(")")
                PosBits = FileLine[Paren+1:Paren2].split("-")
                Start = OverallStart + int(PosBits[0]) - 1
                End = OverallStart + int(PosBits[1])
                Interval = self.IntervalDict.get((Start, End))
                #print "  Exon from %s to %s (local: %s to %s)"%(Start, End, PosBits[0], PosBits[1])
                if Interval:
                    Interval.Occurrences += 1
                    if GINumber not in Interval.IDList:
                        Interval.IDList.append(GINumber)
                else:
                    Interval = IntervalClass(Start, End)
                    Interval.Occurrences = 1
                    Interval.IDList = [GINumber]
                    self.IntervalDict[(Start, End)] = Interval
                    self.Intervals.append(Interval)
                # Link from PreviousInterval to Interval:
                if PreviousInterval:
                    try:
                        Index = PreviousInterval.ForwardLinks.index(Interval)
                        PreviousInterval.ForwardLinkCounts[Index] += 1
                    except:
                        PreviousInterval.ForwardLinks.append(Interval)
                        PreviousInterval.ForwardLinkCounts.append(1)
                    try:
                        Index = Interval.BackwardLinks.index(PreviousInterval)
                        Interval.BackwardLinkCounts[Index] += 1
                    except:
                        Interval.BackwardLinks.append(PreviousInterval)
                        Interval.BackwardLinkCounts.append(1)
                PreviousInterval = Interval
    def WriteGenes(self, OutputFileName):
        """
        'Walk' across the genome, processing one ActivePool of intervals at once.
        Intersect these intervals to form a disjoint collection, then form exons for the
        intervals and write the exons out.
        """
        print "Total of %d intervals found."%len(self.Intervals)
        self.OutputFile = open(OutputFileName, "wb")
        ActivePool = []
        FreshMeat = []
        OldStragglers = []
        Stragglers = []
        while (1):
            print "%s intervals remain to be processed as genes."%len(self.Intervals)
            if len(self.Intervals) == 0:
                break 
            self.Intervals.sort()
            OldStragglers = Stragglers
            if OldStragglers:
                # There are some intervals which we didn't finish in a prior run.  We added them, but not
                # their forward links.  This time we'll add them and their forward links, but not their backward links.
                ActivePool = OldStragglers[:]
                FreshMeat = OldStragglers[:]
                MaxPos = 0
                MinPos = 999999999999
                if VerboseFlag:
                    print "Handling OLDSTRAGGLERS:"
                Limit = MinPos + 100000
                for Interval in ActivePool:
                    MaxPos = max(Interval.End, MaxPos)
                    MinPos = min(Interval.Start, MinPos)
                    if VerboseFlag:
                        print "Interval %s-%s with edges:"%(Interval.Start, Interval.End)
                    for Next in Interval.ForwardLinks:
                        if VerboseFlag:
                            print " to %s"%Next.Start
                        Limit = max(Limit, Next.End)
            else:
                # Accept the next interval, together with everything it overlaps and joins to...up to a point
                Start = self.Intervals[0].Start
                MaxPos = self.Intervals[0].End
                Limit = Start + 100000 # maximum span for any exon group
                FreshMeat = [self.Intervals[0]]
                ActivePool = [self.Intervals[0]]
                del self.Intervals[0]
            Stragglers = []
            # Assimilate all the intervals linked to an interval in ActivePool, stopping the cascade if it
            # walks too far along the chromosome
            while (1):
                print "FM %d, ap %d"%(len(FreshMeat), len(ActivePool))
                if len(FreshMeat)<1:
                    break # We're all done assimilating more intervals
                OldMeat = FreshMeat
                FreshMeat = []
                # Assimilate any intervals linked to by the current active set:
                for Interval in OldMeat:
                    for LinkedInterval in Interval.ForwardLinks:
                        if LinkedInterval.End > Limit:
                            Stragglers.append(Interval)
                        else:
                            if LinkedInterval not in ActivePool:
                                FreshMeat.append(LinkedInterval)
                                ActivePool.append(LinkedInterval)
                                try:
                                    self.Intervals.remove(LinkedInterval)
                                except:
                                    pass
                                MaxPos = max(LinkedInterval.End, MaxPos)
                    if Interval not in OldStragglers:
                        for LinkedInterval in Interval.BackwardLinks:
                            if LinkedInterval.End > Limit:
                                Stragglers.append(Interval)
                            else:
                                if LinkedInterval not in ActivePool:
                                    FreshMeat.append(LinkedInterval)
                                    ActivePool.append(LinkedInterval)
                                    try:
                                        self.Intervals.remove(LinkedInterval)
                                    except:
                                        pass
                                    MaxPos = max(LinkedInterval.End, MaxPos)
                # Assimilate any intervals that are covered by the intervals so far:
                MaxIndex = len(self.Intervals)
                Index = 0
                while Index < MaxIndex:
                    Interval = self.Intervals[Index]
                    if Interval.Start < MaxPos and Interval.End < Limit:
                        #self.Intervals.remove(Interval)
                        FreshMeat.append(Interval)
                        ActivePool.append(Interval)
                        MaxPos = max(MaxPos, Interval.End)
                        del self.Intervals[Index]
                        MaxIndex -= 1
                        continue
                    if Interval.Start >= MaxPos:
                        break
                    Index += 1
            print "Active pool of %s intervals."%len(ActivePool)
            #self.DisplayIntervalGraph(ActivePool, "Before")
            # We've assimilated a collection of intervals into ActivePool.  This collection of intervals will
            # be intersected and sliced and diced to form a collection of exons, at which point we can write out
            # one gene.  Do so now:
            self.ConsolidateIntervalsIntoGene(ActivePool)
    def DisplayIntervalGraph(self, Intervals, NameExtra):
        """
        Construct a plot showing how these intervals cover the genome.
        """
        EndPos = 0
        StartPos = Intervals[0].Start
        for Interval in Intervals:
            EndPos = max(Interval.End, EndPos)
            StartPos = min(Interval.Start, StartPos)
        print "DIGraph!  %s-%s"%(StartPos, EndPos)
        Header = 20
        SideMargin = 50
        Height = 200
        Width = (EndPos - StartPos) + SideMargin
        # Show the 'timeline':
        CoverageImage = Image.new("RGB", (Width, Height), Colors.Background)  # mode, size, [startcolor]
        Draw = ImageDraw.Draw(CoverageImage)
        Draw.line((0, Header, Width - SideMargin, Header), Colors.ResidueNumber)
        Pos = StartPos - (StartPos%100) + 100
        print Pos
        while Pos < EndPos:
            X = Pos-StartPos
            Draw.text((X, 2), "%s"%Pos, Colors.ResidueNumber)
            Draw.line((X, Header, X, Header - 2), Colors.ResidueNumber)
            Pos += 100
        # Plot the intervals:
        for IntervalIndex in range(len(Intervals)):
            Interval = Intervals[IntervalIndex]
            Y = Header + 2 + 2*(IntervalIndex%50)
            X = Interval.Start - StartPos
            X2 = Interval.End - StartPos
            ColorIndex = min(len(Colors.Occurrences) - 1, Interval.Occurrences)
            Color = Colors.Occurrences[ColorIndex]
            Draw.line((X, Y, X2, Y), Color)
            print "%s (%s) to %s (%s)"%(Interval.Start, X, Interval.End, X2)
        OutputFileName = "SpliceGraph\\%s-%s.%s.png"%(StartPos, EndPos, NameExtra)
        CoverageImage.save(OutputFileName, "png")
        print "Wrote %s"%OutputFileName
    def RemoveLocalInterval(self, Interval):
        try:
            del self.LocalIntervalDict[(Interval.Start, Interval.End)]
        except:
            pass
        try:
            self.LocalIntervals.remove(Interval)
        except:
            pass
    def AddLocalInterval(self, Interval):
        if Interval in self.LocalIntervals:
            return
        OldInterval = self.LocalIntervalDict.get((Interval.Start, Interval.End), None)
        if OldInterval:
            self.AssimilateLinks(Interval, OldInterval, 1)
            self.AssimilateLinks(Interval, OldInterval, 0)
            for ID in Interval.IDList:
                if ID not in OldInterval.IDList:
                    OldInterval.IDList.append(ID)
            OldInterval.Occurrences += Interval.Occurrences
            return
        self.LocalIntervals.append(Interval)
        self.LocalIntervalDict[(Interval.Start, Interval.End)] = Interval
    def ConsolidateIntervalsIntoGene(self, Intervals):
        Intervals.sort()
        self.LocalIntervals = Intervals
        print 
        print "="*79
        print "Formed an active pool of %s intervals that stretch from %s to %s"%(len(Intervals), Intervals[0].Start, Intervals[-1].End)
        # Create a LOCAL interval dictionary:
        
        self.LocalIntervalDict = {}
        for Interval in self.LocalIntervals:
            self.LocalIntervalDict[(Interval.Start, Interval.End)] = Interval
        for Interval in self.LocalIntervals:
            print Interval.Start, Interval.End
##        if VerboseFlag:
##            for Interval in self.LocalIntervals:
##                print Interval, Interval.Start, Interval.End
##                for Linky in Interval.ForwardLinks:
##                    print "  ", Linky, Linky.Start, Linky.End
        while (1):
            if VerboseFlag:
                print 
                print "-=-"*10
                print "Process gene #%s"%self.RecordNumber
                print "Interval count: %s"%len(self.LocalIntervals)
            ## HYPERVERBOSE:
##            if self.RecordNumber in (37, 732, 733):
##                for Interval in self.LocalIntervals:
##                    print Interval.Start, Interval.End
##                    for ForwardInterval in Interval.ForwardLinks:
##                        print " -> %s-%s"%(ForwardInterval.Start, ForwardInterval.End)
##                    for BackInterval in Interval.BackwardLinks:
##                        print " %s-%s <-"%(BackInterval.Start, BackInterval.End)
                        
            # First, perform intersections until no more are necessary:
            IntervalA = None
            IntervalB = None
            MergeA = None
            MergeB = None
            for IndexA in range(len(self.LocalIntervals) - 1):
                TryA = self.LocalIntervals[IndexA]
##                TryB = Intervals[IndexA + 1]
##                if TryA.Start == TryB.Start and TryA.End == TryB.End:
##                    IntervalA = TryA
##                    IntervalB = TryB
##                    break
                for IndexB in range(IndexA + 1, len(self.LocalIntervals)):
                    TryB = self.LocalIntervals[IndexB]
                    if TryB.Start < TryA.End:
                        if TryA.Start == TryB.Start and TryA.End == TryB.End:
                            MergeA = TryA
                            MergeB = TryB
                            break
                        if not IntervalA:
                            IntervalA = TryA
                            IntervalB = TryB
                        
                if MergeA:
                    break
            if MergeA:
                IntervalA = MergeA
                IntervalB = MergeB
            if not IntervalA:
                break # we're done carrying out intersections!
            print "Interval overlap: %s-%s and %s-%s"%(IntervalA.Start, IntervalA.End, IntervalB.Start, IntervalB.End)
            print "A has %s links back, %s forward"%(len(IntervalA.BackwardLinks), len(IntervalA.ForwardLinks))
            if IntervalA.BackwardLinks:
                print " Link back to %s-%s"%(IntervalA.BackwardLinks[0].Start,IntervalA.BackwardLinks[0].End)
            if IntervalA.ForwardLinks:
                print " Link fwrd to %s-%s"%(IntervalA.ForwardLinks[0].Start,IntervalA.ForwardLinks[0].End)                
            print "B has %s links back, %s forward"%(len(IntervalB.BackwardLinks), len(IntervalB.ForwardLinks))
            if IntervalB.BackwardLinks:
                print " Link back to %s-%s"%(IntervalB.BackwardLinks[0].Start,IntervalB.BackwardLinks[0].End)
            if IntervalB.ForwardLinks:
                print " Link back to %s-%s"%(IntervalB.ForwardLinks[0].Start,IntervalB.ForwardLinks[0].End)                
            #print " A %s, B %s"%(IntervalA, IntervalB)
            print " size %s, size %s, Overlap size %s"%(IntervalA.End - IntervalA.Start,
                                                                    IntervalB.End - IntervalB.Start,
                                                                    min(IntervalA.End, IntervalB.End) - max(IntervalA.Start, IntervalB.Start))
##            
##            print " gi %s size %s, gi %s size %s, Overlap size %s"%(IntervalA.ID, IntervalA.End - IntervalA.Start,
##                                                                    IntervalB.ID, IntervalB.End - IntervalB.Start,
##                                                                    min(IntervalA.End, IntervalB.End) - max(IntervalA.Start, IntervalB.Start))
##            # Ok, time to perform an intersection:
            
            if IntervalA.Start == IntervalB.Start and IntervalA.End == IntervalB.End:
                # Sometimes, after processing other overlaps, we end up with two identical intervals.
                # Trivial merge ensues:
                # |---------------| A
                # |---------------| B
                self.AssimilateLinks(IntervalA, IntervalB, 1)
                self.AssimilateLinks(IntervalA, IntervalB, 0)
                for ID in IntervalA.IDList:
                    if ID not in IntervalB.IDList:
                        IntervalB.IDList.append(ID)
                IntervalB.Occurrences += IntervalA.Occurrences
                self.LocalIntervals.remove(IntervalA)
                continue
            if IntervalA.Start == IntervalB.Start:
                # |---------------| A
                # |------------------------| B

                # |---------------| A
                #                 |--------| C
                #            or
                # |------------------------| B
                if IntervalA.End > IntervalB.End:
                    Swap = IntervalA
                    IntervalA = IntervalB
                    IntervalB = Swap
                if not IntervalA.ForwardLinkCounts:
                    # We don't gain anything by keeping A around, because the only reasonable way forward
                    # is to keep moving through B.  So, let's get rid of A!
                    # NOTE: It's not totally accurate to increment B's occurrence count by A's
                    # occurrence count.  However, keeping track of the level of coverage
                    # across B would be an extra layer of bookkeeping for little benefit.
                    IntervalB.Occurrences += IntervalA.Occurrences
                    self.AssimilateLinks(IntervalA, IntervalB, 0)
                    self.RemoveLocalInterval(IntervalA)
                    continue # no need to re-sort
                OldCountA = IntervalA.Occurrences
                IntervalA.Occurrences += IntervalB.Occurrences
                self.AssimilateLinks(IntervalB, IntervalA, 0)
                IntervalC = IntervalClass(IntervalA.End, IntervalB.End)
                IntervalC.IDList = IntervalB.IDList
                #IntervalC.ForwardLinks = IntervalB.ForwardLinks
                #IntervalC.ForwardLinkCounts = IntervalB.ForwardLinkCounts
                self.AssimilateLinks(IntervalB, IntervalC, 1)
                IntervalC.BackwardLinks = [IntervalA]
                IntervalC.BackwardLinkCounts = [IntervalB.Occurrences] #[-1] # negative count means "I come from an intersection"
                IntervalC.Occurrences = IntervalB.Occurrences
                IntervalA.ForwardLinks.append(IntervalC)
                IntervalA.ForwardLinkCounts.append(IntervalB.Occurrences)
                self.RemoveLocalInterval(IntervalB)
                self.AddLocalInterval(IntervalC)
                #Intervals.remove(IntervalB)
                #Intervals.append(IntervalC)
                self.LocalIntervals.sort()
                continue
            if IntervalA.End == IntervalB.End:
                #          |---------------| A
                # |------------------------| B

                #          |---------------| A
                # |--------| C
                if IntervalA.Start < IntervalB.Start:
                    Swap = IntervalA
                    IntervalA = IntervalB
                    IntervalB = Swap
                if not IntervalA.BackwardLinkCounts:
                    # Assimilate A into B:
                    IntervalB.Occurrences += IntervalA.Occurrences
                    self.AssimilateLinks(IntervalA, IntervalB, 1)
                    self.RemoveLocalInterval(IntervalA)
                    continue # no need to re-sort
                OldCountA = IntervalA.Occurrences
                IntervalA.Occurrences += IntervalB.Occurrences
                self.AssimilateLinks(IntervalB, IntervalA, 1)
                IntervalC = IntervalClass(IntervalB.Start, IntervalA.Start)
                IntervalC.IDList = IntervalB.IDList
                #IntervalC.BackwardLinks = IntervalB.BackwardLinks
                #IntervalC.BackwardLinkCounts = IntervalB.BackwardLinkCounts
                self.AssimilateLinks(IntervalB, IntervalC, 0)
                IntervalC.ForwardLinks = [IntervalA]
                IntervalC.ForwardLinkCounts = [IntervalB.Occurrences] #[-1] # negative count means "I come from an intersection"
                IntervalC.Occurrences = IntervalB.Occurrences
                IntervalA.BackwardLinks.append(IntervalC)
                IntervalA.BackwardLinkCounts.append(IntervalB.Occurrences)
                self.RemoveLocalInterval(IntervalB)
                self.AddLocalInterval(IntervalC)
                #Intervals.remove(IntervalB)
                #Intervals.append(IntervalC)
                self.LocalIntervals.sort()
                continue
            if IntervalB.Start < IntervalA.Start and IntervalB.End > IntervalA.End:
                # |------------------------| B
                #        |---------|         A
                Swap = IntervalA
                IntervalA = IntervalB
                IntervalB = Swap
            if IntervalA.Start < IntervalB.Start and IntervalA.End > IntervalB.End:
                # |------------------------| A
                #        |---------|         B

                if not IntervalB.ForwardLinks and not IntervalB.BackwardLinks:
                    IntervalA.Occurrences += IntervalB.Occurrences
                    self.RemoveLocalInterval(IntervalB)
                    continue
                # |----------------|         A2
                #                  |-------| B2
                if not IntervalB.BackwardLinks:
                    IntervalA2 = IntervalClass(IntervalA.Start, IntervalB.End)
                    IntervalB2 = IntervalClass(IntervalB.End, IntervalA.End)
                    IntervalA2.ForwardLinks = [IntervalB2]
                    IntervalA2.ForwardLinkCounts = [IntervalA.Occurrences]
                    IntervalB2.BackwardLinks = [IntervalA2]
                    IntervalB2.BackwardLinkCounts = [IntervalA.Occurrences]
                    self.AssimilateLinks(IntervalB, IntervalA2, 1)
                    self.AssimilateLinks(IntervalA, IntervalB2, 1)
                    self.AssimilateLinks(IntervalA, IntervalA2, 0)
                    self.RemoveLocalInterval(IntervalA)
                    self.RemoveLocalInterval(IntervalB)
                    self.AddLocalInterval(IntervalA2)
                    self.AddLocalInterval(IntervalB2)
                    self.LocalIntervals.sort()
                    print "Added subintervals from %s-%s, from %s-%s"%(IntervalA2.Start, IntervalA2.End,
                                                                       IntervalB2.Start, IntervalB2.End)                    
                    continue
                #         |----------------| A2
                # |-------|                  B2
                if not IntervalB.ForwardLinks:
                    IntervalA2 = IntervalClass(IntervalB.Start, IntervalA.End)
                    IntervalB2 = IntervalClass(IntervalA.Start, IntervalB.Start)
                    IntervalB2.ForwardLinks = [IntervalA2]
                    IntervalB2.ForwardLinkCounts = [IntervalA.Occurrences]
                    IntervalA2.BackwardLinks = [IntervalB2]
                    IntervalA2.BackwardLinkCounts = [IntervalA.Occurrences]
                    self.AssimilateLinks(IntervalB, IntervalA2, 0)
                    self.AssimilateLinks(IntervalA, IntervalA2, 1)
                    self.AssimilateLinks(IntervalA, IntervalB2, 0)
                    self.RemoveLocalInterval(IntervalA)
                    self.RemoveLocalInterval(IntervalB)
                    self.AddLocalInterval(IntervalA2)
                    self.AddLocalInterval(IntervalB2)
                    self.LocalIntervals.sort()
                    print "Added subintervals from %s-%s, from %s-%s"%(IntervalB2.Start, IntervalB2.End,
                                                                       IntervalA2.Start, IntervalA2.End)
                    continue
                # |------|    B    |-------| 
                #    C   |---------|   D
                #
                IntervalC = IntervalClass(IntervalA.Start, IntervalB.Start)
                IntervalC.Occurrences = IntervalA.Occurrences
                self.AssimilateLinks(IntervalA, IntervalC, 0)
                #IntervalC.BackwardLinks = IntervalA.BackwardLinks
                #IntervalC.BackwardLinkCounts = IntervalA.BackwardLinkCounts
                IntervalC.ForwardLinks = [IntervalB]
                IntervalC.ForwardLinkCounts = [IntervalA.Occurrences]
                IntervalC.IDList = IntervalA.IDList
                #
                IntervalD = IntervalClass(IntervalB.End, IntervalA.End)
                IntervalD.Occurrences = IntervalA.Occurrences
                self.AssimilateLinks(IntervalA, IntervalD, 1)
                #IntervalD.ForwardLinks = IntervalA.ForwardLinks
                #IntervalD.ForwardLinkCounts = IntervalA.ForwardLinkCounts
                IntervalD.BackwardLinks = [IntervalB]
                IntervalD.BackwardLinkCounts = [IntervalA.Occurrences] # [-1]
                IntervalD.IDList = IntervalA.IDList
                #
                IntervalB.BackwardLinks.append(IntervalC)
                IntervalB.BackwardLinkCounts.append(IntervalA.Occurrences) #(-1)
                IntervalB.ForwardLinks.append(IntervalD)
                IntervalB.ForwardLinkCounts.append(IntervalA.Occurrences)
                IntervalB.Occurrences += IntervalA.Occurrences
                #Intervals.remove(IntervalA)
                #Intervals.append(IntervalC)
                #Intervals.append(IntervalD)
                self.RemoveLocalInterval(IntervalA)
                self.AddLocalInterval(IntervalC)
                self.AddLocalInterval(IntervalD)
                self.LocalIntervals.sort()
                continue
            # And partial overlap case, too:
            if IntervalA.Start > IntervalB.Start and IntervalA.End > IntervalB.End:
                # |-------------| B
                #        |--------------| A
                Swap = IntervalA
                IntervalA = IntervalB
                IntervalB = Swap
                
            if IntervalA.Start < IntervalB.Start and IntervalA.End < IntervalB.End:
                # |-------------| A
                #        |--------------| B

                # |---------------------| A
                if not IntervalB.BackwardLinks and not IntervalA.ForwardLinks:
                    IntervalA.End = IntervalB.End
                    self.AssimilateLinks(IntervalB, IntervalA, 1)
                    self.RemoveLocalInterval(IntervalB)
                    continue
                
                # |------|  B2  |-------| 
                #    C   |------|  D
                IntervalB2 = IntervalClass(IntervalB.Start, IntervalA.End)
                IntervalB2.Occurrences = IntervalA.Occurrences + IntervalB.Occurrences
                self.AssimilateLinks(IntervalB, IntervalB2, 0)
                #IntervalB2.BackwardLinks = IntervalB.BackwardLinks
                #IntervalB2.BackwardLinkCounts = IntervalB.BackwardLinkCounts
                self.AssimilateLinks(IntervalA, IntervalB2, 1)
                #IntervalB2.ForwardLinks = IntervalA.ForwardLinks
                #IntervalB2.ForwardLinkCounts = IntervalA.ForwardLinkCounts
                IntervalB2.IDList = IntervalB.IDList
                #
                IntervalC = IntervalClass(IntervalA.Start, IntervalB.Start)
                IntervalC.Occurrences = IntervalA.Occurrences
                self.AssimilateLinks(IntervalA, IntervalC, 0)
                #IntervalC.BackwardLinks = IntervalA.BackwardLinks
                #IntervalC.BackwardLinkCounts = IntervalA.BackwardLinkCounts
                IntervalC.ForwardLinks = [IntervalB2]
                IntervalC.ForwardLinkCounts = [IntervalA.Occurrences]
                IntervalC.IDList = IntervalA.IDList
                IntervalB2.BackwardLinks.append(IntervalC)
                IntervalB2.BackwardLinkCounts.append(IntervalA.Occurrences)
                #
                IntervalD = IntervalClass(IntervalA.End, IntervalB.End)
                IntervalD.Occurrences = IntervalB.Occurrences
                #IntervalD.ForwardLinks = IntervalB.ForwardLinks
                #IntervalD.ForwardLinkCounts = IntervalB.ForwardLinkCounts
                self.AssimilateLinks(IntervalB, IntervalD, 1)
                IntervalD.BackwardLinks = [IntervalB2]
                IntervalD.BackwardLinkCounts = [IntervalB.Occurrences]
                IntervalD.IDList = IntervalB.IDList
                IntervalB2.ForwardLinks.append(IntervalD)
                IntervalB2.ForwardLinkCounts.append(IntervalB.Occurrences)
                
                #
                self.RemoveLocalInterval(IntervalA)
                self.RemoveLocalInterval(IntervalB)
                self.AddLocalInterval(IntervalC)
                self.AddLocalInterval(IntervalB2)
                self.AddLocalInterval(IntervalD)
##                Intervals.remove(IntervalA)
##                Intervals.remove(IntervalB)
##                Intervals.append(IntervalC)
##                Intervals.append(IntervalB2)
##                Intervals.append(IntervalD)
                self.LocalIntervals.sort()
                continue
        # We have now performed enough intersections.  PERHAPS we've performed too many.
        # If an interval is adjacent to another, and has a link into it, and has no other links,
        # then we can merge it:
        self.MergeRedundantLocalIntervals()
        if VerboseFlag:            
            print "Interval intersections are all complete."
            for Interval in self.LocalIntervals:
                print Interval, Interval.Start, Interval.End
                for Linky in Interval.ForwardLinks:
                    print "  ", Linky, Linky.Start, Linky.End
            print "<<end>>"
        #self.DisplayIntervalGraph(self.LocalIntervals, "After")
        # Read the DNA for each interval:
        self.LocalIntervals.sort()
        if self.ReverseFlag:
            self.LocalIntervals.reverse()
        Exons = self.GenerateExons(self.LocalIntervals)
        # Now we can write out this 'gene' to the file
        if self.ReverseFlag:
            Direction = "-"
        else:
            Direction = "+"
        Name = "%s%s Gene %d of %d exons, %d intervals from %d to %d"%(self.Chromosome, Direction, self.RecordNumber, len(Exons), len(Intervals), Intervals[0].Start, Intervals[-1].End)
        self.RecordNumber += 1
        Str = struct.pack("<256s", Name)
        self.OutputFile.write(Str)
        Str = struct.pack("<256s", Name)
        self.OutputFile.write(Str)
        if VerboseFlag:
            print "Chromosome %s exon count %s"%(self.Chromosome, len(Exons))
        IDList = []
        for Interval in self.LocalIntervals:
            for ID in Interval.IDList:
                if ID not in IDList:
                    IDList.append(ID)
        Str = struct.pack("<ii", self.Chromosome, len(Exons))
        self.OutputFile.write(Str)
        # Write (up to) 10 GI IDs for the GLOB.  Note that if the GLOB spans multiple genes, then we may have ended
        # up with ESTs from an adjacent gene.  
        for X in range(10):
            if X >= len(IDList):
                Str = struct.pack("<i", -1)
            else:
                Str = struct.pack("<i", IDList[X])
            self.OutputFile.write(Str)
        for Interval in self.LocalIntervals:
            for Exon in Interval.Exons:
                self.WriteExon(self.OutputFile, Interval, Exon)
    def MergeRedundantLocalIntervals(self):
        """
        Merge any adjacent intervals with no other informative edges.
        |------------| A
                     |------------| B
        """
        MaxIndex = len(self.LocalIntervals) - 1
        Index = 0
        while Index < MaxIndex:
            IntervalA = self.LocalIntervals[Index]
            IntervalB = self.LocalIntervals[Index+1]
            if IntervalA.End != IntervalB.Start:
                Index += 1
                continue
            if len(IntervalA.ForwardLinks)>1 or IntervalA.ForwardLinks[0] != IntervalB:
                Index += 1
                continue
            if len(IntervalB.BackwardLinks)>1 or IntervalB.BackwardLinks[0] != IntervalA:
                Index += 1
                continue
            # A eats B:
            IntervalA.End = IntervalB.End
            self.AssimilateLinks(IntervalB, IntervalA, 1)
            self.RemoveLocalInterval(IntervalB)
            MaxIndex -= 1
            
    def WriteExon(self, File, Interval, Exon):
        global SequenceLengthHistogram
        global EdgeCountHistogram
        if VerboseFlag:
            print 
            print "%s Exon %s-%s occurrences %s sequencelen %s"%(Exon.Index, Interval.Start, Interval.End,  Interval.Occurrences, len(Exon.Sequence))
        Len = len(Exon.Sequence)
        SequenceLengthHistogram[Len] = SequenceLengthHistogram.get(Len, 0) + 1
        Count = len(Exon.ForwardLinks)
        EdgeCountHistogram[Count] = EdgeCountHistogram.get(Count, 0) + 1
        if VerboseFlag:
            if len(Exon.Sequence) < 50:
                print "  ", Exon.Sequence
            else:
                print "  %s...%s"%(Exon.Sequence[:20], Exon.Sequence[-20:])
        Str = struct.pack("<iiii", Interval.Start, Interval.End, len(Exon.Sequence), Interval.Occurrences)
        File.write(Str)
        File.write(Exon.Sequence)
        if self.ReverseFlag:
##            ForwardLinks = Exon.BackwardLinks
##            BackwardLinks = Exon.ForwardLinks
##            BackwardLinkAA = Exon.ForwardLinkAA
##            BackwardLinkPower = Exon.ForwardLinkPower
            ForwardLinks = Exon.ForwardLinks
            BackwardLinks = Exon.BackwardLinks
            BackwardLinkAA = Exon.BackwardLinkAA
            BackwardLinkPower = Exon.BackwardLinkPower
        else:
            ForwardLinks = Exon.ForwardLinks
            BackwardLinks = Exon.BackwardLinks
            BackwardLinkAA = Exon.BackwardLinkAA
            BackwardLinkPower = Exon.BackwardLinkPower
        if VerboseFlag:
            print "Prefix %s suffix %s backlinks %s forwardlinks %s"%(Exon.Prefix, Exon.Suffix, len(BackwardLinks), len(ForwardLinks))
        Str = struct.pack("<2s2sii", Exon.Prefix, Exon.Suffix, len(BackwardLinks), len(ForwardLinks))
        File.write(Str)
        # Write the exon-index, and power, and amino acid, for each backward edge:
        for Index in range(len(BackwardLinks)):
            OtherExon = BackwardLinks[Index]
            if VerboseFlag:
                print "Back link to exon %s, AA %s, support %s"%(OtherExon.Index, BackwardLinkAA[Index], BackwardLinkPower[Index])
            Str = struct.pack("<iic", OtherExon.Index, BackwardLinkPower[Index], BackwardLinkAA[Index])
            File.write(Str)
##            OtherInterval = BackwardLinks[Index]
##            for OtherExon in OtherInterval.Exons:
##                DNA = OtherExon.Suffix + Exon.Prefix
##                if len(DNA) == 0:
##                    AA = '\0'
##                    print "Back link to exon %s, %s, AA %s, support %s"%(OtherExon.Index, DNA, AA, BackwardLinkCounts[Index])
##                    Str = struct.pack("<iic", OtherExon.Index, BackwardLinkCounts[Index], AA)
##                    File.write(Str)
##                elif len(DNA) == 3:
##                    AA = GeneticCode.get(DNA, 'X')
##                    print "Back link to exon %s, %s, AA %s, support %s"%(OtherExon.Index, DNA, AA, BackwardLinkCounts[Index])
##                    Str = struct.pack("<iic", OtherExon.Index, BackwardLinkCounts[Index], AA)
##                    File.write(Str)
##                
##        for (AA, OtherExon, Flags) in Exon.EdgesBackward:
##            Index = Gene.Exons.index(OtherExon)
##            Str = struct.pack("<ic", Index, AA)
##            File.write(Str)

    def GenerateExons(self, Intervals):
        "Generate exons by interpreting each interval in each reading frame."
        Exons = []
        GenomeFile = open(ChromosomePaths[self.Chromosome], "r")
        for Interval in Intervals:
            Interval.Exons = []
        for Interval in Intervals:
            print "Interval %s-%s DNA"%(Interval.Start, Interval.End)
            GenomeFile.seek(Interval.Start)
            DNA = GenomeFile.read(Interval.End - Interval.Start).upper()
            #print DNA
            if self.ReverseFlag:
                DNA = ReverseComplement(DNA)
            # Three exons, for three frameshifts:
            Exon = ExonClass(Interval.Start, Interval.End)
            Exon.Shift = 0
            Exon.Sequence = Translate(DNA)
            Exon.Sequence = PruneExon(Exon.Sequence)
            Exon.Prefix = ""
            Remainder = len(DNA)%3
            if Remainder:
                Exon.Suffix = DNA[-Remainder:]
            Interval.Exons.append(Exon)
            Exon.Index = len(Exons)
            Exons.append(Exon)
            if VerboseFlag:
                print "Interval %s-%s exon %s pre %s suff %s"%(Interval.Start, Interval.End, Exon.Index, Exon.Prefix, Exon.Suffix)
            #
            Exon = ExonClass(Interval.Start, Interval.End)
            Exon.Shift = 1
            Exon.Sequence = Translate(DNA[1:])
            Exon.Sequence = PruneExon(Exon.Sequence)
            Exon.Prefix = DNA[:1]
            Remainder = (len(DNA)-1)%3
            if Remainder:
                Exon.Suffix = DNA[-Remainder:]
            Interval.Exons.append(Exon)
            Exon.Index = len(Exons)
            Exons.append(Exon)
            if VerboseFlag:
                print "Interval %s-%s exon %s pre %s suff %s"%(Interval.Start, Interval.End, Exon.Index, Exon.Prefix, Exon.Suffix)
            if len(DNA) > 1:            
                Exon = ExonClass(Interval.Start, Interval.End)
                Exon.Shift = 2
                Exon.Sequence = Translate(DNA[2:])
                Exon.Sequence = PruneExon(Exon.Sequence)
                Exon.Prefix = DNA[:2]
                Remainder = (len(DNA)-2)%3
                if Remainder:
                    Exon.Suffix = DNA[-Remainder:]
                Interval.Exons.append(Exon)
                Exon.Index = len(Exons)
                Exons.append(Exon)
                if VerboseFlag:
                    print "Interval %s-%s exon %s pre %s suff %s"%(Interval.Start, Interval.End, Exon.Index, Exon.Prefix, Exon.Suffix)
        if VerboseFlag:                    
            for Interval in Intervals:
                print Interval, "%s-%s: %s exons"%(Interval.Start, Interval.End, len(Interval.Exons))
        self.AddExonLinks(Intervals, Exons)
        return Exons
    def AddExonLinks(self, Intervals, Exons):
        "Construct all legal edges between exons."
        print "AXL"
        for IntervalA in Intervals:
            if VerboseFlag:
                print "Interval A %s %s-%s"%(IntervalA, IntervalA.Start, IntervalA.End)
            for ExonA in IntervalA.Exons:
                SuffixA = len(ExonA.Suffix)
                if self.ReverseFlag:
                    LinksA = IntervalA.BackwardLinks
                    LinkScoresA = IntervalA.BackwardLinkCounts
                else:
                    LinksA = IntervalA.ForwardLinks
                    LinkScoresA = IntervalA.ForwardLinkCounts
                for LinkIndexA in range(len(LinksA)):
                    LinkQualityA = LinkScoresA[LinkIndexA]
                    IntervalB = LinksA[LinkIndexA]
                    # If that interval goes beyond our ken, then skip it:
                    if IntervalB not in Intervals:
                        continue
                    if len(IntervalB.Exons)<2:
                        print "*** Error: Target interval has no exons!"
                        print IntervalB, IntervalB.Start, IntervalB.End
                        print Intervals.index(IntervalB)
                    if SuffixA == 0:
                        AA = "\0"
                        ExonB = IntervalB.Exons[0]
                        ExonA.ForwardLinks.append(ExonB)
                        ExonA.ForwardLinkPower.append(LinkQualityA)
                        ExonA.ForwardLinkAA.append(AA)
                        ExonB.BackwardLinks.append(ExonA)
                        ExonB.BackwardLinkPower.append(LinkQualityA)
                        ExonB.BackwardLinkAA.append(AA)
                        if VerboseFlag:
                            print "Exon %d (%s/%s) joins to exon %d (%s/%s) via '%s'"%(ExonA.Index, ExonA.Prefix, ExonA.Suffix,
                                                                              ExonB.Index, ExonB.Prefix, ExonB.Suffix, AA)
                    elif SuffixA == 2:
                        ExonB = IntervalB.Exons[1]
                        AA = Translate(ExonA.Suffix + ExonB.Prefix)
                        ExonA.ForwardLinks.append(ExonB)
                        ExonA.ForwardLinkPower.append(LinkQualityA)
                        ExonA.ForwardLinkAA.append(AA)
                        ExonB.BackwardLinks.append(ExonA)
                        ExonB.BackwardLinkPower.append(LinkQualityA)
                        ExonB.BackwardLinkAA.append(AA)
                        if VerboseFlag:
                            print "Exon %d (%s/%s) joins to exon %d (%s/%s) via '%s'"%(ExonA.Index, ExonA.Prefix, ExonA.Suffix,
                                                                              ExonB.Index, ExonB.Prefix, ExonB.Suffix, AA)
                        
                    elif SuffixA == 1 and len(IntervalB.Exons) > 2:
                        ExonB = IntervalB.Exons[2]
                        AA = Translate(ExonA.Suffix + ExonB.Prefix)
                        ExonA.ForwardLinks.append(ExonB)
                        ExonA.ForwardLinkPower.append(LinkQualityA)
                        ExonA.ForwardLinkAA.append(AA)
                        ExonB.BackwardLinks.append(ExonA)
                        ExonB.BackwardLinkPower.append(LinkQualityA)
                        ExonB.BackwardLinkAA.append(AA)
                        if VerboseFlag:
                            print "Exon %d (%s/%s) joins to exon %d (%s/%s) via '%s'"%(ExonA.Index, ExonA.Prefix, ExonA.Suffix,
                                                                              ExonB.Index, ExonB.Prefix, ExonB.Suffix, AA)
                        
                    else:
                        # This is the tricky case.  The source exon has one base left,
                        # and the target exon is only one base long.  Link from the ExonA
                        # to all targets of ExonB:
                        ExonB = IntervalB.Exons[0]
                        if self.ReverseFlag:
                            LinksB = IntervalB.BackwardLinks
                            LinkScoresB = IntervalB.BackwardLinkCounts
                        else:
                            LinksB = IntervalB.ForwardLinks
                            LinkScoresB = IntervalB.ForwardLinkCounts
                        for LinkIndexB in range(len(LinksB)):
                            LinkQualityB = LinkScoresB[LinkIndexB]
                            IntervalC = LinksB[LinkIndexB]
                            if IntervalC not in Intervals:
                                continue
                            
                            ExonC = IntervalC.Exons[1]
                            if LinkQualityA < 0:
                                Power = LinkQualityB
                            elif LinkQualityB < 0:
                                Power = LinkQualityA
                            else:
                                Power = min(LinkQualityA, LinkQualityB)
                            AA = Translate(ExonA.Suffix + ExonB.Suffix + ExonC.Prefix)
                            ExonA.ForwardLinks.append(ExonC)
                            ExonA.ForwardLinkPower.append(Power)
                            ExonA.ForwardLinkAA.append(AA)
                            ExonC.BackwardLinks.append(ExonA)
                            ExonC.BackwardLinkPower.append(Power)
                            ExonC.BackwardLinkAA.append(AA)
                            if VerboseFlag:
                                print "** Double-link %s+%s+%s"%(ExonA.Suffix, ExonB.Suffix, ExonC.Prefix)                                
                                print "Exon %d (%s/%s) joins to exon %d (%s/%s) via '%s'"%(ExonA.Index, ExonA.Prefix, ExonA.Suffix,
                                                                                  ExonC.Index, ExonC.Prefix, ExonC.Suffix, AA)
                            #sys.stdin.readline()
                            
    def AssimilateLinks(self, Source, Master, ForwardFlag):
        if ForwardFlag:
            LinkSource = Source.ForwardLinks
            CountSource = Source.ForwardLinkCounts
            LinkMaster = Master.ForwardLinks
            CountMaster = Master.ForwardLinkCounts
        else:
            LinkSource = Source.BackwardLinks
            CountSource = Source.BackwardLinkCounts
            LinkMaster = Master.BackwardLinks
            CountMaster = Master.BackwardLinkCounts
        for Index in range(len(LinkSource)):
            Link = LinkSource[Index]
            try:
                MasterIndex = LinkMaster.index(Link)
                if CountMaster[MasterIndex]>=0 and CountSource[Index] >= 0:
                    CountMaster[MasterIndex] += CountSource[Index]
                else:
                    CountMaster[MasterIndex] = -1 # contiguous!
            except:
                LinkMaster.append(Link)
                CountMaster.append(CountSource[Index])
        # Cut links from other intervals ("third party" intervals) into the source interval:
        for Index in range(len(LinkSource)):
            ThirdParty = LinkSource[Index]
            if ForwardFlag:
                try:
                    ThirdPartySourceIndex = ThirdParty.BackwardLinks.index(Source)
                except:
                    ThirdPartySourceIndex = None
                try:
                    ThirdPartyMasterIndex = ThirdParty.BackwardLinks.index(Master)
                except:
                    ThirdPartyMasterIndex = None
                if ThirdPartySourceIndex!=None:
                    if ThirdPartyMasterIndex!=None:
                        ThirdParty.BackwardLinkCounts[ThirdPartyMasterIndex] += ThirdParty.BackwardLinkCounts[ThirdPartySourceIndex]
                        del ThirdParty.BackwardLinks[ThirdPartySourceIndex]
                        del ThirdParty.BackwardLinkCounts[ThirdPartySourceIndex]
                    else:
                        ThirdParty.BackwardLinks[ThirdPartySourceIndex] = Master
            else:
                try:
                    ThirdPartySourceIndex = ThirdParty.ForwardLinks.index(Source)
                except:
                    ThirdPartySourceIndex = None
                try:
                    ThirdPartyMasterIndex = ThirdParty.ForwardLinks.index(Master)
                except:
                    ThirdPartyMasterIndex = None
                if ThirdPartySourceIndex!=None:
                    if ThirdPartyMasterIndex!=None:
                        ThirdParty.ForwardLinkCounts[ThirdPartyMasterIndex] += ThirdParty.ForwardLinkCounts[ThirdPartySourceIndex]
                        del ThirdParty.ForwardLinks[ThirdPartySourceIndex]
                        del ThirdParty.ForwardLinkCounts[ThirdPartySourceIndex]
                    else:
                        ThirdParty.ForwardLinks[ThirdPartySourceIndex] = Master
##    def AssimilateLinks(self, LinkSource, CountSource, LinkMaster, CountMaster):
##        "Merge all entries from LinkSource, CountSource into LinkMaster, CountMaster"
##        for Index in range(len(LinkSource)):
##            Link = LinkSource[Index]
##            try:
##                MasterIndex = LinkMaster.index(Link)
##                if CountMaster[MasterIndex]>=0 and CountSource[Index] >= 0:
##                    CountMaster[MasterIndex] += CountSource[Index]
##                else:
##                    CountMaster[MasterIndex] = -1 # contiguous!
##            except:
##                LinkMaster.append(Link)
##                CountMaster.append(CountSource[Index])
    def WriteGIIDIndex(self, FileName):
        Keys = GIIDDict.keys()
        Keys.sort()
        print "Write out file positions for %s GIIDs to %s..."%(len(Keys), FileName)
        File = open(FileName, "wb")
        for Key in Keys:
            Value = GIIDDict[Key]
            if type(Value)==type([]):
                Str = "%s\t"%Key
                for Entry in Value:
                    Str += "%s,"%Entry
            else:
                Str = "%s\t%s\t"%(Key, GIIDDict[Key])
            File.write(Str+"\n")
        File.close()
def PrintHistogram(Histogram, Title):
    Count = 0
    print 
    print Title
    Keys = Histogram.keys()
    Keys.sort()
    for Key in Keys:
        Count += Histogram[Key]
    for Key in Keys:

        print "%s: %s (%.2f%%)"%(Key, Histogram[Key], 100*Histogram[Key] / float(Count))


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "<no psyco>"
    # Command-line arguments: chromosome number and + or -
    InputFileName = sys.argv[1]
    OutputFileName = sys.argv[2]
##    ChromosomeNumber = int(sys.argv[1])
##    try:
##        Direction = sys.argv[2][0]
##    except:
##        Direction = "+"
##    if Direction == "-":
##        ReverseFlag = 1
##    else:
##        Direction = "+"
##        ReverseFlag = 0
##    if ChromosomeNumber == 1 and ReverseFlag == 0:
##        GIIDDict = {}
##        print "Track GIIDs."
    Collector = ExonCollector(None, None)
    Stub = os.path.split(InputFileName)[1]
    Stub = os.path.splitext(Stub)[0]
    Stub = Stub[:-1]
    Collector.Chromosome = int(Stub)
    #Collector.GatherExons(r"C:\ftproot\liliana\TestRecord2.txt")
    #Collector.GatherExons(r"C:\ftproot\liliana\humrs1-hg17.polishes-good.noal")
    #Collector.GatherExons(r"C:\ftproot\liliana\polishes-good.noal")
    Collector.GatherExons(InputFileName)
##    if ChromosomeNumber == 1 and ReverseFlag == 0:
##        Collector.WriteGIIDIndex("C:\\ftproot\\liliana\\FullGIIDIndex.txt")
    Collector.WriteGenes(OutputFileName)
    #Collector.WriteGenes("NewSpliceDB\\%s%s.dat"%(ChromosomeNumber, Direction))
    PrintHistogram(SequenceLengthHistogram, "Sequence Lengths")
    PrintHistogram(EdgeCountHistogram, "Edge count histogram")
    
        