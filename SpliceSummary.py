"""
Post-processing of splice-tolerant search results.  Here's the plan:
- Read in all the results that hit the chromosome strand in question
- Select the RefSeq entry which covers the most spectra.
(An entry "covers" a spectrum if the RefSeq entry's intervals cover all
genomic positions for that spectrum's interpretation)  Flag all of these
spectra as "covered"
- Repeat the previous step until the next RefSeq entry can explain few additional spectra
- Report all "overlapping" refseq entries that were iteratively selected; these are
 evidence for alternative splicing
- Report genomic positions that were covered by many good search hits, but never covered by
 a RefSeq entry.  These are evidence for uncharacterized (or mis-annotated) genes.
"""
import sys
import os
import psyco
import scipy
import MakeImage
import cStringIO
from Utils import *
Initialize()
import MSSpectrum

class GeneClass:
    def __init__(self):
        self.Name = None
        self.Intervals = []
        self.CanCoverList = []
        self.CanCoverPower = 0
        self.AcceptedFlag = 0
        self.CoverPower = 0

class SplicedPeptideClass:
    def __init__(self):
        self.Score = None
        self.StartA = None
        self.StartB = None
        self.EndA = None
        self.EndB = None
        self.Count = 0
        self.Cover = None
        
class SpliceReader:
    def __init__(self, ChromosomeNumber, Strand):
        self.ChromosomeNumber = ChromosomeNumber
        self.Strand = Strand
        self.MinMQScore = 1.0
        self.MinLength = 8
        self.Annotations = {}
        self.KnownGenes = []
        # Iterative gene selection - The 'next' gene must cover this many:
        self.MinGeneCoverCount = 5
        # Number of as-yet-uncovered spectra that must remain for a position to be
        # considered 'interesting':
        self.MinimumUncoveredInteresting = 10
    def ReadKnownGenes(self, FileName):
        File = open(FileName, "r")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if self.Strand == 0:
                if Bits[3]!="-":
                    continue
            else:
                if Bits[3]!="+":
                    continue
            Chromosome = ChromosomeMap.get(Bits[2], None)
            if Chromosome != self.ChromosomeNumber:
                continue
            Gene = GeneClass()
            Gene.Name = Bits[0]
            Gene.Name2 = Bits[1]
            StartBits = Bits[9].split(",")
            EndBits = Bits[10].split(",")
            for Index in range(len(StartBits)-1):
                Start = int(StartBits[Index])
                End = int(EndBits[Index])
                Gene.Intervals.append((Start, End))
            Gene.Start = Gene.Intervals[0][0]
            Gene.End = Gene.Intervals[0][1]
            self.KnownGenes.append(Gene)
        File.close()
        print "#Loaded a total of %s known genes."%len(self.KnownGenes)
    def ReadAnnotations(self, FileName):
        try:
            File = open(FileName, "r")
        except:
            return # directory, probably.
        OldSpectrum = None
        # The Annotations dictionary maps INTERVALS to tuples
        # of the form (AnnotationString, BestMQScore, Occurrences)
        # Assume that each peptide has at most one 'true' splice junction.
        # So: keys have the form (StartA, EndA, StartB, EndB), where StartB/EndB can be None.
        self.Annotations = {}
        AnnotationCount = 0
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if len(Bits)<2:
                continue
            Spectrum = (Bits[0], Bits[1])
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            #print "Sc %s chrom %s specfilepos %s"%(Bits[5], Bits[21], Bits[20])
            try:
                Score = float(Bits[5])
            except:
                continue
            if Score < self.MinMQScore:
                continue
            try:
                Chrom = int(Bits[21])
            except:
                continue
            if Chrom != self.ChromosomeNumber:
                continue
            if int(Bits[22]) != self.Strand:
                continue
            try:
                (OverallStart, OverallEnd) = Bits[23].split("-")
            except:
                print FileName
                print FileLine
                print Bits[23]
                raise
            OverallStart = int(OverallStart)
            OverallEnd = int(OverallEnd)
            AnnotationCount += 1
            if len(Bits)>25 and Bits[25]:
                SplicePoints = Bits[25].split("-")
                if len(SplicePoints) > 2:
                    print FileLine
                    continue # For now, drop hits that span two splice boundaries
                A = int(SplicePoints[0])
                B = int(SplicePoints[1])
                StartA = OverallStart
                EndA = min(A, B)
                StartB = max(A, B)
                EndB = OverallEnd
            else:
                StartA = OverallStart
                EndA = OverallEnd
                StartB = None
                EndB = None
                
            Key = (StartA, EndA, StartB, EndB)
            if self.Annotations.has_key(Key):
                Peptide = self.Annotations[Key]
                Peptide.Score = max(Peptide.Score, Score)
                Peptide.Count += 1
            else:
                Peptide = SplicedPeptideClass()
                Peptide.Score = Score
                Peptide.Count = 1
                Peptide.Annotation = Bits[2]
                Peptide.StartA = StartA
                Peptide.EndA = EndA
                Peptide.StartB = StartB
                Peptide.EndB = EndB
                self.Annotations[Key] = Peptide
            #self.Annotations[Key] = (Annotation, Score, Count)
        print "#Loaded a total of %s annotations in %s peptides"%(AnnotationCount, len(self.Annotations.keys()))
        File.close()
    def IterativelySelectGenes(self):
        while (1):
            for Gene in self.KnownGenes:
                Gene.CanCoverPower = 0
                Gene.CanCoverList = []
            for Peptide in self.Annotations.values():
                if Peptide.Cover:
                    continue
                Start = Peptide.StartA
                End = Peptide.EndA
                if Peptide.EndB:
                    End = Peptide.EndB
                else:
                    CoveredB = 1
                for Gene in self.KnownGenes:
                    if Gene.CoverPower or Gene.Start > End or Gene.End < Start:
                        continue
                    # Verify that the gene covers the peptide:
                    CoveredA = 0
                    if Peptide.EndB:
                        CoveredB = 0
                    for (Start, End) in Gene.Intervals:
                        if (Start <= Peptide.StartA and End >= Peptide.EndA):
                            CoveredA = 1
                        if not CoveredB:
                            if (Start <= Peptide.StartB and End >= Peptide.EndB):
                                CoveredB = 1
                    if CoveredA and CoveredB:
                        Gene.CanCoverPower += Peptide.Count
                        Gene.CanCoverList.append(Peptide)
            BestPower = 0
            BestGene = None
            for Gene in self.KnownGenes:
                if Gene.CanCoverPower > BestPower:
                    BestPower = Gene.CanCoverPower
                    BestGene = Gene
            if BestPower < self.MinGeneCoverCount:
                print "#Best power is %s < %s - stop iteration!"%(BestPower, self.MinGeneCoverCount)
                break
            BestGene.CoverPower = BestPower
            for Peptide in BestGene.CanCoverList:
                Peptide.Cover = BestGene
            print "Accept gene %s (%s) to cover %s annotations."%(BestGene.Name, BestGene.Name2, BestPower)
    def SummarizeUncoveredPeptides(self):
        Len = 200000000
        UncoveredSpectra = scipy.zeros(Len, scipy.int8)
        Peptides = self.Annotations.values()
        PeptideCount = len(Peptides)
        for PeptideIndex in range(PeptideCount):
            print "%s of %s..."%(PeptideIndex, PeptideCount)
            Peptide = Peptides[PeptideIndex]
            if Peptide.Cover:
                continue
            End = Peptide.EndA
            if Peptide.EndB:
                End = Peptide.EndB
            if End > Len:
                print "**!!*! Extend UncoveredSpectra... by %s zeroes"%(End - Len)
                return
            for X in range(Peptide.StartA, Peptide.EndA):
                print X
                UncoveredSpectra[X] += 1
            #for Pos in range(Peptide.StartA, Peptide.EndA):
            #    UncoveredSpectra[Pos] += Peptide.Count
            #if Peptide.EndB:
            #    for Pos in range(Peptide.StartB, Peptide.EndB):
            #        UncoveredSpectra[Pos] += Peptide.Count
        IntervalStart = None
        IntervalMax = 0
        print "# Find uncovered positions..."
        for Pos in xrange(Len):
            if Pos%100000 == 0:
                print "#%s..."%Pos
            Count = UncoveredSpectra[Pos]
            if Count > self.MinimumUncoveredInteresting:
                if IntervalStart == None:
                    IntervalStart = Pos
                    IntervalMax = Count
                else:
                    IntervalMax = max(Count, IntervalMax)
            elif (IntervalStart):
                    print "Uncovered\t%s\t%s\t%s\t%s\%s\t"%(self.ChromosomeNumber, self.Strand,
                        IntervalStart, Pos, IntervalMax)
                    IntervalStart = None
                    IntervalMax = 0
    def SummarizeAltSplicing(self):
        """
        Report on all accepted genes that have overlap.
        """
        Len = len(self.KnownGenes)
        for GeneIndex in range(Len):
            Gene = self.KnownGenes[GeneIndex]
            if not Gene.CoverPower:
                continue
            for OtherGeneIndex in range(GeneIndex + 1, Len):
                OtherGene = self.KnownGenes[OtherGeneIndex]
                if not OtherGene.CoverPower:
                    continue
                if OtherGene.End < Gene.Start or OtherGene.Start > Gene.End:
                    continue
                print "AltSplicing\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(self.ChromosomeNumber, self.Strand,
                    Gene.Name, Gene.Start, Gene.End, Gene.CoverPower,
                    OtherGene.Name, OtherGene.Start, OtherGene.End, OtherGene.CoverPower)

if __name__ == "__main__":
    psyco.full()
    Reader = SpliceReader(int(sys.argv[1]), int(sys.argv[2]))
    print "#Load known genes..."
    sys.stdout.flush()
    Reader.ReadKnownGenes("Splice\\RefFlatt.txt")
    print "#Load annotations..."
    sys.stdout.flush()
    Dir = "e:\\ms\\peptideatlas\\GenomicResults\\"
    for FileName in os.listdir(Dir):
        Path = os.path.join(Dir, FileName)
        Reader.ReadAnnotations(Path)
    print "#Iterative gene selection..."
    sys.stdout.flush()
    Reader.IterativelySelectGenes()
    print "#Summarize alternative splicing..."
    sys.stdout.flush()
    Reader.SummarizeAltSplicing()
    print "#Summarize uncovered positions..."
    sys.stdout.flush()
    Reader.SummarizeUncoveredPeptides()