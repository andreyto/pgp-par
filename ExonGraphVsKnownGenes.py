"""
Split exon-graph search results into various categories based on whether
the results correspond to known proteins.

One way we want to analyze the exon graph search results is by checking which exons
and splice junctions come from known genes, which contain SNPs, and which appear to
be new.
Let's categorize all the high-scoring search hits:

Category 0: p-value >0.05 or length <8.  Discard!

Category 1: exon from a known gene, sequence matches
Category 2: exon from a known gene, sequence mismatch
Category 3: single exon, not contained in a known gene

Category 4: spliced peptide, from a known gene, exon matches
Category 5: spliced peptide, from a known gene, sequence mismatch
Category 6: spliced peptide, one endpoint matching a known junction
Category 7: spliced peptide, neither endpoint matching a known junction

The following categories are rarely (or never) produced by new searches:
Category 8: Not contained in a known protein, and doesn't match the genome sequence.
 (May use too many SNPs, or may have mis-reported genomic coordinates in older hits)
Category 9: Unreasonably short internal exon, e.g. XXX;X;XX:XXX

Once we've categorized all high-scoring hits for a scan, we output
the "best" of those hit(s).  Categories 1 and 4 are preferred, then
categories 2 and 5, then categories 3 and 6, then category 7.

"""
import struct
import os
import sys
import traceback
import string
import time
import GeneMapper

NO_ORACLE = 0

RESULTS_EXTENSION = ".txt"
MAX_CATEGORY = 10
PValueLimit = 0.05
MinimumPeptideLength = 8
class ExonClass:
    def __init__(self):
        self.EdgeAA = None

class GeneClass:
    "A gene (a known protein that's been mapped to the genome)"
    def __init__(self):
        self.Exons = [] # instances of ExonClas
        self.Chromosome = None
        self.Strand = 1 # 1 or -1
        self.Start = None
        self.End = None
        self.ProteinNumber = None
        self.ProteinName = None
        self.Sequence = None
    def __cmp__(self, Other):
        if not Other:
            return 1
        if self.Start < Other.Start:
            return -1
        if self.Start > Other.Start:
            return 1
        if self.End < Other.End:
            return -1
        if self.End > Other.End:
            return 1
        return 0


class HitClass:
    "A single search result"
    def __init__(self, Bits):
        self.Donor = None
        self.Acceptor = None
        self.DonorA = None
        self.AcceptorA = None
        self.DonorB = None
        self.AcceptorB = None
        self.Bits = list(Bits)
        while len(self.Bits) < 21:
            self.Bits.append("")
        #self.Spectrum = (Bits[0], Bits[1])
        self.Chromosome = int(Bits[16])
        self.Strand = int(Bits[17])
        if self.Strand == 0:
            self.Strand = -1
        self.Aminos = Bits[2][2:-2]
        self.PValue = float(Bits[10])
        if len(Bits) > 20:
            SpliceInfo = Bits[20].split("-")
        else:
            SpliceInfo = []
        if len(SpliceInfo) > 1:
            self.ExonCount = len(SpliceInfo)
        else:
            self.ExonCount = 1
        (Start, End) = self.Bits[18].split("-")
        self.Start = int(Start)
        self.End = int(End)
    def GetCategoryRank(self):
        """
        Return the RANK of our category.  If our category's rank isn't equal to
        the best category-rank for the scan, then we shall be discarded.
        """
        if NO_ORACLE:
            # All "real" categories get the same rank!
            if self.Category == 0:
                self.CategoryRank = 0
            elif self.Category in (8, 9):
                self.CategoryRank = 1
            else:
                self.CategoryRank = 2
            return self.CategoryRank
        if self.Category == 0:
            self.CategoryRank = 0
        elif self.Category in (8,9):
            self.CategoryRank = 1
        elif self.Category in (7,):
            self.CategoryRank = 2
        elif self.Category in (3, 6):
            self.CategoryRank = 3
        elif self.Category in (2, 5):
            self.CategoryRank = 4
        elif self.Category in (1, 4):
            self.CategoryRank = 5
        return self.CategoryRank
    def Categorize(self):
        """
        Categorize this peptide, according to these criteria
        - Single exon, two exons, or three+ exons?
        - Contained in a known protein?
        - If so, matches the known protein's sequence?
        """
        if self.PValue > PValueLimit:
            self.Category = 0
            self.Bits.append("Cat0")
            return self.Category
        if len(self.Aminos) < MinimumPeptideLength:
            self.Category = 0
            self.Bits.append("Cat0")
            return self.Category
        if self.ExonCount == 1:
            self.CategorizeSingleExonHit()
            return self.Category
        if self.ExonCount == 2:
            self.CategorizeTwoExonHit()
            return self.Category
        if self.ExonCount == 3:
            self.CategorizeThreeExonHit()
            return self.Category
        # More than three exons?  Call it category 9, TOO MANY exons:    
        self.Category = 9
        self.Bits.append("Cat9")
        self.Bits.append("")
        self.Bits.append("")
        self.Bits.append("")
        self.Bits.append("")
        self.Bits.append("")
        return self.Category
    def GetGenomicSNPInfo(self):
        if self.ExonCount == 1:
            GenomePeptide = Analyzer.GetGenomeSequence(self.Chromosome, self.Strand, self.Start, self.End)
        elif self.ExonCount == 2:
            if (self.Donor < self.Start or self.Acceptor < self.Start):
                GenomicSNPInfo = "* Invalid coords"
                return (len(self.Aminos), GenomicSNPInfo)
            if (self.Donor > self.End or self.Acceptor > self.End):
                GenomicSNPInfo = "* Invalid coords"
                return (len(self.Aminos), GenomicSNPInfo)
            GenomePeptide = Analyzer.GetGenomeSequence(self.Chromosome, self.Strand, self.Start, self.End, self.Donor, self.Acceptor)
        elif self.ExonCount == 3:
            if (self.DonorA < self.Start or self.AcceptorA < self.Start or
                self.DonorB < self.Start or self.AcceptorB < self.Start or
                self.DonorA > self.End or self.AcceptorB > self.End or
                self.DonorB > self.End or self.AcceptorA > self.End):
                GenomicSNPInfo = "* Invalid coords"
                return (len(self.Aminos), GenomicSNPInfo)
            GenomePeptide = Analyzer.GetGenomeSequence3(self.Chromosome, self.Strand, self.Start, self.End, self.DonorA, self.AcceptorA,
                                                        self.DonorB, self.AcceptorB)
        GenomicSNPInfo = ""
        if self.Strand == 1:
            GenomePos = self.Start
        else:
            GenomePos = self.End - 1
        MismatchCount = 0
        Len = min(len(self.Aminos), len(GenomePeptide))
        for Index in range(Len):
            if self.Aminos[Index] != GenomePeptide[Index]:
                MismatchCount += 1
                GenomicSNPInfo += "%s (%s/%s),"%(GenomePos, self.Aminos[Index], GenomePeptide[Index])
            if self.Strand == 1:
                GenomePos += 3
                if (self.Donor != None and GenomePos >= self.Donor and GenomePos < self.Acceptor):
                    GenomePos += (self.Acceptor - self.Donor)
                if (self.DonorA != None and GenomePos >= self.DonorA and GenomePos < self.AcceptorA):
                    GenomePos += (self.AcceptorA - self.DonorA)
                if (self.DonorB != None and GenomePos >= self.DonorB and GenomePos < self.AcceptorB):
                    GenomePos += (self.AcceptorB - self.DonorB)
            else:
                GenomePos -= 3
                if (self.Acceptor != None and GenomePos <= self.Donor and GenomePos > self.Acceptor):
                    GenomePos += (self.Acceptor - self.Donor)
                if (self.AcceptorA != None and GenomePos <= self.DonorA and GenomePos > self.AcceptorA):
                    GenomePos += (self.AcceptorA - self.DonorA)
                if (self.AcceptorB != None and GenomePos <= self.DonorB and GenomePos > self.AcceptorB):
                    GenomePos += (self.AcceptorB - self.DonorB)
        if len(GenomePeptide) != len(self.Aminos):
            # The lengths of the peptides are different.  This results from
            # a bug (since fixed) in SetMatchPrefixSuffix, where genomic coordinates
            # could be reported incorrectly.
            # If the genomic peptide includes the peptide, then that's still a match:
            if GenomePeptide.find(self.Aminos) != -1:
                GenomicSNPInfo = ""
                MismatchCount = 0
            else:
                # This *may* confirm a SNP, for now just report it as a (special) mismatch.
                GenomicSNPInfo = "*MISLENGTH (%s %s)"%(GenomePeptide, self.Aminos)
                MismatchCount = len(self.Aminos)
        elif MismatchCount > 3:
            GenomicSNPInfo = "*MISMATCH (%s %s)"%(GenomePeptide, self.Aminos)
        else:
            GenomicSNPInfo = GenomicSNPInfo[:-1] # remove trailing comma        
        return (MismatchCount, GenomicSNPInfo)
    def CategorizeSingleExonHit(self):
        #print "CSEH:", self.Bits[18], self.Bits[2]
        Verbose = 0
        # Parse stuff:
        NearGeneLeft = None
        NearGeneRight = None
        CoveringGenes = []
        SpanningGenes = []
        TerminusFlag = "" # by default, the peptide isn't at either terminus of its protein.
        #######################################################
        # Find the covering gene, and/or near neighbors:
        GeneList = Analyzer.GenesByChromosome[self.Chromosome]
        #if Verbose:
        #    print "%s look through a list of %s genes for chromosome %s strand %s"%(self.Aminos, len(GeneList), self.Chromosome, self.Strand)
        (MismatchCount, GenomicSNPInfo) = self.GetGenomicSNPInfo()
        for Gene in GeneList:
            #if Verbose:
            #    print Gene.Strand, Gene.Start, Gene.End
            if Gene.Strand != self.Strand:
                continue
            if Gene.End < self.Start:
                NearGeneLeft = Gene
                continue                
            if Gene.Start > self.End:
                NearGeneRight = Gene
                break
            CoverFlag = 0
            for Exon in Gene.Exons:
                if Exon.Start <= self.Start and Exon.End >= self.End:
                    # Add this covering-gene to the list:
                    CoveringGenes.append(Gene)
                    CoverFlag = 1
            if not CoverFlag:
                SpanningGenes.append(Gene)
        # If we covered the peptide with known genes, use category 1 or 2:
        if CoveringGenes:
            ##################################################
            # Find the covering-gene which matches us best.  If there's one with
            # a perfect match, note it!
            BestIdentity = -1
            CoveringGeneStr = ""
            BestCoveringGeneStr = ""
            for TestCoveringGene in CoveringGenes:
                CoveringGeneStr += "%s:"%TestCoveringGene.ProteinNumber
                Pos = TestCoveringGene.Sequence.find(self.Aminos)
                if Pos != -1:
                    if BestIdentity < 1:
                        BestCoveringGeneStr = "%s:"%TestCoveringGene.ProteinNumber
                    else:
                        BestCoveringGeneStr = "%s%s:"%(BestCoveringGeneStr, TestCoveringGene.ProteinNumber)
                    BestIdentity = 1
                    BestCoveringGene = TestCoveringGene
                    BestGeneSNPInfo = ""
                    if Pos == 0:
                        TerminusFlag = "N"
                    elif Pos + len(self.Aminos) == len(TestCoveringGene.Sequence):
                        TerminusFlag = "C"
                    else:
                        TerminusFlag = ""
                else:
                    (Identity, GeneSNPInfo) = Analyzer.GetMismatchInfo(self.Bits, TestCoveringGene, self.Aminos, TestCoveringGene.Sequence)
                    if (Identity > BestIdentity) or (Identity == BestIdentity and TestCoveringGene.CoveragePercent > BestCoveringGene.CoveragePercent):
                        if BestIdentity < Identity:
                            BestCoveringGeneStr = "%s:"%TestCoveringGene.ProteinNumber
                        else:
                            BestCoveringGeneStr = "%s%s:"%(BestCoveringGeneStr, TestCoveringGene.ProteinNumber)
                        BestIdentity = Identity
                        BestCoveringGene = TestCoveringGene
                        BestGeneSNPInfo = GeneSNPInfo
            CoveringGeneStr = CoveringGeneStr[:-1] # remove trailing colon
            BestCoveringGeneStr = BestCoveringGeneStr[:-1] # remove trailing colon
            if BestIdentity == 1:
                self.Category = 1
                self.Bits.append("Cat1")
                self.Bits.append(GenomicSNPInfo)
                self.Bits.append(CoveringGeneStr)
                self.Bits.append(BestCoveringGeneStr)
                #self.Bits.append(str(CoveringGene.ProteinName))
                self.Bits.append(BestGeneSNPInfo)
                self.Bits.append(TerminusFlag)
            else:
                self.Category = 2
                self.Bits.append("Cat2")
                self.Bits.append(GenomicSNPInfo)
                self.Bits.append(CoveringGeneStr)
                self.Bits.append(str(BestCoveringGene.ProteinNumber))
                #self.Bits.append(str(CoveringGene.ProteinName))
                self.Bits.append(BestGeneSNPInfo)
                self.Bits.append("")
                #self.Bits.append(Analyzer.GetMismatchInfo(self.Bits, CoveringGene, self.Aminos, CoveringGene.Sequence))
        else:
            # No gene covers this exon.  Report what genes are nearby.
            if MismatchCount > 1:
                # No gene covers this exon, and the exon has more than one mismatch with
                # the genome.  It's probably mis-annotated or uses too many SNPs:
                self.Category = 8
                self.Bits.append("Cat8")
            else:
                # If there's a spanning gene that includes this sequence, then promote it
                # to cat1 with a flag:
                for Gene in SpanningGenes:
                    Pos = Gene.Sequence.find(self.Aminos)
                    if Pos != -1:
                        self.Category = 1
                        self.Bits.append("Cat1")
                        self.Bits.append(GenomicSNPInfo)
                        self.Bits.append("*%s"%Gene.ProteinNumber)
                        self.Bits.append("*%s"%Gene.ProteinNumber)
                        self.Bits.append("")
                        if Pos == 0:
                            TerminusFlag = "N"
                        elif Pos + len(self.Aminos) == len(Gene.Sequence):
                            TerminusFlag = "C"
                        else:
                            TerminusFlag = ""
                        self.Bits.append(TerminusFlag)
                        return self.Category
                self.Category = 3
                self.Bits.append("Cat3")
            self.Bits.append(GenomicSNPInfo)
            Str = ""
            for Gene in SpanningGenes:
                Str += "Span (%s-%s) %s "%(Gene.Start, Gene.End, Gene.ProteinNumber)
            self.Bits.append(Str)
            Str = ""
            if NearGeneLeft:
                Str += "Near (%s-%s): %s "%(NearGeneLeft.Start, NearGeneLeft.End, NearGeneLeft.ProteinNumber)
            if NearGeneRight:
                Str += "Near (%s-%s): %s "%(NearGeneRight.Start, NearGeneRight.End, NearGeneRight.ProteinNumber)
            self.Bits.append(Str)
            self.Bits.append("")
            self.Bits.append("")
        return self.Category
    def CategorizeTwoExonHit(self):
        #print "C2EH:", self.Bits[18], self.Bits[2]
        (DonorStr, AcceptorStr) = self.Bits[20].split("-")
        self.Donor = int(DonorStr)
        self.Acceptor = int(AcceptorStr)
        NearGeneLeft = None
        NearGeneRight = None
        CoveringGenes = []
        DonorGenes = []
        AcceptorGenes = []
        SpanningGenes = []
        TerminusFlag = "" # by default, the peptide isn't at either terminus of its protein.
        ##################################################
        # Get genomic SNP info:
        (MismatchCount, GenomicSNPInfo) = self.GetGenomicSNPInfo()
        #######################################################
        # Find the covering gene(s), and/or near neighbors:
        for Gene in Analyzer.GenesByChromosome[self.Chromosome]:
            if Gene.Strand != self.Strand:
                continue                        
            if Gene.End < self.Start:
                NearGeneLeft = Gene
                continue                
            if Gene.Start > self.End:
                NearGeneRight = Gene
                break
            CoverFlag = 0
            DonorFlag = 0
            AcceptorFlag = 0
            for ExonIndex in range(len(Gene.Exons) - 1):
                Exon = Gene.Exons[ExonIndex]
                NextExon = Gene.Exons[ExonIndex + 1]
                #print self.Chromosome, self.Strand, self.Start, self.Donor, self.Acceptor, self.End
                #print Exon.Start, Exon.End, NextExon.Start, NextExon.End
                # Flag whether this known protein can serve as a splice donor, or acceptor, or
                # both for this protein.  We only set DonorFlag and AcceptorFlag if the gene's
                # exon actually covers this exon.
                #
                #  PEP----TIDE   This example shows a case where a gene gets in DonorGenes and
                #   ex----on     AcceptorGenes, but we do NOT set DonorFlag or AcceptorFlag.    
                if self.Strand == 1:
                    if Exon.End == self.Donor:
                        if Gene not in DonorGenes:
                            DonorGenes.append(Gene)
                        if Exon.Start <= self.Start:
                            DonorFlag = 1
                    if NextExon.Start == self.Acceptor:
                        if Gene not in AcceptorGenes:
                            AcceptorGenes.append(Gene)
                        if NextExon.End >= self.End:
                            AcceptorFlag = 1
                else:
                    if NextExon.Start == self.Donor:
                        if Gene not in DonorGenes:
                            DonorGenes.append(Gene)
                        if NextExon.End >= self.End:
                            DonorFlag = 1
                    if Exon.End == self.Acceptor:
                        if Gene not in AcceptorGenes:
                            AcceptorGenes.append(Gene)
                        if Exon.Start <= self.Start:
                            AcceptorFlag = 1
                #print "D%s A%s"%(DonorFlag, AcceptorFlag)
                if DonorFlag and AcceptorFlag and (Gene not in CoveringGenes):
                    CoveringGenes.append(Gene)
                # Reset donor flag and acceptor flag.  To be cat4, you must be covered by
                # two consecutive exons - you can't skip exons from this isoform!
                DonorFlag = 0
                AcceptorFlag = 0
            SpanningGenes.append(Gene)
        # if we covered the exon, then check to see whether it's an exact match:
        if CoveringGenes:
            ##################################################
            # Loop over covering genes to find the one with best sequence identity:
            BestIdentity = -1
            BestCoveringGene = None
            CoveringGeneStr = ""
            for TestCoveringGene in CoveringGenes:
                CoveringGeneStr += "%s:"%TestCoveringGene.ProteinNumber
                Pos = TestCoveringGene.Sequence.find(self.Aminos)
                if Pos != -1:
                    BestCoveringGene = TestCoveringGene
                    BestIdentity = 1
                    BestGeneSNPInfo = ""
                else:
                    (Identity, GeneSNPInfo) = Analyzer.GetMismatchInfo(self.Bits, TestCoveringGene, self.Aminos, TestCoveringGene.Sequence)
                    if (Identity > BestIdentity) or (Identity == BestIdentity and TestCoveringGene.CoveragePercent > BestCoveringGene.CoveragePercent):
                        BestIdentity = Identity
                        BestCoveringGene = TestCoveringGene
                        BestGeneSNPInfo = GeneSNPInfo
            CoveringGeneStr = CoveringGeneStr[:-1] # remove trailing comma
            if BestIdentity == 1:
                self.Category = 4
                self.Bits.append("Cat4")
                self.Bits.append(GenomicSNPInfo)
                self.Bits.append(CoveringGeneStr)
                self.Bits.append(str(BestCoveringGene.ProteinNumber))
                self.Bits.append(BestGeneSNPInfo)
                Pos = BestCoveringGene.Sequence.find(self.Aminos)
                if Pos == 0:
                    TerminusFlag = "N"
                elif Pos + len(self.Aminos) == len(BestCoveringGene.Sequence):
                    TerminusFlag = "C"
                self.Bits.append(TerminusFlag)
            else:
                self.Category = 5
                self.Bits.append("Cat5")
                self.Bits.append(GenomicSNPInfo)
                self.Bits.append(CoveringGeneStr)
                self.Bits.append(str(BestCoveringGene.ProteinNumber))
                self.Bits.append(BestGeneSNPInfo)
                self.Bits.append("")
        elif (DonorGenes or AcceptorGenes):
            # One side of the splice junction is covered.  If the peptide is an exact
            # match to one of these proteins, promote it to category 4 but star the
            # protein ID to indicate that our splicing is not quite right.
            CheckGeneList = DonorGenes[:]
            CheckGeneList.extend(AcceptorGenes)
            for CheckGene in CheckGeneList:
                if CheckGene.Sequence.find(self.Aminos)!=-1:
                    self.Category = 4
                    self.Bits.append("Cat4")
                    self.Bits.append(GenomicSNPInfo)
                    self.Bits.append(str(CheckGene.ProteinNumber))
                    self.Bits.append("*%s"%CheckGene.ProteinNumber)
                    self.Bits.append("")
                    self.Bits.append("")
                    return self.Category
            # One side of our intron is correct, but the other isn't (or belongs to another protein)
            self.Category = 6
            self.Bits.append("Cat6")
            self.Bits.append(GenomicSNPInfo) #bit22
            DonorStr = ""
            for DonorGene in DonorGenes:
                DonorStr += "%s:"%DonorGene.ProteinNumber
            if DonorStr:
                DonorStr = DonorStr[:-1] # remove trailing comma
            AcceptorStr = ""
            for AcceptorGene in AcceptorGenes:
                AcceptorStr += "%s:"%AcceptorGene.ProteinNumber
            if AcceptorStr:
                AcceptorStr = AcceptorStr[:-1] # remove trailing comma
            self.Bits.append(DonorStr) #bit23
            self.Bits.append(AcceptorStr) #bit24
            self.Bits.append("")
            self.Bits.append("")
        else:
            # NEITHER side of the splice junction is covered:
            if MismatchCount > 1:
                # Doesn't match genome, probably relies upon bogus SNPs (or is mis-annotated)
                self.Category = 8
                self.Bits.append("Cat8")
            else:
                self.Category = 7
                self.Bits.append("Cat7")
            self.Bits.append(GenomicSNPInfo)
            Str = ""
            if SpanningGenes:
                for Gene in SpanningGenes:
                    Str += "Span (%s-%s) %s "%(Gene.Start, Gene.End, Gene.ProteinNumber)
            else:
                if NearGeneLeft:
                    Str += "Near (%s-%s): %s "%(NearGeneLeft.Start, NearGeneLeft.End, NearGeneLeft.ProteinNumber)
                if NearGeneRight:
                    Str += "Near (%s-%s): %s "%(NearGeneRight.Start, NearGeneRight.End, NearGeneRight.ProteinNumber)
            self.Bits.append(Str)
            self.Bits.append("")
            self.Bits.append("")
            self.Bits.append("")
        return self.Category
    def CategorizeThreeExonHit(self):
        "Copied-and-pasted from C2EH due to laziness"
        #print "C3EH:", self.Bits[18], self.Bits[2]
        Junctions = self.Bits[20].split()
        (DonorStrA, AcceptorStrA) = Junctions[0].split("-")
        self.DonorA = int(DonorStrA)
        self.AcceptorA = int(AcceptorStrA)
        (DonorStrB, AcceptorStrB) = Junctions[1].split("-")
        self.DonorB = int(DonorStrB)
        self.AcceptorB = int(AcceptorStrB)
        NearGeneLeft = None
        NearGeneRight = None
        CoveringGenes = []
        DonorGenes = []
        AcceptorGenes = []
        SpanningGenes = []
        TerminusFlag = ""
        ##################################################
        # Get genomic SNP info:
        (MismatchCount, GenomicSNPInfo) = self.GetGenomicSNPInfo()
        #######################################################
        # Find the covering gene(s), and/or near neighbors:
        for Gene in Analyzer.GenesByChromosome[self.Chromosome]:
            if Gene.Strand != self.Strand:
                continue                        
            if Gene.End < self.Start:
                NearGeneLeft = Gene
                continue                
            if Gene.Start > self.End:
                NearGeneRight = Gene
                break
            CoverFlag = 0
            DonorFlagA = 0
            AcceptorFlagA = 0
            DonorFlagB = 0
            AcceptorFlagB = 0
            for ExonIndex in range(len(Gene.Exons) - 1):
                Exon = Gene.Exons[ExonIndex]
                NextExon = Gene.Exons[ExonIndex + 1]
                if self.Strand == 1:
                    if Exon.End == self.DonorA:
                        DonorGenes.append(Gene)
                        if Exon.Start <= self.Start:
                            DonorFlagA = 1
                    if Exon.End == self.DonorB:
                        DonorGenes.append(Gene)
                        DonorFlagB = 1
                    if NextExon.Start == self.AcceptorA:
                        AcceptorGenes.append(Gene)
                        AcceptorFlagA = 1
                    if NextExon.Start == self.AcceptorB:
                        AcceptorGenes.append(Gene)
                        if NextExon.End >= self.End:
                            AcceptorFlagB = 1
                else:
                    if NextExon.Start == self.DonorA:
                        DonorGenes.append(Gene)
                        if NextExon.End >= self.End:
                            DonorFlagA = 1
                    if NextExon.Start == self.DonorB:
                        DonorGenes.append(Gene)
                        DonorFlagB = 1
                    if Exon.End == self.AcceptorA:
                        AcceptorGenes.append(Gene)
                        AcceptorFlagA = 1
                    if Exon.End == self.AcceptorB:
                        AcceptorGenes.append(Gene)
                        if Exon.Start <= self.Start:
                            AcceptorFlagB = 1
                if DonorFlagA and AcceptorFlagA and DonorFlagB and AcceptorFlagB and (Gene not in CoveringGenes):
                    CoveringGenes.append(Gene)
                # Reset donor flag and acceptor flag.  To be cat4, you must be covered by
                # two consecutive exons - you can't skip exons from this isoform!
                DonorFlagA = 0
                DonorFlagB = 0
                AcceptorFlagA = 0
                AcceptorFlagB = 0
                    
            SpanningGenes.append(Gene)
        # if we covered the exon, then check to see whether it's an exact match:
        if CoveringGenes:
            ##################################################
            # Loop over covering genes to find the one with best sequence identity:
            BestIdentity = -1
            BestCoveringGene = None
            CoveringGeneStr = ""
            for TestCoveringGene in CoveringGenes:
                CoveringGeneStr += "%s:"%TestCoveringGene.ProteinNumber
                Pos = TestCoveringGene.Sequence.find(self.Aminos)
                if Pos != -1:
                    BestCoveringGene = TestCoveringGene
                    BestIdentity = 1
                    BestGeneSNPInfo = ""
                    if Pos == 0:
                        TerminusFlag = "N"
                    elif Pos + len(self.Aminos) == len(BestCoveringGene.Sequence):
                        TerminusFlag = "C"
                else:
                    (Identity, GeneSNPInfo) = Analyzer.GetMismatchInfo(self.Bits, TestCoveringGene, self.Aminos, TestCoveringGene.Sequence)
                    if (Identity > BestIdentity) or (Identity == BestIdentity and TestCoveringGene.CoveragePercent > BestCoveringGene.CoveragePercent):
                        BestIdentity = Identity
                        BestCoveringGene = TestCoveringGene
                        BestGeneSNPInfo = GeneSNPInfo
            CoveringGeneStr = CoveringGeneStr[:-1] # remove trailing comma
            if BestIdentity == 1:
                self.Category = 4
                self.Bits.append("Cat4")
                self.Bits.append(GenomicSNPInfo)
                self.Bits.append(CoveringGeneStr)
                self.Bits.append(str(BestCoveringGene.ProteinNumber))
                self.Bits.append(BestGeneSNPInfo)
                self.Bits.append(TerminusFlag)
            else:
                self.Category = 5
                self.Bits.append("Cat5")
                self.Bits.append(GenomicSNPInfo)
                self.Bits.append(CoveringGeneStr)
                self.Bits.append(str(BestCoveringGene.ProteinNumber))
                self.Bits.append(BestGeneSNPInfo)
                self.Bits.append("")
        elif (DonorGenes or AcceptorGenes):
            # One side of the splice junction is covered:
            CheckGeneList = DonorGenes[:]
            CheckGeneList.extend(AcceptorGenes)
            for CheckGene in CheckGeneList:
                if CheckGene.Sequence.find(self.Aminos)!=-1:
                    self.Category = 4
                    self.Bits.append("Cat4")
                    self.Bits.append(GenomicSNPInfo)
                    self.Bits.append(str(CheckGene.ProteinNumber))
                    self.Bits.append("*%s"%CheckGene.ProteinNumber)
                    self.Bits.append("")
                    self.Bits.append(TerminusFlag)
                    return self.Category
            self.Category = 6
            self.Bits.append("Cat6")
            self.Bits.append(GenomicSNPInfo)
            DonorStr = ""
            for DonorGene in DonorGenes:
                DonorStr += "%s:"%DonorGene.ProteinNumber
            if DonorStr:
                DonorStr = DonorStr[:-1] # remove trailing comma
            AcceptorStr = ""
            for AcceptorGene in AcceptorGenes:
                AcceptorStr += "%s:"%AcceptorGene.ProteinNumber
            if AcceptorStr:
                AcceptorStr = AcceptorStr[:-1] # remove trailing comma
            self.Bits.append(DonorStr)
            self.Bits.append(AcceptorStr)
            self.Bits.append("")
            self.Bits.append(TerminusFlag)
        else:
            # NEITHER side of the splice junction is covered:
            if MismatchCount > 1:
                # Doesn't match genome, probably relies upon bogus SNPs (or is mis-annotated)
                self.Category = 8
                self.Bits.append("Cat8")
            else:
                self.Category = 7
                self.Bits.append("Cat7")
            self.Bits.append(GenomicSNPInfo)
            Str = ""
            if SpanningGenes:
                for Gene in SpanningGenes:
                    Str += "Span (%s-%s) %s "%(Gene.Start, Gene.End, Gene.ProteinNumber)
            else:
                if NearGeneLeft:
                    Str += "Near (%s-%s): %s "%(NearGeneLeft.Start, NearGeneLeft.End, NearGeneLeft.ProteinNumber)
                if NearGeneRight:
                    Str += "Near (%s-%s): %s "%(NearGeneRight.Start, NearGeneRight.End, NearGeneRight.ProteinNumber)
            self.Bits.append(Str)
            self.Bits.append("")
            self.Bits.append("")
            self.Bits.append(TerminusFlag)
        return self.Category
    def GetString(self):
        return string.join(self.Bits, "\t")

class AnalyzerClass:
    "Main worker class for comparing exon-graph search results to known proteins"
    def __init__(self):
        # Track which results-files we've parsed already:
        self.AnalysisCheckpoints = {}
        self.PValueLimit = 0.08
        self.MinimumPeptideLength = 8
        self.DBPath = None
        self.IndexFilePath = None
        self.Genes = {} # ProteinNumber -> gene
        self.GenesByChromosome = [] # chromosome -> sorted list of genes
        self.OutputFiles = []
        Header = "File\tScan\tAnnotation\tGeneName\tCharge\tMQScore\tCutScore\tIntenseBY\tBYPresent\tNTT\tp-value\tDeltaCN\tDeltaCNOther\tRecordNumber\tDBFilePos\tSpecFilePos\tChromosome\tStrand\tGenomicPos\tSplicedSequence\tSplices\tCategory\tExtra1\tExtra2\tExtra3\tExtra4\tExtra5\tMultiHit\t\n"
        for Index in range(0, MAX_CATEGORY):
            FileName = "XGVKG%s.txt"%Index
            # Important point: Append, so we needn't re-run everything when new results come in.
            if os.path.exists(FileName):
                File = open(FileName, "a")
            else:
                File = open(FileName, "a")
                File.write(Header) # write the header only ONE time
            self.OutputFiles.append(File)
        self.OverallCounts = [0]*MAX_CATEGORY
        self.FileCounts = [0]*MAX_CATEGORY
        self.XMLDirectories = {}
        self.AllPeptides = {}
        self.KPPeptides = {}
        self.AllSpectrumCount = 0
        SummaryFileName = "SummaryXG.txt"
        if os.path.exists(SummaryFileName):
            self.SummaryFile = open(SummaryFileName, "a")
        else:
            self.SummaryFile = open(SummaryFileName, "a")
            self.SummaryFile.write("FileName\tLines\tSpectra\tCat0\tCat1\tCat2\tCat3\tCat4\tCat5\tCat6\tCat7\tCat8\tCat9\t\n")
        self.DiscoveryCurveFile = open("DiscoveryCurve.txt", "wb")
        self.KPDiscoveryCurveFile = open("KPDiscoveryCurve.txt", "wb")
        
    def PopulateXMLDirectories(self, RootDir):
        """
        Populate a dictionary of the form XMLFileName -> Containing directory
        """
        for FileName in os.listdir(RootDir):
            FilePath = os.path.join(RootDir, FileName)
            if os.path.isdir(FilePath):
                ##print "Recurse to %s"%FilePath
                self.PopulateXMLDirectories(FilePath)
            else:
                (Stub, Extension) = os.path.splitext(FilePath)
                if Extension.lower() == ".mzxml":
                    if self.XMLDirectories.has_key(FileName):
                        print "** Warning: Multiple entries %s and %s found for XML filename %s"%(RootDir, self.XMLDirectories[FileName], FileName)
                    ##else:
                    ##    print "Found %s in %s"%(FilePath, RootDir)
                    self.XMLDirectories[FileName] = RootDir
    def GetFixedFilePath(self, FilePath):
        """
        The search results have all sorts of crazy file-paths.  Correct them to point to
        the local disk.
        """
        Bits = FilePath.replace("/", "\\").split("\\")
        XMLFileName = Bits[-1]
        Directory = self.XMLDirectories.get(XMLFileName, None)
        if not Directory:
            print "* Warning: Directory not known for XML file %s"%XMLFileName
            return FilePath
        return os.path.join(Directory, XMLFileName)
    def SortChromosomeGenes(self):
        for X in range(49):
            self.GenesByChromosome.append([])
        for Gene in self.Genes.values():
            self.GenesByChromosome[Gene.ChromosomeNumber].append(Gene)
        for Index in range(1, len(self.GenesByChromosome)):
            List = self.GenesByChromosome[Index]
            print "Chromosome #%s has %s genes mapped onto it"%(Index, len(List))
            List.sort()
            if len(List) > 1:
                print "first gene %s-%s"%(List[0].Start,List[0].End)
                print "next gene %s-%s"%(List[1].Start,List[1].End)
                print "last gene %s-%s"%(List[-1].Start,List[-1].End)
    def LoadGeneMappings(self, GeneMappingFileName):
        GMFile = open(GeneMappingFileName, "rb")
        for FileLine in GMFile.xreadlines():
            Bits = FileLine.split("\t")
            if Bits[0] == "Protein":
                Gene = GeneClass()
                Gene.ProteinNumber = int(Bits[1])
                Gene.ProteinName = Bits[2] + " " + Bits[4]
            if Bits[0] == "Seed":
                Gene.Strand = int(Bits[4])
                Gene.ChromosomeNumber = int(Bits[3])
                Gene.CoveragePercent = float(Bits[10])
                ExonStr = Bits[7]
                for Chunk in ExonStr.split(","):
                    DashBits = Chunk.split("-")
                    if len(DashBits) < 4:
                        #print "??", DashBits
                        #print Bits
                        continue
                    Exon = ExonClass()
                    Exon.Start = int(DashBits[0])
                    Exon.End = int(DashBits[1])
                    Exon.ProteinStart = int(DashBits[2])
                    Exon.ProteinEnd = int(DashBits[3])
                    if len(DashBits) > 4:
                        Exon.EdgeAA = int(DashBits[4])
                    Gene.Exons.append(Exon)
                    if (Gene.Start == None or Gene.Start > Exon.Start):
                        Gene.Start = Exon.Start
                    if (Gene.End == None or Gene.End < Exon.End):
                        Gene.End = Exon.End
                self.Genes[Gene.ProteinNumber] = Gene
        GMFile.close()
    def LoadOldFormatGeneMappings(self, GeneMappingFile):
        File = open(GeneMappingFile, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                ProteinNumber = int(Bits[0])
                ChromosomeNumber = int(Bits[2])
                ExonCount = int(Bits[7])
            except:
                #traceback.print_exc()
                continue # header/footer line
            Gene = GeneClass()
            Gene.Strand = int(Bits[3])
            Gene.ProteinNumber = ProteinNumber
            Gene.ProteinName = Bits[1]
            Gene.ChromosomeNumber = ChromosomeNumber
            # Note how well the gene is covered - if coverage% is low, that may indicate
            # that the reference exons are bad:
            try:
                Gene.CoveragePercent = float(Bits[6]) / float(Bits[5])
            except:
                print FileLine
                raise
            ExonStr = Bits[8]
            if ExonStr and ExonStr[0] == '"':
                ExonStr = ExonStr[1:-1]
            for Chunk in ExonStr.split(","):
                DashBits = Chunk.split("-")
                if len(DashBits)<4:
                    #print "??", DashBits
                    #print Bits
                    continue
                Exon = ExonClass()
                Exon.Start = int(DashBits[0])
                Exon.End = int(DashBits[1])
                Exon.ProteinStart = int(DashBits[2])
                Exon.ProteinEnd = int(DashBits[3])
                if len(DashBits) > 4:
                    Exon.EdgeAA = int(DashBits[4])
                Gene.Exons.append(Exon)
                if (Gene.Start == None or Gene.Start > Exon.Start):
                    Gene.Start = Exon.Start
                if (Gene.End == None or Gene.End < Exon.End):
                    Gene.End = Exon.End
            self.Genes[Gene.ProteinNumber] = Gene
        File.close()
    def LoadSequences(self, DBPath, IndexPath):
        LastPos = -1
        #IndexFile = open(IndexPath, "rb")
        DBFile = open(DBPath, "rb")
        ProteinNumber = 0
        Data = DBFile.read()
        while (1):
            Pos = Data.find("*", LastPos + 1)
            if Pos == -1:
                if self.Genes.has_key(ProteinNumber):
                    self.Genes[ProteinNumber].Sequence = Data[LastPos + 1:]
                break
            if self.Genes.has_key(ProteinNumber):
                self.Genes[ProteinNumber].Sequence = Data[LastPos + 1:Pos]
            ProteinNumber += 1
            LastPos = Pos
        #IndexFile.close()
        DBFile.close()
    def AnalyzeResultsMain(self, DirList):
        FileNameList = []
        for Dir in DirList:
            for FileName in os.listdir(Dir):
                (Stub, Extension) = os.path.splitext(FileName)
                if Extension != RESULTS_EXTENSION:
                    continue
                Path = os.path.join(Dir, FileName)
                
                if not os.path.isdir(Path):
                    FileNameList.append(Path)
        FileNameList.sort()
        ResultFileCount = len(FileNameList)
        Index = 0
        for Path in FileNameList:
            Index += 1
            print
            print "** Read results from file %s of %s, %s"%(Index + 1, ResultFileCount, Path)
            self.AnalyzeResults(Path)
    def AnalyzeResults(self, FilePath):
        #print "\n\n\n%s"%FileName
        if os.path.isdir(FilePath):
            SubFileNames = os.listdir(FilePath)
            SubFileNames.sort()
            for SubFileName in SubFileNames:
                Path = os.path.join(FilePath, SubFileName)
                print "Analyze results from:", Path
                self.AnalyzeResults(Path) # recurse
        else:
            FileName = os.path.split(FilePath)[1]
            FileModTime = os.stat(FilePath).st_mtime
            XGInfo = self.AnalysisCheckpoints.get(FileName, None)
            if XGInfo:
                # We've tried parsing a results-file of this name before.
                # Check whether we are already up-to-date in analyzing this file:
                TimeDiff = FileModTime - XGInfo.ResultsTime
                if TimeDiff < 1:
                    print "File %s up-to-date, skipping"%(FileName)
                    return # We're up-to-date.
            else:
                # We've never seen this file before.  Create a new XGInfo object.
                XGInfo = XGFileInfoClass(FileName)
                XGInfo.ParseTime = time.time()
                XGInfo.ResultsTime = FileModTime
                self.AnalysisCheckpoints[FileName] = XGInfo
            print FilePath, FileName
            self.AnalyzeResultsFromFile(FilePath, XGInfo)
            self.SaveCheckpoint() # Save the checkpoint, just in case!
            #sys.exit(1)
    def AnalyzeResultsFromFile(self, FilePath, XGInfo):
        try:
            File = open(FilePath, "rb")
        except:
            print "** Error: Couldn't open the file '%s' for reading.  (Perhaps it's copying home from grid?)  Skipping it!"%FilePath
            return
        FileName = os.path.split(FilePath)[1]
        PendingHits = []
        PendingPeptides = {}
        OldSpectrum = None
        
        FixedFileName = None
        LineNumber = 0
        self.FileCounts = [0]*MAX_CATEGORY
        self.FileSpectrumCount = 0
        # Write a spacer-line to our discovery curve files:
        self.DiscoveryCurveFile.write("#Results file %s\n"%FileName)
        self.KPDiscoveryCurveFile.write("#Results file %s\n"%FileName)
        if NO_ORACLE:
            #DeltaScoreA = -0.25
            #DeltaScoreB = -0.5
            DeltaScoreA = -1.5
            DeltaScoreB = -2.5
        else:
            DeltaScoreA = -1.5
            DeltaScoreB = -2.5
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 100 == 0:
                print "%s %s..."%(FileName, LineNumber)
            Bits = FileLine.strip().split("\t")
            try:
                Spectrum = (Bits[0], Bits[1])
                DeltaScore = float(Bits[11])
                Peptide = Bits[2][2:-2]
                Score = float(Bits[5])
                PValue = float(Bits[10])
            except:
                continue
            if Spectrum != OldSpectrum:
                if OldSpectrum != None:
                    XGInfo.ScanCount += 1
                    if len(PendingHits):
                        self.OutputPendingHits(PendingHits)
                PendingHits = []
                PendingPeptides = {}
                OldSpectrum = Spectrum
                FixedFileName = self.GetFixedFilePath(Bits[0])
            if (Score < 10 and DeltaScore < DeltaScoreA) or (Score >= 10 and DeltaScore < DeltaScoreB):
                if not PendingPeptides.has_key(Peptide):
                    # Skip over this hit, because:
                    # - score is substantially lower than the top hit
                    # - and, the peptide sequence isn't the same as an 'acceptable' hit
                    continue
            # TEMP:
            # Some searches were run without a p-value calibration curve, and so the
            # p-value is *always* reported as zero.  We don't want to report those
            # results unless the MQScore is good!
            if (Score < 1 and PValue == 0):
                print "** Warning: P-value is 0, but score is %s"%(Score)
                continue
            try:                
                Hit = HitClass(Bits)
            except:
                continue # not a valid hit
            PendingPeptides[Peptide] = 1
            Hit.Bits[0] = FixedFileName
            Hit.Categorize()
            PendingHits.append(Hit)
        XGInfo.ScanCount += 1
        self.OutputPendingHits(PendingHits)
        #################
        # Summarize the hits of each type that went into this file:
        Str = "%s\t%s\t%s\t"%(FileName, LineNumber, self.FileSpectrumCount)
        for X in range(MAX_CATEGORY):
            Str += "%s\t"%self.FileCounts[X]
        self.SummaryFile.write(Str+"\n")
        self.SummaryFile.flush()
        #################
        print "Overall counts:", self.OverallCounts
        print "File counts:", self.FileCounts
    def GetGenomeSequence3(self, Chromosome, Strand, Start, End, DonorA, AcceptorA,
                           DonorB, AcceptorB):
        GenomeFileName = "e:\\chromosome\\chr%s.trie"%Chromosome
        try:
            GenomeFile = open(GenomeFileName, "rb")
            GenomeFile.seek(Start)
            if Strand == 1:
                DNA = GenomeFile.read(DonorA - Start)
                GenomeFile.seek(AcceptorA)
                DNA += GenomeFile.read(DonorB - AcceptorA)
                GenomeFile.seek(AcceptorB)
                DNA += GenomeFile.read(End - AcceptorB)
            else:            
                DNA = GenomeFile.read(AcceptorB - Start)
                GenomeFile.seek(DonorB)
                DNA += GenomeFile.read(AcceptorA - DonorB)
                GenomeFile.seek(DonorA)
                DNA += GenomeFile.read(End - DonorA)
            GenomeFile.close()
            if Strand == -1:
                DNA = GeneMapper.ReverseComplement(DNA)
            GenomicPeptide = GeneMapper.Translate(DNA)
        except:
            traceback.print_exc()
            return "??????"
        return GenomicPeptide
    def GetGenomeSequence(self, Chromosome, Strand, Start, End,
            Donor = None, Acceptor = None):
        GenomeFileName = "e:\\chromosome\\chr%s.trie"%Chromosome
        try:
            GenomeFile = open(GenomeFileName, "rb")
            GenomeFile.seek(Start)
            if Donor != None:
                DNA = GenomeFile.read(min(Donor, Acceptor) - Start)
                GenomeFile.seek(max(Donor, Acceptor))
                DNA += GenomeFile.read(End - max(Donor, Acceptor))
            else:            
                DNA = GenomeFile.read(End - Start)
            GenomeFile.close()
            if Strand == -1:
                DNA = GeneMapper.ReverseComplement(DNA)
            GenomicPeptide = GeneMapper.Translate(DNA)
        except:
            traceback.print_exc()
            return "??????"
        return GenomicPeptide
    def GetMismatchInfo(self, Bits, Gene, Peptide, Protein):
        (Start, End) = Bits[18].split("-")
        Start = int(Start)
        End = int(End)
        DashCount = Bits[20].count("-")
        Strand = int(Bits[17])
        if DashCount > 1:
            # UGH.  Let's not read the protein sequence across all these splice junctions,
            # let's just find the closest match:
            MaxStartPos = len(Protein) - len(Peptide)
            BestMatchCount = 0
            BestSNPInfo = ""
            for StartPos in range(MaxStartPos):
                MatchCount = 0
                SNPInfo = ""
                for X in range(len(Peptide)):
                    if Protein[StartPos + X] == Peptide[X]:
                        MatchCount += 1
                        SNPInfo += Protein[StartPos + X]
                    else:
                        SNPInfo += Protein[StartPos + X].lower()
                    if MatchCount > BestMatchCount:
                        BestMatchCount = MatchCount
                        BestSNPInfo = SNPInfo
            Identity = BestMatchCount / float(len(Peptide))
            return (Identity, BestSNPInfo)
        elif DashCount == 1:
            (Donor, Acceptor) = Bits[20].split("-")
            Donor = int(Donor)
            Acceptor = int(Acceptor)
        else:
            Donor = None
            Acceptor = None
        ##################################
        # Get the IPI reference sequence:
        IPISequence = "???"
        if Donor:
            IPISequenceA = ""
            IPISequenceB = ""
            EdgeAA = ""
            for Exon in Gene.Exons:
                if Exon.Start <= Start and Exon.End >= min(Donor, Acceptor):
                    IPIEdge = ""
                    if Gene.Strand == 1:
                        ProteinStart = Exon.ProteinStart + (Start - Exon.Start) / 3
                        ProteinEnd = ProteinStart + (min(Donor, Acceptor) - Start) / 3
                        IPISequenceA = Gene.Sequence[ProteinStart:ProteinEnd]
                        if Exon.EdgeAA and (Exon.EdgeAA < ProteinStart or Exon.EdgeAA >= ProteinEnd):
                            EdgeAA = Gene.Sequence[Exon.EdgeAA]
                    else:
                        ProteinStart = Exon.ProteinStart + (Exon.End - min(Donor, Acceptor)) / 3
                        ProteinEnd = ProteinStart + (min(Donor, Acceptor) - Start) / 3
                        IPISequenceB = Gene.Sequence[ProteinStart:ProteinEnd]
                        if Exon.EdgeAA and (Exon.EdgeAA < ProteinStart or Exon.EdgeAA >= ProteinEnd):
                            EdgeAA = Gene.Sequence[Exon.EdgeAA]
                if Exon.Start <= max(Donor, Acceptor) and Exon.End >= End:
                    if Gene.Strand == 1:
                        ProteinStart = Exon.ProteinStart + (max(Donor, Acceptor) - Exon.Start) / 3
                        ProteinEnd = ProteinStart + (End - max(Donor, Acceptor)) / 3
                        IPISequenceB = Gene.Sequence[ProteinStart:ProteinEnd]
                    else:
                        ProteinStart = Exon.ProteinStart + (Exon.End - End) / 3
                        ProteinEnd = ProteinStart + (End - max(Donor, Acceptor)) / 3
                        IPISequenceA = Gene.Sequence[ProteinStart:ProteinEnd]
            IPISequence = IPISequenceA + EdgeAA + IPISequenceB
            #print "%s-%s-%s-%s"%(Start, min(Donor, Acceptor), max(Donor, Acceptor), End)
            #print "IPISequence: %s + %s + %s -> %s"%(IPISequenceA, EdgeAA, IPISequenceB, IPISequence)
        else:
            for Exon in Gene.Exons:
                if Exon.Start <= Start and Exon.End >= End:
                    if Gene.Strand == 1:
                        ProteinStart = Exon.ProteinStart + (Start - Exon.Start) / 3
                        ProteinEnd = Exon.ProteinStart + (End - Exon.Start) / 3
                    else:
                        ProteinStart = Exon.ProteinStart + (Exon.End - End) / 3
                        ProteinEnd = Exon.ProteinStart + (Exon.End - Start) / 3
                    IPISequence = Gene.Sequence[ProteinStart:ProteinEnd]
        MatchCount = 0
        MismatchCount = 0
        NiceIPISequence = ""
        Len = min(len(IPISequence), len(Peptide))
        SNPPepString = ""
        SNPString = ""
        if Strand == 1:
            GenomePos = Start
        else:
            GenomePos = End - 1
        for Index in range(Len):
            if IPISequence[Index] == Peptide[Index]:
                MatchCount += 1
                SNPPepString += IPISequence[Index]
            else:
                MismatchCount += 1
                SNPPepString += IPISequence[Index].lower()
                SNPString += "%s (%s/%s),"%(GenomePos, Peptide[Index], IPISequence[Index])
            if Strand == 1:
                GenomePos += 3
                if (Donor != None and GenomePos >= Donor and GenomePos < Acceptor):
                    GenomePos += (Acceptor - Donor)
            else:
                GenomePos -= 3
                if (Acceptor != None and GenomePos <= Donor and GenomePos > Acceptor):
                    GenomePos += (Acceptor - Donor)
        #print "%s %s-%s %s %s"%(Bits[3], Start, End, Donor, Acceptor)
        #print "%s vs %s: %s match %s mismatch"%(Peptide, IPISequence, MatchCount, MismatchCount)
        if len(IPISequence) != len(Peptide):
            print "*Length (%s %s)"%(IPISequence, Peptide)
            #sys.stdin.readline()
        if (MatchCount + MismatchCount == 0):
            Identity = 0
            SNPString = "ERROR"
        else:
            Identity = MatchCount / float(MatchCount + MismatchCount)
            if MismatchCount == 0:
                SNPString = ""
            elif Identity < 0.75:
                if len(IPISequence) != len(Peptide):
                    SNPString = "*len (%s %s)"%(IPISequence, Peptide) # out of frame
                else:
                    SNPString = "*oof (%s %s)"%(IPISequence, Peptide) # out of frame
                #sys.stdin.readline()
            else:
                # There's one mismatch, but we're in frame:
                SNPString += " %s"%SNPPepString
        return (Identity, SNPString)
        #return "g:%s p:%s ipi:%s"%(GenomicPeptide, NicePeptide, NiceIPISequence)
    def OutputPendingHits(self, PendingHits):
        self.FileSpectrumCount += 1
        self.AllSpectrumCount += 1
        BestRank = 0
        HitsForBestRank = 0
        for Hit in PendingHits:
            Rank = Hit.GetCategoryRank()
            if Rank > BestRank:
                HitsForBestRank = 1
                BestRank = Rank
            elif Rank == BestRank:
                HitsForBestRank += 1
        if BestRank == 0:
            return # nothing worth reporting!
        FirstFlag = 1
        for Hit in PendingHits:
            if Hit.CategoryRank < BestRank:
                continue
            if FirstFlag:
                self.OverallCounts[Hit.Category] += 1
                self.FileCounts[Hit.Category] += 1
                self.FileCounts[0] += 1
                self.OverallCounts[0] += 1
                FirstFlag = 0
                if not self.AllPeptides.has_key(Hit.Aminos):
                    self.AllPeptides[Hit.Aminos] = 1
                    Str = "%s\t%s\t%s\t%s\t\n"%(self.AllSpectrumCount, self.OverallCounts[0], len(self.AllPeptides.keys()), Hit.Aminos)
                    self.DiscoveryCurveFile.write(Str)
                if Hit.Category in (1, 4) and not self.KPPeptides.has_key(Hit.Aminos):
                    self.KPPeptides[Hit.Aminos] = 1
                    Str = "%s\t%s\t%s\t%s\t\n"%(self.AllSpectrumCount, self.OverallCounts[1] + self.OverallCounts[4],
                        len(self.KPPeptides.keys()), Hit.Aminos)
                    self.KPDiscoveryCurveFile.write(Str)
            if HitsForBestRank > 1:
                Hit.Bits.append("Multi:%s"%HitsForBestRank)
            # Write to the output file for this category:
            OutputFile = self.OutputFiles[Hit.Category]
            OutputFile.write(Hit.GetString())
            OutputFile.write("\n")
            # Also write to the overall merged output file:
            OutputFile = self.OutputFiles[0]
            OutputFile.write(Hit.GetString())
            OutputFile.write("\n")
            
    def SaveCheckpoint(self):
        """
        Read our analysis-status: What mzxml files have we read results for?  When were
        those results modified, and when did we generate the results?
        """
        print "Save XGCheckpoint..."
        File = open("XGCheckpoint.txt", "wb")
        Keys = self.AnalysisCheckpoints.keys()
        Keys.sort()
        for Key in Keys:
            XGInfo = self.AnalysisCheckpoints[Key]
            Str = "%s\t%s\t%s\t%s\t"%(XGInfo.FileName, XGInfo.ResultsTime, XGInfo.ParseTime, XGInfo.ScanCount)
            File.write(Str + "\n")
        File.close()
    def LoadCheckpoint(self):
        print "Load XGCheckpoint..."
        try:
            File = open("XGCheckpoint.txt", "rb")
        except:
            print "XGCheckpoint not found - starting analysis fresh, I hope you're happy."
            return
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if len(Bits) < 3:
                # Assume it's a blank or header line.
                continue
            Info = XGFileInfoClass(Bits[0])
            Info.ResultsTime = float(Bits[1])
            Info.ParseTime = float(Bits[2])
            Info.ScanCount = int(Bits[3])
            self.AnalysisCheckpoints[Bits[0]] = Info
        File.close()
        
class XGFileInfoClass:
    def __init__(self, FileName):
        self.FileName = FileName # results file name, no path
        self.ResultsTime = 0
        self.ParseTime = 0
        self.ScanCount = 0
        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running non-optimized)"
    Analyzer = AnalyzerClass()
    Analyzer.LoadCheckpoint()
    Analyzer.DBPath = os.path.join("Database", "IPIv315.trie")
    Analyzer.IndexFilePath = os.path.join("Database", "IPIv315.index")
    print "PopulateXMLDirectories:"
    Analyzer.PopulateXMLDirectories("e:\\ms\\PeptideAtlas")
    Analyzer.PopulateXMLDirectories("e:\\ms\\Briggs")
    print "Load gene mappings..."
    #Analyzer.LoadGeneMappings("GeneMappings.IPI315.txt")
    #Analyzer.LoadGeneMappings("GeneMapperOutput.10.txt")
    #Analyzer.LoadGeneMappings("GeneMapperOutput\\GeneMappings.F.txt")
    Analyzer.LoadGeneMappings("GeneMappings.Best.txt")
    print "Sort chromosome genes..."
    Analyzer.SortChromosomeGenes()
    print "Load protein sequences..."
    Analyzer.LoadSequences(Analyzer.DBPath, Analyzer.IndexFilePath)
    print "Analyze results..."
    #Analyzer.AnalyzeResults(r"E:\ms\PeptideAtlas\ResultsX\xp091003_08.txt")
    #Dirs = ["e:\\ms\\Briggs\\ResultsXX", r"E:\ms\PeptideAtlas\ResultsXX"]
    Dirs = ["e:\\ms\\Briggs\\ResultsXFixed", r"E:\ms\PeptideAtlas\ResultsXFixed"]
    Analyzer.AnalyzeResultsMain(Dirs)
    Analyzer.SaveCheckpoint()
    #Analyzer.AnalyzeResults("patest.out")