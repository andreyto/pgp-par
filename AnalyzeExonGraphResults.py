"""
Analysis of the exon graph search results:
- First, run ExonGraphVsKnownGenes, to categorize spectra according to whether
they are explainable by "easy" exons (known genes, matching sequence) or "hard"
exons (unknown genomic intervals and splice junctions).
- Next, run this script (AnalyzeExonGraphResults) to examine the results in more
detail.  Run for categories 1+4, 2+5, or others.
"""
import os
import sys
import traceback
import GeneBrowserPlot
from Utils import *
Initialize()

XGVKG0 = "XGWithIPIFilter\\XGVKG0.txt"
#XGVKG0 = "XGNoIPIFilter\\XGVKG0.txt"
#XGVKG0 = "ENCODEXG0.txt"
#XGVKG0 = "Tempy.txt"

# We want to see the coverage of these guys:
JuicyProteins = [29924,16715,8259,1967,9222,34935,40142,6760,17683,32091]

ShowNovelPeptidesForGenes = [18448, 17173, 17096, 13956, 2615, 239420, 252232,
                             249900, 249900, 249928, 254305, 254339, 253530,
                             278118, 276531, 278831, 276029, 273093, 273093,
                             270577, 289702, 289702, 292891, 302832, 313732,
                             346115, 348455, 365825, 371269, 370090, 366734,
                             366070, 43961, 39163, 49614, 394278, 395338,
                             396335, 403554, 403226, 71106, 83716, 103017,
                             102850, 126795, 129632, 125373, 129129, 125202, 
                             
                             130416, 124214, 148140, 148140, 148140, 175374, 
                             179298, 167564, 174670, 166968, 193780,
                             193154, 211913,

                            # And now, the reverse strand!
                            21039, 30117, 29565, 247235, 267490, 263340, 265997, 262805, 282281,
                             279548, 286670, 279829, 297143, 305768, 305768, 334884,
                             334076, 349921, 356464, 356464, 349736, 356464, 356285, 351138,
                             349736, 373284, 373284, 373284, 374824, 376017, 377680, 374366, 
                             376786, 377680, 377680, 374942, 372977, 70252, 60294, 59942, 69335, 
                             65403, 69040, 54826, 399694, 400098, 400106, 415583, 93340, 
                             88661, 85846, 112396, 112396, 112396, 112396, 110226, 422818, 137221,
                             140666, 134541, 159205, 159205, 159205, 159205, 159205, 
                             159205, 159205, 159205, 159205, 159205, 159205, 159205, 159205, 159205,
                             159205, 159205, 159205, 159205, 159205, 159205, 166740,
                             166095, 161007, 157272, 161007, 190656, 190656, 181392, 187970, 187970,
                             187970, 205518, 205518, 205518, 210424, 208755
                             ]

                             
                             

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
try:
    TheFont = ImageFont.truetype("Times.ttf", 12)
except:
    TheFont = ImageFont.load_default()

MAX_PM_LINE_NUMBER = 0 # for debugging: bail out during parse for speed
MAX_LINE_NUMBER = 0 # for debugging: bail out during parse for speed

ReadingFrameColors = [(200, 55, 55), (55, 200, 55), (55, 55, 200)]

class Colors:
    "Bag of constants specifying a color scheme"
    Background = (255, 255, 255)
    Exon = (155, 155, 155)
    SNP = (55, 55, 55)
    Mismatch = (255, 100, 100)
    Peptide = (200, 0, 0)
    PeptideSpliceGhost = (220, 155, 155)
    ESTWrong = (100, 255, 100)
    LinkCorrect = (200, 0, 0)
    LinkWrong = (125, 125, 255)
    GF = (0, 0, 100)
    Label = (0, 0, 0)
    Border = (200, 200, 200)
    SpanningEdge = (0, 155, 155)

class ExonClass:
    def __init__(self):
        self.HitCount = 0
        self.DonorHitCount = 0
        self.AcceptorHitCount = 0

class ProteinClass:
    "A known protein."
    def __init__(self):
        self.Exons = [] # instances of ExonClass
        self.ChromosomeNumber = None
        self.Strand = 1 # 1 or -1
        self.Start = None
        self.End = None
        self.ProteinNumber = None
        self.ProteinName = None
        self.Sequence = None
        self.HitCount = 0
        self.SharedHitCount = 0
        self.AcceptedFlag = 0
        self.Peptides = {}
        self.SharedPeptides = {}
        self.HitDetails = {}
        self.DetailPeptides = {}
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



class GeneClass:
    "An exon-graph gene."
    def __init__(self):
        self.Name = ""
        self.ChromosomeNumber = None
        self.Strand = None
        self.HitDetails = {}
        self.DetailPeptides = {}
        self.PeptideCategories = {}
        self.GeneNumber = None
        self.PeptideBestScores = {}
        
class AnalyzerClass:
    def __init__(self):
        self.Proteins = {} # ProteinNumber -> ProteinClass instance
        self.Genes = {} # exon graph recordnumber -> GeneClass instance
        self.ProteinsByChromosome = [] # chromosome -> sorted list of Proteins
        # Demand at least this many category 1 and 4 hits to a known protein
        # before we make the call that it's present:
        self.MinSpectraForKnownProtein = 2
        self.MinSpectraForXGene = 15
        # Count the peptides, and spectra, accounted for by
        # known proteins.
        self.SpectraForProtein = {}
        self.PeptidesForProtein = {} 
        # Count the peptides, and spectra, accounted for by genes
        # in the exon graph.  
        self.PeptidesForGene = {}
        self.SpectraForGene = {}
        self.UnexplainedGeneHits = {}
        self.InterestingGenes = []
        self.InterestingGeneDict = {}
        self.DrawPresentProteinImages = 0
    def SortChromosomeProteins(self):
        for X in range(49):
            self.ProteinsByChromosome.append([])
        for Protein in self.Proteins.values():
            self.ProteinsByChromosome[Protein.ChromosomeNumber].append(Protein)
        for Index in range(1, len(self.ProteinsByChromosome)):
            List = self.ProteinsByChromosome[Index]
            print "Chromosome #%s has %s Proteins mapped onto it"%(Index, len(List))
            List.sort()
            if len(List) > 1:
                print "first Protein %s-%s"%(List[0].Start,List[0].End)
                print "next Protein %s-%s"%(List[1].Start,List[1].End)
                print "last Protein %s-%s"%(List[-1].Start,List[-1].End)
    def LoadProteinMappings(self, ProteinMappingFile):
        File = open(ProteinMappingFile, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if MAX_PM_LINE_NUMBER > 0 and LineNumber > MAX_PM_LINE_NUMBER:
                break #bail out quickly
            Bits = FileLine.split("\t")
            if Bits[0] == "Protein":
                Protein = ProteinClass()
                Protein.ProteinNumber = int(Bits[1])
                Protein.ProteinName = Bits[2]
                Protein.Description = Bits[4]
            if Bits[0] == "Seed":
                Protein.ChromosomeNumber = int(Bits[3])
                Protein.CoveragePercent = float(Bits[10])
                Strand = int(Bits[4])
                # The file holds a "reverse-flag", we want 1 for forward and -1 for reverse.
                if Strand == 1:
                    Protein.Strand = -1
                else:
                    Protein.Strand = 1
                ExonStr = Bits[7]
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
                    Protein.Exons.append(Exon)
                    if (Protein.Start == None or Protein.Start > Exon.Start):
                        Protein.Start = Exon.Start
                    if (Protein.End == None or Protein.End < Exon.End):
                        Protein.End = Exon.End
                # we'll only keep proteins with at least a chromosome#
                self.Proteins[Protein.ProteinNumber] = Protein
        File.close()
    def LoadOldProteinMappings(self, ProteinMappingFile):
        File = open(ProteinMappingFile, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if MAX_PM_LINE_NUMBER > 0 and LineNumber > MAX_PM_LINE_NUMBER:
                break #bail out quickly
            Bits = FileLine.split("\t")
            try:
                ProteinNumber = int(Bits[0])
                ChromosomeNumber = int(Bits[2])
                ExonCount = int(Bits[7])
            except:
                #traceback.print_exc()
                continue # header/footer line
            Protein = ProteinClass()
            Protein.Strand = int(Bits[3])
            Protein.ProteinNumber = ProteinNumber
            Protein.ProteinName = Bits[1]
            if Protein.ProteinName[0] == '"':
                Protein.ProteinName = Protein.ProteinName[1:-1]
            Protein.ChromosomeNumber = ChromosomeNumber
            # Note how well the Protein is covered - if coverage% is low, that may indicate
            # that the reference exons are bad:
            try:
                Protein.CoveragePercent = float(Bits[6]) / float(Bits[5])
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
                Protein.Exons.append(Exon)
                if (Protein.Start == None or Protein.Start > Exon.Start):
                    Protein.Start = Exon.Start
                if (Protein.End == None or Protein.End < Exon.End):
                    Protein.End = Exon.End
            self.Proteins[Protein.ProteinNumber] = Protein
        File.close()
    def ReportPresentProteins(self):
        """
        Called after SelectPresentProteins.
        For each known protein that was judged "present", report the number
        of spectra hits and shared-hits, and the number of peptide hits and shared-hits.
        Goal:
        Populate members of Protein objects: HitCount, SharedHitCount, Peptides, SharedPeptides
        (optional) generate "coverage" image for the protein
        """
        for Protein in self.Proteins.values():
            Protein.HitCount = 0
            Protein.SharedHitCount = 0
            Protein.CoverageFlags = [0]*len(Protein.Sequence)
        File = open(XGVKG0, "rb")
        OldSpectrum = None
        LineNumber = 0
        CurrentSpectrumHits = []
        # Iterate over file lines.  Once we get all the hits for a spectrum,
        # update the hits (and maybe shared-hits) of any acceptable proteins.
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "%s..."%(LineNumber)
                if MAX_LINE_NUMBER > 0 and LineNumber > MAX_LINE_NUMBER:
                    break
            Bits = FileLine.strip().split("\t")
            Spectrum = (Bits[0], Bits[1])
            if Spectrum != OldSpectrum:
                if len(CurrentSpectrumHits):
                    self.CountSpectrumHits(CurrentSpectrumHits)
                OldSpectrum = Spectrum
                CurrentSpectrumHits = []
            # Skip over any lines not in category 1 or 4:
            if len(Bits)<22 or Bits[21] not in ("Cat1", "Cat4"):
                continue
            if Bits[24][0] == "*":
                continue
            try:
                ProteinNumber = int(Bits[23])
                ProteinNumbers = [ProteinNumber]
            except:
                try:
                    ProteinNumbers = Bits[23].split(":")
                    ProteinNumbers = map(int, ProteinNumbers)
                except:
                    continue # header line
            for ProteinNumber in ProteinNumbers:
                # Remember any "acceptable" proteins, for sending to CountSpectrumHits.
                if ProteinNumber in self.AcceptedProteins:
                    Info = (self.AcceptedProteins.index(ProteinNumber), ProteinNumber, Bits)
                    CurrentSpectrumHits.append(Info)
        if len(CurrentSpectrumHits):
            self.CountSpectrumHits(CurrentSpectrumHits)
        File.close()
        ####################################################################
        # Output the list of proteins judged present due to cat1 and cat4 results
        PeptideFrequencyDistribution = {}
        File = open("ReportPresentProteins.txt", "wb")
        for AcceptedProteinIndex in range(len(self.AcceptedProteins)):
            ProteinNumber = self.AcceptedProteins[AcceptedProteinIndex]
            Protein = self.Proteins[ProteinNumber]
            TotalHitCount = 0
            PeptideCount = 0
            SortedPeptideList = []
            for (Key, Value) in Protein.HitDetails.items():
                TotalHitCount += Value
                PeptideCount += 1
                SortedPeptideList.append((Value, Key))
                PeptideFrequencyDistribution[Value] = PeptideFrequencyDistribution.get(Value, 0) + 1
            Str = "Protein\t%s\t%s\t%s\t"%(ProteinNumber, TotalHitCount, PeptideCount)
            #########################
            # Report how many INTRONS and how many EXONS were hit:
            IntronHitCount = 0
            ExonHitCount = 0
            for Exon in Protein.Exons:
                if Exon.HitCount:
                    ExonHitCount += 1
                if Exon.DonorHitCount:
                    IntronHitCount += 1
            ExonCount = len(Protein.Exons)
            IntronCount = max(0, len(Protein.Exons) - 1)
            Str += "%s\t%s\t%s\t%s\t"%(ExonHitCount, IntronHitCount, ExonCount, IntronCount)
            #########################
            # Report our overlap with other proteins that have already been accepted:
            if TotalHitCount > 10:
                for PriorProteinIndex in range(AcceptedProteinIndex):
                    OtherProtein = self.Proteins[self.AcceptedProteins[PriorProteinIndex]]
                    if Protein.ChromosomeNumber == OtherProtein.ChromosomeNumber and Protein.Strand == OtherProtein.Strand:
                        if Protein.Exons[0].Start < OtherProtein.Exons[-1].End and Protein.Exons[-1].End > OtherProtein.Exons[0].Start:
                            Str += "%s\t"%OtherProtein.ProteinNumber
            File.write(Str + "\n")
            # Temp: don't report the peptide list, we don't want to overflow excel row-count with
            # this report:
##            SortedPeptideList.sort()
##            SortedPeptideList.reverse()
##            for (Count, Key) in SortedPeptideList:
##                Str = "%s\t%s\t%s\t\n"%(Count, Protein.DetailPeptides[Key], Key)
##                File.write(Str)
        File.close()
        ########################################################################
        # Report all peptides, see if they follow a Pareto distribution:
        print "Report peptide frequency distribution..."
        File = open("PeptideFrequencyDistribution.txt", "wb")
        Keys = PeptideFrequencyDistribution.keys()
        Keys.sort()
        for Key in Keys:
            File.write("%s\t%s\t\n"%(Key, PeptideFrequencyDistribution[Key]))
        File.close()
        #######################################################
        # Protein coverage image, if requested.
        if self.DrawPresentProteinImages:
            for ProteinNumber in self.AcceptedProteins:
                Protein = self.Proteins[ProteinNumber]
                self.GenerateProteinImage(Protein)
        for ProteinNumber in self.AcceptedProteins:
            if ProteinNumber in JuicyProteins:
                Protein = self.Proteins[ProteinNumber]
                self.GenerateProteinImage(Protein)                
    def CountUnexplainedSpectrumHits(self, SpectrumHits):
        """
        Helper for ReportUnexplainedGenes: Count the hits for the
        "most interesting" gene.
        """
        SpectrumHits.sort()
        GeneNumber = self.InterestingGenes[SpectrumHits[0][0]]
        Bits = SpectrumHits[0][1]
        Gene = self.Genes[GeneNumber]
        Peptide = Bits[2][2:-2]
        (Start, End) = Bits[18].split("-")
        Start = int(Start)
        End = int(End)
        if Bits[20]:
            # spliced!
            DashBillions = Bits[20].split("-")
            if len(DashBillions) > 2:
                return # skip three-exon hit!
            (Donor, Acceptor) = DashBillions
            Donor = int(Donor)
            Acceptor = int(Acceptor)
            Key = (Start, min(Donor, Acceptor), max(Donor, Acceptor), End)
        else:
            Key = (Start, End)
        Gene.HitDetails[Key] = Gene.HitDetails.get(Key, 0) + 1
        Gene.DetailPeptides[Key] = "%s.%s.%s"%(Bits[2][0], Bits[19], Bits[2][-1])  #Bits[2]
        Gene.PeptideCategories[Key] = Bits[21]
        Gene.PeptideBestScores[Key] = max(float(Bits[5]), Gene.PeptideBestScores.get(Key, -99))
        
        #################################################################
        # Append these bits to the spectrum-hit file for the gene:
        GeneFile = open("NovelPeptides\\%s.txt"%GeneNumber, "a")
        Str = string.join(Bits, "\t")
        GeneFile.write(Str)
        GeneFile.write("\n")
        GeneFile.close()
    def ReportProteinCoverage(self):
        ProteinCoverageFile = open("ProteinCoverage.txt", "wb")
        for (ProteinID, Protein) in self.Proteins.items():
            Length = len(Protein.Sequence)
            Coverage = 0
            Coverage2 = 0
            for Int in Protein.CoverageFlags:
                if Int:
                    Coverage += 1
                if Int > 1:
                    Coverage2 += 1
            CoveragePercent = Coverage / float(max(1, Length))
            CoveragePercent2 = Coverage2 / float(max(1, Length))
            Str = "%s\t%s\t%s\t%s\t%s\t%s\t"%(ProteinID, Length, Coverage, CoveragePercent, Coverage2, CoveragePercent2)
            ProteinCoverageFile.write(Str + "\n")
        ProteinCoverageFile.close() 
    def CountSpectrumHits(self, SpectrumHits):
        """
        Helper for ReportPresentProteins: Count the hits and shared-hits for the
        known protein, and for the exon graph gene.
        """
        SpectrumHits.sort()
        # The FIRST protein on the list gets hits, any OTHER proteins in this list
        # get shared-hits.
        MatchedFlag = 0
        for MatchIndex in range(len(SpectrumHits)):
            Bits = SpectrumHits[MatchIndex][2]
            ProteinNumber = SpectrumHits[MatchIndex][1]
            GeneNumber = int(Bits[13])
            Protein = self.Proteins[ProteinNumber]
            Peptide = Bits[2][2:-2]
            (Start, End) = Bits[18].split("-")
            Start = int(Start)
            End = int(End)
            if (MatchedFlag):
                # Another protein containing the same peptide:
                Protein.SharedHitCount += 1
                Protein.SharedPeptides[Peptide] = Protein.SharedPeptides.get(Peptide, 0) + 1
                self.SpectraForGene[GeneNumber] = self.SpectraForGene.get(GeneNumber, 0) + 1
                if not self.PeptidesForGene.has_key(GeneNumber):
                    self.PeptidesForGene[GeneNumber] = {}
                self.PeptidesForGene[GeneNumber][Peptide] = self.PeptidesForGene[GeneNumber].get(Peptide, 0) + 1
                if Peptide == TopHitPeptide:
                    # It's the same peptide, so let's flag coverage:
                    Pos = Protein.Sequence.find(Peptide)
                    if Pos != -1:
                        for ResidueNumber in range(Pos, min(len(Protein.CoverageFlags), Pos + len(Peptide))):
                            Protein.CoverageFlags[ResidueNumber] += 1
                            #Protein.CoverageFlags[ResidueNumber] = "x"
            else:
                Pos = Protein.Sequence.find(Peptide)                
                if Pos != -1:
                    MatchedFlag = 1
                    for ResidueNumber in range(Pos, min(len(Protein.CoverageFlags), Pos + len(Peptide))):
                        Protein.CoverageFlags[ResidueNumber] += 1
                        #if Protein.CoverageFlags[ResidueNumber] != "X":
                        #    Protein.CoverageFlags = Protein.CoverageFlags[:ResidueNumber] + "X" + Protein.CoverageFlags[ResidueNumber + 1:]
                        #Protein.CoverageFlags[ResidueNumber] = "x"
                    # The CORRECT protein:
                    Protein.HitCount += 1
                    Protein.Peptides[Peptide] = Protein.Peptides.get(Peptide, 0) + 1
                    TrueStart = Start
                    TrueEnd = End
                    TopHitPeptide = Peptide
                    Pos = Protein.Sequence.find(Peptide)
                    if Bits[20]:
                        # spliced!
                        DashBillions = Bits[20].split("-")
                        if len(DashBillions) > 2:
                            continue # three-exons, don't bother with it!
                        (Donor, Acceptor) = DashBillions
                        Donor = int(Donor)
                        Acceptor = int(Acceptor)
                        for Exon in Protein.Exons:
                            if Donor in (Exon.Start, Exon.End):
                                Exon.HitCount += 1
                                Exon.DonorHitCount += 1
                            if Acceptor in (Exon.Start, Exon.End):
                                Exon.HitCount += 1
                                Exon.AcceptorHitCount += 1
                        Key = (Start, min(Donor, Acceptor), max(Donor, Acceptor), End)
                        Protein.HitDetails[Key] = Protein.HitDetails.get(Key, 0) + 1
                        Protein.DetailPeptides[Key] = Bits[2]
                    else:
                        # Not spliced!
                        Key = (Start, End)
                        Protein.HitDetails[Key] = Protein.HitDetails.get(Key, 0) + 1
                        Protein.DetailPeptides[Key] = Bits[2]
                        for Exon in Protein.Exons:
                            if Exon.Start <= Start and Exon.End >= End:
                                Exon.HitCount += 1        
                else:
                    print "WARNING: Peptide '%s' not found in protein %s %s"%(Peptide, Protein.ProteinName, Protein.ProteinNumber)
    def ReportUnexplainedGenes(self):
        """
        Called after FindUnexplainedGenes has populated self.InterestingGenes.
        Read all the non-perfect-match hits.  For each spectrum, remember the
        hit to the first "interesting" gene.  Once we've processed all the
        spectrum hits, we're ready to report the unexplained peptides
        for each gene record.
        """
        InterestingGeneFile = open("PutativeNovelExons.txt", "wb")
        Header = "Gene\tName\tKnownSpectra\tKnownPeptides\ttNewSpectra\tNewPeptides\tCategory\tNotes\n"
        InterestingGeneFile.write(Header)
        File = open(XGVKG0, "rb")
        OldSpectrum = None
        LineNumber = 0
        CurrentSpectrumHits = []
        # Iterate over file lines.  Once we get all the hits for a spectrum,
        # update the hits (and maybe shared-hits) of any acceptable proteins.
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "%s..."%(LineNumber)
                if MAX_LINE_NUMBER > 0 and LineNumber > MAX_LINE_NUMBER:
                    break
            Bits = FileLine.strip().split("\t")
            try:
                Spectrum = (Bits[0], Bits[1])
            except:
                continue # bad line
            if Spectrum != OldSpectrum:
                if len(CurrentSpectrumHits):
                    self.CountUnexplainedSpectrumHits(CurrentSpectrumHits)
                OldSpectrum = Spectrum
                CurrentSpectrumHits = []
            if len(Bits)<22 or Bits[21] in ("Cat1", "Cat4", "Cat8", "Cat9"):
                continue
            try:
                GeneNumber = int(Bits[13])
            except:
                continue # header line
            if self.InterestingGeneDict.has_key(GeneNumber):
                Info = (self.InterestingGenes.index(GeneNumber), Bits)
                CurrentSpectrumHits.append(Info)
        if len(CurrentSpectrumHits):
            self.CountUnexplainedSpectrumHits(CurrentSpectrumHits)
        File.close()
        ######################################################################
        # Iterate over the "interesting" genes.  Write out an annotated image
        # and a spreadsheet for each one.
        UnexplainedPeptideFile = open("UnexplainedPeptides.txt", "wb")
        for GeneNumber in self.InterestingGenes:
            Gene = self.Genes[GeneNumber]
            SpectrumCount = self.SpectraForGene.get(GeneNumber, 0)
            PeptideCount = len(self.PeptidesForGene.get(GeneNumber, {}).keys())
            #print "  Known hits: %s in %s peptides"%(SpectrumCount, PeptideCount)
            Str = "%s\t%s\t%s\t%s\t"%(GeneNumber, Gene.Name, SpectrumCount, PeptideCount)
            PepCount = 0
            SpecCount = 0
            for (Key, Count) in Gene.HitDetails.items():
                PepCount += 1
                SpecCount += Count
                if Count >= 10: 
                    # Let's report this peptide to the UnexplainedPeptides file.  We'd like to know
                    # the gene-id, the number of hits, and the protein(s) containing the peptide.
                    Annotation = Gene.DetailPeptides[Key]
                    PeptideString = Annotation[2:-2]
                    PeptideAminos = PeptideString.replace(":", "").replace(";","")
                    PeptideCategory = Gene.PeptideCategories[Key]
                    UPStr = "%s\t%s\t%s\t%s\t%s\t"%(GeneNumber, PeptideString, PeptideCategory, Key, Count)
                    Chunks = PeptideString.replace(";",":").split(":")
                    AllChunksFound = 1
                    for Chunk in Chunks:
                        ChunkFound = 0
                        if len(Chunk)<5:
                            continue
                        ProteinCoverageString = "%s:"%Chunk
                        for Protein in self.Proteins.values():
                            Pos = Protein.Sequence.find(Chunk)
                            if Pos!=-1:
                                ProteinCoverageString += "%s(%s) "%(Protein.ProteinNumber, Pos)
                                ChunkFound = 1
                        UPStr += "%s\t"%ProteinCoverageString
                        if not ChunkFound:
                            AllChunksFound = 0
                    if len(Chunks) == 1:
                        UPStr += "\t"
                    # Add a flag, so we can prioritize the peptides with an
                    # unexplained "chunk" of info
                    UPStr += "%s\t"%AllChunksFound
                    # Add a flag for checking tryptic-ness:
                    TrypticFlag = 0
                    if PeptideString[-1] in ("R","K") or Annotation[-1] == '-':
                        if Annotation[0] in ("-", "R", "K"):
                            TrypticFlag = 1
                    UPStr += "%s\t"%TrypticFlag
                    UPStr += "%s\t"%Gene.PeptideBestScores[Key]
                    UnexplainedPeptideFile.write(UPStr+"\n")
            Str += "%s\t%s\t"%(PepCount, SpecCount)
            Str += "\n"
            InterestingGeneFile.write(Str)
            InterestingGeneFile.flush()
            try:
                self.ReportUnexplainedPeptidesForGene(Gene)
            except:
                traceback.print_exc() # image generation failure!
        File.close()        
    def ReportUnexplainedPeptidesForGene(self, Gene):
        """
        Generate an image summarizing the unexplained peptides for this gene.
        """
        Plotter = GeneBrowserPlot.PlotterClass()
        SpectrumCount = 0
        PeptideCount = 0
        for (Key, Count) in Gene.HitDetails.items():
            SpectrumCount += Count
            PeptideCount += 1
        Plotter.AddHeaderLine("Unexplained peptides for %s %s: %s spectra (%s peptides)"%(Gene.GeneNumber, Gene.Name, SpectrumCount, PeptideCount))
        KnownPeptides = Plotter.AddTrack("KnownPeptides")
        KnownPeptides.Color = GeneBrowserPlot.Colors.Peptide
        ########################################################################
        # Add tracks for known proteins.  As we do so, we also add
        # features to the KnownPeptides track.
        ProteinHitCounts = {}
        for Protein in self.ProteinsByChromosome[Gene.ChromosomeNumber]:
            HitCount = 0
            if Protein.Start <= Gene.End and Protein.End >= Gene.Start:
                for (Key, Count) in Protein.HitDetails.items():
                    HitCount += Count
                    if len(Key) == 2:
                        Feature = KnownPeptides.AddFeature(Key[0], Key[1])
                        if Protein.Strand == 1:
                            ReadingFrame = Key[0] % 3
                        else:
                            ReadingFrame = (Key[1] - 1) % 3
                        Feature.Color = ReadingFrameColors[ReadingFrame]    
                    else:
                        Feature = KnownPeptides.AddFeature(Key[0], min(Key[1],Key[2]))
                        if Protein.Strand == 1:
                            ReadingFrame = Key[0] % 3
                        else:
                            ReadingFrame = (min(Key[1], Key[2]) - 1) % 3
                        Feature.Color = ReadingFrameColors[ReadingFrame]
                        Feature = KnownPeptides.AddFeature(min(Key[1],Key[2]), max(Key[1],Key[2]))
                        Feature.Type = GeneBrowserPlot.FeatureTypes.Dashed
                        Feature = KnownPeptides.AddFeature(max(Key[1], Key[2]), Key[3])
                        if Protein.Strand == 1:
                            ReadingFrame = max(Key[1], Key[2]) % 3
                        else:
                            ReadingFrame = (Key[3] - 1) % 3
                        Feature.Color = ReadingFrameColors[ReadingFrame]
            ProteinHitCounts[Protein.ProteinNumber] = HitCount
        SortedList = []
        for (ProteinNumber, HitCount) in ProteinHitCounts.items():
            SortedList.append((HitCount, ProteinNumber))
        SortedList.sort()
        SortedList.reverse()
        HeaderLineCount = 1
        for (HitCount, ProteinNumber) in SortedList:
            Protein = self.Proteins[ProteinNumber]
            if HeaderLineCount < 7:
                Plotter.AddHeaderLine("Near protein %s (%s hits): %s %s"%(ProteinNumber, HitCount, Protein.ProteinName, Protein.Description))
                HeaderLineCount += 1
            ProteinTrack = Plotter.AddTrack("Prot:%s"%Protein.ProteinNumber)
            for Exon in Protein.Exons:
                ProteinTrack.AddFeature(Exon.Start, Exon.End)
            
        ########################################################################
        # Add tracks for unknown peptides, sorted by spectrum-count.
        SortedList = []
        for (Key, Count) in Gene.HitDetails.items():
            SortedList.append((Count, Key))
        SortedList.sort()
        SortedList.reverse()
        for (Count, Key) in SortedList:
            Peptide = Gene.DetailPeptides[Key]
            Track = Plotter.AddTrack(Gene.PeptideCategories[Key])
            Track.RequiredFlag = 1
            Track.Color = GeneBrowserPlot.Colors.Peptide2
            Track.LabelLineFlag = 0
            if len(Key) == 2:
                Feature = Track.AddFeature(Key[0], Key[1])
                Feature.Label = "%s (%s)"%(Peptide, Count)
                if Gene.Strand == 1:
                    ReadingFrame = Key[0] % 3
                else:
                    ReadingFrame = (Key[1] - 1) % 3
                Feature.Color = ReadingFrameColors[ReadingFrame]                
            else:
                Feature = Track.AddFeature(Key[0], min(Key[1],Key[2]))
                if Gene.Strand == 1:
                    ReadingFrame = Key[0] % 3
                else:
                    ReadingFrame = (min(Key[1], Key[2]) - 1) % 3
                Feature.Color = ReadingFrameColors[ReadingFrame]
                #Feature.Label = "%s (%s)"%(Peptide, Count)
                Feature = Track.AddFeature(min(Key[1],Key[2]), max(Key[1],Key[2]))
                Feature.Type = GeneBrowserPlot.FeatureTypes.Dashed
                Feature = Track.AddFeature(max(Key[1], Key[2]), Key[3])
                if Gene.Strand == 1:
                    ReadingFrame = max(Key[1], Key[2]) % 3
                else:
                    ReadingFrame = (Key[3] - 1) % 3
                Feature.Color = ReadingFrameColors[ReadingFrame]
                Feature.Label = "%s (%s)"%(Peptide, Count)
        FileName = os.path.join("NovelPeptides", "%s.png"%Gene.GeneNumber)
        Plotter.Plot(FileName)
    def FindUnexplainedGenes(self):
        """
        Goal: Populate self.InterestingGenes, self.InterestingGeneDict

        Read all matches from category 2,3,5,6,7.  Iteratively accept
        exon-graph genes that cover as many of these matches as possible,
        until the next gene would explain fewer than self.MinSpectraForXGene spectra
        """
        self.InterestingGenes = []
        self.InterestingGeneDict = {}
        SharedPeptides = {} 
        CurrentSpectrumGenes = [] # list of gene-numbers for current spectrum
        File = open(XGVKG0, "rb")
        OldSpectrum = None
        LineNumber = 0
        #########################################################
        # Read the non-perfect-match annotations:
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "%s..."%(LineNumber)
                if MAX_LINE_NUMBER > 0 and LineNumber > MAX_LINE_NUMBER:
                    break
                
            Bits = FileLine.strip().split("\t")
            # If this spectrum is already accounted for by a known protein, or if it's
            # a "bad" category, then don't hypothesize any new gene:
            if len(Bits)<22 or Bits[21] in ("Cat1", "Cat4", "Cat8", "Cat9"):
                continue
            try:
                MQScore = float(Bits[5])
                PValue = float(Bits[10])
                GeneNumber = int(Bits[13])
            except:
                traceback.print_exc()
                continue # header line
            if MQScore < 0:
                continue # filter!
            print Bits[21], MQScore, GeneNumber
            # Note this gene:
            Gene = self.Genes.get(GeneNumber, None)
            if not Gene:
                Gene = GeneClass()
                Gene.GeneNumber = GeneNumber
                Gene.ChromosomeNumber = int(Bits[16])
                Gene.Strand = int(Bits[17])
                Gene.Name = Bits[3]
                CommaBits = Bits[3].split(",")
                (Start, End) = CommaBits[-1].split("-")
                Gene.Start = int(Start)
                Gene.End = int(End)
                self.Genes[GeneNumber] = Gene
            Spectrum = (Bits[0], Bits[1])
            if Spectrum != OldSpectrum:
                OldSpectrum = Spectrum
                CurrentSpectrumGenes = []
            else:
                # Update shared-counts:
                for OldNumber in CurrentSpectrumGenes:
                    if not SharedPeptides.has_key(OldNumber):
                        SharedPeptides[OldNumber] = {}
                    SharedPeptides[OldNumber][GeneNumber] = SharedPeptides[OldNumber].get(GeneNumber, 0) + 1
                    if not SharedPeptides.has_key(GeneNumber):
                        SharedPeptides[GeneNumber] = {}
                    SharedPeptides[GeneNumber][OldNumber] = SharedPeptides[GeneNumber].get(OldNumber, 0) + 1                
            CurrentSpectrumGenes.append(GeneNumber)
            self.UnexplainedGeneHits[GeneNumber] = self.UnexplainedGeneHits.get(GeneNumber, 0) + 1
        File.close()
        #########################################################
        # Iteratively pick the genes we judge "interesting", on the basis
        # of the number of unexplained hits (self.UnexplainedGeneHits)
        print "Select 'interesting' genes with unexplained hits:"
        while (1):
            BestHitCount = 0
            BestHitGene = None
            for (GeneNumber, HitCount) in self.UnexplainedGeneHits.items():
                if self.InterestingGeneDict.has_key(GeneNumber):
                    continue
                if HitCount > BestHitCount:
                    BestHitCount = HitCount
                    BestHitGene = GeneNumber
            print "BestHitCount: %s hits to XG record %s"%(BestHitCount, BestHitGene)
            if BestHitCount < self.MinSpectraForXGene:
                break # done!
            Gene = self.Genes[BestHitGene]
            print "Gene %s (%s): %s unexplained hits."%(BestHitGene, Gene.Name, BestHitCount)
            SpectrumCount = self.SpectraForGene.get(BestHitGene, 0)
            PeptideCount = len(self.PeptidesForGene.get(BestHitGene, {}).keys())
            print "  Known hits: %s in %s peptides"%(SpectrumCount, PeptideCount)
            self.InterestingGenes.append(BestHitGene)
            self.InterestingGeneDict[BestHitGene] = 1
            # Decrement the hit count of any other genes that share peptides with this one:
            Dict = SharedPeptides.get(BestHitGene, {})
            for (OtherGeneID, SharedCount) in Dict.items():
                if not self.InterestingGeneDict.has_key(OtherGeneID):
                    self.UnexplainedGeneHits[OtherGeneID] -= SharedCount
        # Manually add in some genes which were deemed interesting by an EST-based analysis:
        print "Check for TARGET genes in list:"
        for GeneNumber in ShowNovelPeptidesForGenes:
            if GeneNumber not in self.InterestingGenes:
                if self.Genes.has_key(GeneNumber):
                    self.InterestingGenes.append(GeneNumber)
                    self.InterestingGeneDict[GeneNumber] = 1
                    print "Manually added gene %s to 'interesting' list!"%GeneNumber
                else:
                    print "Can't add gene '%s', we never saw it!"%GeneNumber
            else:
                print "Target gene '%s' is on the interesting list."%GeneNumber
    def SelectPresentProteins(self):
        """
        Goal: Populate self.AcceptedProteins
        
        Read all matches from categories 1 and 4: Peptides that map to a known gene,
        without or with a splice boundary.  Iteratively accept known proteins that
        explain all these hits, until the next addition would explain fewer than
        self.MinSpectraForKnownProtein new spectra.  Then loop over the file again,
        this time counting the spectra and peptides for each protein, and for each gene
        in the exon-graph.
        """
        # List of accepted protein-IDs, in the order in which they were picked:
        self.AcceptedProteins = [] 
        # SharedPeptides[ProteinA][ProteinB] = Count
        # We keep track of peptides shared between proteins, so that a
        # more common protein can "explain away" the shared peptides of a rarer protein.
        SharedPeptides = {} 
        CurrentSpectrumProteins = [] # list of protein-numbers for current spectrum
        ######################################################################
        # READ the perfect-match annotations:
        File = open(XGVKG0, "rb")
        OldSpectrum = None
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "%s..."%(LineNumber)
                if MAX_LINE_NUMBER > 0 and LineNumber > MAX_LINE_NUMBER:
                    break
            Bits = FileLine.strip().split("\t")
            if len(Bits)<22 or Bits[21] not in ("Cat1", "Cat4"):
                continue
            if Bits[24][0] == "*":
                continue
            try:
                ProteinNumber = int(Bits[23])
                ProteinNumbers = [ProteinNumber]
            except:
                try:
                    ProteinNumbers = Bits[23].split(":")
                    ProteinNumbers = map(int, ProteinNumbers)
                except:
                    continue # header line
            for ProteinNumber in ProteinNumbers:
                Spectrum = (Bits[0], Bits[1])
                if Spectrum != OldSpectrum:
                    OldSpectrum = Spectrum
                    CurrentSpectrumProteins = []
                if ProteinNumber not in CurrentSpectrumProteins:
                    # Update shared peptide flags:
                    for OldNumber in CurrentSpectrumProteins:
                        if not SharedPeptides.has_key(OldNumber):
                            SharedPeptides[OldNumber] = {}
                        SharedPeptides[OldNumber][ProteinNumber] = SharedPeptides[OldNumber].get(ProteinNumber, 0) + 1
                        if not SharedPeptides.has_key(ProteinNumber):
                            SharedPeptides[ProteinNumber] = {}
                        SharedPeptides[ProteinNumber][OldNumber] = SharedPeptides[ProteinNumber].get(OldNumber, 0) + 1
                    CurrentSpectrumProteins.append(ProteinNumber)
                Protein = self.Proteins.get(ProteinNumber, None)
                if MAX_LINE_NUMBER <=0 and not Protein:
                    print "ERROR: Unknown protein %s"%ProteinNumber
                if Protein:
                    Protein.HitCount += 1
##                #print "Protein %s hitcount %s"%(Protein.ProteinNumber, Protein.HitCount)
##                (Start, End) = Bits[18].split("-")
##                Start = int(Start)
##                End = int(End)
##                if Bits[20]:
##                    # spliced!
##                    (Donor, Acceptor) = Bits[20].split("-")
##                    Donor = int(Donor)
##                    Acceptor = int(Acceptor)
##                    for Exon in Protein.Exons:
##                        if Donor in (Exon.Start, Exon.End):
##                            Exon.HitCount += 1
##                            Exon.DonorHitCount += 1
##                        if Acceptor in (Exon.Start, Exon.End):
##                            Exon.HitCount += 1
##                            Exon.AcceptorHitCount += 1
##                    Key = (Start, min(Donor, Acceptor), max(Donor, Acceptor), End)
##                    Protein.HitDetails[Key] = Protein.HitDetails.get(Key, 0) + 1
##                    Protein.Peptides[Key] = Bits[2]
##                else:
##                    Key = (Start, End)
##                    Protein.HitDetails[Key] = Protein.HitDetails.get(Key, 0) + 1
##                    Protein.Peptides[Key] = Bits[2]
##                    # Not spliced!
##                    for Exon in Protein.Exons:
##                        if Exon.Start <= Start and Exon.End >= End:
##                            Exon.HitCount += 1
        File.close()
        ######################################################################
        print "Read perfect-match results.  Select proteins:"
        # ITERATE: Pick a Protein that picks as many of the remaining annotations as possible.
        CoveredProteinCount = [0, 0, 0, 0] # 1, 2, 5, 10
        AllPeptides = {}
        AllSplicedPeptides = {}
        AllExonCount = 0
        AllJunctionCount = 0
        CoveredExonCount = 0
        CoveredJunctionCount = 0
        while (1):
            BestHitCount = 0
            BestHitProtein = None
            for Protein in self.Proteins.values():
                if Protein.AcceptedFlag:
                    continue
                if Protein.HitCount > BestHitCount:
                    BestHitCount = Protein.HitCount
                    BestHitProtein = Protein
            if BestHitProtein:
                print "BestHitCount:", BestHitCount, BestHitProtein.ProteinNumber
            if BestHitCount < self.MinSpectraForKnownProtein:
                break # done!
            self.AcceptedProteins.append(BestHitProtein.ProteinNumber)
            BestHitProtein.AcceptedFlag = 1
            
##            Protein = BestHitProtein
##            print
##            print "%s(%s) hits to protein %s (genome-protein coverage %s)\n%s"%(BestHitCount, len(Protein.HitDetails.keys()), Protein.ProteinNumber, 100*Protein.CoveragePercent, Protein.ProteinName)
##            for ExonIndex in range(len(Protein.Exons)):
##                Exon = Protein.Exons[ExonIndex]
##                print "Exon %s (%s-%s): %s hits, %s donor hits, %s acceptor hits"%(ExonIndex,
##                    Exon.Start, Exon.End, Exon.HitCount, Exon.DonorHitCount,
##                    Exon.AcceptorHitCount)
##                AllExonCount += 1
##                if Exon.HitCount:
##                    CoveredExonCount += 1
##                if Exon.AcceptorHitCount:
##                    CoveredJunctionCount += 1
##                if ExonIndex:
##                    AllJunctionCount += 1
            # Now that Protein is accepted, decrease the hit-count of proteins that share
            # peptides with Protein.
            Dict = SharedPeptides.get(BestHitProtein.ProteinNumber, {})
            for (OtherProteinID, SharedCount) in Dict.items():
                OtherProtein = self.Proteins.get(OtherProteinID, None)
                if OtherProtein and not OtherProtein.AcceptedFlag:
                    OtherProtein.HitCount -= SharedCount
##            for Key in Protein.HitDetails.keys():
##                AllPeptides[Key] = AllPeptides.get(Key, 0) + 1
##                if len(Key) > 2:
##                    AllSplicedPeptides[Key] = AllSplicedPeptides.get(Key, 0) + 1
##            PeptideCount = len(Protein.HitDetails.keys())
##            if PeptideCount > 0:
##                CoveredProteinCount[0] += 1
##            if PeptideCount > 1:
##                CoveredProteinCount[1] += 1
##            if PeptideCount > 4:
##                CoveredProteinCount[2] += 1
##            if PeptideCount > 9:
##                CoveredProteinCount[3] += 1
##            #self.GenerateGeneImage(Gene)
##        print "CoveredProteinCount:", CoveredProteinCount
##        print "Peptides:", len(AllPeptides.keys())
##        Count = 0
##        for Value in AllPeptides.values():
##            if Value >= 1:
##                Count += 1
##        print "Peptides 2+", Count
##        print "SplicedPeptides:", len(AllSplicedPeptides.keys())
##        Count = 0
##        for Value in AllSplicedPeptides.values():
##            if Value >= 1:
##                Count += 1
##        print "SplicedPeptides 2+", Count
##        print "AllExonCount", AllExonCount 
##        print "AllJunctionCount", AllJunctionCount 
##        print "CoveredExonCount", CoveredExonCount 
##        print "CoveredJunctionCount", CoveredJunctionCount 
        
    def GetX(self, GenomePosition):
        X = (GenomePosition - self.GroupGenomeStart) # could scale here
        X = max(0, min(X, self.GroupWidth))
        X += self.GroupLeft
        return X
    def GenerateProteinImage(self, Protein):
        TotalHitCount = 0
        PeptideCount = 0
        for (Key, Value) in Protein.HitDetails.items():
            TotalHitCount += Value
            PeptideCount += 1
        PeptideCount = min(50, PeptideCount)
        # First, let's decide how to group the exons, so that we know how
        # wide the image will be:
        ExonGroups = []
        NextExon = 0
        CurrentGroupEnd = None
        while (1):
            if NextExon >= len(Protein.Exons):
                break
            Exon = Protein.Exons[NextExon]
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
        Height = 90 + PeptideCount*10
        self.GenomeY = 55
        self.Image = Image.new("RGB", (Width, Height), Colors.Background)
        print "->-> Created image %s x %s"%(Width, Height)
        self.Draw = ImageDraw.Draw(self.Image)
        # Caption
        if Protein.Strand == 1:
            StrandName = "forward"
        else:
            StrandName = "reverse"
        Str = "Protein %s maps to %s exons on chromosome %s (%s strand)"%(Protein.ProteinNumber, len(Protein.Exons), Protein.ChromosomeNumber, StrandName)
        self.Draw.text((5, 2), Str, Colors.Label)
        self.Draw.text((5, 10), Protein.ProteinName, Colors.Label)
        Str = "%s total hits"%TotalHitCount
        self.Draw.text((5, 18), Str, Colors.Label)
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
            LabelY = 30 + (GroupIndex % 2)*10
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
##                for Mismatch in Exon.Mismatches:
##                    X = self.GetX(Mismatch[1])
##                    self.Draw.line((X, self.GenomeY - 4, X, self.GenomeY + 5), Colors.Mismatch)
##                    self.Draw.line((X-1, self.GenomeY - 4, X-1, self.GenomeY + 5), Colors.Mismatch)
##                    self.Draw.line((X+1, self.GenomeY - 4, X+1, self.GenomeY + 5), Colors.Mismatch)
##                for SNP in Exon.SNPs:
##                    X = self.GetX(SNP[1])
##                    self.Draw.line((X, self.GenomeY - 4, X, self.GenomeY + 5), Colors.SNP)
##                    self.Draw.line((X-1, self.GenomeY - 4, X-1, self.GenomeY + 5), Colors.SNP)
##                    self.Draw.line((X+1, self.GenomeY - 4, X+1, self.GenomeY + 5), Colors.SNP)
                Str = "%s-%s"%(Exon.ProteinStart, Exon.ProteinEnd)
                self.Draw.text((LeftX, self.GenomeY + 8), Str, Colors.Label)
                Sequence = Protein.Sequence[Exon.ProteinStart:Exon.ProteinEnd]
                MaxPrintLen = max(8, (RightX - LeftX) / 7)
                if len(Sequence) < MaxPrintLen:
                    Str = Sequence
                else:
                    Half = (MaxPrintLen / 2) - 1
                    Str = "%s..%s"%(Sequence[:Half], Sequence[-Half:])
                self.Draw.text((LeftX, self.GenomeY + 16), Str, Colors.Label)
            # Draw border with the next exon-group:
            if (GroupIndex < len(ExonGroups)):
                self.Draw.line((Right, 0, Right, Height), Colors.Border)
            self.GroupLeft += self.GroupWidth + 1
        ########################################################
        # Draw the single-exon and two-exon peptides:
        SortedList = []
        for (Key, Count) in Protein.HitDetails.items():
            SortedList.append((Count, Key))
        SortedList.sort()
        SortedList.reverse()
        Y = 90
        for (Count, Key) in SortedList[:50]:
            # Find the exon that covers this one:
            FoundFlag = 0
            for Exon in Protein.Exons:
                if Exon.Start <= Key[0] and Exon.End >= Key[1]:
                    X1 = Exon.LeftX + (Key[0] - Exon.Start)
                    X2 = Exon.LeftX + (Key[1] - Exon.Start)
                    FoundFlag = 1
            if not FoundFlag:
                print "** Error: couldn't find exon covering key!", Key
                Y += 10
                continue
            else:
                self.Draw.line((X1, Y, X2, Y), Colors.Peptide)
                self.Draw.line((X1, Y+1, X2, Y+1), Colors.Peptide)
                print "Found:", Key
            Str = "%s (%s)"%(Protein.DetailPeptides[Key], Count)
            if len(Key) == 4:
                FoundFlag = 0
                for Exon in Protein.Exons:
                    if Exon.Start <= Key[2] and Exon.End >= Key[3]:
                        X3 = Exon.LeftX + (Key[2] - Exon.Start)
                        X4 = Exon.LeftX + (Key[3] - Exon.Start)
                        FoundFlag = 1
                if not FoundFlag:
                    print "** Error: couldn't find exon covering key!", Key
                else:
                    self.Draw.line((X3, Y, X4, Y), Colors.Peptide)
                    self.Draw.line((X3, Y+1, X4, Y+1), Colors.Peptide)
                    for XX in range(X2, X3, 3):
                        self.Draw.line((XX, Y, XX+1, Y), Colors.PeptideSpliceGhost)
                        self.Draw.line((XX, Y+1, XX+1, Y+1), Colors.PeptideSpliceGhost)
                    print "Found:", Key
                    self.Draw.text((X4, Y-4), Str, Colors.Label)
                    
            else:
                self.Draw.text((X2, Y-4), Str, Colors.Label)
            Y += 10
        OutputFileName = "GeneMap\\Peptides%s.png"%Protein.ProteinNumber
        self.Image.save(OutputFileName, "png")        
    def SummarizeMisses(self):
        pass
    def LoadSequences(self, DBPath, IndexPath):
        LastPos = -1
        #IndexFile = open(IndexPath, "rb")
        DBFile = open(DBPath, "rb")
        ProteinNumber = 0
        Data = DBFile.read()
        while (1):
            Pos = Data.find("*", LastPos + 1)
            if Pos == -1:
                if self.Proteins.has_key(ProteinNumber):
                    self.Proteins[ProteinNumber].Sequence = Data[LastPos + 1:]
                break
            if self.Proteins.has_key(ProteinNumber):
                self.Proteins[ProteinNumber].Sequence = Data[LastPos + 1:Pos]
            ProteinNumber += 1
            LastPos = Pos
        #IndexFile.close()
        DBFile.close()
    def ReportSNPs(self):
        class SNPClass:
            def __init__(self):
                self.ChromosomeNumber = None
                self.Position = None
                # Count how many times we spanned this genomic pos'n without the polymorphism:
                self.WTOccurrenceCount = 0
                self.OccurrenceCount = 0
                self.WTBestScore = (-999, -999, "", "", "")
                self.OtherBestScore = (-999, -999, "", "", "")
                # Type 1: Match protein but not genome (category 1/4); default
                # Type 2: Match genome but not protein (category 2/5)
                # Type 3: Match neither genome nor protein (category 2/5)
                self.Type = 0
            def __cmp__(self, Other):
                "Sort by chromosome number, then position on chromosome"
                if not isinstance(Other, SNPClass):
                    return 1
                if self.ChromosomeNumber != Other.ChromosomeNumber:
                    return cmp(self.ChromosomeNumber, Other.ChromosomeNumber)
                if self.Position < Other.Position:
                    return -1
                if self.Position > Other.Position:
                    return 1
                return cmp(id(self), id(Other))
        self.SNPDict = {}
        self.SNPsByProtein = {}
        self.SNPsByChromosome = {}
        File = open(XGVKG0, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.split("\t")
            try:
                ChromosomeNumber = int(Bits[16])
                Strand = int(Bits[17])
                Category = int(Bits[21][3:]) # cat1 -> 1
            except:
                continue
            if Category not in (1,2,4,5):
                continue
            if len(Bits)<26:
                continue # temp: for running report while XGVKG is still running!
            # Skip over hits that don't have valid genomic coords, or that are out of frame with the known protein:
            if Bits[22] and Bits[22][0] == "*":
                continue
            if Bits[25] and Bits[25][0] == "*":
                continue
            MQScore = float(Bits[5])
            DeltaScore = float(Bits[12])
            for BitIndex in (22, 25):
                for SNPCommaBit in Bits[BitIndex].split(","):
                    SNPCommaBit = SNPCommaBit.strip()
                    if not SNPCommaBit:
                        continue
                    # 105279024 (N/K)
                    SNPSpaceBits = SNPCommaBit.split()
                    #print SNPSpaceBits
                    if len(SNPSpaceBits) != 2:
                        continue
                    Position = int(SNPSpaceBits[0])
                    Poly = SNPSpaceBits[1]
                    Key = (ChromosomeNumber, Position, Poly[1])
                    print SNPCommaBit, Key
                    if Bits[24][0] == "*":
                        ProteinID = int(Bits[24][1:])
                    elif Bits[24].find(":"):
                        ProteinID = int(Bits[24].split(":")[0])
                    else:
                        ProteinID = int(Bits[24])
                    SNP = self.SNPDict.get(Key, None)
                    if not SNP:
                        SNP = SNPClass()
                        SNP.ChromosomeNumber = ChromosomeNumber
                        SNP.Position = Position
                        SNP.Poly = Poly
                        SNP.ProteinID = ProteinID
                        self.SNPDict[Key] = SNP
                        if not self.SNPsByProtein.has_key(ProteinID):
                            self.SNPsByProtein[ProteinID] = []
                        self.SNPsByProtein[ProteinID].append(SNP)
                        if not self.SNPsByChromosome.has_key(ChromosomeNumber):
                            self.SNPsByChromosome[ChromosomeNumber] = []
                        self.SNPsByChromosome[ChromosomeNumber].append(SNP)
                    SNP.OccurrenceCount += 1
                    Score = (MQScore, DeltaScore, Bits[0], Bits[1], Bits[2])
                    if Score > SNP.OtherBestScore:
                        SNP.OtherBestScore = Score
                    if BitIndex == 22: # genome mismatch
                        if SNP.Type == 0:
                            SNP.Type = 1
                        elif SNP.Type == 2:
                            SNP.Type = 3
                            SNP.Poly += Poly
                    else: # known-protein mismatch
                        if SNP.Type == 0:
                            SNP.Type = 2
                        elif SNP.Type == 1:
                            SNP.Type = 3
                            SNP.Poly += Poly
        File.close()
        print "#TEMP: Report initial SNP dictionary."
        self.ReportSNPDictionary(self.SNPDict) #% TEMP!
        print "#Sort SNP lists..."
        for List in self.SNPsByChromosome.values():
            List.sort()
        print "#Read results again, to get non-snippy coverage of SNPs..."
        File = open(XGVKG0, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                ChromosomeNumber = int(Bits[16])
                Strand = int(Bits[17])
                Category = int(Bits[21][3:]) # cat1 -> 1
            except:
                continue
            # Skip over all hits except SNPless category 1 or 4 hits:
            if Category not in (1,4):
                continue
            if Bits[22] or Bits[25]:
                continue
            (Start, End) = Bits[18].split("-")
            Start = int(Start)
            End = int(End)
            MQScore = float(Bits[5])
            DeltaScore = float(Bits[12])
            Score = (MQScore, DeltaScore, Bits[0], Bits[1], Bits[2])
            if Category == 1:
                self.CoverSNPs(ChromosomeNumber, Start, End, Score)
            else:
                DashBits = Bits[20].split("-")
                if len(DashBits) != 2:
                    continue # Ugh, skip three-exon guys for now
                (Donor, Acceptor) = DashBits
                Donor = int(Donor)
                Acceptor = int(Acceptor)
                if Strand == 1:
                    self.CoverSNPs(ChromosomeNumber, Start, Donor, Score)
                    self.CoverSNPs(ChromosomeNumber, Acceptor, End, Score)
                else:
                    self.CoverSNPs(ChromosomeNumber, Start, Acceptor, Score)
                    self.CoverSNPs(ChromosomeNumber, Donor, End, Score)
        File.close()
        self.ReportSNPDictionary(self.SNPDict)
    def ReportSNPDictionary(self, SNPDict):
        Keys = SNPDict.keys()
        Keys.sort()
        File = open("XGSNPReport.txt", "wb")
        Header = "Chromosome\tPosition\tPolymorphism\tType\tSNPCount\tWTCount\tProtein\t"
        Header += "WT MQScore\tWT DeltaScore\tWT Run\tWT Scan#\t"
        Header += "MQScore\tDeltaScore\tRun\tScan#\t"
        File.write(Header)
        for Key in Keys:
            SNP = SNPDict[Key]
            Str = "%s\t%s\t%s\t%s\t"%(SNP.ChromosomeNumber, SNP.Position, SNP.Poly, SNP.Type)
            ProteinName = self.Proteins[SNP.ProteinID].ProteinName
            Str += "%s\t%s\t%s\t"%(SNP.OccurrenceCount, SNP.WTOccurrenceCount, ProteinName)
            FileName = SNP.WTBestScore[2].replace("/","\\").split("\\")[-1]
            Str += "%s\t%s\t%s\t%s\t%s\t"%(SNP.WTBestScore[0],SNP.WTBestScore[1], FileName, SNP.WTBestScore[3], SNP.WTBestScore[4])
            FileName = SNP.OtherBestScore[2].replace("/","\\").split("\\")[-1]
            Str += "%s\t%s\t%s\t%s\t%s\t"%(SNP.OtherBestScore[0],SNP.OtherBestScore[1], FileName, SNP.OtherBestScore[3], SNP.OtherBestScore[4])
            print Str
            File.write(Str + "\n")
        File.close()
    def CoverSNPs(self, Chromosome, Start, End, Score):
        SNPList = self.SNPsByChromosome.get(Chromosome, [])
        for SNP in SNPList:
            if SNP.Position >= Start and SNP.Position < End:
                SNP.WTOccurrenceCount += 1
                if Score > SNP.WTBestScore:
                    SNP.WTBestScore = Score
                
    def ReportAlternativeSplicing(self):
        """
        Alternative splicing report:
        Find all splice donors and splice acceptors that have multiple putative partners.
        For each point, keep track of its partners...and their categories!
        """
        AltSplicingAnchorThreshold = 7 # base pairs
        class SplicingPoint:
            def __init__(self, ChromosomeNumber, Strand, Pos):
                self.ChromosomeNumber = ChromosomeNumber
                self.Strand = Strand
                self.Pos = Pos
                self.CategoryCounts = [0]*15
                self.Partners = []
                self.ReportFlag = 0
                self.HappyFlag = 0
                self.PartnerCounts = {} # partner pos -> count
                self.PartnerCounts4 = {} # partner pos -> count, cat4 only
                self.PartnerProteinID = {}
                self.PartnerPeptides = {}
            def GetOccurrenceInfo(self):
                Str = "All %s"%self.CategoryCounts[0]
                if self.CategoryCounts[4]:
                    Str += " Match %s"%self.CategoryCounts[4]
                if self.CategoryCounts[5]:
                    Str += " Mismatch %s"%self.CategoryCounts[5]
                if self.CategoryCounts[6]:
                    Str += " OneSide %s"%self.CategoryCounts[6]
                if self.CategoryCounts[7]:
                    Str += " Unknown %s"%self.CategoryCounts[7]
                if self.CategoryCounts[8]:
                    Str += " GenomicSNP %s"%self.CategoryCounts[8]
                if self.CategoryCounts[9]:
                    Str += " ShortX %s"%self.CategoryCounts[9]
                return Str
        # (chrom, strand, pos) -> SplicingPoint
        self.SplicePoints = {}
        SpecialTestFile = open("SpecialSplicingReportFile", "wb")
        File = open(XGVKG0, "rb")
        LineNumber = 0
        SpliceLineCount = 0
        self.PeptideBestScores = {}
        # Iterate over all high-quality annotations and remember the splice points in
        # categories 4, 5, 6.  
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 10000 == 0:
                print "%s..."%LineNumber
                if MAX_LINE_NUMBER > 0 and LineNumber > MAX_LINE_NUMBER:
                    break
            Bits = FileLine.split("\t")
            if len(Bits) < 20:
                continue
            try:
                ChromosomeNumber = int(Bits[16])
                Strand = int(Bits[17])
                Category = int(Bits[21][3:]) # cat1 -> 1
            except:
                continue
            if Category == 4 and Bits[24][0] == "*":
                continue # SKIP this guy, splicing is not reported correctly.
            JunctionStr = Bits[20].strip()
            if not JunctionStr:
                continue
            if Category in (8, 9):
                continue # Bogus match, skip it entirely!
            ###################################################
            # Decide whether we'll consider this match "happy" for purposes
            # of highlighting alternative splicing.  
            if Category in (1, 4):
                HappyRowFlag = 1
            else:
                HappyRowFlag = 0
            # We consider category-6 hits to be happy, provided the same
            # protein serves as an acceptor and a donor.  (That would be
            # explainable with exon-skipping)
            if Category == 6:
                Donors = Bits[24].strip()
                if Donors and Donors[0] == "*":
                    Donors = Donors[1:]
                DonorList = Donors.split(":")
                Acceptors = Bits[25]
                if Acceptors and Acceptors[0] == "*":
                    Acceptors = Acceptors[1:]
                AcceptorList = Acceptors.split(":")
                # MORE GENEROUS flagging:
                if Donors and Acceptors:
                    HappyRowFlag = 1
##                for ProteinID in DonorList:
##                    if not ProteinID:
##                        continue
##                    if ProteinID in AcceptorList:
##                        HappyRowFlag = 1
            # Genomic mismatches make us unhappy:
            if Bits[23] and Bits[23][0] == "*":
                HappyRowFlag = 0
            # Non-tryptic peptides should *maybe* make us unhappy:
            # Inferior matches make us unhappy:
            MQScore = float(Bits[5])
            DeltaScore = float(Bits[12])
            if float(Bits[11]) < -0.1:
                HappyRowFlag = 0
            # Poor-quality matches make us unhappy:
            if float(Bits[5]) < -1:
                HappyRowFlag = 0
            if HappyRowFlag and Category != 4:
                print "H", Category, DonorList, AcceptorList
            ###################################################
            SpliceLineCount += 1
            (Start, End) = Bits[18].split("-")
            Start = int(Start)
            End = int(End)
            Junctions = JunctionStr.split()
            PrevDonor = None
            PrevAcceptor = None
            SpliceList = []
            for Junction in Junctions:
                Junction = Junction.strip()
                if not Junction:
                    continue
                (Donor, Acceptor) = Junction.split("-")
                Donor = int(Donor)
                Acceptor = int(Acceptor)
                SpliceList.append((Donor, Acceptor))
            for Index in range(len(SpliceList)):
                (Donor, Acceptor) = SpliceList[Index]
                # Decide whether this portion of the peptide is long enough to "anchor" us
                if Strand == 1:
                    if Index == 0:
                        LeftChunkLength = Donor - Start
                    else:
                        LeftChunkLength = Donor - SpliceList[Index - 1][1]
                    if Index == len(SpliceList) - 1:
                        RightChunkLength = End - Acceptor
                    else:
                        RightChunkLength = SpliceList[Index + 1][0] - Acceptor
                else:
                    # 6418339-6418694
                    # 6418686-6418536 6418477-6418344
                    if Index == 0:
                        LeftChunkLength = End - Donor
                    else:
                        LeftChunkLength = SpliceList[Index - 1][1] - Donor
                    if Index == len(SpliceList) - 1:
                        RightChunkLength = Acceptor - Start
                    else:
                        RightChunkLength = Acceptor - SpliceList[Index + 1][0]
                if LeftChunkLength > AltSplicingAnchorThreshold:
                    # Find (or create) the donor point object:
                    Key = (ChromosomeNumber, Strand, Donor)
                    DonorPoint = self.SplicePoints.get(Key, None)
                    if not DonorPoint:
                        DonorPoint = SplicingPoint(ChromosomeNumber, Strand, Donor)
                        self.SplicePoints[Key] = DonorPoint
                    if HappyRowFlag:
                        DonorPoint.HappyFlag += 1
                else:
                    DonorPoint = None
                if RightChunkLength > AltSplicingAnchorThreshold:
                    # Find (or create) the acceptor point object:
                    Key = (ChromosomeNumber, Strand, Acceptor)
                    AcceptorPoint = self.SplicePoints.get(Key, None)
                    if not AcceptorPoint:
                        AcceptorPoint = SplicingPoint(ChromosomeNumber, Strand, Acceptor)
                        self.SplicePoints[Key] = AcceptorPoint
                    if HappyRowFlag:
                        AcceptorPoint.HappyFlag += 1
                else:
                    AcceptorPoint = None
                if DonorPoint and AcceptorPoint:
                    # Remember partners:
                    if AcceptorPoint not in DonorPoint.Partners:
                        DonorPoint.Partners.append(AcceptorPoint)
                    if DonorPoint not in AcceptorPoint.Partners:
                        AcceptorPoint.Partners.append(DonorPoint)
                    AcceptorPoint.PartnerCounts[DonorPoint.Pos] = AcceptorPoint.PartnerCounts.get(DonorPoint.Pos, 0) + 1
                    DonorPoint.PartnerCounts[AcceptorPoint.Pos] = DonorPoint.PartnerCounts.get(AcceptorPoint.Pos, 0) + 1
                    if Category == 4:
                        AcceptorPoint.PartnerCounts4[DonorPoint.Pos] = AcceptorPoint.PartnerCounts4.get(DonorPoint.Pos, 0) + 1
                        DonorPoint.PartnerCounts4[AcceptorPoint.Pos] = DonorPoint.PartnerCounts4.get(AcceptorPoint.Pos, 0) + 1
                        
                    #################################################################################
                    # Make note of this protein ID for the pair:
                    IDListA = Bits[23].split(":")
                    IDListB = Bits[24]
                    if IDListB and IDListB[0] == "*":
                        IDListB = IDListB[1:]
                    IDListB = IDListB.split(":")
                    ProteinIDs = []
                    IDListA.extend(IDListB)
                    for ID in IDListA:
                        try:
                            IntID = int(ID)
                            ProteinIDs.append(IntID)
                        except:
                            continue
                    if not AcceptorPoint.PartnerProteinID.has_key(DonorPoint.Pos):
                        AcceptorPoint.PartnerProteinID[DonorPoint.Pos] = {}
                    for ID in ProteinIDs:
                        AcceptorPoint.PartnerProteinID[DonorPoint.Pos][ID] = AcceptorPoint.PartnerProteinID[DonorPoint.Pos].get(ID, 0) + 1
                    if not DonorPoint.PartnerProteinID.has_key(AcceptorPoint.Pos):
                        DonorPoint.PartnerProteinID[AcceptorPoint.Pos] = {}
                    for ID in ProteinIDs:
                        DonorPoint.PartnerProteinID[AcceptorPoint.Pos][ID] = DonorPoint.PartnerProteinID[AcceptorPoint.Pos].get(ID, 0) + 1
                    #################################################################################
                    # Make note of this peptide for the pair:
                    PeptideString = "%s.%s.%s"%(Bits[2][0], Bits[19], Bits[2][-1])
                    (OldScore, OldDelta, OldSpec, OldScan) = self.PeptideBestScores.get(PeptideString, (-999, -999, None, None))
                    if MQScore > OldScore:
                        self.PeptideBestScores[PeptideString] = (MQScore, DeltaScore, Bits[0], Bits[1])
                    if not AcceptorPoint.PartnerPeptides.has_key(DonorPoint.Pos):
                        AcceptorPoint.PartnerPeptides[DonorPoint.Pos] = {}
                    if not AcceptorPoint.PartnerPeptides[DonorPoint.Pos].has_key(PeptideString):
                        AcceptorPoint.PartnerPeptides[DonorPoint.Pos][PeptideString] = (float(Bits[10]), 1)
                    else:
                        (PValue, Count) = AcceptorPoint.PartnerPeptides[DonorPoint.Pos][PeptideString]
                        PValue = min(PValue, float(Bits[10]))
                        Count += 1
                        AcceptorPoint.PartnerPeptides[DonorPoint.Pos][PeptideString] = (PValue, Count)
                    if not DonorPoint.PartnerPeptides.has_key(AcceptorPoint.Pos):
                        DonorPoint.PartnerPeptides[AcceptorPoint.Pos] = {}
                    if not DonorPoint.PartnerPeptides[AcceptorPoint.Pos].has_key(PeptideString):
                        DonorPoint.PartnerPeptides[AcceptorPoint.Pos][PeptideString] = (float(Bits[10]), 1)
                    else:
                        (PValue, Count) = DonorPoint.PartnerPeptides[AcceptorPoint.Pos][PeptideString]
                        PValue = min(PValue, float(Bits[10]))
                        Count += 1
                        DonorPoint.PartnerPeptides[AcceptorPoint.Pos][PeptideString] = (PValue, Count)                                    
                                                                    
                # Count number-of-hits:
                if DonorPoint:
                    DonorPoint.CategoryCounts[0] += 1
                    DonorPoint.CategoryCounts[Category] += 1
                if AcceptorPoint:
                    AcceptorPoint.CategoryCounts[0] += 1
                    AcceptorPoint.CategoryCounts[Category] += 1
        File.close()
        print "Parsed %s lines, %s with splicing; found a total of %s splice points"%(LineNumber, SpliceLineCount, len(self.SplicePoints.keys()))
        #############################################################################
        # Iterate over all splicing points, and decide whether they are worth reporting.
        # We consider a point "happy" if it has at least one hit in category 1/4.
        # Exception: A point is not made happy by a hit with genomic mismatch.
        # Exception: A cat6 hit makes a point happy, if a protein acts as donor and acceptor 
        # for the point (i.e. exon-skipping!)
        # We'll report a splicing point if it is happy, and has two or more happy partners.
        # We'll sort points by the occurrence-counts of their happy partners.
        Keys = self.SplicePoints.keys()
        Keys.sort()
        for Key in Keys:
            Point = self.SplicePoints[Key]
            if not Point.HappyFlag:
                continue
            Point.HappyPartnerCount = 0
            for Partner in Point.Partners:
                if Partner.HappyFlag:
                    Point.HappyPartnerCount += 1
            if Point.HappyPartnerCount > 1:
                Point.ReportFlag = 1
        # Sort the "interesting" points by the second-highest number of times it partners
        # with a happy partner.
        SortedList = []
        for SplicePoint in self.SplicePoints.values():
            if SplicePoint.ReportFlag:
                HappyPartnerCounts = []
                for (Pos, Count) in SplicePoint.PartnerCounts.items():
                    # Decide whether this is a happy partner:
                    HappyFlag = 0
                    for Partner in SplicePoint.Partners:
                        if Partner.Pos == Pos and Partner.HappyFlag:
                            HappyFlag = 1
                    if HappyFlag:
                        HappyPartnerCounts.append(Count)
                HappyPartnerCounts.sort()
                HappyPartnerCounts.reverse()
                if len(HappyPartnerCounts) > 1:
                    SortedList.append((HappyPartnerCounts[1], SplicePoint))
                else:
                    SortedList.append((-1, SplicePoint))
        SortedList.sort()
        SortedList.reverse()
        SplicePointList = []
        for (Count, Point) in SortedList:
            SplicePointList.append(Point)
        SummaryReportFile = open("AltSplicingSummaryTable.txt", "wb")
        Header = "Index\tCount\tPosition\tStrand\tChromosomt\tHappyPartners\tCat4\tCat5\tCat6\tCat7\tCat8\tCat9\tAP4\tAP5\tAP6\tAP7\tAP8\tAP9\tPartnerA\tPA4\tPAAll\tPartnerB\tPB4\tPBAll\tPartnerC\tPC4\tPCAll\t"
        print Header
        for PointIndex in range(len(SortedList)):
            (Count, Point) = SortedList[PointIndex]
            Str = "%s\t%s\t%s\t%s\t%s\t%s\t"%(PointIndex, Count, Point.Pos, Point.Strand, Point.ChromosomeNumber, Point.HappyPartnerCount)
            Str += "%s\t%s\t%s\t%s\t%s\t%s\t"%(Point.CategoryCounts[4], Point.CategoryCounts[5],
                                       Point.CategoryCounts[6], Point.CategoryCounts[7],
                                       Point.CategoryCounts[8], Point.CategoryCounts[9])
            TotalCounts = [0]*10
            for Partner in Point.Partners:
                for X in range(10):
                    TotalCounts[X] += Partner.CategoryCounts[X]
            Str += "%s\t%s\t%s\t%s\t%s\t%s\t"%(TotalCounts[4], TotalCounts[5],
                                       TotalCounts[6], TotalCounts[7],
                                       TotalCounts[8], TotalCounts[9])
            SortedPartners = []
            for Partner in Point.Partners:
                SortedPartners.append((Point.PartnerCounts.get(Partner.Pos, 0), Partner.CategoryCounts[4], Partner))
            SortedPartners.sort()
            SortedPartners.reverse()
            print Str
            #print "Sorted partners:", SortedPartners
            ###########################################################################
            # SUMMARY REPORT:
            for (PCount, Cat4Count, Partner) in SortedPartners:
                if not Partner.HappyFlag:
                    continue
                # Note the best protein for this partnership:
                ProList = []
                for (ProteinID, Count) in Point.PartnerProteinID[Partner.Pos].items():
                    ProList.append((Count, ProteinID))
                ProList.sort()
                if not ProList:
                    print "*** ERROR: UNKNOWN PROTEIN ID!"
                    ProteinID = 0
                else:
                    ProteinID = ProList[-1][1]
                ###ProteinID = self.PartnerProteinID[Partner.Pos]
                Protein = self.Proteins[ProteinID]
                Str = "%s\t"%PCount
                if Point.PartnerCounts4.get(Partner.Pos, 0):
                    Str += "%s\t"%(Protein.ProteinName)
                else:
                    Str += "*%s\t"%(Protein.ProteinName) # star for exon skipping
                # Note the best peptide for this partnership:
                PPList = []
                for (Key, Tuple) in Point.PartnerPeptides[Partner.Pos].items():
                    PPList.append((Tuple[0], -Tuple[1], Key))
                PPList.sort()
                PeptideString = PPList[0][-1]
                Str += "%s\t"%PeptideString
                # Add the best MQScore / delta for the peptide:
                (MQScore, DeltaScore, FileName, ScanNumber) = self.PeptideBestScores.get(PeptideString, ("???", "???", "",""))
                FileName = FileName.replace("/","\\").split("\\")[-1]
                Str += "%s\t%s\t%s\t%s\t"%(MQScore, DeltaScore, FileName, ScanNumber)
                # Possibly-used fields:
                Str += "chr%s\t"%Point.ChromosomeNumber
                Min = min(Point.Pos, Partner.Pos)
                Max = max(Point.Pos, Partner.Pos)
                if Point.Strand == 1:
                    Str += "forward\t"
                    Str += "%s-%s\t"%(Min, Max)
                else:
                    Str += "reverse\t"
                    Str += "%s-%s\t"%(Max, Min)
                Str += "%s\t"%Protein.Description
                # Extra fields, not for printing, useful for us:
                Str += "%s\t%s\t"%(PointIndex, ProteinID)
                print Str
                SummaryReportFile.write(Str+"\n")
                #PCount = Point.PartnerCounts.get(Partner.Pos, 0)
                # One row per partner...
                #print "%s\t%s\t%s\t%s\t%s\t%s\t"%(PointIndex, Point.Pos, Point.HappyFlag, 
                #    Partner.Pos, Partner.HappyFlag, Point.PartnerCounts.get(Partner.Pos, 0))
                #Str += "%s\t%s\t%s\t"%(Partner.Pos, Cat4Count, Partner.CategoryCounts[0])
            SummaryReportFile.write("\n") # padding to separate groups
            ###########################################################################
##            print Str
##            print S
##            print "Splice boundary at %s on strand %s chromosome %s"%(Point.Pos, Point.Strand, Point.ChromosomeNumber)
##            print "Occurrences:", Point.GetOccurrenceInfo()
##            for Partner in Point.Partners:
##                print "  %s %s %s (%s)"%(Partner.ChromosomeNumber, Partner.Strand, Partner.Pos, Partner.GetOccurrenceInfo())
        #############################################################################
        # Iterate over results again, and print out the lines for all interesting
        # alt-splicing points.  For each line, check each junction point, and if the
        # junction point is one of our "interesting" points (if ReportFlag is true),
        # the print the file line.
        # Note: We print rows whether they're "happy" or not!)
        # Note: We might print a row two or more times, if it corresponds to two or
        # more "interesting" junction points.
        #
        # Header for second half of report:
        print "\nAltSpliceIndex\tFile\tScan\tAnnotation\tGeneName\tCharge\tMQScore\tCutScore\tIntenseBY\tBYPresent\tNTT\tp-value\tDeltaCN\tDeltaCNOther\tRecordNumber\tDBFilePos\tSpecFilePos\tChromosome\tStrand\tGenomicPos\tSplicedSequence\tSplices\tCategory\tExtra1\tExtra2\tExtra3\tExtra4\tExtra5\tMultiHit\t"
        
        File = open(XGVKG0, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if len(Bits) < 20:
                continue
            try:
                ChromosomeNumber = int(Bits[16])
                Strand = int(Bits[17])
                Category = int(Bits[21][3:]) # cat1 -> 1
            except:
                continue
            if Category in (8, 9):
                continue
            if not Bits[20]:
                continue
            ReportThisLine = []
            Junctions = Bits[20].split()
            for Junction in Junctions:
                Junction = Junction.strip()
                if not Junction:
                    continue
                (Donor, Acceptor) = Junction.split("-")
                Donor = int(Donor)
                Acceptor = int(Acceptor)
                Key = (ChromosomeNumber, Strand, Donor)
                DonorPoint = self.SplicePoints.get(Key, None)
                if DonorPoint and DonorPoint.ReportFlag:
                    #ReportThisLine = min(ReportThisLine, SplicePointList.index(DonorPoint))
                    ReportThisLine.append(SplicePointList.index(DonorPoint))
                Key = (ChromosomeNumber, Strand, Acceptor)
                AcceptorPoint = self.SplicePoints.get(Key, None)
                if AcceptorPoint and AcceptorPoint.ReportFlag:
                    #ReportThisLine = ReportThisLine = min(ReportThisLine, SplicePointList.index(AcceptorPoint))
                    ReportThisLine.append(SplicePointList.index(AcceptorPoint))
            for ReportThisLineIndex in ReportThisLine:
                print "%s\t%s"%(ReportThisLineIndex, FileLine.strip())
            # SORT this report by columns A, T, U!
        File.close()
    def ReportHypotheticalProteinCoverage(self):
        MinimumPeptidesForHypoProtein = 1
        MinimumSpectraForHypoProtein = 2
        OutputFile = open("HypoProteinCoverage.txt", "wb")
        Header = "ProteinID\tIPIID\tAnnotation\tAnnotationLevel\tFullyTrypticPeptideCount\tPeptideCount\tSpectrumCount\t"
        for Letter in ("A", "B", "C", "D"):
            Header += "Peptide%s\tSpectra%s\tMQScore%s\tDeltaScore%s\tSpecFile\tScan#\t"%(Letter, Letter, Letter, Letter)
        Header += "\n"
        OutputFile.write(Header)
        # First, decide the "annotation level" of each protein.
        LevelCounts = [0, 0, 0,]
        for Protein in self.Proteins.values():
            if Protein.Description.find("HYPOTHETICAL")!=-1:
                Protein.AnnotationLevel = 0
            elif Protein.Description.find("PUTATIVE")!=-1:
                Protein.AnnotationLevel = 1
            else:
                Protein.AnnotationLevel = 2
            LevelCounts[Protein.AnnotationLevel] += 1
            Protein.Peptides = {} # peptide -> count
            Protein.HitDetails = {} # peptide -> info
        print "annotation level counts:", LevelCounts
        ##########################################################################################
        # Now, iterate over cat1 / cat4 hits, and track peptides from proteins with low annotation level.
        LineNumber = 0
        for FileName in ("XGVKG1.txt", "XGVKG4.txt"):
            print "Parsing lines from %s..."%FileName
            File = open(FileName, "rb")
            OldSpectrum = None
            PendingHits = []
            for FileLine in File.xreadlines():
                LineNumber += 1
                Bits = FileLine.split("\t")
                Spectrum = (Bits[0], Bits[1])
                if Spectrum != OldSpectrum:
                    self.StoreHypoProteinHits(OldSpectrum, PendingHits)
                    PendingHits = []
                    OldSpectrum = Spectrum
                ProteinBits = Bits[23].split(":")
                for ProteinBit in ProteinBits:
                    try:
                        ProteinNumber = int(ProteinBit)
                    except:
                        continue
                    Protein = self.Proteins.get(ProteinNumber, None)
                    if Protein:
                        PendingHit = (Protein.AnnotationLevel, float(Bits[5]), Bits, Protein, Bits[2], Bits[18])
                        PendingHits.append(PendingHit)
                    else:
                        print "UNKNOWN protein number:", ProteinNumber
            self.StoreHypoProteinHits(OldSpectrum, PendingHits)
        ##########################################################################################
        # Now, iterate the "novel" proteins with many peptides and many spectra:
        SortedList = []
        for Protein in self.Proteins.values():
            if Protein.AnnotationLevel < 2:
                SpectrumCount = 0
                PeptideCount = 0
                for (Peptide, Count) in Protein.Peptides.items():
                    SpectrumCount += Count
                    PeptideCount += 1
                Tuple = (-Protein.AnnotationLevel, SpectrumCount, Protein)
                SortedList.append(Tuple)
        SortedList.sort()
        SortedList.reverse()
        for (Dummy, DummyB, Protein) in SortedList:
            SpectrumCount = 0
            PeptideCount = 0
            TrypPeptideCount = 0
            CommonPepList = []
            for (Peptide, Count) in Protein.Peptides.items():
                SpectrumCount += Count
                PeptideCount += 1
                CommonPepList.append((Count, Peptide))
                if Peptide[0] in ("*","-", "K", "R"):
                    if Peptide[-1] in ("*","-") or Peptide[-3] in ("K","R"):
                        TrypPeptideCount += 1
            # Apply a cutoff:
            if SpectrumCount < MinimumSpectraForHypoProtein:
                continue
            if PeptideCount < MinimumPeptidesForHypoProtein:
                continue
            # Report:
            CommonPepList.sort()
            CommonPepList.reverse()
            Str = "%s\t%s\t%s\t"%(Protein.ProteinNumber, Protein.ProteinName, Protein.Description)
            Str += "%s\t%s\t%s\t%s\t"%(Protein.AnnotationLevel, TrypPeptideCount, PeptideCount, SpectrumCount)
            for CommonPeptideIndex in range(min(4, len(CommonPepList))):
                (Count, Peptide) = CommonPepList[CommonPeptideIndex]
                Str += "%s\t%s\t"%(Peptide, Count)
                BestScore = Protein.BestPeptideScores.get(Peptide, ("", "", "", ""))
                SpecFileName = BestScore[2].replace("/", "\\").split("\\")[-1]
                Str += "%s\t%s\t%s\t%s\t"%(BestScore[0], BestScore[1], SpecFileName, BestScore[3])
            print Str
            OutputFile.write(Str + "\n")
    def StoreHypoProteinHits(self, Spectrum, PendingHits):
        if not len(PendingHits):
            return
        #print "StoreHypoProteinHits: %s"%len(PendingHits)
        PendingHits.sort()
        PendingHits.reverse()
        MaxAnnotationLevel = PendingHits[0][0]
        for PendingHit in PendingHits:
            if PendingHit[0] < MaxAnnotationLevel:
                break # don't bother to remember these!
            Bits = PendingHit[2]
            Protein = PendingHit[3]
            Peptide = PendingHit[4]
            Protein.Peptides[Peptide] = Protein.Peptides.get(Peptide, 0) + 1
            Protein.HitDetails[Peptide] = (PendingHit[5])
            if not hasattr(Protein, "BestPeptideScores"):
                Protein.BestPeptideScores = {}
            OldScore = Protein.BestPeptideScores.get(Peptide, (-999, -999, None, None))
            NewScore = (PendingHit[1], float(Bits[12]), Bits[0], Bits[1])
            if NewScore > OldScore:
                Protein.BestPeptideScores[Peptide] = NewScore
        
    def ReportKnownProteinIsoforms(self):
        GenesByChromosome = []
        RGenesByChromosome = []
        for X in range(49):
            GenesByChromosome.append([])
            RGenesByChromosome.append([])
        for Protein in self.Proteins.values():
            if Protein.CoveragePercent < 0.9:
                continue
            if Protein.Strand == 1:
                GenesByChromosome[Protein.ChromosomeNumber].append((Protein.Start, Protein.End, Protein))
            else:
                RGenesByChromosome[Protein.ChromosomeNumber].append((Protein.Start, Protein.End, Protein))
        for List in GenesByChromosome:
            List.sort()
        for List in RGenesByChromosome:
            List.sort()
        for ListOfLists in (GenesByChromosome, RGenesByChromosome):
            for ProteinList in ListOfLists:
                for ProteinIndex in range(len(ProteinList)):
                    Protein = ProteinList[ProteinIndex][2]
                    Protein.Master = None
                    # Consider adding this protein as an isoform of a preceding protein:
                    OverlapFlag = 0
                    for OldProteinIndex in range(ProteinIndex - 1, -1, -1):
                        OldProtein = ProteinList[OldProteinIndex][2]
                        if OldProtein.End < Protein.Start:
                            break
                        for Exon in Protein.Exons:
                            for OldExon in OldProtein.Exons:
                                if Exon.End > OldExon.Start and Exon.Start < OldExon.End:
                                    OverlapFlag = 1
                                    break
                        if OverlapFlag:
                            if OldProtein.Master:
                                OldProtein.Master.IsoformCount += 1
                                Protein.Master = OldProtein.Master
                            else:
                                OldProtein.IsoformCount += 1
                                Protein.Master = OldProtein
                            Protein.IsoformCount = 0
                            break
                    if not OverlapFlag:
                        Protein.IsoformCount = 1
        IsoCountHistogram = {}
        ExonCountHistogram = {}
        for ListOfLists in (GenesByChromosome, RGenesByChromosome):
            for ProteinList in ListOfLists:
                for Tuple in ProteinList:
                    Protein = Tuple[2]
                    ExonCount = len(Protein.Exons)
                    ExonCountHistogram[ExonCount] = ExonCountHistogram.get(ExonCount, 0) + 1
                    if Protein.IsoformCount:
                        IsoCountHistogram[Protein.IsoformCount] = IsoCountHistogram.get(Protein.IsoformCount, 0) + 1
        Keys = ExonCountHistogram.keys()
        print "\nExon count histogram:"
        for Key in Keys:
            print "%s\t%s\t"%(Key, ExonCountHistogram[Key])
        Keys = IsoCountHistogram.keys()
        print "\nIsoform count histogram:"
        for Key in Keys:
            print "%s\t%s\t"%(Key, IsoCountHistogram[Key])
        
                        
    def SaveHits(self, GeneID, FileName):
        """
        Find all the hits to the gene of interest.
        """
        print "Find hits to gene %s..."%GeneID
        HitsFile = open(FileName, "wb")
        File = open(XGVKG0, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                HitGeneID = int(Bits[13])
                PValue = float(Bits[10])
            except:
                continue
            if HitGeneID != GeneID:
                continue
            # Save the hits, for reference:
            HitsFile.write(FileLine)
        HitsFile.close()
    def PlotHits(self, GeneID, ProteinIDList, HitFileName):
        #Protein = self.Proteins[ProteinID]
        File = open(HitFileName, "rb")
        Gene = None
        FullNotes = ""
        ##############################
        # Read the hits from the file:        
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                GeneNumber = float(Bits[13])
            except:
                continue
            #################
            # Create gene object, if necessary:
            if not Gene:
                Gene = GeneClass()
                Gene.GeneNumber = GeneID
                Gene.ChromosomeNumber = int(Bits[16])
                Gene.Strand = int(Bits[17])
                Gene.Name = Bits[3]
                if Gene.Name[0] == '"':
                    Gene.Name = Gene.Name[1:-1]
                CommaBits = Gene.Name.split(",")
                (Start, End) = CommaBits[-1].split("-")
                Gene.Start = int(Start)
                Gene.End = int(End)
                self.Genes[GeneNumber] = Gene
            Peptide = Bits[2][2:-2]
            ##########################
            # Remember this peptide:
            (Start, End) = Bits[18].split("-")
            Start = int(Start)
            End = int(End)
            if Bits[20]:
                # spliced!
                DashBillions = Bits[20].split("-")
                if len(DashBillions) > 2:
                    continue # skip three-exon hit!
                (Donor, Acceptor) = DashBillions
                Donor = int(Donor)
                Acceptor = int(Acceptor)
                Key = (Start, min(Donor, Acceptor), max(Donor, Acceptor), End)
            else:
                Key = (Start, End)
            Gene.HitDetails[Key] = Gene.HitDetails.get(Key, 0) + 1
            Gene.DetailPeptides[Key] = "%s.%s.%s"%(Bits[2][0], Bits[19], Bits[2][-1])  #Bits[2]
            Gene.PeptideCategories[Key] = Bits[21]
            Gene.PeptideBestScores[Key] = max(float(Bits[5]), Gene.PeptideBestScores.get(Key, -99))
        ################################
        # Plot the known protein and the hits.
        Plotter = GeneBrowserPlot.PlotterClass()
        Plotter.LabelPositionsFlag = 0
        SpectrumCount = 0
        PeptideCount = 0
        for (Key, Count) in Gene.HitDetails.items():
            SpectrumCount += Count
            PeptideCount += 1
        #Plotter.AddHeaderLine("%s %s: %s spectra (%s peptides)"%(Gene.GeneNumber, Gene.Name, SpectrumCount, PeptideCount))
        #
        for ProteinID in ProteinIDList:
            Protein = self.Proteins[ProteinID]
            ProteinTrack = Plotter.AddTrack("%s"%Protein.ProteinName)
            ProteinTrack.Color = GeneBrowserPlot.Colors.Black
            ProteinTrack.RequiredFlag = 1
            for ExonIndex in range(len(Protein.Exons)):
                Exon = Protein.Exons[ExonIndex]
                Feature = ProteinTrack.AddFeature(Exon.Start, Exon.End)
                Feature.Type = GeneBrowserPlot.FeatureTypes.Fat
                if ExonIndex < len(Protein.Exons) - 1:
                    Feature = ProteinTrack.AddFeature(Exon.End, Protein.Exons[ExonIndex + 1].Start)
                    Feature.Type = GeneBrowserPlot.FeatureTypes.FatDashed
                    Feature.Color = GeneBrowserPlot.Colors.Grey
        ########################################################################
        # Add tracks for unknown peptides, sorted by spectrum-count.
        SortedList = []
        for (Key, Count) in Gene.HitDetails.items():
            SortedList.append((Count, Key))
        SortedList.sort()
        SortedList.reverse()
        TrackLabels = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
        for SortedListPos in range(len(SortedList)):
            (Count, Key) = SortedList[SortedListPos]
            Peptide = Gene.DetailPeptides[Key]
            Track = Plotter.AddTrack("") #Gene.PeptideCategories[Key])
            Track.RequiredFlag = 1
            Track.Color = GeneBrowserPlot.Colors.Label
            Track.LabelLineFlag = 0
            if len(Key) == 2:
                Feature = Track.AddFeature(Key[0], Key[1])
                Feature.Label = "%s (%s)"%(Peptide, Count)
                if Gene.Strand == 1:
                    ReadingFrame = Key[0] % 3
                else:
                    ReadingFrame = (Key[1] - 1) % 3
                #Feature.Color = ReadingFrameColors[ReadingFrame]                
            else:
                # Intron included!
                Feature = Track.AddFeature(Key[0], min(Key[1],Key[2]))
                if Gene.Strand == 1:
                    ReadingFrame = Key[0] % 3
                else:
                    ReadingFrame = (min(Key[1], Key[2]) - 1) % 3
                #Feature.Color = ReadingFrameColors[ReadingFrame]
                
                Feature = Track.AddFeature(min(Key[1],Key[2]), max(Key[1],Key[2]))
                Feature.Type = GeneBrowserPlot.FeatureTypes.Dashed
                Feature.Color = GeneBrowserPlot.Colors.Grey
                Feature = Track.AddFeature(max(Key[1], Key[2]), Key[3])
                Feature.Label = "%s (%s)"%(Peptide, Count)
                if Gene.Strand == 1:
                    ReadingFrame = max(Key[1], Key[2]) % 3
                else:
                    ReadingFrame = (Key[3] - 1) % 3
                #Feature.Color = ReadingFrameColors[ReadingFrame]
            FullNotes += "%s: %s (%s)\n"%(TrackLabels[SortedListPos], Peptide, Count)
        FileName = "Gene.%s.png"%Gene.GeneNumber
        Plotter.Plot(FileName)
        print FullNotes
    def SpewProteinNames(self):
        File = open("ReportPresentProteins.txt", "rb")
        for FileLine in File.xreadlines():
            Bits = list(FileLine.strip().split("\t"))
            try:
                ProteinID = int(Bits[1])
            except:
                continue
            Bits = list(Bits)
            Bits[0] = self.Proteins[ProteinID].ProteinName
            print string.join(Bits, "\t")
        File.close()
def CountDistinctPeptides():
    """
    Count the unique spectra + peptides in each file.
    (VERY high-level analysis)
    """
    for X in range(1, 10):
        File = open("XGVKG%s.txt"%X, "rb")
        OldSpectrum = None
        SpectrumCount = 0 
        PepDict = {}
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if len(Bits)<3:
                continue
            Spectrum = (Bits[0], Bits[1])
            # (optional) Skip over additional hits for the same spectrum
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            Peptide = Bits[2][2:-2]
            if len(Peptide) < 8:
                continue
            PepDict[Peptide] = 1
            SpectrumCount += 1
        print X, len(PepDict.keys()), SpectrumCount
        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running non-optimized)"
    if len(sys.argv)<2:
        print "Please provide an analysis command ('splice' or 'full' are nice)"
        sys.exit(-1)
    Command = sys.argv[1].lower()
    if Command == "count":
        CountDistinctPeptides()
        sys.exit()
    Analyzer = AnalyzerClass()
    if Command == "splice":
        # Alt splicing report:
        Analyzer.LoadProteinMappings("GeneMappings.Best.txt")
        Analyzer.ReportAlternativeSplicing()
        sys.exit()
    if Command == "hypo":
        # hypothetical protein annotation report:
        Analyzer.LoadProteinMappings("GeneMappings.Best.txt")
        Analyzer.ReportHypotheticalProteinCoverage()
        sys.exit()
    if Command == "snp":
        # SNP report:
        Analyzer.LoadProteinMappings("GeneMappings.Best.txt")
        Analyzer.ReportSNPs()
        sys.exit()
    if Command == "iso":
        # Report on the exon counts and isoform counts of known proteins:
        Analyzer.LoadProteinMappings("GeneMappings.Best.txt")
        Analyzer.ReportKnownProteinIsoforms()
        sys.exit()
    if Command == "plot":
        # Given a PROTEIN ID and a GENE ID, we want to plot the protein and plot all the peptides
        # from that gene ID (whether they hit the protein or not!)
        GeneID = int(sys.argv[2])
        ProteinIDList = []
        for X in range(3, len(sys.argv)):
            ProteinID = int(sys.argv[X])
            ProteinIDList.append(ProteinID)
        HitFileName = "Hits.%s.txt"%GeneID
        # Don't troll for hits if we already saved them:
        if not os.path.exists(HitFileName):
            Analyzer.SaveHits(GeneID, HitFileName)
        Analyzer.LoadProteinMappings("GeneMappings.Best.txt")
        Analyzer.PlotHits(GeneID, ProteinIDList, HitFileName)
        #Analyzer.PlotNovelExons(ProteinID, GeneID)
        sys.exit()
    if Command == "spewproteinnames":
        Analyzer.LoadProteinMappings("GeneMappings.Best.txt")
        Analyzer.SpewProteinNames()
    if Command == "full":
        Analyzer.DBPath = os.path.join("Database", "IPIv315.trie")
        Analyzer.IndexFilePath = os.path.join("Database", "IPIv315.index")
        print "Load gene mappings:"
        #Analyzer.LoadProteinMappings("GeneMapperOutput\\GeneMappings.E.txt")
        Analyzer.LoadProteinMappings("GeneMappings.Best.txt")
        print "Sort chromosome proteins:"
        Analyzer.SortChromosomeProteins()
        print "Load protein sequences:"
        Analyzer.LoadSequences(Analyzer.DBPath, Analyzer.IndexFilePath)
        print "Analyze..."
        Analyzer.SelectPresentProteins()
        Analyzer.ReportPresentProteins()
        Analyzer.ReportProteinCoverage()
        # Delete old 'unexplained gene' stuff so we can start a fresh analysis:
        os.system("del /q NovelPeptides\\*")
        Analyzer.FindUnexplainedGenes()
        Analyzer.ReportUnexplainedGenes()
        sys.exit(1)
    print "Unknown command:", Command