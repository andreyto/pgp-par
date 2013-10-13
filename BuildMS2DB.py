"""
Code to support input of various data into .ms2db file format.
Simplest case: Import a FASTA file
Fancy case: Import swiss-prot format (also used by IPI)
Fancy case: Import genbank format
Fanciest case: Import data from various intron and exon sources (e.g. ESTs,
 de novo gene predictions).
"""
import os
import sys
import struct
import time
import getopt
import DNA

UsageInfo = """
BuildMS2DB.py - Generate a .ms2db format database from input.
Arguments:
-d [FILENAME]: Database to convert
-f [FORMAT]: Source database format.  Options are:
  fasta - Import a fasta file.  Due to the limitations of FASTA format,
    no exon information or cross-references can be imported.
  knowngene - Import genes from a tab-delimited format, based on the
    "KnownGene.txt" file from the UCSC genome browser.  See code
    documentation for format details.
  merge - Merge <Gene> records from all the .ms2db files in a directory.
-w [FILENAME]: Output filename.
-q: Quick test flag (for debugging)

Example:
BuildMS2db.py -d Database\Lens.fasta -f FASTA -w Database\Lens.ms2db
BuildMS2db.py -d F:\HGenome\KnownGene.txt -f KnownGene -w Database\KnownGene.ms2db -q
"""

class KGBits:
    Name = 0
    Chromosome = 1
    Strand = 2
    TranscriptionStart = 3
    TranscriptionEnd = 4
    CodingStart = 5
    CodingEnd = 6
    ExonCount = 7
    ExonStarts = 8
    ExonEnds = 9
    SprotID = 10
    UniqueID = 11

class GroupingGeneClass:
    def __init__(self):
        self.Name = None
        self.Exons = []
        self.CrossReferences = []
    def BuildFromExon(self, Exon):
        self.Exons.append(Exon)
        Exon.OwningGene = self
        for BackExon in Exon.BackwardExons:
            if BackExon.OwningGene:
                if BackExon.OwningGene == self:
                    continue
                else:
                    print "** Error: Trying to assimilate an exon with an owning gene!"
                    print self
                    print BackExon.OwningGene
                    raise ValueError
            self.AccumulateExons(BackExon)
        for ForwardExon in Exon.ForwardExons:
            if ForwardExon.OwningGene:
                if ForwardExon.OwningGene == self:
                    continue
                else:
                    print "** Error: Trying to assimilate an exon with an owning gene!"
                    print self
                    print BackExon.OwningGene
                    raise ValueError
            self.AccumulateExons(ForwardExon)
    def AccumulateExons(self, Exon):
        if Exon not in self.Exons:
            self.Exons.append(Exon)
            Exon.OwningGene = self
        for BackExon in Exon.BackwardExons:
            if BackExon.OwningGene:
                if BackExon.OwningGene == self:
                    continue
                else:
                    print "** Error: Trying to assimilate an exon with an owning gene!"
                    print self
                    print BackExon.OwningGene
                    raise ValueError
            self.AccumulateExons(BackExon)
        for ForwardExon in Exon.ForwardExons:
            if ForwardExon.OwningGene:
                if ForwardExon.OwningGene == self:
                    continue
                else:
                    print "** Error: Trying to assimilate an exon with an owning gene!"
                    print self
                    print BackExon.OwningGene
                    raise ValueError
            self.AccumulateExons(ForwardExon)

class ExonClass:
    def __init__(self, Start, End, ReadingFrame):
        self.Start = Start
        self.End = End
        self.ReadingFrame = ReadingFrame
        self.ForwardExons = [] # "forward" in the "genome residue number is larger" sense!
        self.BackwardExons = []
        self.OwningGene = None
    def LinkForward(self, Exon):
        # SanityChecking:
        if Exon.Start <= self.Start:
            print "*Error: Exon %s-%s can't link forward to exon %s-%s"%(self.Start, self.End, Exon.Start, Exon.End)
            raise ValueError
        if Exon not in self.ForwardExons:
            self.ForwardExons.append(Exon)
        if self not in Exon.BackwardExons:
            Exon.BackwardExons.append(self)
    def InheritForwardEdges(self, DeadExon):
        # Sanity check:
        if self.End != DeadExon.End:
            print "** Error: Exon %s-%s can't inherit forward edges from %s-%s"%(self.Start, self.End, DeadExon.Start, DeadExon.End)
            raise ValueError
        for OtherForwardExon in DeadExon.ForwardExons:
            if OtherForwardExon.Start < self.End:
                print "** Error: Exon %s-%s can inherit link to %s-%s from %s-%s"%(self.Start, self.End,
                   OtherForwardExon.Start, OtherForwardExon.End, DeadExon.Start, DeadExon.End)
                raise ValueError
            # We get their link:
            if OtherForwardExon not in self.ForwardExons:
                self.ForwardExons.append(OtherForwardExon)
            # And OtherForwardExon stops linking to DeadExon:
            OtherForwardExon.BackwardExons.remove(DeadExon)
            if self not in OtherForwardExon.BackwardExons:
                OtherForwardExon.BackwardExons.append(self)
    def InheritBackwardEdges(self, DeadExon):
        if self.Start != DeadExon.Start:
            print "** Error: Exon %s-%s can't inherit forward edges from %s-%s"%(self.Start, self.End, DeadExon.Start, DeadExon.End)
            raise ValueError
        for OtherBackwardExon in DeadExon.BackwardExons:
            if OtherBackwardExon.End > self.Start:
                print "** Error: Exon %s-%s can inherit link back to %s-%s from %s-%s"%(self.Start, self.End,
                   OtherBackwardExon.Start, OtherBackwardExon.End, DeadExon.Start, DeadExon.End)
                raise ValueError
            # We get the link:
            if OtherBackwardExon not in self.BackwardExons:
                self.BackwardExons.append(OtherBackwardExon)
            # And OtherBackwardExon stops linking to DeadExon:
            OtherBackwardExon.ForwardExons.remove(DeadExon)
            if self not in OtherBackwardExon.ForwardExons:
                OtherBackwardExon.ForwardExons.append(self)
    def __str__(self):
        return "<exon %s-%s>"%(self.Start, self.End)

class KGCrossReference:
    def __init__(self):
        self.Name = None
        self.ExonIndexes = []
        
class StrandClass:
    def __init__(self):
        self.Chromosome = None
        self.ForwardFlag = 1
        self.ExonsByInterval = {}
        self.ExonsByStart = {}
        self.ExonsByEnd = {}
        self.Genes = []
        self.GroupingGenes = []
    def __str__(self):
        return "<Strand %s-%s>"%(self.Chromosome, self.ForwardFlag)
    
class ReferenceGeneClass:
    def __init__(self, Bits):
        self.ParseBits(Bits)
    def __str__(self):
        return "<gene %s>"%self.Name
    def ParseBits(self, Bits):
        self.Name = Bits[KGBits.Name]
        self.SprotID = Bits[KGBits.SprotID]
        self.UniqueID = Bits[KGBits.UniqueID]
        if Bits[KGBits.Strand] == "+":
            self.ForwardFlag = 1
        elif Bits[KGBits.Strand] == "-":
            self.ForwardFlag = 0
        
        # Prepare for exons:
        CDStart = int(Bits[KGBits.CodingStart])
        CDEnd = int(Bits[KGBits.CodingEnd])
        ExonCount = int(Bits[KGBits.ExonCount])
        ExonStarts = Bits[KGBits.ExonStarts]
        if ExonStarts[0] == '"':
            ExonStarts = ExonStarts[1:-1]
        ExonStarts = list(ExonStarts.split(","))
        ExonEnds = Bits[KGBits.ExonEnds]
        if ExonEnds[0] == '"':
            ExonEnds = ExonEnds[1:-1]
        ExonEnds = list(ExonEnds.split(","))
        # Iterate over exons in gene order, keeping track of reading frame:
        ExtraAminoAcids = 0
        self.Exons = []
        if self.ForwardFlag:
            ExonIndexes = range(ExonCount)
        else:
            ExonIndexes = range(ExonCount - 1, -1, -1)
        for ExonIndex in ExonIndexes:
            Start = int(ExonStarts[ExonIndex])
            End = int(ExonEnds[ExonIndex])
            #print "Exon %s: %s-%s"%(ExonIndex, Start, End)
            # Ignore exons outside the translation portion:
            if (Start >= CDEnd) or (End <= CDStart):
                continue
            OldExtraAA = ExtraAminoAcids
            if self.ForwardFlag:
                ExonStart = max(Start, CDStart)
                ExonEnd = min(End, CDEnd)
                CodonStart = ExonStart
                if ExtraAminoAcids == 1:
                    CodonStart += 2
                elif ExtraAminoAcids == 2:
                    CodonStart += 1
                ExonKey = (ExonStart, ExonEnd, CodonStart % 3)
                self.Exons.append(ExonKey)
                ExtraAminoAcids = (ExonEnd - CodonStart) % 3
            else:
                ExonStart = max(Start, CDStart)
                ExonEnd = min(End, CDEnd)
                CodonStart = ExonEnd - 1
                if ExtraAminoAcids == 1:
                    CodonStart -= 2
                elif ExtraAminoAcids == 2:
                    CodonStart -= 1
                ExonKey = (ExonStart, ExonEnd, CodonStart % 3)
                self.Exons.append(ExonKey)
                ExtraAminoAcids = (CodonStart - ExonStart + 1) % 3
            #print "  ->Exon %s-%s"%(ExonStart, ExonEnd)
            #print "%s (For%s) exon %s:"%(self.Name, self.ForwardFlag, ExonIndex)
            #print "  %s-%s (ExtraAA %s) becomes %s-%s RF%s ExtraAA %s"%(Start, End, OldExtraAA, ExonStart, ExonEnd, ExonKey[2], ExtraAminoAcids)
            #print "  ->CodonStart %s"%CodonStart
            #print "ForwardFlag%s: (%s, %s) with CD (%s, %s) becomes (%s, %s)"%(self.ForwardFlag, Start, End, CDStart, CDEnd, ExonStart, ExonEnd)
            #print "RF %s (Start %% 3 == %s"%(ExonKey[2], ExonKey[0] % 3)
    def DebugPrint(self, GenomePath):
        GenomeFile = open(GenomePath, "rb")
        print ">>>Debug print Gene '%s' For%s"%(self.Name, self.ForwardFlag)
        for ExonIndex in range(len(self.Exons)):
            (Start, End, ReadingFrame) = self.Exons[ExonIndex]
            GenomeFile.seek(Start)
            DNASequence = GenomeFile.read(End - Start)
            if self.ForwardFlag:
                if Start % 3 == ReadingFrame:
                    pass
                elif Start % 3 == (ReadingFrame + 1) % 3:
                    DNASequence = DNASequence[2:]
                else:
                    DNASequence = DNASequence[1:]
            else:
                DNASequence = DNA.ReverseComplement(DNASequence)
                if (End - 1) % 3 == ReadingFrame:
                    pass
                elif (End - 1) % 3 == (ReadingFrame - 1) % 3:
                    DNASequence = DNASequence[2:]
                else:
                    DNASequence = DNASequence[1:]
            Protein = DNA.Translate(DNASequence)
            print "  Exon %s (%s-%s) RF%s"%(ExonIndex, Start, End, ReadingFrame)
            print "    Sequence %s...%s"%(Protein[:20], Protein[-20:])
                                                    
class MS2DBGenerator:
    def __init__(self):
        self.OutputFileName = os.path.join("Database", "Database.ms2db")
        self.InputFileName = None
        self.GeneCount = 0
        self.ExonCount = 0
        self.GenomeDir = None
        self.QuickParseLineLimit = None
    def GetXMLTimeStamp(self):
        """
        Return an XML datetime, like this: 2002-05-30T09:00:00
        """
        TimeTuple = time.localtime(time.time())
        Str = "%s-%s-%sT%02d:%02d:%02d"%(TimeTuple[0], TimeTuple[1], TimeTuple[2],
            TimeTuple[3], TimeTuple[4], TimeTuple[5])
        return Str
    def GetDatabaseTag(self):
        Str = """<Database CreationTime="%s" CreatedBy="BuildMS2DB.py" SourceDBFile="%s">\n"""%(self.GetXMLTimeStamp(), self.InputFileName)
        return Str
    def OpenDBFiles(self):
        if not self.InputFileName:
            print UsageInfo
            sys.exit(-1)
        if not os.path.isdir(self.InputFileName):
                self.InputFile = open(self.InputFileName, "rb")
        self.OutputFile = open(self.OutputFileName, "wb")
        self.IndexFileName = os.path.splitext(self.OutputFileName)[0] + ".ms2index"
        self.IndexFile = open(self.IndexFileName, "wb")
    def GenerateByMerge(self):
        self.OutputFile.write("<Database>\n")
        for FileName in os.listdir(self.InputFileName):
            Path = os.path.join(self.InputFileName, FileName)
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension.lower() != ".ms2db":
                continue
            File = open(Path, "rb")
            for FileLine in File.xreadlines():
                if FileLine.find("<Database")!=-1:
                    continue
                if FileLine.find("</Database>")!=-1:
                    continue
                self.OutputFile.write(FileLine)
            File.close()
        self.OutputFile.write("\n</Database>\n")
    def GenerateDatabase(self):
        Format = self.InputFormat.lower()
        if Format == "fasta":
            self.GenerateFromFASTA()
        elif Format == "knowngene":
            self.GenerateFromKnownGenes()
        elif Format == "merge":
            self.GenerateByMerge()
        else:
            print "* Warning: Unhandled format '%s'"%self.InputFormat
            print UsageInfo
            sys.exit(-1)
    def FinishFASTARecord(self, SourceFilePos, Name, Sequence):
        "Helper for GenerateFromFASTA"
        FilePos = self.OutputFile.tell()
        self.OutputFile.write("""<Gene Name="%s" ExonCount="1">\n"""%Name)
        self.OutputFile.write("""  <Exon Index="0">\n""")
        self.OutputFile.write("""    <ExonSequence Length="%s">"""%len(Sequence))
        self.OutputFile.write(Sequence)
        self.OutputFile.write("</ExonSequence>\n")
        self.OutputFile.write("  </Exon>\n")
        self.OutputFile.write("</Gene>\n")
        IndexBlock = struct.pack("<qi80s", SourceFilePos, FilePos, Name)
        self.IndexFile.write(IndexBlock)
        self.ExonCount += 1
        self.GeneCount += 1
    def GenerateFromFASTA(self):
        """
        The simplest possible case - grab FASTA input, spit out a minimalistic
        .ms2db database.
        """
        Tag = self.GetDatabaseTag()
        self.OutputFile.write(Tag)
        CurrentName = None
        CurrentSequence = None
        LineNumber = 0
        SourceFilePos = 0
        for FileLine in self.InputFile.xreadlines():
            LineNumber += 1
            if FileLine[0] == ">":
                if CurrentSequence:
                    self.FinishFASTARecord(CurrentRecordFilePos, CurrentName, CurrentSequence)
                CurrentSequence = ""
                CurrentName = FileLine[1:].strip()
                CurrentRecordFilePos = SourceFilePos
            else:
                CurrentSequence += FileLine.strip()
            SourceFilePos += len(FileLine)
        if CurrentSequence:
            self.FinishFASTARecord(CurrentRecordFilePos, CurrentName, CurrentSequence)
        self.OutputFile.write("</Database>\n")
        print "--> Wrote %s genes and %s exons, from %s input lines."%(self.GeneCount, self.ExonCount, LineNumber)
    def ParseCommandLine(self):
        (Options, Args) = getopt.getopt(sys.argv[1:], "d:f:w:q")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                # -d database filename
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(-1)
                self.InputFileName = Value
            elif Option == "-f":
                self.InputFormat = Value
            elif Option == "-w":
                self.OutputFileName = Value
            elif Option == "-q":
                self.QuickParseLineLimit = 25
        if not self.InputFileName:
            print UsageInfo
            sys.exit(-1)
        # Default: Genome sequence is in the same directory as the input file.
        if not self.GenomeDir:
            self.GenomeDir = os.path.split(self.InputFileName)[0]
    def GenerateFromKnownGenes(self):
        """
        Parse KnownGene.txt (or an equivalently-formatted file) to accumulate
        exons, introns, and cross-references to protein databases.
        Subfunctions of this function have the "KG" prefix.
        Outline:
        - Read input file, accumulate strands and exons
        - Split strand exons to get a minimal disjoint set
        - Read input file, accumulate edges
        - Group exons into genes using edges
        - Read input file, add cross-references to gene records
        """
        self.Strands = {}
        # Read input files, accumulate strands and genes and exons:
        print "KGParseExons:"
        self.KGParseExons()
        self.StrandKeys = self.Strands.keys()
        self.StrandKeys.sort()
        #self.KGDebugPrintStrandExons(0)
        # Split strand exons, to get a minimal disjoint set:
        print "KGSplitStrandExons:"
        self.KGSplitStrandExons()
        #self.KGDebugPrintStrandExons(1)
        # Assemble exons into genes, by linking together their edges:
        print "KGGroupExonsIntoGenes:"
        self.KGGroupExonsIntoGenes()
        # Build cross-references for each gene; specify the exons it covers
        print "KGBuildCrossReferences:"
        self.KGBuildCrossReferences()
        # Write out XML:
        print "KGGenerateXML:"
        self.KGGenerateXML()
    def KGDebugPrintStrandExons(self, SplitCompleteFlag):
        """
        For debugging the assimilation of known genes:
        Print all the exons on each strand.
        """
        for Strand in self.Strands.values():
            print ">>>Strand %s exons:"%Strand
            if SplitCompleteFlag:
                ExonList = Strand.ExonsByEnd.values()
            else:
                ExonList = Strand.ExonsByInterval.values()
            for Exon in ExonList:
                print " Exon %s-%s RF %s"%(Exon.Start, Exon.End, Exon.ReadingFrame)
                for LinkExon in Exon.ForwardExons:
                    print "  >>Forward to exon %s-%s RF %s"%(LinkExon.Start, LinkExon.End, LinkExon.ReadingFrame)
                for LinkExon in Exon.BackwardExons:
                    print "  <<Backward to exon %s-%s RF %s"%(LinkExon.Start, LinkExon.End, LinkExon.ReadingFrame)
    def KGBuildCrossReferences(self):
        """
        For each known gene, attach it to a GroupingGene, and find the finalized exon(s) which cover
        each exon of the known gene.
        """
        for StrandKey in self.StrandKeys:
            Strand = self.Strands[StrandKey]
            for Gene in Strand.Genes:
                # Find the GroupingGene that covers this gene's exons:
                if not Gene.Exons:
                    continue
                (Start, End, ReadingFrame) = Gene.Exons[0]
                if Strand.ForwardFlag:
                    Key = (End, ReadingFrame)
                    Exon = Strand.ExonsByEnd.get(Key, None)
                    if not Exon:
                        print "** Error: %s ExonsByStart has no entry for %s for %s"%(Strand, Key, Gene)
                        continue
                else:
                    Key = (Start, ReadingFrame)
                    Exon = Strand.ExonsByStart.get(Key, None)
                    if not Exon:
                        print "** Error: %s ExonsByStart has no entry for %s for %s"%(Strand, Key, Gene)
                        continue
                ContainingGene = Exon.OwningGene
                # Add a cross reference to the GroupingGene:
                CR = KGCrossReference()
                CR.Name = Gene.Name
                CR.ID = Gene.SprotID
                ContainingGene.CrossReferences.append(CR)
                # Iterate over the exons in the gene.  Find the exons from the GroupingGene which
                # cover them.
                for (Start, End, ReadingFrame) in Gene.Exons:
                    while 1:
                        (NewStart, NewEnd) = self.KGAddCrossReferenceExons(Strand, CR, Start, End, ReadingFrame, ContainingGene)
                        if NewStart >= End or NewEnd <= Start:
                            break
                        Start = NewStart
                        End = NewEnd
        # And now that known genes have been assimilated, we'll hazard a "name" for these
        # new genes:
        for Strand in self.Strands.values():
            for GeneIndex in range(len(Strand.GroupingGenes)):
                Gene = Strand.GroupingGenes[GeneIndex]
                if not Gene.CrossReferences:
                    print "* Warning: Gene %d %s has no crossreferences!"%(GeneIndex, Gene)
                    Gene.Name = "Gene %s on %s"%(GeneIndex, Strand)
                else:
                    Gene.Name = Gene.CrossReferences[0].Name # arbitrary, we could do better!
    def KGAddCrossReferenceExons(self, Strand, CR, Start, End, ReadingFrame, ContainingGene):
        """
        Look up the gene that covers part of [Start, End).  On forward strand, work from
        start to end; on reverse strand, work from end to start.
        Return the new interval-to-be-covered.
        """
        if Strand.ForwardFlag:
            Key = (Start, ReadingFrame)
            Exon = Strand.ExonsByStart[Key]
            CR.ExonIndexes.append(Exon.Index)
            return (Exon.End, End)
        else:
            Key = (End, ReadingFrame)
            Exon = Strand.ExonsByEnd[Key]
            CR.ExonIndexes.append(Exon.Index)
            return (Start, Exon.Start)
    def KGGenerateXML(self):
        Tag = self.GetDatabaseTag()
        self.OutputFile.write(Tag)
        for StrandKey in self.StrandKeys:
            Strand = self.Strands[StrandKey]
            ChromosomePath = os.path.join(self.GenomeDir, "%s.trie"%Strand.Chromosome)
            self.GenomeFile = open(ChromosomePath, "rb")
            for Gene in Strand.GroupingGenes:
                self.KGGenerateGeneXML(Strand, Gene)
            self.GenomeFile.close()
        self.OutputFile.write("</Database>\n\n")
    def KGGenerateGeneXML(self, Strand, Gene):
        # Write to the index file:
        FilePos = self.OutputFile.tell()
        IndexStr = struct.pack("<qi80s", 0, FilePos, Gene.Name)
        self.IndexFile.write(IndexStr)
        # Write the xml for this gene:
        GeneTag = """<Gene Name="%s" ExonCount="%s" Chromosome="%s" ForwardFlag="%s">\n"""%(Gene.Name, len(Gene.Exons), Strand.Chromosome, Strand.ForwardFlag)
        self.OutputFile.write(GeneTag)
        # Grab exon sequences (and prefixes and suffixes):
        for Exon in Gene.Exons:
            self.KGGetExonSequence(self.GenomeFile, Exon, Gene.ForwardFlag)
        for ExonIndex in range(len(Gene.Exons)):
            Exon = Gene.Exons[ExonIndex]
            ExonTag = """  <Exon Index="%s" Start="%s" End="%s" Length="%s">\n"""%(ExonIndex, Exon.Start, Exon.End, len(Exon.Sequence))
            self.OutputFile.write(ExonTag)
            self.OutputFile.write("    <ExonSequence>%s</ExonSequence>\n"%Exon.Sequence)
            if Gene.ForwardFlag:
                LinkExonList = Exon.BackwardExons
            else:
                LinkExonList = Exon.ForwardExons
            for LinkExon in LinkExonList:
                if LinkExon.End == Exon.Start or LinkExon.Start == Exon.End:
                    TagName = "ExtendsExon"
                else:
                    TagName = "LinkFrom"
                # Get bridging AA:
                if Exon.Prefix:
                    Codon = LinkExon.Suffix + Exon.Prefix
                    LinkAA = DNA.Translate(Codon)
                else:
                    LinkAA = ""
                if LinkAA:
                    AAStr = "AA=\"%s\""%LinkAA
                else:
                    AAStr = ""
                LinkTag = """    <%s Index="%s" %s>\n"""%(TagName, LinkExon.Index, AAStr)
                self.OutputFile.write(LinkTag)
            self.OutputFile.write("  </Exon>\n")
        for CrossReference in Gene.CrossReferences:
            CRTag = """  <CrossReference Database="Swiss-Prot" URL="http://ca.expasy.org/uniprot/%s" ID="%s">\n"""%(CrossReference.ID, CrossReference.ID)
            self.OutputFile.write(CRTag)
            IDString = ""
            for CRExonIndex in CrossReference.ExonIndexes:
                IDString += "%s, "%CRExonIndex
            IDString = IDString[:-2] # strip trailing comma+space
            self.OutputFile.write("""    <CRExons Index="%s" />\n"""%IDString)
            self.OutputFile.write("  </CrossReference>\n")
        self.OutputFile.write("</Gene>\n\n")
        # Drop exon sequences, to save space:
        for Exon in Gene.Exons:
            Exon.Sequence = None
    def KGGroupExonsIntoGenes(self):
        StrandCount = 0
        for Strand in self.Strands.values():
            ExonCount = 0
            GeneCount = 0
            # Iterate over the exons in the strand.  When you find an exon that's
            # not already part of a gene, create a new gene based around the exon.
            for Exon in Strand.ExonsByEnd.values():
                ExonCount += 1
                if not Exon.OwningGene:
                    NewGene = GroupingGeneClass()
                    NewGene.ForwardFlag = Strand.ForwardFlag
                    NewGene.BuildFromExon(Exon)
                    Strand.GroupingGenes.append(NewGene)
                    GeneCount += 1
            print "%s: Built %s genes for a total of %s exons"%(Strand, GeneCount, ExonCount)
            # Now let's sort the exons within each gene, so that they move FORWARD (along the protein)
            for Gene in Strand.GroupingGenes:
                SortedList = []
                for Exon in Gene.Exons:
                    SortedList.append((Exon.Start, Exon.End, Exon))
                SortedList.sort()
                if not Strand.ForwardFlag:
                    SortedList.reverse()
                Gene.Exons = []
                for Tuple in SortedList:
                    ExonIndex = len(Gene.Exons)
                    Gene.Exons.append(Tuple[-1])
                    Gene.Exons[-1].Index = ExonIndex
            StrandCount += 1
        print "Processed %s strands."%StrandCount
    def KGParseExons(self):
        File = open(self.InputFileName, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 100 == 0:
                print "Line %s..."%LineNumber
            if self.QuickParseLineLimit!=None and LineNumber > self.QuickParseLineLimit:
                break
            if FileLine[0] == "#":
                continue
            Bits = FileLine.split("\t")
            if len(Bits) < KGBits.UniqueID:
                continue
            if Bits[KGBits.Strand] == "+":
                ForwardFlag = 1
            elif Bits[KGBits.Strand] == "-":
                ForwardFlag = 0
            else:
                print "* Error: Line %s has illegal forward flag '%s'"%(LineNumber, Bits[KGBits.Strand])
            # Get the STRAND:
            Chromosome = Bits[KGBits.Chromosome]
            StrandKey = (Chromosome, ForwardFlag)
            if not self.Strands.has_key(StrandKey):
                 Strand = StrandClass()
                 Strand.Chromosome = Chromosome
                 Strand.ForwardFlag = ForwardFlag
                 self.Strands[StrandKey] = Strand
            else:
                Strand = self.Strands[StrandKey]
            Gene = ReferenceGeneClass(Bits)
            Strand.Genes.append(Gene)
            ##### Temp debugging:
            ####GenomeFilePath = os.path.join(self.GenomeDir, "%s.trie"%Chromosome)
            ####Gene.DebugPrint(GenomeFilePath)
            # Add the gene's exons to our exon-dictionary:
            PrevExon = None
            for ExonIndex in range(len(Gene.Exons)):
                ExonKey = Gene.Exons[ExonIndex]
                # Add an exon to Strand.ExonsByInterval, if necessary
                if Strand.ExonsByInterval.has_key(ExonKey):
                    Exon = Strand.ExonsByInterval[ExonKey]
                else:
                    Exon = ExonClass(ExonKey[0], ExonKey[1], ExonKey[2])
                    Strand.ExonsByInterval[ExonKey] = Exon
                #print "Exon %s: %s, prev %s"%(ExonIndex, Exon, PrevExon)
                # Add edges between exons, if necessary:
                if PrevExon:
                    if Strand.ForwardFlag:
                        PrevExon.LinkForward(Exon)
                        #if Exon not in PrevExon.ForwardExons:
                        #    PrevExon.ForwardExons.append(Exon)
                        #if PrevExon not in Exon.BackwardExons:
                        #    Exon.BackwardExons.append(PrevExon)
                    else:
                        Exon.LinkForward(PrevExon)
                        #if Exon not in PrevExon.BackwardExons:
                        #    PrevExon.BackwardExons.append(Exon)
                        #if PrevExon not in Exon.ForwardExons:
                        #    Exon.ForwardExons.append(PrevExon)
                # Remember this exon, so that we link to it the next cycle through the loop:
                PrevExon = Exon
    def KGSplitStrandExons(self):
        """
        For each strand: Find all pairs of exons (with compatible reading frame) which overlap.
        When overlap is found, split the overlapping exons into two (or three) sub-exons.
        """
        for StrandKey in self.StrandKeys:
            Strand = self.Strands[StrandKey]
            print "Split exons on %s..."%Strand
            ExonList = []
            for (Key, Exon) in Strand.ExonsByInterval.items():
                (Start, End, ReadingFrame) = Key
                ExonList.append(Exon)
            ExonIndexA = 0
            while ExonIndexA < len(ExonList):
                ExonA = ExonList[ExonIndexA]
                ExonADeletedFlag = 0
                ExonIndexB = ExonIndexA + 1
                while ExonIndexB < len(ExonList):
                    ExonB = ExonList[ExonIndexB]
                    if ExonA.ReadingFrame != ExonB.ReadingFrame:
                        ExonIndexB += 1
                        continue
                    if ExonB.Start >= ExonA.End:
                        ExonIndexB += 1
                        continue
                    if ExonA.Start >= ExonB.End:
                        ExonIndexB += 1
                        continue
                    #############################################################
                    # Delete redundant exons:
                    if ExonA.Start == ExonB.Start and ExonA.End == ExonB.End:
                        ExonA.InheritForwardEdges(ExonB)
                        ExonA.InheritBackwardEdges(ExonB)
                        del ExonList[ExonIndexB]
                        # DON'T increment index b, just continue
                        continue
                    if ExonA.Start == ExonB.Start:
                        if ExonA.End > ExonB.End:
                            # A-----
                            # B---
                            #     **
                            NewExon = ExonClass(ExonB.End, ExonA.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon)
                            NewExon.InheritForwardEdges(ExonA)
                            ExonB.InheritBackwardEdges(ExonA)
                            del ExonList[ExonIndexA]
                            ExonADeletedFlag = 1
                            break
                        else:
                            # A---
                            # B-----
                            #     **
                            NewExon = ExonClass(ExonA.End, ExonB.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon)
                            NewExon.InheritForwardEdges(ExonB)
                            ExonA.InheritBackwardEdges(ExonB)
                            del ExonList[ExonIndexB]
                            continue
                    elif ExonA.End == ExonB.End:
                        if ExonA.Start < ExonB.Start:
                            # A-----
                            # B  ---
                            #  11
                            NewExon = ExonClass(ExonA.Start, ExonB.Start, ExonA.ReadingFrame)
                            NewExon.InheritBackwardEdges(ExonA)
                            ExonList.append(NewExon)
                            ExonB.InheritForwardEdges(ExonA)
                            del ExonList[ExonIndexA]
                            ExonADeletedFlag = 1
                            break
                        else:
                            # A  ---
                            # B-----
                            #  11
                            NewExon = ExonClass(ExonB.Start, ExonA.Start, ExonA.ReadingFrame)
                            ExonList.append(NewExon)
                            NewExon.InheritBackwardEdges(ExonB)
                            ExonA.InheritForwardEdges(ExonB)
                            del ExonList[ExonIndexB]
                            continue
                    else:
                        if ExonA.Start < ExonB.Start and ExonA.End < ExonB.End:
                            # A------
                            # B   ------
                            #  111222333
                            NewExon1 = ExonClass(ExonA.Start, ExonB.Start, ExonA.ReadingFrame)
                            ExonList.append(NewExon1)
                            NewExon2 = ExonClass(ExonB.Start, ExonA.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon2)
                            NewExon3 = ExonClass(ExonA.End, ExonB.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon3)
                            NewExon1.InheritBackwardEdges(ExonA)
                            NewExon2.InheritForwardEdges(ExonA)
                            NewExon2.InheritBackwardEdges(ExonB)
                            NewExon3.InheritForwardEdges(ExonB)
                            # Delete B first, since index B > index A:
                            del ExonList[ExonIndexB]
                            del ExonList[ExonIndexA]
                            ExonADeletedFlag = 1
                            break
                        if ExonA.Start < ExonB.Start and ExonA.End > ExonB.End:
                            # A---------
                            # B   ---
                            #  111   222
                            NewExon1 = ExonClass(ExonA.Start, ExonB.Start, ExonA.ReadingFrame)
                            ExonList.append(NewExon1)
                            NewExon2 = ExonClass(ExonB.End, ExonA.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon2)
                            NewExon1.InheritBackwardEdges(ExonA)
                            NewExon2.InheritForwardEdges(ExonA)
                            del ExonList[ExonIndexA]
                            ExonADeletedFlag = 1
                            break
                        if ExonA.Start > ExonB.Start and ExonA.End > ExonB.End:
                            # A    ------
                            # B ------
                            #   111222333
                            NewExon1 = ExonClass(ExonB.Start, ExonA.Start, ExonA.ReadingFrame)
                            ExonList.append(NewExon1)
                            NewExon2 = ExonClass(ExonA.Start, ExonB.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon2)
                            NewExon3 = ExonClass(ExonB.End, ExonA.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon3)
                            NewExon1.InheritBackwardEdges(ExonB)
                            NewExon2.InheritForwardEdges(ExonB)
                            NewExon2.InheritBackwardEdges(ExonA)
                            NewExon3.InheritForwardEdges(ExonA)
                            # Delete B first, since index B > index A:
                            del ExonList[ExonIndexB]
                            del ExonList[ExonIndexA]
                            ExonADeletedFlag = 1
                            break
                        if ExonA.Start > ExonB.Start and ExonA.End < ExonB.End:
                            # A    ---
                            # B ---------
                            #   111   222
                            NewExon1 = ExonClass(ExonB.Start, ExonA.Start, ExonA.ReadingFrame)
                            ExonList.append(NewExon1)
                            NewExon2 = ExonClass(ExonA.End, ExonB.End, ExonA.ReadingFrame)
                            ExonList.append(NewExon2)
                            NewExon1.InheritBackwardEdges(ExonB)
                            NewExon2.InheritForwardEdges(ExonB)
                            del ExonList[ExonIndexB]
                            continue
                # We finihed our B loop:
                if not ExonADeletedFlag:
                    ExonIndexA += 1
            # We finished A loop; ExonList is final.
            # (re)generate the exon dictionaries for the strand!
            for Exon in ExonList:
                Key = (Exon.Start, Exon.ReadingFrame)
                if Strand.ExonsByStart.has_key(Key):
                    print "* Error: Strand %s has redundant start-point exon %s!"%(Strand, Key)
                Strand.ExonsByStart[Key] = Exon
                Key = (Exon.End, Exon.ReadingFrame)
                if Strand.ExonsByEnd.has_key(Key):
                    print "* Error: Strand %s has redundant start-point exon %s!"%(Strand, Key)
                Strand.ExonsByEnd[Key] = Exon
            # Add "adjacent-edges" wherever necessary:
            for Exon in ExonList:
                CheckForwardKey = (Exon.End, Exon.ReadingFrame)
                LinkExon = Strand.ExonsByStart.get(CheckForwardKey, None)
                if LinkExon:
                    Exon.LinkForward(LinkExon)
                CheckBackwardKey = (Exon.Start, Exon.ReadingFrame)
                LinkExon = Strand.ExonsByEnd.get(CheckBackwardKey, None)
                if LinkExon:
                    LinkExon.LinkForward(Exon)
    
    def KGGetExonSequence(self, GenomeFile, Exon, ForwardFlag):
        """
        Given an exon, read its sequence (DNA prefix, coded protein, and DNA suffix)
        from the genome file.
        """
        GenomeFile.seek(Exon.Start)
        DNASequence = GenomeFile.read(Exon.End - Exon.Start)
        if ForwardFlag:
            if Exon.Start % 3 == Exon.ReadingFrame:
                Exon.Prefix = ""
            elif Exon.Start % 3 == (Exon.ReadingFrame + 1) % 3:
                Exon.Prefix = DNASequence[:2]
                DNASequence = DNASequence[2:]
            else:
                Exon.Prefix = DNASequence[:1]
                DNASequence = DNASequence[1:]
        else:
            DNASequence = DNA.ReverseComplement(DNASequence)
            if (Exon.End - 1) % 3 == Exon.ReadingFrame:
                Exon.Prefix = ""
            elif (Exon.End - 1) % 3 == (Exon.ReadingFrame - 1) % 3:
                Exon.Prefix = DNASequence[:2]
                DNASequence = DNASequence[2:]
            else:
                Exon.Prefix = DNASequence[:1]
                DNASequence = DNASequence[1:]
        # Set suffix:
        DNALen = len(DNASequence) - len(Exon.Prefix)
        Leftovers = DNALen % 3
        if Leftovers:
            Exon.Suffix = DNASequence[-Leftovers:]
        else:
            Exon.Suffix = ""
        Protein = DNA.Translate(DNASequence)
        # Remove stop codon, if present:
        if Protein and Protein[-1] == "X":
            Protein = Protein[:-1]
        Exon.Sequence = Protein
        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    MS2DB = MS2DBGenerator()
    MS2DB.ParseCommandLine()
    MS2DB.OpenDBFiles()
    MS2DB.GenerateDatabase()