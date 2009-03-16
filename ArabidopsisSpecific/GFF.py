"""GFF.py
Methods for parsing and dealing with the GFF format, or the GFF-like gbrowse format
Be careful!!!!!! GFFL and GFF are different, but have nearly identical methods/classes

Notes for the GFF format are on:
http://www.arabidopsis.org/cgi-bin/gbrowse/arabidopsis/?help=annotation

reference = B0511
EST	yk260e10.5	15569-15724
EST	yk672a12.5	537-618,3187-3294

NOTE: GFF Format lists CDS indicies in numerical order.  Meaning that things on the reverse strand
are still written like:
Chr1	TAIR7	CDS	12245807	12246014	.	-	0	Parent=AT1G33780.1
the opposite is true for the GFFL.

"""

#[1,2,3,4,5,C,M]
ArabChrList = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"]
TAIR7ChrSize = [30432563, 19705359, 23470805, 18585042, 26992728, 154478, 366924]

class ObservedSequenceClass:
    """container for GFFL lines, possibly spliced items.  This can get pretty messy, so be careful
    Observed Sequences are written in the (aweful) GFFL differently for items belonging to
    different strands.  For Example, on the reverse strand we write
    525..521,424..388
    See how the numbers are totally backwards.  both within a unit (525..521) and between?
    The variable self.Sequences is mean to be a numberical ordering.  Always lowest to highest.
    """
    def __init__(self):
        self.Sequences = [] # [(start, stop), (start, stop), ..] values are 'int'
        self.Name = None
        self.Chromosome = None
        self.Frame = [] # one for each sequence array element
        self.Start = None #first base of this sequence, type = int
        self.Stop = None #last base of this sequence, type = int
    def PrintMe(self):
        print "**I'm ObservedSequenceClass object %s"%self.Name
        print "\tI'm on chromosome %s, start %d, stop %d"%(self.Chromosome, self.Start, self.Stop)
        print ""

class ExonClass:
    "container for genomic information"
    def __init__(self):
        self.Start = None #int  Exon.Start < Exon.Stop NO EXCEPTIONS.
        self.Stop = None #int
        self.TotalDNABases = None # int
        self.Phase = None #int
        self.Strand = None
        self.ParentTranscript = None
    def PrintMe(self):
        print "Exon: Start: %s, Stop %s"%(self.Start, self.Stop)
        print "\tTotal bases covered: %s"%self.TotalDNABases
        print "\tParentTranscript: %s"%self.ParentTranscript

class TranscriptClass:
    def __init__(self):
        self.Exons = [] #should be a sorted list of ExonClass objects
        self.Chromosome= None
        self.Strand = None
        self.Name = None #ID from the mRNA line
        self.ParentLocus = None
        self.CodingLength = None # length of unspliced DNA in the coding region, #DNA length from ATG to STOP
    def PrintMe(self):
        print "Print Transcript, %s"%self.Name
        print "Chr: %s, Strand %s, ParentLocus %s"%(self.Chromosome, self.Strand, self.ParentLocus)
        print "My coding length is %d bases from ATG to Stop"%self.CodingLength
        print "Exons in the transcript"
        for Exon in self.Exons:
            Exon.PrintMe()
        print "\n\n"

class LocusClass:
    def __init__(self):
        self.Name = None
        self.Chromosome = None
        self.FirstBase = None # Order on the chromosome, NOT start of transcription, type = int
        self.LastBase = None # Order on the chromosome, NOT start of transcription, type = int
        self.Transcripts = []
        self.Introns = {} # (start, stop) = observedCount
        self.IntronsByTranscript = {} #transcript -> dictionary of (start, stop) ->dummy
        ## used to compare introns so that we can distinguish alternative transcripts
    def PrintMe(self):
        print "Locus %s on %s"%(self.Name, self.Chromosome)
        #print "Intron set:%s"%self.Introns
        for Transcript in self.Transcripts:
            Transcript.PrintMe()

class GFFLClass:
    "the gbrowse GFF-Like format.  NOT real gff"
    def __init__(self):
        "constructor"
        self.Columns = GFFLColumns
        self.AllObservedSequences = []

    def ParseGFFL(self, FileName):
        """At this moment, we can't really get the phase of the spliced exons
        so we postpone setting the sequences and phases arrays until later.
        """
        Handle = open(FileName, "rb")
        Chromosome = None
        self.AllObservedSequences = [] #clear it out evert tune we pasrse a new file
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            Line = Line.strip()
            if Line[:9] == "reference":
                Chromosome = Line[10:]
                continue 
            Bits = list(Line.split("\t"))
            if len(Bits) < 3:
                #other meta data
                continue
            Sequence = ObservedSequenceClass()
            Sequence.Chromosome = Chromosome
            Sequence.Name = Bits[self.Columns.Name]
            IndexBits = Bits[self.Columns.Sequences].split(",")
            FrameString = Bits[self.Columns.Comment]
            if not self.IsOnReverseStrand(IndexBits):
                continue
            self.FigureOutGFFLSequences(Sequence, IndexBits, FrameString)
            self.AllObservedSequences.append(Sequence)
            if 0:# Sequence.Name == "Prot37823":
                print "ParseGFFL, found sequence of interest"
                Sequence.PrintMe()
        Handle.close()
        return self.AllObservedSequences

    def IsOnReverseStrand(self,IndexBits):
        "determines whether this is on the reverse strand or not"
        (Bit1First,Bit1Last) = IndexBits[0].split("..")
        Reverse = 0 #whether I should reverse the arrays
        if Bit1First < Bit1Last:
            ## forward strand, numerical ordering
            Reverse = 0
        elif Bit1First == Bit1Last:
            ## rare case of only a single nuc included in this exon
            (Bit2First,Bit2Last) = IndexBits[1].split("..")
            if Bit2First < Bit2Last:
                ## vanishingly small possibility that both the first and second exons have only a single base
                ## so I don't have to have an option here for ==
                Reverse =0
            else:
                Reverse = 1
        else:
            Reverse = 1
        return Reverse

    def FigureOutGFFLSequences(self, Sequence, IndexBits, FrameString):
        """This can get ugly.  Read about ObservedSequenceClass objects before you touch this code, capiche?
        IndexBits originate from the sequence specified in GFFL
        forward strand (normal easy):  123..145,234..344
        reverse strand (messy):  344..234,145..123
        We want to create the Sequence.Sequences = [(123,145),(234,344)]  the numerical ordering of sequence.
        """
        ## 1. determine whether stuff is in numerical order or not (forward or reverse strand)
        Reverse = self.IsOnReverseStrand(IndexBits)
        if Reverse:
            IndexBits.reverse()
        ## now do the actual work.  stupid double helix wasting my time
        for Bit in IndexBits:
            (Start, Stop) = Bit.split("..")
            Start = int(Start)
            Stop = int(Stop)
            if Reverse:
                Swap = Stop
                Stop = Start
                Start = Swap
            Sequence.Sequences.append((Start, Stop))
            if not Sequence.Start:
                Sequence.Start = Start
        #after running through it all
        Sequence.Stop = Stop #trusting from the last iteration.
        if FrameString.find("Frame") == -1:
            #NO Frame Information Given
            Sequence.Frame.append(None)
        else:
            FrameString = FrameString.replace("#Frame:", "")
            FrameBits = FrameString.split(",")
            if Reverse:
                FrameBits.reverse()
            for Bit in FrameBits:
                if Bit == '':
                    continue
                Sequence.Frame.append(int(Bit))
class GFFLColumns:
    Type = 0
    Name = 1
    Sequences = 2
    Comment= 3


class GFFColumns:
    Chr = 0
    Type = 2
    Start = 3
    Stop = 4
    Phase = 7
    Strand = 6
    Group = 8


class GFFClass:
    """the real GFF format
        Possible Types include: gene, pseudogene, ncRNA, exon, CDS, tRNA, three_prime_UTR, five_prime_UTR
        cDNA_Match
    """

    def __init__(self):
        self.Columns = GFFColumns
        self.AllExons = []
        self.Transcripts = {} #ID -> TranscriptClass object
        self.LociBoundaries = {} #ID -> (start, stop) 'int'
        self.AllExonCount = 0

    def ParseGeneralFeatureGFF(self, FileName, LocusTypeList):
        """this parses the GFF file looking for general locus features, and then
        determines what to do with those features based on some parameters
        """
        print "Parsing GFF File %s"%FileName
        Handle = open(FileName, "rb")
        ## set up an array for returning data of each type.
        ## array will be tuples (name, chr, start, stop), which the
        ## calling function will do something with
        ReturnLocations = {}
        for LocusType in LocusTypeList:
            #print "Locus Type ", LocusType
            ReturnLocations[LocusType] = []
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            Bits = Line.strip().split("\t")
            LocusType = Bits[self.Columns.Type]
            if not LocusType in LocusTypeList:
                continue # we don't care about it
            Start = int(Bits[self.Columns.Start])
            Stop = int(Bits[self.Columns.Stop])
            if Stop < Start:
                print "SNAFU !!: ParseGeneralFeatureGFF.  stop < stop"
            Chr = Bits[self.Columns.Chr]
            GeneID = self.ParseIDFromGFF(Bits)
            #print "Appending %s, %s,%s, %s"%(GeneID, Chr, Start, Stop)
            ReturnLocations[LocusType].append((GeneID, Chr, Start, Stop))
        return ReturnLocations

    def ParseGFF(self, FileName):
        """parse out the GFF into exons, and coordinates
            1. start looking for the mRNA tag, make a gene out of it's name/ID
            2. put all exons into the gene to which they belong (parent=ID)
        """
        print "Parsing GFF File %s"%FileName
        Handle = open(FileName, "rb")
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            Bits = Line.strip().split("\t")
            if Bits[self.Columns.Type] == "mRNA":
                GeneID = self.ParseIDFromGFF(Bits)
                if not GeneID:
                    #SNAFU parsing
                    continue
                self.Transcripts[GeneID] = TranscriptClass()
                self.Transcripts[GeneID].Name = GeneID
                self.Transcripts[GeneID].ParentLocus = self.ParseParentFromGFF(Bits)
                self.Transcripts[GeneID].Chromosome = Bits[self.Columns.Chr]
                self.Transcripts[GeneID].Strand = Bits[self.Columns.Strand]
            if Bits[self.Columns.Type] == "CDS":
                self.ParseCDSFromGFF(Bits)
            if Bits[self.Columns.Type] == "gene":
                self.ParseGeneFromGFF(Bits)
        # Go through and assign parentage
        TempTranscripts = {}
        while(len(self.AllExons) > 0):
            Exon = self.AllExons.pop()
            if not TempTranscripts.has_key(Exon.ParentTranscript):
                TempTranscripts[Exon.ParentTranscript] = []
            TempTranscripts[Exon.ParentTranscript].append(Exon)
        # now we must sort the exons
        for TranscriptName in TempTranscripts.keys():
            List = TempTranscripts[TranscriptName]
            self.Transcripts[TranscriptName].Exons = self.SortExonList(List, self.Transcripts[TranscriptName].Strand)
        return self.Transcripts

    def ParsecDNAGFF(self, FileName):
        """parse out the GFF into exons, and coordinates
            1. start looking for the mRNA tag, make a gene out of it's name/ID
            2. put all exons into the gene to which they belong (parent=ID)
        """
        print "Parsing GFF File %s"%FileName
        Handle = open(FileName, "rb")
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            Bits = Line.strip().split("\t")
            if Bits[self.Columns.Type] in ["cDNA_match", "EST_match"]:
                GeneID = self.ParseIDFromGFF(Bits)
                if not GeneID:
                    #SNAFU parsing
                    continue
                self.Transcripts[GeneID] = TranscriptClass()
                self.Transcripts[GeneID].Name = GeneID
                self.Transcripts[GeneID].ParentLocus = self.ParseParentFromGFF(Bits)
                self.Transcripts[GeneID].Chromosome = Bits[self.Columns.Chr]
                self.Transcripts[GeneID].Strand = Bits[self.Columns.Strand]
                self.ParsecDNAExonFromGFF(Bits)
        # Go through and assign parentage
        TempTranscripts = {}
        while(len(self.AllExons) > 0):
            Exon = self.AllExons.pop()
            if not TempTranscripts.has_key(Exon.ParentTranscript):
                TempTranscripts[Exon.ParentTranscript] = []
            TempTranscripts[Exon.ParentTranscript].append(Exon)
        # now we must sort the exons
        for TranscriptName in TempTranscripts.keys():
            List = TempTranscripts[TranscriptName]
            self.Transcripts[TranscriptName].Exons = self.SortExonList(List, self.Transcripts[TranscriptName].Strand)
        return self.Transcripts


    def ParseGeneFromGFF(self, Bits):
        ID = self.ParseIDFromGFF(Bits)
        Start = int(Bits[self.Columns.Start])
        Stop = int(Bits[self.Columns.Stop])
        self.LociBoundaries[ID] = (Start, Stop)

    def ParseIDFromGFF(self, Bits):
        GroupInfo = Bits[self.Columns.Group]
        Bytes = GroupInfo.split(";")
        ID = None
        for Byte in Bytes:
            if Byte[:2] == "ID":
                ID = Byte[3:]
        return ID

    def ParseParentFromGFF(self, Bits):
        GroupInfo = Bits[self.Columns.Group]
        Bytes = GroupInfo.split(";")
        for Byte in Bytes:
            if Byte[:6] == "Parent":
                #print Byte[7:], Byte
                return Byte[7:]
        return None


    def ParseCDSFromGFF(self, Bits):
        "make an exon object and put it into the GFFExons.  Here Exon == CODING exon"
        self.AllExonCount += 1
        Exon = ExonClass()
        Exon.Start = int(Bits[self.Columns.Start])
        Exon.Stop = int(Bits[self.Columns.Stop])
        Exon.TotalDNABases = Exon.Stop - Exon.Start
        Exon.Phase = int(Bits[self.Columns.Phase])
        Exon.Strand = Bits[self.Columns.Strand]
        Exon.ParentTranscript = self.ParseParentFromGFF(Bits)
        if Exon.ParentTranscript:
            self.AllExons.append(Exon)

    def ParsecDNAExonFromGFF(self, Bits):
        "make an exon object and put it into the GFFExons.  Here Exon == CODING exon"
        self.AllExonCount += 1
        Exon = ExonClass()
        Exon.Start = int(Bits[self.Columns.Start])
        Exon.Stop = int(Bits[self.Columns.Stop])
        Exon.TotalDNABases = Exon.Stop - Exon.Start
        #Exon.Phase = int(Bits[self.Columns.Phase])
        Exon.Strand = Bits[self.Columns.Strand]
        Exon.ParentTranscript = self.ParseIDFromGFF(Bits)
        if Exon.ParentTranscript:
            self.AllExons.append(Exon)

    def SortExonList(self, List, Strand):
        "assumption that there are no overlapping exons in this gene.  it's a logical gene"
        ReturnList = []
        Hash = {}
        for Exon in List:
            Hash[Exon.Start] = Exon
        Keys = Hash.keys()
        Keys.sort()
        if Strand == "-":  # this is to make the first exon in the list, the first coding exon, not first lexographically
            Keys.reverse()
        for Key in Keys:
            ReturnList.append(Hash[Key])
        return ReturnList

    def GetGenomicLocation(self, GFFKey, Residue, Offset):
        "look for the Ith residue in the exon structure of the GFF gene, possible offset for END"
        ExonList = self.Transcripts[GFFKey].Exons
        DNAIndex = Residue * 3
        #print "GetGenomicLocation: %s, %s, %s"%(GFFKey, Residue, Offset)
        if Offset == "END":
            DNAIndex += 2
        #print "DNA Bases needed for %s, %s is %s"%(Residue, Offset, DNAIndex)
        # now find the exon with that index
        DNAExonStart = -1
        DNAExonStop = -1
        for Exon in ExonList:
            DNAExonStart = DNAExonStop + 1
            DNAExonStop = DNAExonStart + Exon.TotalDNABases
            #print "This exon covers DNA bases %d to %d"%(DNAExonStart, DNAExonStop)
            if DNAExonStop >= DNAIndex:
                #Location Found in this exon.
                BasesIntoThisExon = DNAIndex - DNAExonStart
                GenomicLocation = BasesIntoThisExon + Exon.Start
                #print "Getting the genomic location for %d, %s is %d"%(Residue, Offset, GenomicLocation)
                return (GenomicLocation, Exon.Start)
        print "ERROR: No genomic location (GFF) found %s, %s, %s"%(GFFKey, Residue, Offset)
        return (None, None)

    def GetGenomicLocationMinusStrand(self, GFFKey, Residue, Offset):
        """look for the Ith residue in the exon structure of the GFF gene, possible offset for END
        This gets a bit tricky for the minus stranded stuff, so we blow it out into a separate function.
        """
        ExonList = self.Transcripts[GFFKey].Exons
        DNAIndex = Residue * 3
        if Offset == "END":
            DNAIndex += 2
        #print "DNA Bases needed for %s, %s is %s"%(Residue, Offset, DNAIndex)
        # now find the exon with that index
        DNAExonStart = -1
        DNAExonStop = -1
        for Exon in ExonList:
            DNAExonStart = DNAExonStop + 1
            DNAExonStop = DNAExonStart + Exon.TotalDNABases
            #print "This exon covers DNA bases %d to %d"%(DNAExonStart, DNAExonStop)
            if DNAExonStop >= DNAIndex:
                #Location Found in this exon.
                BasesIntoThisExon = DNAIndex - DNAExonStart
                GenomicLocation =  Exon.Stop - BasesIntoThisExon 
                #print "Getting the genomic location for %d, %s is %d"%(Residue, Offset, GenomicLocation)
                return (GenomicLocation, Exon.Start)
        print "ERROR: NO genomic location found"
        print GFFKey, Residue, Offset
        return (None, None)
    