"""ProkaryoticGFF3.py
Methods for reading or writing the GFF3 format

Notes for the GFF format are on: http://www.sequenceontology.org/gff3.shtml

##gff-version   3
##sequence-region   ctg123 1 1497228       
ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN
ctg123 . TF_binding_site 1000  1012  .  +  .  ID=tfbs00001;Parent=gene00001
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123 . mRNA            1300  9000  .  +  .  ID=mRNA00003;Parent=gene00001;Name=EDEN.3

ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mRNA00003
ctg123 . exon            1050  1500  .  +  .  ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123 . exon            3000  3902  .  +  .  ID=exon00003;Parent=mRNA00001,mRNA00003
ctg123 . exon            5000  5500  .  +  .  ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123 . exon            7000  9000  .  +  .  ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003

ctg123 . CDS             1201  1500  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123 . CDS             3000  3902  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123 . CDS             5000  5500  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123 . CDS             7000  7600  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1

ctg123 . CDS             1201  1500  .  +  0  ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123 . CDS             5000  5500  .  +  0  ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123 . CDS         7000  7600     .  +  0  ID=cds00002;Parent=mRNA00002;Name=edenprotein.2

Column 1: "seqid"

The ID of the landmark used to establish the coordinate system for the
current feature. e.g. Chromosome1, or Contig34

Column 2: "source"

The source is a free text qualifier intended to describe the algorithm
or operating procedure that generated this feature. here "Proteomics"

Column 3: "type"

The type of the feature is CDS for us, because these are coding sequences.

Columns 4 & 5: "start" and "end"

The start and end of the feature, in 1-based integer coordinates,
relative to the landmark given in column 1.  Start is always less than
or equal to end.

Column 6: "score"

The score of the feature, a floating point number.  We will use
the pvalue (best the localFDR) of the best spectrum for this peptide

Column 7: "strand"

The strand of the feature.  + for positive strand (relative to the
landmark), - for minus strand.

Column 8: "phase"

For features of type "CDS", the phase indicates where the feature
begins with reference to the reading frame.  The phase is one of the
integers 0, 1, or 2, indicating the number of bases that should be
removed from the beginning of this feature to reach the first base of
the next codon. In other words, a phase of "0" indicates that the next
codon begins at the first base of the region described by the current
line, a phase of "1" indicates that the next codon begins at the
second base of this region, and a phase of "2" indicates that the
codon begins at the third base of this region. This is NOT to be
confused with the frame, which is simply start modulo 3.

For forward strand features, phase is counted from the start
field. For reverse strand features, phase is counted from the end
field.

The phase is REQUIRED for all CDS features.

If we're dealing with unspliced peptides, the phase should always be zero.

Column 9: "attributes"

A list of feature attributes in the format tag=value.  Multiple
tag=value pairs are separated by semicolons.  

These tags have predefined meanings:

    ID       Indicates the name of the feature.  IDs must be unique
       within the scope of the GFF file. e.g. Peptide00001

    Name   Display name for the feature.  This is the name to be
           displayed to the user. Here use the peptide sequence
"""

class GFFColumns:  # see big comment above, these are zero indexed.
    Chr = 0
    Source = 1 #"Proteomics" or "Inspect" or "6frame"
    Type = 2
    Start = 3
    Stop = 4
    Score = 5
    Strand = 6
    Phase = 7
    Attributes = 8


class GFFClass:
    """The GFFClass is methods used for reading or writing GFF lines.
    """

    def __init__(self):
        self.Columns = GFFColumns
        self.NumColumns = 9
        self.Values = [] # an array of values to represent data for the various columns.  ALWAYS clear out before use!
        self.PeptideCount = 0
        
    def ResetPeptideCount(self):
        """If for some reason you want to reset the peptide ID counter to zero.
        """
        self.PeptideCount = 0

    def WriteProteomicEvidenceLine(self, Peptide, ProteinName, StartLocation, PValue):
        """Given the peptide and some other information, we return a line sufficient for gff3.
        For the moment we will assume that each peptide is its own entity, and not 
        give it parents or anything like that. Peptide is the amino acid sequence
        ProteinName comes from the fasta header and should have the format of "Protein0.Chr:Chr1.Frame2.StartNucleotide2.Strand+"
        WARNING: This does not do conversion! it takes data ready to write.
        """
        #self.Values must be filled with STRINGS
        self.Values = [0] * self.NumColumns #clean out before use
        self.Values[self.Columns.Source] = "Proteomics"
        self.Values[self.Columns.Type] = "CDS"
        self.Values[self.Columns.Phase] = "0" # all unspliced peptides are phase 0
        self.Values[self.Columns.Score] = "%s"%PValue
        self.PeptideCount += 1
        PeptideAttributeID = "Peptide%s"%self.PeptideCount
        self.Values[self.Columns.Attributes] = "ID=%s;Name=%s"%(PeptideAttributeID, Peptide)
        ## start parsing out things from the peptide
        ## format of Protein0.Chr:Chr1.Frame2.StartNucleotide2.Strand+
        ChrStringStart = ProteinName.find("Chr:") + 4
        ChrStringEnd = ProteinName.find(".Frame")
        self.Values[self.Columns.Chr] = ProteinName[ChrStringStart:ChrStringEnd]
        self.Values[self.Columns.Strand] = "+"
        if (ProteinName.find("Strand-") > -1):
            self.Values[self.Columns.Strand] = "-"
            (PeptideStartNucleotide, PeptideStopNucleotide) = self.GetCoordsRevStrand(Peptide, ProteinName, StartLocation)
        else:
            #now the position.  We first get the protein start nuc, and then offset to the peptide start
            StrandStartPos = ProteinName.find(".Strand")
            StartNucPos = ProteinName.find("StartNuc") + 8 # offset for the string "StartNucleotide"
            ProteinStartNucleotide = int(ProteinName[StartNucPos:StrandStartPos])
            #print "From %s I parsed out %s"%(ProteinName, ProteinStartNucleotide)
            # now offset from the protein residue to the start nucleotide
            PeptideStartOffset = StartLocation * 3 # because StartLocation is the start in amino acid space
            PeptideStartNucleotide = PeptideStartOffset + ProteinStartNucleotide
            PeptideStopNucleotide = PeptideStartNucleotide + (len(Peptide) * 3) - 1 # we do a minus one because the bases are inclusive
        self.Values[self.Columns.Start] = "%s"%PeptideStartNucleotide
        self.Values[self.Columns.Stop] = "%s"%PeptideStopNucleotide
        GFF3FormatedLine = "\t".join(self.Values)
        return GFF3FormatedLine

    def GetCoordsRevStrand(self, Peptide, ProteinName, StartLocation):
        """Separate method because doing math on the reverse strand is difficult.
        For the visually minded:  codons on the reverse listed below
        1   4   7   10  13  16  19
        AAA TTT CCC GGG AAA TTT CCC
        TTT AAA GGG CCC TTT AAA GGG AGT
         F   K   G   P   F   K   G   *
        Peptide GKFPG starts at position 21 (last C) and includes all bases up to 7, so we should do 
        start = 7, stop = 21
        Peptide KFPG starts at position 18.  So start protein (21) minus AA offset (3)
        """
        StrandStartPos = ProteinName.find(".Strand")
        StartNucPos = ProteinName.find("StartNuc") + 8 # offset for the string "StartNucleotide"
        ProteinStartNucleotide = int(ProteinName[StartNucPos:StrandStartPos])
        PeptideStartOffset = StartLocation *3
        PeptideStartNucleotide = ProteinStartNucleotide - PeptideStartOffset
        PeptideStopNucleotide = PeptideStartNucleotide - (len(Peptide)*3) + 1
        #now because this is on the reverse, and GFF3 requires start< stop, we send
        #them back in the reverse order
        return (PeptideStopNucleotide, PeptideStartNucleotide)
    
        
        
         
        
        
        
        
        
        