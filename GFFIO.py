"""GFFIO.py
This auxiliary set of classes is to be used for input/output
related needs with the GFF3 format.

NOTE: this is a utility, and not executable from the command line

"""

from bioseq import FlatFileIO

class File(FlatFileIO):
    def __iter__(self):
        for line in self.io.xreadlines():
            if line[0] in ["#", "\n", ""]:
                continue
            yield Record(line)

    def write(self,record):
        attributes = '.'
        if len(record.attributes) > 0:
            # Attributes is a dictionary, which we join into key=value
            # pairs with each separated by a ;
            pairs = [(p[0],str(p[1])) for p in record.attributes.items()]
            attributes = ";".join(["=".join(x) for x in pairs]) 

        self.io.write("\t".join([
            record.seqid,
            record.source,
            record.type,
            str(record.start),
            str(record.end),
            str(record.score),
            record.strand,
            str(record.phase),
            attributes
            ]) + "\n")

class Record(object):
    def __init__(self, gffline="\t".join(list('.'*9)) ):
        """Parameters: String line of GFF
        Return: GFFIO.Record object filled with values
        Description: Take a GFF line and parse it out into the object. You
        can do whatever you want with the object.  The format should be invariant
        """    
        cols = gffline.strip().split("\t")
        self.seqid = cols[0]
        self.source= cols[1]
        self.type  = cols[2]
        self.start = cols[3] == '.' and '.' or int(cols[3])
        self.end   = cols[4] == '.' and '.' or int(cols[4])
        self.score = cols[5] == '.' and '.' or float(cols[5])
        self.strand= cols[6]
        self.phase = cols[7] == '.' and '.' or int(cols[7])
        self.attributes = {}
        for Pair in cols[8].split(";"):
            if '.' == Pair:
                continue
            (Key, Value) = Pair.split("=")
            self.attributes[Key] = Value


"""
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
