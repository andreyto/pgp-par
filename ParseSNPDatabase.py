"""
Given the file snp.txt from the ucsc genome browser, convert it to a more convenient
binary format.  
"""
import os
import sys
import traceback
import struct
from Utils import *

class PolymorphismTypes:
    SNP = 0

"""
SNP Map (combined table of SNPs, hg17 and later)

table snp
"Polymorphism data from dbSnp XML files or genotyping arrays"
    (
    string  chrom;      "Chromosome"
    uint    chromStart; "Start position in chrom"
    uint    chromEnd;   "End position in chrom"
    string  name;       "Reference SNP identifier or Affy SNP name"
    uint    score;      "Not used"
    char[1] strand;     "Which DNA strand contains the observed alleles"
    string  observed;   "The sequences of the observed alleles"
    string  molType;    "Sample type from exemplar ss"
    string  class;      "The class of variant"
    string  valid;      "The validation status of the SNP"
    float   avHet;      "The average heterozygosity from all observations"
    float   avHetSE;    "The Standard Error for the average heterozygosity"
    string  func;       "The functional category of the SNP"
    string  locType;    "How the variant affects the reference sequence"
    string  source;     "Source of the data - dbSnp, Affymetrix, ..."
    string  exception;  "List of exceptionIds for 'invariant' conditions"
    )
"""

RCDict = {"A":"T", "T":"A", "G":"C", "C":"G"}
def GetRC(String):
    "return the reverse complement of a nucleotide sequence"
    Len = len(String)
    RCStr = ""
    for X in range(Len-1, -1, -1):
        RCStr += RCDict.get(String[X], String[X])
    return RCStr

def ParseDatabase(SNPFileName):
    SNPBinaryFileDir = "SNP"
    SNPLists = []
    for ChromosomeNumber in range(0, 49):
        SNPLists.append([])
    ################################################
    # Read the SNP records, distill the binary info we care about, and write out
    # to the right chromosome file:
    File = open(SNPFileName, "r")
    LineNumber = 0
    for FileLine in File.xreadlines():
        LineNumber += 1
        if LineNumber % 1000 == 0:
            print LineNumber,
        Bits = FileLine.split("\t")
        if len(Bits)<11:
            continue
        # There are several types of polymorphism reported in snp.txt
        # Most common are:
        #  snp (two forms, each 1base)
        #  microsat (short repeated sequence, variable count)
        #  mixed (three forms or more, possibly of different lengths)
        #  in-del (one empty form, one of length 1 or more)
        #  named (large deletions and such)
        # Handling SNPs is easy.  Handling in-dels is harder: The sequence may or may not
        # be contained in the reference genome.  Determining whether it's present probably
        # requires us to compare the surrounding sequence, since - especially if it's a
        # single nucleotide - the short form might still have the right sequence in the
        # right spot.  Also, in-dels do make the whole picture much trickier when we look at
        # our interval graph...  
        SNPType = Bits[9]
        # For now, we ONLY consider SNPs. That way our exons' genomic coordinates will be
        # accurate, and we don't have to do hairy stuff to handle deletions that
        # span/touch interval boundaries.
        if SNPType != "snp":
            continue
        # We may want to skip SNPs that have exceptions listed.
        if len(Bits)>16:
            Exceptions = Bits[16].strip().split(",")
            if "21" in Exceptions:
                # Neither allele matches the genomic sequence!  We could possibly try
                # all three possibilities (two from the snp, one from the genome), but
                # that seems a bit silly.
                continue
            if len(Exceptions[0]):
                continue
        Type = PolymorphismTypes.SNP # snp
        ChromosomeName = Bits[1]
        ChromosomeNumber = ChromosomeMap.get(ChromosomeName, None)
        if not ChromosomeNumber:
            print "Weird chromosome number:", ChromosomeName
            continue # comes from chr1_random or other jargle
        Start = int(Bits[2])
        End = int(Bits[3])
        if End > Start + 1:
            # We're ONLY handling SNPs, so the end should always be equal to start + 1.  Sometimes
            # it's not (and that's reported as an exception); let's skip those.
            continue
        Forms = Bits[7].split("/")
        if Bits[6] == "+":
            Forms = list(Forms)
        elif Bits[6] == "-":
            FixedForms = []
            for Form in Forms:
                FixedForms.append(GetRC(Form))
            Forms = FixedForms
        else:
            print "** Error: Unknown strand!"
            print FileLine
        for X in range(len(Forms)):
            Forms[X] = Forms[X].replace("-","")
        Str = "%s%s"%(Forms[0], Forms[1])
        SNPLists[ChromosomeNumber].append((Start, Str))
##        SNPFile = ChromosomeFiles[ChromosomeNumber]
##        # Record format:
##        # Position (int), type (char), <forms>
##        # for type 0 (snp), <forms> is 2 characters.
##        SNPFile.write(struct.pack("<i", Start))
##        SNPFile.write(chr(Type))
##        if Type == PolymorphismTypes.SNP:
##            SNPFile.write(Forms[0])
##            SNPFile.write(Forms[1])
    ################################################
    # Open 24 files, one for each chromosome:
    for ChromosomeNumber in range(1, 49):
        Path = os.path.join(SNPBinaryFileDir, "%s.snp"%ChromosomeNumber)
        File = open(Path, "wb")
        List = SNPLists[ChromosomeNumber]
        List.sort()
        ListLen = len(List)
        OldPos = None
        OldForms = ""
        for (Pos, Forms) in List:
            if Pos == OldPos:
                # Two SNPs at the same locus.  Merge them.
                for Form in Forms:
                    if Form not in OldForms:
                        OldForms += Form
            else:
                if OldForms:
                    File.write(struct.pack("<i", OldPos))
                    File.write(chr(len(OldForms) - 2))
                    if len(OldForms)>2:
                        print "** TRIPLE!", ChromosomeNumber, OldPos, OldForms
                    File.write(OldForms)
                OldPos = Pos
                OldForms = Forms
        if OldForms:
            File.write(struct.pack("<i", OldPos))
            File.write(chr(len(OldForms) - 2))
            File.write(OldForms)
                
        File.close()

def DebugPrintDatabase(DBPath, FirstRecord, LastRecord):
    print DBPath
    File = open(DBPath, "rb")
    RecordNumber = 0
    while (1):
        Data = File.read(4)
        if not Data:
            break # eof
        GenomeStart = struct.unpack("<i", Data)[0]
        PolyType = ord(File.read(1))
        FormCount = PolyType + 2
        Forms = File.read(FormCount)
        print "SNP record %d: Pos %d forms %s"%(RecordNumber, GenomeStart, Forms)
        RecordNumber += 1
        if LastRecord!=None and RecordNumber > LastRecord:
            break
        
if __name__ == "__main__":
    #import psyco
    #psyco.full()
    if len(sys.argv)<2:
        Command = "parse"
    else:
        Command = sys.argv[1].lower()
    if Command == "parse":
        #ParseDatabase("SNPHead.txt")
        ParseDatabase("e:\\chromosome\\snp.txt")
    elif Command == "print":
        #DebugPrintDatabase("SNP\\21.snp", 0, None)
        #DebugPrintDatabase("SNP\\6.snp", 0, None)
        DebugPrintDatabase("SNP\\9.snp", 0, None)
    else:
        print "What does '%s' mean?  I know 'parse' and 'print'."%(Command)