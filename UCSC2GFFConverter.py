"""
Convert the KnownGene annotation track (from the UCSC genome browser) to a GFF file.
"""
import os
import sys
import traceback

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

class ConverterClass:
    def __init__(self):
        pass
    def ConvertFile(self, InputFileName, Dir):
        """
        Parse KnownGenes.txt.  Write out a GFF-format file containing one line
        for each exon in each gene.
        """
        InputFile = open(InputFileName, "rb")
        OutputFiles = {}
        LineNumber = 0
        for FileLine in InputFile.xreadlines():
            LineNumber += 1
            if LineNumber % 100 == 0:
                print "Line %s..."%LineNumber
            if FileLine[0] == "#":
                continue
            Bits = list(FileLine.split("\t"))
            if len(Bits) < KGBits.UniqueID:
                continue
            for Index in range(len(Bits)):
                if Bits[Index][0] == '"':
                    Bits[Index] = Bits[Index][1:-1]
            Chromosome = Bits[KGBits.Chromosome]
            OutputFile = OutputFiles.get(Chromosome, None)
            if not OutputFile:
                OutputPath = os.path.join(Dir, "%s.gff"%Chromosome)
                OutputFile = open(OutputPath, "wb")
                OutputFiles[Chromosome] = OutputFile
            SequenceName = Bits[KGBits.UniqueID].strip()
            StartStr = "%s\t%s\texon\t"%(SequenceName, Chromosome)
            if Bits[KGBits.Strand] == "+":
                ForwardFlag = 1
            elif Bits[KGBits.Strand] == "-":
                ForwardFlag = 0
            else:
                print "* Error: Line %s has illegal forward flag '%s'"%(LineNumber, Bits[KGBits.Strand])
                continue
            ExtraAA = 0
            ExonCount = int(Bits[KGBits.ExonCount])
            ExonStarts = Bits[KGBits.ExonStarts].split(",")
            ExonEnds = Bits[KGBits.ExonEnds].split(",")
            CodingStart = int(Bits[KGBits.CodingStart])
            CodingEnd = int(Bits[KGBits.CodingEnd])
            if ForwardFlag:
                ExonIndexes = range(ExonCount)
            else:
                ExonIndexes = range(ExonCount - 1, -1, -1)
            for ExonIndex in ExonIndexes:
                Start = int(ExonStarts[ExonIndex])
                End = int(ExonEnds[ExonIndex])
                if Start > CodingEnd or End < CodingStart:
                    continue
                Start = max(Start, CodingStart)
                End = min(End, CodingEnd)
                # convert from knowngene (0-based numbering) to GFF (1-based numbering).
                # Later, MS2DB will convert us back to 0-based numbering.
                Start += 1  
                Str = StartStr + "%s\t%s\t1.0\t%s\t"%(Start, End, Bits[KGBits.Strand])
                # Frame:
                if ExtraAA == 0:
                    Frame = 0
                elif ExtraAA == 1:
                    Frame = 2
                else:
                    Frame = 1
                Str += "%s\t"%Frame
                Length = ExtraAA + (End - Start + 1) # add 1 to correct length
                ExtraAA = Length % 3
                Str += """ID "%s" SprotID "%s"\t"""%(Bits[KGBits.Name], Bits[KGBits.SprotID])
                OutputFile.write(Str + "\n")
if __name__ == "__main__":
    Converter = ConverterClass()
    #Converter.ConvertFile("e:\\chromosome\\knowngene.txt", "gff\\KnownGene")
    Converter.ConvertFile("f:\\HGenome\\knowngene.txt", "gff\\KnownGene")
    