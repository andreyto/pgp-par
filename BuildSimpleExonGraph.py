"""
Build simple exon graphs, for testing exon-graph-to-sequence and
exon-graph-to-exon-graph alignment.
"""
import os
import sys
import struct

class ExonClass:
    def __init__(self):
        self.Sequence = ""
        self.BackExons = []
        self.BackAA = []
        self.ForwardEdgeCount = 0
    def WriteExon(self, OutFile):
        OutFile.write(struct.pack("<i", 1)) # start (on genome)
        OutFile.write(struct.pack("<i", 500)) # end (on genome)
        OutFile.write(struct.pack("<i", len(self.Sequence))) # seqlen
        OutFile.write(struct.pack("<i", 1)) # OccurrenceCount
        if self.Sequence:
            OutFile.write(self.Sequence)
        OutFile.write("XX") # prefix
        OutFile.write("XX") # suffix
        OutFile.write(struct.pack("<i", len(self.BackExons))) # back link count
        OutFile.write(struct.pack("<i", self.ForwardEdgeCount))
        for EdgeIndex in range(len(self.BackExons)):
            OutFile.write(struct.pack("<i", self.BackExons[EdgeIndex].Index))
            OutFile.write(struct.pack("<i", 5))
            if self.BackAA[EdgeIndex]:
                OutFile.write(self.BackAA[EdgeIndex])
            else:
                OutFile.write(chr(0))
class GeneClass:
    def __init__(self):
        self.Exons = []
    def WriteGene(self, OutFileName):
        OutFile = open(OutFileName, "wb")
        Name = "TestGene"
        OutFile.write(struct.pack("<256s", Name))
        OutFile.write(struct.pack("<256s", Name))
        OutFile.write(struct.pack("<i", 1)) # chromnumber
        OutFile.write(struct.pack("<i", len(self.Exons))) #ExonCount
        ExonIndex = 0
        for Exon in self.Exons:
            Exon.Index = ExonIndex
            Exon.WriteExon(OutFile)
            ExonIndex += 1
            
def BuildSimpleGene(Sequence, OutFileName):
    Gene = GeneClass()
    Pos = len(Sequence) / 2
    ExonA = ExonClass()
    ExonA.Sequence = Sequence[:Pos]
    ExonA.ForwardEdgeCount = 1
    Gene.Exons.append(ExonA)
    ExonB = ExonClass()
    ExonB.Sequence = Sequence[Pos:]
    Gene.Exons.append(ExonB)
    ExonB.BackAA.append(None)
    ExonB.BackExons.append(ExonA)
    Gene.WriteGene(OutFileName)
    print ExonA.Sequence
    print ExonB.Sequence

if __name__ == "__main__":
    BuildSimpleGene(sys.argv[1], sys.argv[2]) # sequence, output file name