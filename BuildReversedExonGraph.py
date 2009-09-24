"""
Write out the REVERSED version of a exon graph database.

We build this database to estimate our false positive rate.  ALMOST all
matches against the reversed database are invalid.  Therefore, so the
distribution of top MQScores provides an empirical p-value distribution
for the search against the real database (which has the same size, same
number of nodes of various degree, etc).
"""
import struct
import os
import sys
import string
import traceback

VERBOSE_FLAG = 0

def Main(InputFileName, OutputFileName):
    InputFile = open(InputFileName, "rb")
    OutputFile = open(OutputFileName, "wb")
    while (1):
        print "File position:", InputFile.tell()
        Name = InputFile.read(256)
        if not Name:
            break
        ReverseOneGene(InputFile, OutputFile, Name)

class ExonClass:
    def __init__(self):
        self.BackwardLinks = []
        self.ForwardLinks = []
        
def ReverseOneGene(InputFile, OutputFile, Name):
    OutputFile.write(Name)
    InputFile.read(256) # junk
    OutputFile.write(Name)
    NullPos = Name.find("\0")
    if NullPos:
        Name = Name[:NullPos]
    print Name
    # Chromosome number:
    Data = InputFile.read(4)
    OutputFile.write(Data)
    if VERBOSE_FLAG:
        print "Chromosome:", struct.unpack("<i", Data)[0]
    # Forward flag:
    Data = InputFile.read(1)
    if Data == chr(0):
        OutputFile.write(chr(1))
    else:
        OutputFile.write(chr(0))
    # Exon count:
    Data = InputFile.read(4)
    ExonCount = struct.unpack("<i", Data)[0]
    if VERBOSE_FLAG:
        print "Exon count:", ExonCount
    OutputFile.write(Data)
    # Iterate over exons:
    ExonList = []
    for ExonIndex in range(ExonCount):
        Exon = ExonClass()
        Exon.Start = struct.unpack("<i", InputFile.read(4))[0]
        Exon.End = struct.unpack("<i", InputFile.read(4))[0]
        Exon.Length = struct.unpack("<i", InputFile.read(4))[0]
        Exon.Occurrences = struct.unpack("<i", InputFile.read(4))[0]
        if VERBOSE_FLAG:
            print "Exon %s: %s-%s (%s)"%(ExonIndex, Exon.Start, Exon.End, Exon.Length)
        if Exon.Length > 10000:
            print "ERROR: Exon length is TOO HUGE!", Exon.Length
            return
        if Exon.Length:
            Exon.Sequence = InputFile.read(Exon.Length)
            if VERBOSE_FLAG:
                print Exon.Sequence[:250]
            Exon.Sequence = list(Exon.Sequence)
            Exon.Sequence.reverse()
            Exon.Sequence = string.join(Exon.Sequence, "")
            if VERBOSE_FLAG:
                print Exon.Sequence[:250]
        else:
            Exon.Sequence = ""
        Exon.Index = len(ExonList)
        ExonList.append(Exon)
        Exon.Prefix = InputFile.read(2)
        Exon.Suffix = InputFile.read(2)
        if VERBOSE_FLAG:
            print "Prefix %s Suffix %s"%(Exon.Prefix, Exon.Suffix)
        Exon.BackwardLinkCount = struct.unpack("<i", InputFile.read(4))[0]
        Exon.ForwardLinkCount = struct.unpack("<i", InputFile.read(4))[0]
        if VERBOSE_FLAG:
            print "Back %s forward %s"%(Exon.BackwardLinkCount, Exon.ForwardLinkCount)
        for LinkIndex in range(Exon.BackwardLinkCount):
            BackExonIndex = struct.unpack("<i", InputFile.read(4))[0]
            BackLinkPower = struct.unpack("<i", InputFile.read(4))[0]
            BackLinkAA = InputFile.read(1)
            if VERBOSE_FLAG:
                print "Back link power %d via '%c' to exon %d"%(BackLinkPower, BackLinkAA, BackExonIndex)
            BackExon = ExonList[BackExonIndex]
            Exon.BackwardLinks.append((BackExon, BackLinkPower, BackLinkAA))
            BackExon.ForwardLinks.append((Exon, BackLinkPower, BackLinkAA))
    ####################################
    ExonList.reverse()
    for Index in range(len(ExonList)):
        ExonList[Index].Index = Index
    for Index in range(len(ExonList)):
        Exon = ExonList[Index]
        OutputFile.write(struct.pack("<i", Exon.Start))
        OutputFile.write(struct.pack("<i", Exon.End))
        OutputFile.write(struct.pack("<i", Exon.Length))
        OutputFile.write(struct.pack("<i", Exon.Occurrences))
        if Exon.Sequence:
            OutputFile.write(Exon.Sequence)
        OutputFile.write(Exon.Prefix)
        OutputFile.write(Exon.Suffix)
        # Forward and backward links have been swapped:
        OutputFile.write(struct.pack("<i", Exon.ForwardLinkCount))
        OutputFile.write(struct.pack("<i", Exon.BackwardLinkCount))
        
        for (LinkedExon, Power, AA) in Exon.ForwardLinks:
            OutputFile.write(struct.pack("<i", LinkedExon.Index))
            OutputFile.write(struct.pack("<i", Power))
            OutputFile.write(AA)

if __name__ == "__main__":
    import psyco
    psyco.full()
    #Main("ESTSpliceDB\\Genome.dat", "ESTSpliceDB\\GenomeReversed.dat")
    #Main("ESTSpliceDB\\1+.dat", "TestReversed.dat")
    #Main("TestReversed.dat", "TestReversedReversed.dat")
    Main(sys.argv[1], sys.argv[2])
