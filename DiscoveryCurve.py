"""
"""


def SplitFile(InputFileName, ColumnA, ColumnB, OutputFileName):
    N = 38
    InputFile = open(InputFileName, "rb")
    OutputFile = open(OutputFileName, "wb")
    LineNumber = 0
    for FileLine in InputFile.xreadlines():
        Bits = FileLine.strip().split("\t")
        if len(Bits) < max(ColumnA, ColumnB) + 1:
            OutputFile.write(FileLine) # include "parsing..." messages
            continue
        LineNumber += 1
        if LineNumber % N == 0:
            OutputFile.write("%s\t%s\t\n"%(Bits[ColumnA], Bits[ColumnB]))
    InputFile.close()
    OutputFile.close()
        


#SplitFile("DiscoveryCurve.txt", 1, 2, "DCAllPeptides.txt")
#SplitFile("KPDiscoveryCurve.txt", 1, 2, "DCKnownPeptides.txt")
#SplitFile("RDiscoveryCurve.txt", 2, 3, "DCRandomKP.txt")
#SplitFile("RDiscoveryCurve.txt", 0, 1, "RAllPeptides.txt")
SplitFile("RDiscoveryCurve.txt", 2, 3, "RKnownPeptides.txt")
#SplitFile("SDisco0005.txt", 1, 2, "ShewDisco.txt")
