import struct

File = open("TrainingFiles\\SVMToPValueTable.txt", "r")
Values = []
for FileLine in File.xreadlines():
    Bits = FileLine.split("\t")
    try:
        SVMScore = float(Bits[0])
        CumulativeFraction = float(Bits[1]) / 100.0
    except:
        continue
    Values.append((SVMScore, CumulativeFraction))

# PValues.dat has the following format:
# Bin count (as an int)
# Minimum bin SVM score
# then, a list of p-values for each SVM score bin, starting with the lowest.
OutputFile = open("PValues.dat", "wb")
OutputFile.write(struct.pack("<i", len(Values)))
OutputFile.write(struct.pack("<f", Values[0][0]))
BinIndex = 0
for (SVMScore, CumulativeFraction) in Values:
    PValue = 1.0 - CumulativeFraction
    PValue = max(PValue, 0.00001)
    print BinIndex, PValue
    BinIndex += 1
    OutputFile.write(struct.pack("<f", PValue))