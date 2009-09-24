"""
Generate equal number of false and true elements:
"""
import sys
import random
import os

random.seed(1) # for reproducibility

File = open(sys.argv[1], "r")
ProportionalFlag = 0
if len(sys.argv)>2:
    MaxRecordCount = int(sys.argv[2])
    if sys.argv[2][0] == "+":
        ProportionalFlag = 1
else:
    MaxRecordCount = -1

if len(sys.argv) > 3:
    OutFileName = sys.argv[3]
else:
    OutFileName = "TrainingSet.txt"
OutFile = open(OutFileName, "w")
Positives = []
Negatives = []
for FileLine in File.xreadlines():
    if FileLine[0] in ("+1"):
        Positives.append(FileLine)
    else:
        Negatives.append(FileLine)
random.shuffle(Positives)
random.shuffle(Negatives)
if ProportionalFlag:
    PositiveCount = min(MaxRecordCount, len(Positives))
    NegativeCount = int(round(PositiveCount * len(Negatives) / float(len(Positives))))
else:
    if MaxRecordCount > 0:
        MaxRecordCount = min(MaxRecordCount, len(Positives), len(Negatives))
    else:
        MaxRecordCount = min(len(Positives), len(Negatives))
    PositiveCount = MaxRecordCount
    NegativeCount = MaxRecordCount
        
print "%s positives, %s negatives...report %s + %s -."%(len(Positives), len(Negatives), PositiveCount, NegativeCount)
Positives = Positives[:MaxRecordCount]
Negatives = Negatives[:MaxRecordCount]
for Line in Positives:
    OutFile.write(Line)
for Line in Negatives:
    OutFile.write(Line)    