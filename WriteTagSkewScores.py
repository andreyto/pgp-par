"""
Parses TagSkewScores.txt (written by DriveScorpion) and writes out TagSkewScores.dat (binary
file, read by Inspect for penalizing tags)
"""
import os
import struct

File = open("TagSkewScores.txt", "rb")
ScoresByAbsTotalSkew = []
ScoresByTotalAbsSkew = []
for FileLine in File.xreadlines():
    Bits = FileLine.split("\t")
    try:
        ScoreA = float(Bits[0])
        ScoreB = float(Bits[1])
    except:
        continue
    ScoresByAbsTotalSkew.append(ScoreA)
    ScoresByTotalAbsSkew.append(ScoreB)
File.close()


File = open("TagSkewScores.dat", "wb")
ScoreCount = len(ScoresByAbsTotalSkew)
Str = struct.pack("<i", ScoreCount)
File.write(Str)
for Index in range(ScoreCount):
    Str = struct.pack("<f", ScoresByAbsTotalSkew[Index])
    File.write(Str)
for Index in range(ScoreCount):
    Str = struct.pack("<f", ScoresByTotalAbsSkew[Index])
    File.write(Str)
File.close()
    