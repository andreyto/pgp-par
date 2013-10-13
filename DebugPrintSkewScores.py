"""
DebugPrintSkewScores: Ensure that tag edge skew scores are computed properly
"""
import os
import sys
import struct

File = open(sys.argv[1], "rb")
BinCount = struct.unpack("<i", File.read(struct.calcsize("i")))[0]
print "Bin count:", BinCount
FloatSize = struct.calcsize("<f")
for BlobIndex in range(2):
    for BinIndex in range(BinCount):
        ScoreString = File.read(FloatSize)
        Score = struct.unpack("<f", ScoreString)
        print "  %d: %s"%(BinIndex, Score)
        