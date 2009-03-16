"""
Convert a tab-delimited file to an SVM training file
"""
import sys
import string
import os

# Arguments:
#  Tab-delimited file [min column] [max column]
File = open(sys.argv[1], "r")

for FileLine in File.xreadlines():
    Bits = list(FileLine.strip().split("\t"))
    if len(Bits)<2:
        continue
    if int(Bits[0]):
        Bits[0] = "+1"
    else:
        Bits[0] = "-1"
    for Index in range(1, len(Bits)):
        Bits[Index] = "%d:%s"%(Index, Bits[Index])
    print string.join(Bits, " ")
    
        