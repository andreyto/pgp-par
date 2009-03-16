"""
Split the exon-ID file by chromosome
"""
from Utils import *

ChromosomeFiles = [None]
for Index in range(1, 49):
    File = open("E:\\chromosome\\GeneIDOutput\\Exons%d.txt"%Index, "wb")
    ChromosomeFiles.append(File)
    
File = open(r"E:\Chromosome\GeneIDOutput\all.exons.geneid","rb")
for FileLine in File.xreadlines():
    Bits = FileLine.split("\t")
    ChromNumber = ChromosomeMap.get(Bits[0], 0)
    if ChromNumber:
        Score = float(Bits[5])
        if Score >= -1:
            ChromosomeFiles[ChromNumber].write(FileLine)
    else:
        print "Weird chromosome encountered:", Bits[0]