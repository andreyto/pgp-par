"""
Given a spliced-protein DB, write out an index.
"""
import struct
import os
import sys
import traceback

def Main():
    File = open(sys.argv[1], "rb")
    IndexFileName = os.path.splitext(sys.argv[1])[0] + ".index"
    IndexFile = open(IndexFileName, "wb")
    OneIntSize = struct.calcsize("i")
    TwoIntSize = struct.calcsize("ii")
    LateBlockSize = struct.calcsize("2s2sii")
    BackLinkSize = struct.calcsize("iic")
    GIIDBlockSize = struct.calcsize("i")*10
    GeneCount = 0
    while (1):
        FilePos = File.tell()
        Name = File.read(256)
        if not Name:
            break
        Pos = Name.find(chr(0))
        if Pos!=-1:
            Name = Name[:Pos]
        Dummy = File.read(256)
        ChromosomeNumber = struct.unpack("<i", File.read(4))[0]
        ForwardFlag = File.read(1)
        ExonCount = struct.unpack("<i", File.read(4))[0]
        #IDBlock = File.read(GIIDBlockSize)
        GeneCount += 1
        #if GeneCount % 1000 == 0:
        #    print "#%d: Chromosome %s, exon count %s, name '%s'"%(GeneCount, Chrom, ExonCount, Name)
        # Scan past the exons:
        for ExonIndex in range(ExonCount):
            File.read(TwoIntSize) # Start and end
            Str = File.read(OneIntSize) # length
            SequenceLength = struct.unpack("<i", Str)[0]
            File.read(OneIntSize) # occurrences
            File.read(SequenceLength) # Sequence
            Str = File.read(LateBlockSize) # prefixchars, suffixchars, linkcounts
            (Prefix, Suffix, BackLinkCount, ForwardLinkCount) = struct.unpack("<2s2sii", Str)
            for LinkIndex in range(BackLinkCount):
                File.read(BackLinkSize) # linked exon index, link power, amino acid
        # Now we've finished scanning past the record...write to the index file:
        IndexStr = struct.pack("<qi80s", 0, FilePos, Name)
        IndexFile.write(IndexStr)
        #IndexFile.write(IDBlock)
        print Name, ChromosomeNumber, ord(ForwardFlag)
    print "%s\t%s\t"%(sys.argv[1], GeneCount)    
            

if __name__ == "__main__":
    import psyco
    psyco.full()
    Main()
