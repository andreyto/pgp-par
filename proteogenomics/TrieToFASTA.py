"""
Translate a .trie file into FASTA format.  Inverse of PrepDB.py
"""
import os
import sys
import struct

def SplitLongLines(Str):
    Start = 0
    LineLen = 78
    End = LineLen
    Output = ""
    while (Start < len(Str)):
        Output += Str[Start:End] + "\n"
        Start += LineLen
        End += LineLen
    return Output
        

def Translate(SourceFileName, IndexFileName, TargetFileName, MaximumRecordLength = None):
    SourceFile = open(SourceFileName, "rb")
    Text = SourceFile.read()
    SourceFile.close()
    IndexFile = open(IndexFileName, "rb")
    TargetFile = open(TargetFileName, "wb")
    Pos = 0
    LastSequence = 0
    RecordNumber = 0
    while (1):
        RecordEnd = Text.find("*", Pos)
        if RecordEnd == -1:
            Aminos = Text[Pos:]
            LastSequence = 1
        else:
            Aminos = Text[Pos:RecordEnd]
        IndexFile.seek(92*RecordNumber + 12)
        Name = IndexFile.read(80)
        NullPos = Name.find('\0')
        if NullPos>-1:
            Name = Name[:NullPos]
        #Name = FixName(struct.unpack("<40s", Str))
        if MaximumRecordLength != None:
            ChunkCount = len(Aminos) / MaximumRecordLength
            if len(Aminos) % MaximumRecordLength:
                ChunkCount += 1
            for ChunkIndex in range(ChunkCount):
                ChunkStart = MaximumRecordLength * ChunkIndex
                ChunkEnd = MaximumRecordLength * (ChunkIndex + 1)
                #print "Chunk %s: %s...%s"%(ChunkIndex, ChunkStart, ChunkEnd)
                if ChunkCount > 1:
                    ChunkName = "%s.%s"%(Name, ChunkIndex)
                else:
                    ChunkName = Name
                TargetFile.write(">%s\n"%ChunkName)
                AminoString = SplitLongLines(Aminos[ChunkStart:ChunkEnd])
                TargetFile.write(AminoString)
        else:
            TargetFile.write(">%s\n"%Name)
            AminoString = SplitLongLines(Aminos)
            TargetFile.write(AminoString)
        Pos = RecordEnd + 1
        RecordNumber += 1
        if (LastSequence):
            break

TriePath = sys.argv[1]
(Stub, Extension) = os.path.splitext(TriePath)
IndexPath = "%s.index"%Stub
FASTAPath = "%s.fasta"%Stub
Translate(TriePath, IndexPath, FASTAPath)

