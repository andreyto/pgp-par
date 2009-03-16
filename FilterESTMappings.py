"""
The script ESTDatabaseSplit converts the raw EST-to-genome alignments to separate binary
files for each chromosome strand.  This script reads these binary files, and examines all
putative splice junctions.  Spurious splice junctions are removed.  This script then writes
out intervals and their splice junctions (in order) as a new file.

We accept a splice junction if:
- It's supported by 2 or more ESTs, or
- Its consensus splice signal score is good

The signal consists of 13 bases (9 from donor, 4 from acceptor); the most important and most
heavily weighted part is the "GU...AG" consensus signal.
A    |    x
C A G|G U x A G U . . . . . . . . C A G|G 
3 4 5 6 7 8 9 1011                15161718
"""
import os
import struct
import sys
import traceback
import math
import string

ESTSubDirectory = "ESTREF" # EST or ESTREF

RC = {"A":"T","T":"A","C":"G","G":"C","N":"N","X":"X"}

def GetRC(Str):
    NewStr = ""
    for X in range(len(Str) - 1, -1, -1):
        NewStr += RC[Str[X]]
    return NewStr

Profile = [
      {"A":0.3501  ,"C":0.3420  ,"G":0.1873  ,"T":0.1207},
      {"A":0.6085  ,"C":0.1094  ,"G":0.1220  ,"T":0.1601},
      {"A":0.0956  ,"C":0.0381  ,"G":0.7933  ,"T":0.0730},
      {"A":0.0028  ,"C":0.0010  ,"G":0.9949  ,"T":0.0014},
      {"A":0.0028  ,"C":0.0124  ,"G":0.0012  ,"T":0.9836},
      {"A":0.5512  ,"C":0.0304  ,"G":0.3846  ,"T":0.0338},
      {"A":0.7011  ,"C":0.0764  ,"G":0.1142  ,"T":0.1084},
      {"A":0.0787  ,"C":0.0545  ,"G":0.8026  ,"T":0.0642},
      {"A":0.1570  ,"C":0.1569  ,"G":0.1926  ,"T":0.4935},
      {"A":0.0560  ,"C":0.6400  ,"G":0.0094  ,"T":0.2947},
      {"A":0.9937  ,"C":0.0011  ,"G":0.0013  ,"T":0.0039},
      {"A":0.0011  ,"C":0.0014  ,"G":0.9943  ,"T":0.0032},
      {"A":0.2453  ,"C":0.1438  ,"G":0.4948  ,"T":0.1160},
    ]

class SpliceJunction:
    def __init__(self):
        self.Start = None
        self.End = None
        self.Sequence = ""
        self.Count = 0
        self.Score = 0 
    def GetScore(self):
        Score = 0
        for X in range(len(self.Sequence)):
            Score += math.log(Profile[X].get(self.Sequence[X], 0.25))
        self.Score = Score
        return Score
    
class FilterMeister:
    def __init__(self):
        # Key = end (higher-numbered residue) of junction
        # Value = list of junction objects
        self.JunctionEnds = {}
        self.JunctionStarts = {}
        # Key (start, end), value OccurrenceCount
        self.Intervals = {} 
    def Read(self, ChromosomeNumber, ReverseFlag, Task):
        if ReverseFlag:
            ReverseChar = "-"
        else:
            ReverseChar = "+"
        ChromosomePath = "e:\\Chromosome\\Chr%s.trie"%ChromosomeNumber
        try:
            ChromFile = open(ChromosomePath, "rb")
        except:
            print "Couldn't open chromosome path '%s' - exiting"%ChromosomePath
            return
        ESTFile = open("%s\\%s%s.sorted"%(ESTSubDirectory, ChromosomeNumber, ReverseChar), "rb")
        #############################################
        RecordCount = 0
        while 1:
            Intervals = []
            IntervalCount = ESTFile.read(4)
            if not IntervalCount:
                break
            RecordCount += 1
            #if (RecordCount%10000 == 0):
            #    print RecordCount
            IntervalCount = struct.unpack("<i", IntervalCount)[0]
            ESTFilePos = struct.unpack("<I", ESTFile.read(4))[0]
            for Index in range(IntervalCount):
                (Start, End) = struct.unpack("<ii", ESTFile.read(8))
                Intervals.append((Start, End))
                self.Intervals[(Start, End)] = self.Intervals.get((Start, End), 0) + 1
            # Count these intervals, and count these splicejunctions:
            for Index in range(len(Intervals) - 1):
                End = Intervals[Index+1][0]
                Start = Intervals[Index][1]
                FoundFlag = 0
                if not self.JunctionEnds.has_key(End):
                    self.JunctionEnds[End] = []
                if not self.JunctionStarts.has_key(Start):
                    self.JunctionStarts[Start] = []
                for Junction in self.JunctionEnds[End]:
                    if Junction.Start == Start:
                        Junction.Count += 1
                        FoundFlag = 1
                if not FoundFlag:
                    self.AddNewJunction(Start, End, ChromFile, ReverseFlag)
        ESTFile.close()
        ChromFile.close()
        if Task == "reportsize":
            # report the total size of all intervals seen.
            self.ReportTotalIntervalSizes()
        elif Task == "writedb":
            # Now, write everything out:
            DBName = "%s\\%s%s.filtered"%(ESTSubDirectory, ChromosomeNumber, ReverseChar)
            self.WriteIntervals(DBName)
        elif Task == "findaltsplice":
            self.FindAltSplicing(ChromosomeNumber, ReverseChar)
        elif Task == "ipi":
            self.FindAltSplicingIPI(ChromosomeNumber, ReverseChar)
        elif Task == "debugprint":
            self.DebugPrintIntervals(8997000, 9012000)
        else:
            print "** Error: Unknown task %s"%Task
    def ReportTotalIntervalSizes(self):
        TotalSize = 0
        for (Key, Count) in self.Intervals.items():
            (Start, End) = Key
            TotalSize += (End - Start) * Count
        print "Total EST interval size:", TotalSize
    def FindAltSplicingIPI(self, ChromosomeNumber, ReverseChar):
        """
        Given the database IPI.txt, write out the best evidence for alternative splicing
        within each gene.  For now, only show genes which also had an entry in KnownGene.txt.
        For the other genes, we don't know the strand, we're not sure which intervals on the
        chromosome we hit, and we may not even know the chromosome!
        """
        # COPIED AND MODIFIED from FindAltSplicing
        File = open("IPI.txt", "r")
        JSKeys = self.JunctionStarts.keys()
        JSKeys.sort()
        JEKeys = self.JunctionEnds.keys()
        JEKeys.sort()
        for FileLine in File.xreadlines():
            Bits = list(FileLine.strip().split("\t"))
            if len(Bits)<9:
                continue # blank or incomplete line
            if len(Bits[5])<4:
                continue # chromosome not known!
            Chrom = Bits[5][3:]
            if Chrom == "X":
                Chrom = 23
            if Chrom == "Y":
                Chrom = 24
            else:
                try:
                    Chrom = int(Chrom)
                except:
                    continue
            if Chrom!=ChromosomeNumber or Bits[6]!=ReverseChar:
                continue
            BestAltSpliceCount = 0
            Bits.append(0) # number of ESTs supporting the second-best splice point
            if Bits[7][0] == '"':
                Bits[7] = Bits[7][1:-1]
            if Bits[8][0] == '"':
                Bits[8] = Bits[8][1:-1]
            CommaBits = Bits[7].split(",")
            Start = int(CommaBits[0])
            CommaBits = Bits[8].split(",")
            End = int(CommaBits[-2])
            for Key in JSKeys:
                if Key > End:
                    break
                if Key < Start:
                    continue
                Junctions = self.JunctionStarts[Key]
                GoodCount = 0
                Str = ""
                Counts = []
                for Junction in Junctions:
                    if Junction.Count > 1 or Junction.Score > -15:
                        GoodCount += 1
                        Str += "to %s (%d, %.2f) "%(Junction.End, Junction.Count, Junction.Score)
                        Counts.append(Junction.Count)
                if GoodCount<2:
                    continue
                Bits.append(str(Key))
                Bits.append(Str)
                Counts.sort()
                BestAltSpliceCount = max(BestAltSpliceCount, Counts[-2])
            for Key in JEKeys:
                if Key > End:
                    break
                if Key < Start:
                    continue
                Junctions = self.JunctionEnds[Key]
                GoodCount = 0
                Str = ""
                Counts = []
                for Junction in Junctions:
                    if Junction.Count > 1 or Junction.Score > -15:
                        GoodCount += 1
                        Str += "from %s (%d, %.2f) "%(Junction.Start, Junction.Count, Junction.Score)
                        Counts.append(Junction.Count)
                if GoodCount<2:
                    continue
                Bits.append(str(Key))
                Bits.append(Str)
                Counts.sort()
                BestAltSpliceCount = max(BestAltSpliceCount, Counts[-2])
            Bits[9] = str(BestAltSpliceCount)
            print string.join(Bits, "\t")
        File.close()
    def FindAltSplicing(self, ChromosomeNumber, ReverseChar):
        """
        Read knownGene.txt, and for each gene, write out the best evidence for
        alternative splicing within that gene.
        """
        File = open("Splice\\knownGene.txt", "r")
        JSKeys = self.JunctionStarts.keys()
        JSKeys.sort()
        JEKeys = self.JunctionEnds.keys()
        JEKeys.sort()
        for FileLine in File.xreadlines():
            Bits = list(FileLine.strip().split("\t"))
            Chrom = Bits[1][3:]
            if Chrom == "X":
                Chrom = 23
            if Chrom == "Y":
                Chrom = 24
            else:
                try:
                    Chrom = int(Chrom)
                except:
                    continue
            if Chrom!=ChromosomeNumber or Bits[2]!=ReverseChar:
                continue
            BestAltSpliceCount = 0
            Bits.append(0) # number of ESTs supporting the second-best splice point
            CommaBits = Bits[8].split(",")
            Start = int(CommaBits[0])
            CommaBits = Bits[9].split(",")
            End = int(CommaBits[-2])
            for Key in JSKeys:
                if Key > End:
                    break
                if Key < Start:
                    continue
                Junctions = self.JunctionStarts[Key]
                GoodCount = 0
                Str = ""
                Counts = []
                for Junction in Junctions:
                    if Junction.Count > 1 or Junction.Score > -15:
                        GoodCount += 1
                        Str += "to %s (%d, %.2f) "%(Junction.End, Junction.Count, Junction.Score)
                        Counts.append(Junction.Count)
                if GoodCount<2:
                    continue
                Bits.append(str(Key))
                Bits.append(Str)
                Counts.sort()
                BestAltSpliceCount = max(BestAltSpliceCount, Counts[-2])
            for Key in JEKeys:
                if Key > End:
                    break
                if Key < Start:
                    continue
                Junctions = self.JunctionEnds[Key]
                GoodCount = 0
                Str = ""
                Counts = []
                for Junction in Junctions:
                    if Junction.Count > 1 or Junction.Score > -15:
                        GoodCount += 1
                        Str += "from %s (%d, %.2f) "%(Junction.Start, Junction.Count, Junction.Score)
                        Counts.append(Junction.Count)
                if GoodCount<2:
                    continue
                Bits.append(str(Key))
                Bits.append(Str)
                Counts.sort()
                BestAltSpliceCount = max(BestAltSpliceCount, Counts[-2])
            Bits[12] = str(BestAltSpliceCount)
            print string.join(Bits, "\t")
        File.close()
        
    def WriteIntervals(self, DBName):
        File = open(DBName, "wb")
        # Iterate over intervals.  When you write an interval whose
        # left point is the end of a junction, and that left endpoint
        # hasn't been processed yet, write the splice junctions leading into
        # that left point.
        Keys = self.Intervals.keys()
        Keys.sort()
        TotalIntervals = 0
        JunctionsKept = 0
        JunctionsTossed = 0
        LastEdgePosition = None
        for (Start, End) in Keys:
            TotalIntervals += 1
            Count = self.Intervals[(Start, End)]
            File.write(struct.pack("<iii", Start, End, Count))
            if Start == LastEdgePosition:
                File.write(struct.pack("<i", 0))
                continue
            JunctionList = self.JunctionEnds.get(Start, [])
            # Decide now which junctions to keep:
            KeepList = []
            for Junction in JunctionList:
                if Junction.Count > 1 or Junction.Score > -12.75:
                    JunctionsKept += 1
                    KeepList.append(Junction)
                else:
                    JunctionsTossed += 1
            File.write(struct.pack("<i", len(KeepList)))
            for Junction in KeepList:
                File.write(struct.pack("<iif", Junction.Start, Junction.Count, Junction.Score))
            #File.write(struct.pack("<i", len(KeepList)))
            #for Junction in KeepList:
            #    File.write(struct.pack("<ii", Junction.Start, Junction.Count))
        File.close()
        print "Wrote to '%s'\t%s\t intervals, \t%s\t junctions (\t%s\t bad junctions)"%(DBName, TotalIntervals, JunctionsKept, JunctionsTossed)
    def AddNewJunction(self, Start, End, ChromFile, ReverseFlag):
        NewJunction = SpliceJunction()
        NewJunction.Start = Start
        NewJunction.End = End
        NewJunction.Count = 1
        if ReverseFlag:
            ChromFile.seek(Start - 1)
            NewJunction.Sequence = ChromFile.read(4)
            ChromFile.seek(End - 6)
            NewJunction.Sequence += ChromFile.read(9)
            NewJunction.Sequence = NewJunction.Sequence.upper()
            NewJunction.Sequence = GetRC(NewJunction.Sequence)
        else:
            ChromFile.seek(Start - 3)
            NewJunction.Sequence = ChromFile.read(9)
            ChromFile.seek(End - 3)
            NewJunction.Sequence += ChromFile.read(4)
            NewJunction.Sequence = NewJunction.Sequence.upper()
        NewJunction.GetScore()
        # TEMP:
        #print NewJunction.Sequence, NewJunction.Score
        self.JunctionEnds[End].append(NewJunction)
        self.JunctionStarts[Start].append(NewJunction)
            
    def ReadKnownGeneDict(self):
        File = open("splice\\KnownGene.txt", "r")
        self.KnownGeneDict = {}
        for FileLine in File.xreadlines():
            FileLine = FileLine.strip()
            Bits = FileLine.split("\t")
            if Bits[1].find("random")!=-1:
                continue # garbage line, I think.
            self.KnownGeneDict[Bits[0].upper()] = Bits
            self.KnownGeneDict[Bits[10].upper()] = Bits
        print "Read %s knowngene keys"%len(self.KnownGeneDict.keys())
        File.close()
    def ProduceIPISpreadsheet(self):
        File = open("Database\\ipi.HUMAN.v3.11.dat", "rb")
        for FileLine in File.xreadlines():
            FileLine = FileLine.strip()
            Code = FileLine[:2]
            if Code == "ID" and len(FileLine)>2 and FileLine[2] == " ":
                CurrentID = FileLine[2:].strip().split(".")[0]
                CurrentChrom = ""
                ChromStart = ""
                ChromEnd = ""
                CurrentSprotID = ""
                CurrentKnownBits = None
                continue
            if Code == "CC":
                Pos = FileLine.find("CHROMOSOME:")
                if Pos!=-1:
                    CurrentChrom = FileLine[Pos + len("CHROMOSOME:"):]
                    if CurrentChrom[-1] == ".":
                        CurrentChrom = CurrentChrom[:-1]
                    continue
                Pos = FileLine.find("START CO-ORDINATE:")
                if Pos!=-1:
                    ChromStart = FileLine[Pos + len("START CO-ORDINATE:"):]
                    if ChromStart[-1] == ".":
                        ChromStart = ChromStart[:-1]
                    continue
                Pos = FileLine.find("END CO-ORDINATE:")
                if Pos!=-1:
                    ChromEnd = FileLine[Pos + len("END CO-ORDINATE:"):]
                    if ChromEnd[-1] == ".":
                        ChromEnd = ChromEnd[:-1]
                    continue
            if Code == "DR":
                # DR   UniProtKB/Swiss-Prot; O95793-1; STAU_HUMAN; M.
                #                                      ^^^^^^^^^^
                # DR   UniProtKB/TrEMBL; Q9NPW0; Q9NPW0_HUMAN; -.
                #                                ^^^^^^^^^^^^
                Bits = FileLine.strip().split(";")
                if len(Bits)>2:
                    TryID = Bits[2].strip()
                    #print TryID
                    if self.KnownGeneDict.has_key(TryID):
                        CurrentSprotID = TryID
                        CurrentKnownBits = self.KnownGeneDict[TryID]
                continue
            if Code == "//": # end of record
                Str = "%s\t%s\t%s\t%s\t"%(CurrentID, CurrentChrom, ChromStart, ChromEnd)
                if CurrentSprotID:
                    Str += "%s\t%s\t%s\t%s\t%s\t"%(CurrentSprotID, CurrentKnownBits[1], CurrentKnownBits[2],
                                               CurrentKnownBits[8], CurrentKnownBits[9])
                print Str
    def DebugPrintIntervals(self, CoverageStart, CoverageEnd):
        "Print all intervals that overlap the specified genomic coordinates."
        Keys = self.Intervals.keys()
        Keys.sort()
        for (Start, End) in Keys:
            if CoverageEnd!=None and Start > CoverageEnd:
                continue
            if CoverageStart!=None and End < CoverageStart:
                continue
            print "Interval (%s, %s) has count %s"%(Start, End, self.Intervals[(Start, End)])
            for Junction in self.JunctionEnds.get(Start, []):
                print "Incoming junction %s-%s %s %s"%(Junction.Start, Junction.End, Junction.Count, Junction.Score)
            for Junction in self.JunctionStarts.get(End, []):
                print "Outgoing junction %s-%s %s %s"%(Junction.Start, Junction.End, Junction.Count, Junction.Score)
        
if __name__ == "__main__":
    import psyco
    psyco.full()
    FM = FilterMeister()
    #FM.Read(3, 1, "debugprint")
    #sys.exit(1)
    ##############
    #FM.ReadKnownGeneDict() # for IPI only
    #FM.ProduceIPISpreadsheet()
    #sys.exit(1)
    ##############
    
    for ChromosomeNumber in range(1, 49):
        for ReverseFlag in (0, 1):
            FM = FilterMeister()
            FM.Read(ChromosomeNumber, ReverseFlag, "reportsize")
            #FM.Read(ChromosomeNumber, ReverseFlag, "writedb")
            #FM.Read(ChromosomeNumber, ReverseFlag, "findaltsplice")
            #FM.Read(ChromosomeNumber, ReverseFlag, "ipi")