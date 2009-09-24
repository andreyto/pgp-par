"""
Generate an image, and summary text file, summaring the EST coverage of a genomic interval.
"""
import os
import sys
import struct
from PIL import Image
from PIL import ImageDraw

class Colors:
    "Color scheme for coverage plots"
    White = (255,255,255)
    Green = (0,255,0)
    Blue = (0,0,255)
    PaleBlue = (10,10,80)
    Red = (255,0,0)
    Grey = (155,155,155)
    ResidueNumber = (155, 155, 155)
    #Grey = (0,0,0)
    Background = (255, 255, 255)
    # Color schemes for different occurrence counts:
    Occurrences = [(255,0,0), (155, 155, 155), (145, 145, 145),
                   (135, 135, 135), (125, 125, 125), (115, 115, 115),
                   (105, 105, 105), (95, 95, 95), (85, 85, 85),
                   (75, 75, 75), (65, 65, 65), (55, 55, 55),
                   (45, 45, 45), (35, 35, 35), (25, 25, 25),
                   (15, 15, 15), (5, 5, 5), (0,0,0)]

def ReadESTs():
    File = open("NewSpliceDB\\%s%s.sorted"%(Chromosome, ReverseChar), "rb")
    ESTs = []
    while (1):
        RecordCount = File.read(4)
        if not RecordCount:
            break
        RecordCount = struct.unpack("<i", RecordCount)[0]
        File.read(4) # original file pos - discard
        ESTs.append([])
        for Index in range(RecordCount):
            (EStart, EEnd) = struct.unpack("<ii", File.read(8))
            #print "Int %s of %s: %s to %s"%(Index, RecordCount, EStart, EEnd)
            ESTs[-1].append((EStart, EEnd))
        if ESTs[-1][0][0] > End or ESTs[-1][-1][1] < Start:
            del ESTs[-1] # doesn't hit interesting area.
    File.close()
    return ESTs
##
##    
##    ESTs = [] # list entries are ESTs.  ESTs have the form [(Start,End), (Start, End),...]
##    File = open(r"C:\ftproot\liliana\polishes-good.noal", "rb")
##    RecordLineNumber = 0
##    OverallLineNumber = 0
##    while 1:
##        FileLine = File.readline()
##        if not FileLine:
##            break
##        OverallLineNumber += 1
##        if OverallLineNumber % 10000 == 0:
##            #print OverallLineNumber, OverallLineNumber / 51849830.0
##            pass
##        FileLine = FileLine.strip()
##        if FileLine == "sim4end":
##            # Remove the last EST, if it's no good:
##            if RecordOK:
##                IntervalsHit = 0
##                for (IntervalStart, IntervalEnd) in ESTs[-1]:
##                    if IntervalStart > Start and IntervalStart < End:
##                        IntervalsHit = 1
##                    if IntervalEnd > Start and IntervalEnd < End:
##                        IntervalsHit = 1
##                if not IntervalsHit:
##                    del ESTs[-1]
##                else:
##                    print "EST found:", ESTs[-1]
##                    print "FilePos %s, line number %s"%(File.tell(), OverallLineNumber)
##            continue
##        if FileLine == "sim4begin":
##            RecordLineNumber = 0
##            RecordOK = 1
##            continue
##        RecordLineNumber += 1
##        if RecordLineNumber == 1:
##            if FileLine.find("-forward-")!=-1:
##                RecordReverseFlag = 0
##            else:
##                RecordReverseFlag = 1
##            if RecordReverseFlag != ReverseFlag:
##                RecordOK = 0
##            BracketPos = FileLine.find("[")
##            BracketPos = FileLine.find("[", BracketPos + 1)
##            BrackEndPos = FileLine.find("]", BracketPos + 1)
##            (StartStr, EndStr) = FileLine[BracketPos + 1:BrackEndPos].split("-")
##            OverallStart = int(StartStr)
##            OverallEnd = int(EndStr)
##            if int(StartStr) > Start or int(EndStr) < End:
##                RecordOK = 0
##            # 843[322-0-0] 0[227952347-227956642] <311-0-96-forward-unknown>
##            if RecordOK:
##                ESTs.append([])
##        if RecordLineNumber > 3 and RecordOK:
##            # 1-197 (1996-2192) <191-0-96>
##            Paren = FileLine.find("(")
##            Paren2 = FileLine.find(")")
##            PosBits = FileLine[Paren+1:Paren2].split("-")
##            IntervalStart = OverallStart + int(PosBits[0]) - 1
##            IntervalEnd = OverallStart + int(PosBits[1])
##            ESTs[-1].append((IntervalStart, IntervalEnd))
##    return ESTs

def DrawImage(ESTs):
    global Start
    global End
    if not ESTs:
        return    
    # Adjust starting position close to the first one seen:
    ImageStart = Start
    if Start < ESTs[0][0][0]:
        ImageStart = min(End - 500, ESTs[0][0][0])
    ImageEnd = min(Start + 2000, End)
    #print Start, End
    Height = min(len(ESTs) + 100, 1024+100)
    ESTs = ESTs[:1024]
    Width = ImageEnd - ImageStart
    ImageIndex = 0
    while 1:
        CoverImage = Image.new("RGB", (Width, Height), Colors.Background)  # mode, size, [startcolor]
        Draw = ImageDraw.Draw(CoverImage)
        Header = 20
        Draw.line((0, Header, Width , Header), Colors.ResidueNumber)
        Pos = ImageStart - (ImageStart % 100) + 100
        #print Pos
        # Draw x axis:
        while Pos < ImageEnd:
            X = Pos - ImageStart
            Draw.text((X, 2), "%s"%Pos, Colors.ResidueNumber)
            Draw.line((X, Header, X, Header - 2), Colors.ResidueNumber)
            Pos += 100
        # Draw true intervals:
        for (IntervalStart, IntervalEnd) in TrueIntervals:
            for Y in range(21,24):
                Draw.line((IntervalStart - ImageStart, Y, IntervalEnd - ImageStart, Y), Colors.Red)
                #print "TRUE:", IntervalStart, IntervalEnd
        for ESTIndex in range(len(ESTs)):
            Y = ESTIndex + 25
            for (ESTStart, ESTEnd) in ESTs[ESTIndex]:
                #print "X:", ESTStart-Start, ESTEnd-Start
                Draw.line((ESTStart - ImageStart, Y, ESTEnd - ImageStart, Y), Colors.Green)
                #print "EST:", ESTStart, ESTEnd
        CoverImage.save("ESTCoverage%s.png"%ImageIndex, "png")
        if ImageEnd >= End:
            break
        ImageIndex += 1
        ImageStart += 1500
        ImageEnd += 1500

def SummarizeESTs(ESTs):
    print "%s ESTs map to the region of interest."%len(ESTs)
    for EST in ESTs:
        Str = ""
        for Interval in EST:
            Str += "%s-%s, "%Interval
        print Str
    print ">>>the end<<<"

if __name__ == "__main__":
    if len(sys.argv) == 1:
        Chromosome = 5
        ReverseFlag = 1
        ReverseChar = "-"
        Start = 17590093 - 100
        End = 17627806 + 100
        TrueIntervals = [(17603658, 17638175)]
        
##        Chromosome = 21
##        ReverseFlag = 0
##        ReverseChar = "+"
##        Start = 43462209 - 100
##        End = 43465982 + 100
##        TrueIntervals = [(43462209, 43462467),
##                         (43463695, 43463818),
##                         (43465249, 43465982)]
    else:
        TrueIntervals = []
        GeneName = sys.argv[1]
        File = open("KnownGeneAltSplicing.txt", "r") #%%
        #File = open("Splice\\KnownGene.txt", "r")
        Start = None
        for FileLine in File.xreadlines():
            Bits = FileLine.strip().split("\t")
            if Bits[10] == GeneName:
                StartBits = Bits[8].split(",")
                EndBits = Bits[9].split(",")
                for Index in range(len(StartBits) - 1):
                    TrueIntervals.append((int(StartBits[Index]), int(EndBits[Index])))
                Chromosome = Bits[1][3:]
                if Chromosome == "X":
                    Chromosome = 23
                elif Chromosome == "Y":
                    Chromosome = 24
                else:
                    Chromosome = int(Chromosome)
                ReverseChar = Bits[2]
                if ReverseChar == "+":
                    ReverseFlag = 0
                else:
                    ReverseFlag = 1
                Start = TrueIntervals[0][0] - 100
                End = TrueIntervals[-1][1] - 100
                break
        if Start == None:
            print "** Unknown gene: '%s'"%GeneName
            sys.exit(1)
        print "Region of interest: %s-%s on chromosome %s"%(Start, End, Chromosome)
    
    ESTs = ReadESTs()
    ########################################################
##    # Weird sort:
##    BozoSort = []
##    for EST in ESTs:
##        BozoSort.append((EST[-1][1], EST))
##    BozoSort.sort()
##    BozoSort.reverse()
##    ESTs = []
##    for (Dummy, EST) in BozoSort:
##        ESTs.append(EST)
##        print EST[-1][1] #%%
    ########################################################    
    #ESTs.sort()
    SummarizeESTs(ESTs)
    DrawImage(ESTs)