"""
For processing EST-to-genome alignments into a splice-tolerant database.
Produce a simple list of intervals-groups, SORTED by startposition.
"""
import sys
import os

def SortIntervals(FileName, OutputFileName):
    File = open(FileName, "rb")
    OutputFile = open(OutputFileName, "wb")
    RecordLineNumber = 0
    LineNumber = 0
    IntervalGroups = []
    while (1):
        FilePos = File.tell()
        FileLine = File.readline()
        if not FileLine:
            break
        FileLine = FileLine.strip()
        LineNumber += 1
        if LineNumber%1000 == 0:
            print "Line %s..."%LineNumber
##            if LineNumber > 50000:
##                break #%%%
        if FileLine == "sim4begin":
            RecordLineNumber = 0
            PreviousInterval = None
            RecordStartFilePos = FilePos
            State = 1 # innocent, until proven guilty :)
            IntervalGroups.append([])
            continue
        if FileLine == "sim4end":
            #print "*" # delimiter
            continue
        RecordLineNumber += 1
        if RecordLineNumber == 1: # Line 1: Chromosome position and orientation
            # Line of the form "2[339-0-0] 0[227615646-227622915] <279-0-96-forward-unknown>"
            Bits = FileLine.split("-")
            if len(Bits)<7:
                # a bogus-looking line.
                State = 0
                continue
            Bracket = FileLine.find("[")
            Bracket = FileLine.find("[", Bracket + 1)
            EndBracket = FileLine.find("]", Bracket)
            RangeBits = FileLine[Bracket+1:EndBracket].split("-")
            OverallStart = int(RangeBits[0])
            OverallEnd = int(RangeBits[1])
        elif RecordLineNumber == 2: #edef, the GI number
            GINumber = int(FileLine.split("|")[1])
        elif RecordLineNumber == 3: #ddef, the chromosome number
            Chrom = FileLine.split()[0]
            Chrom = Chrom[6:]
        else: # an exon line
            # Line format: 377-622 (117076-117321) <246-0-100> ->
            Paren = FileLine.find("(")
            Paren2 = FileLine.find(")")
            PosBits = FileLine[Paren+1:Paren2].split("-")
            Start = OverallStart + int(PosBits[0]) - 1
            End = OverallStart + int(PosBits[1])
            IntervalGroups[-1].append((Start, End))
            #print "%s %s"%(Start, End)
    print "Sort..."
    IntervalGroups.sort()
    print "Write out..."
    for Group in IntervalGroups:
        #print "(group of %d intervals)"%len(Group)
        OutputFile.write("*\n")
        for Interval in Group:
            OutputFile.write("%s %s\n"%Interval)
            #print " ",Interval

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "<no psyco>"
    SortIntervals(sys.argv[1], sys.argv[2])