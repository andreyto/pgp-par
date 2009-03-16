"""
Generate a histogram of intron lengths.
"""

def GetIntronLengthHistogram():
    IntronLengths = {}
    File = open("KnownGeneAltSplicing.txt", "rb")
    AllSpans = []
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        if len(Bits)<10:
            continue
        Starts = Bits[8].split(",")
        Ends = Bits[9].split(",")
        for Index in range(len(Starts) - 2):
            IntronLength = int(Starts[Index + 1]) - int(Ends[Index])
            IntronLengths[IntronLength] = IntronLengths.get(IntronLength, 0) + 1
    Keys = IntronLengths.keys()
    Keys.sort()
    for Key in Keys:
        print "%s\t%s"%(Key, IntronLengths[Key])
            
def GetSpanHistogram():
    """
    Let's look at known mappings to decide how large of a genomic neighborhood to
    consider in each case.
    Here are the largest spans and their genes:
    (2304258, 'CNTP2_HUMAN')
    (2220382, 'P11532-4')
    (2169338, 'Q68CQ8_HUMAN')
    (2092329, 'DMD_HUMAN')
    (2090782, 'DMD_HUMAN')
    (2009201, 'P11532-4')
    (2009201, 'P11532-4')
    (1900275, 'LRP1B_HUMAN')    
    """
    File = open("KnownGeneAltSplicing.txt", "rb")
    AllSpans = []
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        if len(Bits)<10:
            continue
        Starts = Bits[8].split(",")
        Ends = Bits[9].split(",")
        Span = int(Ends[-2]) - int(Starts[0])
        AllSpans.append((Span, Bits[10]))
    AllSpans.sort()
    AllSpans.reverse()
    for X in range(1000):
        print AllSpans[X]


GetIntronLengthHistogram()        