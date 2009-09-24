import math

def ComputeROCCurve(List):
    """
    Compute the ROC curve for a set of tuples of the form (reading, truthflag)
    """
    List.sort()
    List.reverse()
    AllPositive = 0 
    AllNegative = 0
    for (Score, Truth) in List:
        if (Truth):
            AllPositive += 1
        else:
            AllNegative += 1
    Area = 0
    TPCount = 0
    FPCount = 0
    for (Score, Truth) in List:
        if (Truth):
            TPCount += 1
        else:
            FPCount += 1
            TPRate = TPCount / float(AllPositive)
            Area += TPRate
    Area /= float(AllNegative)
    return Area

def GetMedian(List):
    LocalCopy = List[:]
    LocalCopy.sort()
    Len = len(LocalCopy)
    if Len % 2:
        return LocalCopy[Len / 2]
    return (LocalCopy[Len / 2] + LocalCopy[(Len / 2) - 1]) / 2.0
        
def Sum(List):
    Total = 0
    for Entry in List:
        Total += Entry
    return Total

def GetMedian(List):
    SortedList = List[:]
    SortedList.sort()
    Len = len(SortedList)
    if Len % 2 == 1:
        return SortedList[Len / 2]
    Score = (SortedList[Len / 2] + SortedList[(Len / 2) - 1]) / 2.0
    return Score

def GetMean(List):
    if not len(List):
        return None
    Mean = 0
    for Entry in List:
        Mean += Entry
    Mean /= float(len(List))
    return Mean

def GetMeanStdDev(List):
    "Computes mean, standard deviation for a list of numbers"
    if not len(List):
        return (0, 0)
    Mean = 0
    for Entry in List:
        Mean += Entry
    Mean /= float(len(List))
    StdDev = 0
    for Entry in List:
        StdDev += (Entry - Mean) ** 2
    StdDev = math.sqrt(StdDev / float(len(List)))
    return (Mean, StdDev)


def GetStdDev(List):
    "Computes standard deviation for a list of numbers"
    if not len(List):
        return 0.0
    Mean = 0
    for Entry in List:
        Mean += Entry
    Mean /= float(len(List))
    StdDev = 0
    for Entry in List:
        StdDev += (Entry - Mean) ** 2
    StdDev = math.sqrt(StdDev / float(len(List)))
    return StdDev

def CreateHistogram(List, FirstBin, BinSize):
    """Given a list of NUMERICAL values, bin size and first bin, we 
    make a dictionary representing the histogram.  It's your job
    to figure out the bin size and first bin.  We are only computers
    we have no reasoning.
    """
    Histogram = {} #bin->count
    List.sort()
    BiggestValue = List[-1]
    #print "THis is the BiggestValue that I found %s"%BiggestValue
    Bin = FirstBin
    #instantiate bins
    while(Bin <= BiggestValue):
        Histogram[Bin] =0
        Bin += BinSize
    #now fill in bins
    Bin = FirstBin
    BinMax = FirstBin + BinSize
    for Value in List:
        if Value < BinMax:
            Histogram[Bin] += 1
        else:
            while( not(Value < BinMax)):
                Bin += BinSize
                BinMax += BinSize
            Histogram[Bin] += 1
    return Histogram
    
def PrettyPrintHistogram(Histogram, FileHandle = None):
    """You should call this with a histogram that is created
    by BasicStats::CreateHistogram.  Other stuff might not work.
    Also, the FileHandle object must be an actual handle, not a file
    name that you want me to make a handle out of.
    """
    #print header
    if FileHandle:
        FileHandle.write("Histogram\n")
    else:
        print "Histogram\n"
    Keys = Histogram.keys()
    Keys.sort()
    for Key in Keys:
        if FileHandle:
            FileHandle.write("%s\t%s\n"%(Key, Histogram[Key]))
        else:
            print "%s\t%s"%(Key, Histogram[Key])
        
def PrettyPrintCDF(CDF, FileHandle = None):
    """You should call this with a histogram that is created
    by BasicStats::CreateHistogram.  Other stuff might not work.
    Also, the FileHandle object must be an actual handle, not a file
    name that you want me to make a handle out of.
    """
    #print header
    if FileHandle:
        FileHandle.write("Cumulative Distribution Function\n")
    else:
        print "Cumulative Distribution Function\n"
    Keys = CDF.keys()
    Keys.sort()
    for Key in Keys:
        if FileHandle:
            FileHandle.write("%s\t%s\n"%(Key, CDF[Key]))
        else:
            print "%s\t%s"%(Key, CDF[Key])

def CreateCDF(List, FirstBin, BinSize, Granularity):
    """This creates a cumulative distribution function of the data
    List is the list of data
    First Bin and bin size are the same as CreateHistogram.  Granularity
    is the level of granularity that you want things reported in, e.g. 5% increments
    """
    NumElements = len(List)
    CDF = {} # bin=> coverage
    List.sort()
    CurrAccumulationFraction = 0.0 #0 => 1 by steps = granularity
    CurrAccumulationCount =0 # 0 -> NumElements
    Histogram = CreateHistogram(List, FirstBin, BinSize) # I do this just to make counting easier
    HKeys = Histogram.keys()
    HKeys.sort() #this key is really the histogram bin
    for Key in HKeys:
        Count = Histogram[Key]
        CurrAccumulationCount += Count
        #now check to see what fraction I am at with this new data point
        NewFraction = CurrAccumulationCount/ float (NumElements)
        NewPoint = 0 #save something?
        while NewFraction > (CurrAccumulationFraction + Granularity):
            NewPoint = 1
            CurrAccumulationFraction += Granularity
        if NewPoint:
            CDF[Key] = CurrAccumulationFraction
    return CDF
        
        
    
