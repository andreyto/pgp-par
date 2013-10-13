"""
Computing the output of an SVM takes longer than we'd like.  Let's try condensing
the list of support vectors, to save time:
- Find the closest two positive vectors.  Combine them into a vector whose coordinates
are averaged, and whose strength is the sum of the old strengths.  Do the same for
the closest two negative vectors.
- Continue this process until the number of vectors remaining is small, or the distances
between the next two pairs are large.
"""
import sys
import os
import traceback
import math
import psyco

class VectorClass:
    def __init__(self, Weight, Size):
        self.NearestNeighbor = None
        self.NearestDistance = None
        self.Features = [0]*Size
        self.Weight = abs(Weight)
        if Weight > 0:
            self.Category = 1
        else:
            self.Category = -1
    def Assimilate(self, OtherVector):
        NewFeatures = []
        for Z in range(len(self.Features)):
            Value = (self.Features[Z] * self.Weight + OtherVector.Features[Z] * OtherVector.Weight) / float(self.Weight + OtherVector.Weight)
            NewFeatures.append(Value)
        self.Features = NewFeatures
        self.Weight += OtherVector.Weight
    def __str__(self):
        Str = "%.4f "%(self.Weight * self.Category)
        for Index in range(len(self.Features)):
            Str += "%d:%.4f "%(Index+1, self.Features[Index])
        return Str
    def GetDistance(self, OtherVector):
        Distance = 0
        for Z in range(len(self.Features)):
            Distance += (self.Features[Z] - OtherVector.Features[Z])**2
        Distance = math.sqrt(Distance)
        return Distance

def PrintVectors(VectorList):
    for Vector in VectorList:
        print Vector
        
PositiveVectors = []
NegativeVectors = []


def ResetNeighbor(Vector, VectorList):
    Vector.NearestDistance = 9999
    for OtherVector in VectorList:
        if OtherVector == Vector:
            continue
        Distance = Vector.GetDistance(OtherVector)
        if Distance < Vector.NearestDistance:
            Vector.NearestDistance = Distance
            Vector.NearestNeighbor = OtherVector

def ReadVectors(FileName):
    global NegativeVectors
    global PositiveVectors
    File = open(FileName, "r")
    for FileLine in File.xreadlines():
        Bits = FileLine.split()
        try:
            Category = float(Bits[0])
        except:
            continue
        Vector = VectorClass(Category, 6)
        for Bit in Bits[1:]:
            (Coord, Value) = Bit.split(":")
            Vector.Features[int(Coord) - 1] = float(Value)
        if Category>0:
            PositiveVectors.append(Vector)
        else:
            NegativeVectors.append(Vector)
    print "Set positive vector neighbors:"
    Len = len(PositiveVectors)
    for X in range(Len):
        Vector = PositiveVectors[X]
        if X%100 == 0:
            print "%s/%s"%(X, Len)        
        ResetNeighbor(Vector, PositiveVectors)
    print "Set negative vector neighbors:"
    Len = len(NegativeVectors)
    for X in range(Len):
        Vector = NegativeVectors[X]
        if X%100 == 0:
            print "%s/%s"%(X, Len)
        ResetNeighbor(Vector, NegativeVectors)

def GetClosestPair(VectorList):
    BestDistance = 9999
    for Vector in VectorList:
        if Vector.NearestDistance < BestDistance:
            BestDistance = Vector.NearestDistance
            BestVector = Vector
    return (BestVector, BestVector.NearestNeighbor, BestDistance)

if __name__ == "__main__":
    psyco.full()
    ReadVectors(sys.argv[1])
    MergeCutoff = 0.1
    DoneMergingNeg = 0
    DoneMergingPos = 0
    while (1):
        if len(NegativeVectors)<5 or len(PositiveVectors)<5:
            break
        if not DoneMergingNeg:
            (NA, NB, NDist) = GetClosestPair(NegativeVectors)
            if NDist > MergeCutoff:
                DoneMergingNeg = 1
            else:
                NA.Assimilate(NB)
                NegativeVectors.remove(NB)
                for Vector in NegativeVectors:
                    if Vector.NearestNeighbor in (NA, NB):
                        ResetNeighbor(Vector, NegativeVectors)
                    else:
                        Distance = Vector.GetDistance(NA)
                        if Distance < Vector.NearestDistance:
                            Vector.NearestDistance = Distance
                            Vector.NearestNeighbor = NA
                print "%s Merged negative vectors (dist %.4f)"%(len(NegativeVectors), NDist)
        if not DoneMergingPos:
            (NA, NB, NDist) = GetClosestPair(PositiveVectors)
            if NDist > MergeCutoff:
                DoneMergingPos = 1
            else:
                NA.Assimilate(NB)
                PositiveVectors.remove(NB)
                for Vector in PositiveVectors:
                    if Vector.NearestNeighbor in (NA, NB):
                        ResetNeighbor(Vector, PositiveVectors)
                    else:
                        Distance = Vector.GetDistance(NA)
                        if Distance < Vector.NearestDistance:
                            Vector.NearestDistance = Distance
                            Vector.NearestNeighbor = NA
                print "%s Merged positive vectors (dist %.4f)"%(len(PositiveVectors), NDist)
        if DoneMergingNeg and DoneMergingPos:
            break

       
    PrintVectors(PositiveVectors)
    PrintVectors(NegativeVectors)