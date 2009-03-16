"""
Nelder-Mead simplex optimization to find the best weights for tag-scoring
parameters.  Works about as well as an SVM.
"""
import os
import math
import sys
import traceback
import random


# Simplex optimization loaf:
class SimplexPoint:
    def __init__(self, Dimensions):
        self.Coords = []
        for X in range(Dimensions):
            self.Coords.append(10*random.random())
    def xGetScore(self):
        Val = 100 * (self.Coords[1] - self.Coords[0]**2)**2 + (1-self.Coords[0])**2
        print "  SCORED:",self.Coords, Val
        self.Score = Val
        return Val
    def GetScore(self):
        global GoodCoords
        global EvilCoords
        ErrorCount = 0
        OverallScore = 0
        BlipSize = 1
        BlipSize2 = 0.1
        BlipSize3 = 0.01
        for Coords in GoodCoords:
            Score = 0
            for Index in range(len(self.Coords)):
                Score += self.Coords[Index] * Coords[Index]
            #Score = self.X * Coords[0] + self.Y * Coords[1] + self.Z * Coords[2]
            if Score < BlipSize:
                ErrorCount += (BlipSize-Score)**2 * GoodScalor
                #ErrorCount += 1
##            if Score < -BlipSize:
##                ErrorCount += 1 # doubleplusbad!
            
        for Coords in EvilCoords:
            #Score = self.X * Coords[0] + self.Y * Coords[1] + self.Z * Coords[2]
            #print "Bad:", Score
            Score = 0
            for Index in range(len(self.Coords)):
                Score += self.Coords[Index] * Coords[Index]            
            if Score > -BlipSize:
                #ErrorCount += 1
                ErrorCount += (Score + BlipSize)**2
        #print "Overall:", OverallScore
        OverallScore = ErrorCount
        self.Score = OverallScore
        #self.MeanGood = MeanGood / len(GoodCoords)
        #self.MeanEvil = MeanEvil / len(EvilCoords)        
        return OverallScore
        #self.Score = -ErrorCount
        #return self.Score
        #self.Score = (MeanGood - MeanEvil)
        #return self.Score
    def GetMeans(self):
        MeanGood = 0
        MeanEvil = 0
        for Coords in GoodCoords:
            Score = 0
            for Index in range(len(self.Coords)):
                Score += self.Coords[Index] * Coords[Index]            
            MeanGood += Score
        for Coords in EvilCoords:
            Score = 0
            for Index in range(len(self.Coords)):
                Score += self.Coords[Index] * Coords[Index]            
            MeanEvil += Score
        self.MeanGood = MeanGood / len(GoodCoords)
        self.MeanEvil = MeanEvil / len(EvilCoords)        
    def __cmp__(self, Other):
        if self.Score < Other.Score:
            return -1
        if self.Score > Other.Score:
            return 1
        return 0
    def __str__(self):
        Str = "<"
        for X in self.Coords:
            Str += "%.3f,"%X
        return Str+":%.4f>"%self.Score
        #return "<%.3f, %.3f, %.3f sc %.3f>"%(self.X, self.Y, self.Z, self.Score)
    def PrintROCCurve(self, Goods, Evils):
        Cutoff = -10
        print "\n\nROC curve:"
        while Cutoff < 10:
            TP = 0
            FP = 0
            for Coords in Goods:
                Score = 0
                for X in range(len(self.Coords)):
                    Score += self.Coords[X]*Coords[X]
                if Score > Cutoff:
                    TP += 1
            TPRate = TP / float(len(Goods))
            for Coords in Evils:
                Score = 0
                for X in range(len(self.Coords)):
                    Score += self.Coords[X]*Coords[X]
                if Score > Cutoff:
                    FP += 1
            FPRate = FP / float(len(Evils))
            print "%.2f\t%s\t%s\t"%(Cutoff, FPRate, TPRate)
            Cutoff += 0.1



def GetGoodCoords(GoodCoords, EvilCoords, FileName,
                  Columns, TrueFlagColumn):
    File = open(FileName, "r")
    Means = [0]*len(Columns)
    StdDevs = [0]*len(Columns)
    LineCount = 0
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        if len(Bits)<2 or not Bits[0]:
            continue
        LineCount += 1
        for X in range(len(Columns)):
            Means[X] += float(Bits[Columns[X]])
    for X in range(len(Columns)):
        Means[X] /= float(LineCount)
    File.seek(0)
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        if len(Bits)<2 or not Bits[0]:
            continue
        for X in range(len(Columns)):
            StdDevs[X] += (float(Bits[Columns[X]]) - Means[X])**2
    for X in range(len(Columns)):
        StdDevs[X] = math.sqrt(StdDevs[X] / LineCount)
    for X in range(len(Columns)):
        print "Column %s: Mean %.2f, stddev %.2f"%(Columns[X], Means[X], StdDevs[X])
    File.seek(0)
    for FileLine in File.xreadlines():
        FileLine = FileLine.strip()
        Bits = FileLine.split("\t")
        if len(Bits)<2 or not Bits[0] or Bits[0][0]=="#":
            continue
        if Bits[0][0]=="*" and not Means:
            for X in range(len(Columns)):
                Means.append(float(Bits[Columns[X]]))
            continue
        if Bits[0][0]=="%" and not StdDevs:
            for X in range(len(Columns)):
                StdDevs.append(float(Bits[Columns[X]]))
            continue
        Coords = []
        for X in range(len(Columns)):
            Coords.append( (float(Bits[Columns[X]]) - Means[X]) / StdDevs[X])
        if Bits[TrueFlagColumn]=="1":
            GoodCoords.append(Coords)
        else:
            EvilCoords.append(Coords)
    print "Mean:", Means
    print "StdDev:", StdDevs
    #return (Means, StdDev)

def SimplexOptimize(CoordCount):
    SimplexPoints = []
    random.seed()
    for PointIndex in range(CoordCount + 1):
        Point = SimplexPoint(CoordCount)
        Point.GetScore()
        SimplexPoints.append(Point)
    OldError = 9999999
    for Cycle in range(500):
        SimplexPoints.sort()
        #if abs(SimplexPoints[-1].Score - OldError)<0.00005:
        #    break
        OldError = SimplexPoints[-1].Score
        print "Points:", map(str, SimplexPoints)
        #print SimplexPoints[-1] #%%%
        OldCoords = SimplexPoints[-1].Coords[:]
        CentroidCoords = [0]*CoordCount
        for Point in SimplexPoints[:-1]:
            for X in range(CoordCount):
                CentroidCoords[X] += Point.Coords[X]
        for X in range(CoordCount):
            CentroidCoords[X] /= len(SimplexPoints)-1
        
        # Reflection point:
        ReflectionPoint = SimplexPoint(CoordCount)
        for X in range(CoordCount):
            ReflectionPoint.Coords[X] = OldCoords[X] + (CentroidCoords[X] - OldCoords[X])*2
        ReflectionScore = ReflectionPoint.GetScore()
        #print "Centroid:", CentroidCoords, "Reflection:", ReflectionPoint
        # If reflection is average:
        if SimplexPoints[0].Score <= ReflectionScore and ReflectionScore < SimplexPoints[-2].Score:
            SimplexPoints[-1] = ReflectionPoint
            #print "REFLECTED."
            continue
        # If reflection is best:
        if ReflectionScore < SimplexPoints[0].Score:
            ExpansionPoint = SimplexPoint(CoordCount)
            for X in range(CoordCount):
                ExpansionPoint.Coords[X] = OldCoords[X] + (CentroidCoords[X] - OldCoords[X])*3
            ExpansionPoint.GetScore()
            if ExpansionPoint.Score < ReflectionPoint.Score:
                SimplexPoints[-1] = ExpansionPoint
                #print "EXPANDED."
            else:
                #print "REFLECTED."
                SimplexPoints[-1] = ReflectionPoint
            continue
        # If reflection score beats the old point's score:
        if ReflectionScore < SimplexPoints[-1].Score:
            # Outside:
            ContractionPoint = SimplexPoint(CoordCount)
            for X in range(CoordCount):
                ContractionPoint.Coords[X] = OldCoords[X] + (CentroidCoords[X] - OldCoords[X])*1.5
            if ContractionPoint.GetScore() < ReflectionScore:
                SimplexPoints[-1] = ContractionPoint
                #print "OUTER CONTRACTION."
                continue
        else:
            # Inside:
            ContractionPoint = SimplexPoint(CoordCount)
            for X in range(CoordCount):
                ContractionPoint.Coords[X] = OldCoords[X] + (CentroidCoords[X] - OldCoords[X])*0.5
            if ContractionPoint.GetScore() < SimplexPoints[-2].Score:
                SimplexPoints[-1] = ContractionPoint
                #print "INNER CONTRACTION."
                continue
        # Shrikage:
        for Index in range(1, len(SimplexPoints)):
            for X in range(CoordCount):
                SimplexPoints[Index].Coords[X] += (SimplexPoints[0].Coords[X] - SimplexPoints[Index].Coords[X]) / 2.0
        #print "SHRINKAGE."
    return SimplexPoints[0]

def FixDataSet():
    File = open("ISB_Sequences.fa", "r")
    GoodStuff = File.read()
    File.close()
    File = open("xPValueTestSet.txt", "r")
    Shuffler = []
    HeaderDone = 0
    for FileLine in File.xreadlines():
        FileLine = FileLine.strip()
        Bits = FileLine.split("\t")
        if Bits[0]=="Spectrum":
            if not HeaderDone:
                print FileLine
            HeaderDone = 1
            continue
        if Bits[25] == "0":
            print FileLine
        if Bits[2] and GoodStuff.find(Bits[2])==-1:
            Shuffler.append(FileLine)
    random.shuffle(Shuffler)
    for Line in Shuffler[:]:
        print Line
        
#FixDataSet()
#sys.exit()

if __name__ == "__main__":
    TruthColumn = 0
    #KeyColumns = [3,5,6,7,8,9,10]
    #KeyColumns = [1,3,4,5]
    KeyColumns = [1,2,3,4,5]
    GoodCoords = []
    EvilCoords = []
    GetGoodCoords(GoodCoords, EvilCoords, "PRMTrainingSet.txt", KeyColumns, TruthColumn) #"PValues2.txt")
    GoodScalor = (len(EvilCoords) / float(len(GoodCoords)))
    print "GoodScalor:", GoodScalor
    #GoodScalor = 1
######    Bob = SimplexPoint(2)
######    Bob.Coords = [1,1]
######    print Bob.GetScore()
######    sys.exit(1)
    for X in range(1):
        Point = SimplexOptimize(len(KeyColumns))
        #Point = SimplexOptimize(2)
    GoodTestCoords = []
    EvilTestCoords = []
    GetGoodCoords(GoodTestCoords, EvilTestCoords, "PRMTestSet.txt", KeyColumns, TruthColumn)
    print Point
    Point.GetMeans()
    print Point.MeanGood
    print Point.MeanEvil
    print "TRAINING SET roc curve:"
    Point.PrintROCCurve(GoodCoords, EvilCoords)
    print "TEST SET roc curve:"
    Point.PrintROCCurve(GoodTestCoords, EvilTestCoords)
