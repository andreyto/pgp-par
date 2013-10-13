"""
Train an SVM using specified combination of features and kernel.  Wraps EasySVM.py.
Useful for trying out different feature-sets.
"""
import os
import sys
import EasySVM
import random
random.seed(1)

class AllFeatures:
    Charge = "Charge"
    Score = "Score"
    Intensity = "Intensity"
    IntensityFiltered = "IntensityFiltered"
    Peaks = "Peaks"
    PeaksFiltered = "PeaksFiltered"
    Candidates = "Candidates"
    DeltaScore = "DeltaScore"
    DeltaScoreOther = "DeltaScoreOther"
    BYPeaks = "BYPeaks"
    SignalFiltered = "SignalFiltered"
    Signal = "Signal"
    Length = "Length"
    Mass = "Mass"
    MassDelta = "MassDelta"
    

TrainFileName = "PValueTrain.txt"
TestFileName = "PValueTest.txt"
    
def PreparePValueFiles(Features):
    VerboseFileNames = [("+1", "PValueISBTrue.txt", 766),
                     #("+1", "PValueIKKBTrue.txt", 200),
                     ("-1", "PValueISBFalse.txt", 766),
                     #("-1", "PValueIKKBFalse.txt", 1000),
                     ]
    TrainingFile = open(TrainFileName, "w")
    TestingFile = open(TestFileName, "w")
    
    for (IsTrue, FileName, LinesToProcess) in VerboseFileNames:
        File = open(FileName, "r")
        LineCount = 0
        FileLines = File.readlines()
        random.shuffle(FileLines)
        for FileLine in FileLines:
            LineCount += 1
            FileLine = FileLine.strip()
            Bits = FileLine.split("\t")
            if len(Bits)<10 or Bits[0]=="Spectrum":
                continue
            if not Bits[2]:
                continue # no hit at all
            Line = "%s "%IsTrue
            Feature = 1
            #######################################################
            # ADD FEATURES:
            if Features.has_key(AllFeatures.Charge):
                Charge = int(Bits[13])
                IsCharge1 = 0
                IsCharge2 = 0
                IsCharge3 = 0
                if Charge == 1:
                    IsCharge1 = 1
                elif Charge == 2:
                    IsCharge2 = 1
                else:
                    IsCharge3 = 1
                Line += "%s:%s "%(Feature,IsCharge1)
                Feature += 1
                Line += "%s:%s "%(Feature,IsCharge2)
                Feature += 1
                Line += "%s:%s "%(Feature,IsCharge3)
                Feature += 1
            if Features.has_key(AllFeatures.Score):
                Value = float(Bits[5])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.Intensity):
                Intensity = float(Bits[33])
                Line += "%s:%s "%(Feature, Intensity)
                Feature += 1
            if Features.has_key(AllFeatures.IntensityFiltered):
                Intensity = float(Bits[45])
                Line += "%s:%s "%(Feature, Intensity)
                Feature += 1
            if Features.has_key(AllFeatures.Peaks):
                Value = float(Bits[35])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.PeaksFiltered):
                Value = float(Bits[47])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.Candidates):
                Value = float(Bits[24])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.DeltaScore):
                try:
                    Value = float(Bits[28])
                except:
                    Value = 0
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.DeltaScoreOther):
                try:
                    Value = float(Bits[27])
                except:
                    Value = 0
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.BYPeaks):
                Value = float(Bits[39])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.Signal):
                Value = float(Bits[42])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.SignalFiltered):
                Value = float(Bits[43])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.Length):
                Value = len(Bits[2])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.Mass):
                Value = float(Bits[11])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
            if Features.has_key(AllFeatures.MassDelta):
                Value = float(Bits[12])
                Line += "%s:%s "%(Feature, Value)
                Feature += 1
                
            #######################################################
            # Dump this line:
            if (LineCount >= LinesToProcess/2):
                TestingFile.write(Line + "\n")
            else:
                TrainingFile.write(Line + "\n")
            LineCount += 1
            if LineCount >= LinesToProcess:
                break
    TrainingFile.close()
    TestingFile.close()

def PrepSVMFile(InFileName, OutFileName):
    InFile = open(InFileName, "r")
    OutFile = open(OutFileName, "w")
    for FileLine in InFile.xreadlines():
        Bits = FileLine.split("\t")
        if Bits[0] == "1":
            Str = "+1 "
        else:
            Str = "-1 "
        for Index in range(1, len(Bits)):
            Bit = Bits[Index].strip()
            if Bit:
                Str += "%s:%s "%(Index, Bits[Index])
        OutFile.write(Str+"\n")
    InFile.close()
    OutFile.close()
    
if __name__ == "__main__":
##    Features = {#AllFeatures.Charge:1,
##                AllFeatures.Score:1,
##                #AllFeatures.Intensity:1,
##                AllFeatures.IntensityFiltered:1,
##                #AllFeatures.Peaks:1,
##                AllFeatures.PeaksFiltered:1,
##                AllFeatures.Candidates:1,
##                #AllFeatures.DeltaScore:1,
##                AllFeatures.DeltaScoreOther:1,
##                AllFeatures.BYPeaks:1,
##                AllFeatures.SignalFiltered:1,
##                #AllFeatures.Signal:1,
##                #AllFeatures.Length:1,
##                AllFeatures.MassDelta:1,
##                }
##    PreparePValueFiles(Features)
    PrepSVMFile("PRMTrainingSet.txt", "Train.txt")
    PrepSVMFile("PRMTestSet.txt", "Test.txt")
    Result = EasySVM.Main("Train.txt", "Test.txt")
    print Result
    
    
