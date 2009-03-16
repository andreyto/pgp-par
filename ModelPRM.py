"""
Train and test the PRM scoring model, for blind search.  (Also for tag generation?)
"""
import os
import sys
import time
from Utils import *
Global.FixedMods = {"C":57.0518}
Initialize()

InspectPath = "inspect_d.exe"
InspectScriptPath = "c:\\ms\\TrainingScripts\\TrainingSetBlind.txt"

def Main(TrainingOracle, TrainingDir, TestingOracle, TestingDir):
    # Outline:
    # - Run inspect train prm
    # - Run ScorpionTrainer.py to generate BN tables
    # - Run inspect test prm
    # - Summarize the histogram 
    TrainStartTime = time.clock()
    #Command = """%s train cc "%s" "%s" """%(InspectPath, TrainingOracle, TrainingDir)
    #print Command
    #os.system(Command)
    TrainEndTime = time.clock()
    print 
    print "Training runtime: %s seconds"%(TrainEndTime - TrainStartTime)
    TestStartTime = time.clock()
    Command = """%s -i "%s" -o PRMTestingOutput.txt"""%(InspectPath, InspectScriptPath)
    print Command
    os.system(Command)
    TestEndTime = time.clock()
    print 
    print "Testing runtime: %s seconds"%(TestEndTime - TestStartTime)
    Oracle = ReadTrainingOracle(TrainingOracle)
    PrintResultHistogram("PRMTestingOutput.txt", Oracle)

def PrintResultHistogram(ResultsFileName, Oracle):
    CurrentSpectrum = None
    MatchFlag = 0
    DeltaMatchFlag = 0
    MatchIndex = 0
    File = open(ResultsFileName, "r")
    ExactHistogram = {}
    DeltaHistogram = {}
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        FileName = Bits[0].split("\\")[-1]
        if FileName != CurrentSpectrum:
            if not MatchFlag:
                ExactHistogram[999] = ExactHistogram.get(999, 0) + 1
            if not DeltaMatchFlag:
                DeltaHistogram[999] = DeltaHistogram.get(999, 0) + 1
            CurrentSpectrum = FileName
            MatchFlag = 0
            DeltaMatchFlag = 0
            MatchIndex = 0
        MatchIndex += 1
        TruePeptide = Oracle.get(FileName, None)
        if not TruePeptide:
            print "** Warning: file '%s' has no match in the oracle"%FileName
        else:
            Peptide = GetPeptideFromModdedName(Bits[9][2:-2])
            if not DeltaMatchFlag:
                if Bits[1] == TruePeptide.Aminos or Bits[1][1:] == TruePeptide.Aminos or Bits[1][:-1] == TruePeptide.Aminos \
                   or Bits[1] == TruePeptide.Aminos[1:] or Bits[1] == TruePeptide.Aminos[:-1]:
                    DeltaMatchFlag = 1
                    DeltaHistogram[MatchIndex] = DeltaHistogram.get(MatchIndex, 0) + 1
            if not MatchFlag and Peptide.IsSame(TruePeptide):
                MatchFlag = 1
                ExactHistogram[MatchIndex] = ExactHistogram.get(MatchIndex, 0) + 1
    if not MatchFlag:
        ExactHistogram[999] = ExactHistogram.get(999, 0) + 1
    if not DeltaMatchFlag:
        DeltaHistogram[999] = DeltaHistogram.get(999, 0) + 1

    for Histogram in (DeltaHistogram, ExactHistogram):
        if Histogram == DeltaHistogram:
            print "Delta-correct match histogram:"
        else:
            print "Exact match histogram:"
        Keys = Histogram.keys()
        Keys.sort()
        TotalRows = 0
        for Key in Keys:
            TotalRows += Histogram[Key]
        TotalRows = max(1, TotalRows)
        Cumulation = 0
        for Key in Keys:
            Cumulation += Histogram[Key]
            print "%s\t%s\t%.2f\t%.2f\t"%(Key, Histogram[Key], Histogram[Key] / float(TotalRows), Cumulation / float(TotalRows))

if __name__ == "__main__":
    TrainingOracle = sys.argv[1]
    TrainingDir = sys.argv[2]
    TestingOracle = sys.argv[3]
    TestingDir = sys.argv[4]
    Main(TrainingOracle, TrainingDir, TestingOracle, TestingDir)
    #PrintResultHistogram("TrainingSetOldPRMOut.txt", ReadTrainingOracle(TrainingOracle))
    #PrintResultHistogram("TrainingSetNewPRMOut.txt", ReadTrainingOracle(TrainingOracle))
    #PrintResultHistogram("TrainingSetNewPRMOut.txt", ReadTrainingOracle(TrainingOracle))
    