"""
Script for supporting the construction of a great scoring model.
Given some features - nodes in a bayesian network - we'd like to find
an optimal way to link them up. We'll make some assumptions, based on
causal intuitions, about who is allowed to be a parent of whom. And
we'll limit the total number of parents for any single node to a
reasonable number - say, three. Our goal is to pick the set of parents
which has maximum predictive value.
"""

import os
import sys
import traceback
import math
import struct

MORE_INTENSE_IS_ALWAYS_BETTER = 0

class ScorpionTrainer:
    def ReadFeatures(self, FileName = "ScorpionFeatures.txt"):
        File = open(FileName, "r")
        self.FeatureLines = []
        self.Headers = "DynamicRange\ty\tb\tyi\tbi\ty2\tb2\tyh2o\ta\tbh2o\tynh3\tbnh3\tcharge\tflankb\tflanky\tsector\ta\tb\tc\td\te\tf".split("\t")
        for FileLine in File.xreadlines():
            if FileLine[0] == "#":
                self.Headers = []
                Bits = FileLine[1:].strip().split()
                for Index in range(len(Bits)):
                    self.Headers.append(Bits[Index])
                continue # skip comment line
            Bits = FileLine.strip().split()
            Line = []
            for Bit in Bits:
                try:
                    Line.append(int(Bit))
                except:
                    Line.append(Bit)
            self.FeatureLines.append(Line)
        File.close()
    def FindGoodParent(self, Index, PossibleParentList):
        ParentListLen = len(PossibleParentList)
        BestScore = None
        BestSubset = []
        SubList = []
        for FirstParentIndex in range(ParentListLen):
            SubList = [PossibleParentList[FirstParentIndex]]
            Score = self.GetParentingScore(Index, SubList)
            print "Parents %s score %s"%(SubList, Score)
            if (BestScore==None or Score > BestScore):
                BestScore = Score
                BestSubset = SubList
        print "Best single parent is %s with score %s"%(BestSubset, BestScore)
    def FindGoodParent2(self, Index, PossibleParentList):
        ParentListLen = len(PossibleParentList)
        BestScore = None
        BestSubset = []
        SubList = []
        for FirstParentIndex in range(ParentListLen):
            for SecondParentIndex in range(FirstParentIndex+1, ParentListLen):
                SubList = [PossibleParentList[FirstParentIndex], PossibleParentList[SecondParentIndex]]
                Score = self.GetParentingScore(Index, SubList)
                print "Parents %s score %s"%(SubList, Score)
                if (BestScore==None or Score > BestScore):
                    BestScore = Score
                    BestSubset = SubList
        print "Best parents is %s with score %s"%(BestSubset, BestScore)
    def FindGoodParent3(self, Index, PossibleParentList):
        ParentListLen = len(PossibleParentList)
        BestScore = None
        BestSubset = []
        SubList = []
        for FirstParentIndex in range(ParentListLen):
            for SecondParentIndex in range(FirstParentIndex+1, ParentListLen):
                for ThirdParentIndex in range(SecondParentIndex+1, ParentListLen):
                    SubList = [PossibleParentList[FirstParentIndex],
                               PossibleParentList[SecondParentIndex],
                               PossibleParentList[ThirdParentIndex],
                               ]
                    Score = self.GetParentingScore(Index, SubList)
                    print "Parents %s score %s"%(SubList, Score)
                    if (BestScore==None or Score > BestScore):
                        BestScore = Score
                        BestSubset = SubList
        print "Best parents is %s with score %s"%(BestSubset, BestScore)
    def GetParentingScore(self, Feature, ParentFeatures):
        Dict = {}
        for FeatureLine in self.FeatureLines:
            # For now, skip rows that aren't charge 2, within dynamic range
            if FeatureLine[0]!=3 or FeatureLine[12]!=2:
                continue
            ChildValue = FeatureLine[Feature]
            Parentals = []
            for Parent in ParentFeatures:
                Parentals.append(FeatureLine[Parent])
            Parentals = tuple(Parentals)
            if not Dict.has_key(Parentals):
                Dict[Parentals] = {}
            Dict[Parentals][ChildValue] = Dict[Parentals].get(ChildValue, 0) + 1
            Dict[Parentals][None] = Dict[Parentals].get(None, 0) + 1
        Score = 0
        for Parentals in Dict.keys():
            for ChildValue in Dict[Parentals].keys():
                if ChildValue == None:
                    continue
                Probability = Dict[Parentals][ChildValue] / float(Dict[Parentals][None])
                Score += Dict[Parentals][ChildValue] * math.log(Probability)
        return Score
    def FindInformativeFlanks(self, Feature):
        IntensityByLFlank = {None:[0,0,0,0]}
        LFlankCount = {}
        IntensityByRFlank = {None:[0,0,0,0]}
        RFlankCount = {}
        for Line in self.FeatureLines:
            if Line[0]!=3:
                continue # either b or y (or both) is outside the dynamic range
            #if Line[22]>2:
            #    continue # this is the very first or very last break
            if Line[12]!=2:
                continue # charge isn't 2
            if Line[15] not in (0, 1):
                continue # Weird sector.
            L = Line[16]
            if not IntensityByLFlank.has_key(L):
                IntensityByLFlank[L] = [0,0,0,0]
            IntensityByLFlank[L][Line[Feature]] += 1
            IntensityByLFlank[None][Line[Feature]] += 1
            LFlankCount[L] = LFlankCount.get(L, 0) + 1
            LFlankCount[None] = LFlankCount.get(None, 0) + 1
            R = Line[17]
            if not IntensityByRFlank.has_key(R):
                IntensityByRFlank[R] = [0,0,0,0]
            IntensityByRFlank[R][Line[Feature]] += 1
            IntensityByRFlank[None][Line[Feature]] += 1
            RFlankCount[R] = RFlankCount.get(R, 0) + 1
            RFlankCount[None] = RFlankCount.get(None, 0) + 1
            
        Keys = LFlankCount.keys()
        Keys.sort()
        print 
        print "Feature %s"%Feature
        print "Left flanking AA effects:"
        for Key in Keys:
            Count = LFlankCount[Key]
            print "%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t"%(Key, Count,
                                                      IntensityByLFlank[Key][0]*100/float(Count),
                                                      IntensityByLFlank[Key][1]*100/float(Count),
                                                      IntensityByLFlank[Key][2]*100/float(Count),
                                                      IntensityByLFlank[Key][3]*100/float(Count))
        for Key in Keys:
            Count = LFlankCount[Key]
##            print "0: %s of %s for %s, vs %s of %s overall"%(IntensityByLFlank[Key][0], Count, Char, IntensityByLFlank[None][0], LFlankCount[None])
##            print "1: %s of %s for %s, vs %s of %s overall"%(IntensityByLFlank[Key][1], Count, Char, IntensityByLFlank[None][1], LFlankCount[None])
##            print "2: %s of %s for %s, vs %s of %s overall"%(IntensityByLFlank[Key][2], Count, Char, IntensityByLFlank[None][2], LFlankCount[None])
##            print "3: %s of %s for %s, vs %s of %s overall"%(IntensityByLFlank[Key][3], Count, Char, IntensityByLFlank[None][3], LFlankCount[None])
            print "%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t"%(Key, Count,
                                                      (IntensityByLFlank[Key][0]/float(Count)) * (LFlankCount[None]/float(IntensityByLFlank[None][0])),
                                                      (IntensityByLFlank[Key][1]/float(Count)) * (LFlankCount[None]/float(IntensityByLFlank[None][1])),
                                                      (IntensityByLFlank[Key][2]/float(Count)) * (LFlankCount[None]/float(IntensityByLFlank[None][2])),
                                                      (IntensityByLFlank[Key][3]/float(Count)) * (LFlankCount[None]/float(IntensityByLFlank[None][3])),
                                                        )
        print "Right flanking AA effects:"
        for Key in Keys:
            Count = RFlankCount[Key]
            print "%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t"%(Key, Count,
                                                      IntensityByRFlank[Key][0]*100/float(Count),
                                                      IntensityByRFlank[Key][1]*100/float(Count),
                                                      IntensityByRFlank[Key][2]*100/float(Count),
                                                      IntensityByRFlank[Key][3]*100/float(Count))
        for Key in Keys:
            Count = RFlankCount[Key]
            print "%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t"%(Key, Count,
                                                      (IntensityByRFlank[Key][0]/float(Count)) * (RFlankCount[None]/float(IntensityByRFlank[None][0])),
                                                      (IntensityByRFlank[Key][1]/float(Count)) * (RFlankCount[None]/float(IntensityByRFlank[None][1])),
                                                      (IntensityByRFlank[Key][2]/float(Count)) * (RFlankCount[None]/float(IntensityByRFlank[None][2])),
                                                      (IntensityByRFlank[Key][3]/float(Count)) * (RFlankCount[None]/float(IntensityByRFlank[None][3])),
                                                        )
    def ProducePRMProbabilityTables(self, Charge, OutputFileName):
        FeatureCount = len(self.FeatureLines[0])
        MaxValues = [0]*FeatureCount
        for Line in self.FeatureLines:
            try:
                Feature = int(Line[0])
            except:
                continue # header line
            for Index in range(FeatureCount):
                MaxValues[Index] = max(MaxValues[Index], Line[Index])
        ValueCounts = []
        for Index in range(FeatureCount):
            ValueCounts.append(MaxValues[Index] + 1)
        OutputFile = open(OutputFileName, "wb")
        # Write out the network connectivity.  Hard-coded magicks!
        OutputFile.write(struct.pack("<i", 16)) # network size
        # 0: Dynamic range
        OutputFile.write(struct.pack("<ii64s", 0, ValueCounts[0], "Dynamic range"))
        # 1: y fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[1], "y ion"))
        ParentBlock = struct.pack("<iiiiiiii", 15, -1, -1, -1, ValueCounts[1], 0, 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable1(OutputFile, Charge, ValueCounts, 1, [15], "")
        # 2: b fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[2], "b ion"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, -1, -1, ValueCounts[15]*ValueCounts[2], ValueCounts[2], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 2, [1,15], "")
        # 3: y isotopic fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[3], "y+1 isotopic peak"))
        ParentBlock = struct.pack("<iiiiiiii", 1, -1, -1, -1, ValueCounts[3], 0, 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable1(OutputFile, Charge, ValueCounts, 3, [1,], "")
        # 4: b isotopic fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[4], "b+1 isotopic peak"))
        ParentBlock = struct.pack("<iiiiiiii", 2, -1, -1, -1, ValueCounts[4], 0, 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable1(OutputFile, Charge, ValueCounts, 4, [2,], "")
        # 5: doubly-charged y
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[5], "doubly-charged y"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, -1, -1, ValueCounts[5]*ValueCounts[15], ValueCounts[5], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 5, [1,15], "")
        # 6: doubly-charged b
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[6], "doubly-charged b"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[6]*ValueCounts[15], ValueCounts[6], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 6, [2,15], "")
        # 7: y-H2O
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[7], "y-H2O"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, -1, -1, ValueCounts[7]*ValueCounts[15], ValueCounts[7], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 7, [1,15], "")
        # 8: a
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[8], "a ion"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[8]*ValueCounts[15], ValueCounts[8], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 8, [2,15], "")
        # 9: b-h2o
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[9], "b-H2O"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[9]*ValueCounts[15], ValueCounts[9], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 9, [2,15], "")
        # 10: y-NH3
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[10], "y-NH3"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, -1, -1, ValueCounts[10]*ValueCounts[15], ValueCounts[10], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 10, [1,15], "")
        # 11: b-nh3
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[11], "b-NH3"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[11]*ValueCounts[15], ValueCounts[11], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 9, [2,15], "")
        # 12: charge
        OutputFile.write(struct.pack("<ii64s", 0, 5, "Charge"))
        # 13, 14: FlankB and FlankY, not set.
        OutputFile.write(struct.pack("<ii64s", 0, 0, "Dummy"))
        OutputFile.write(struct.pack("<ii64s", 0, 0, "Dummy"))
        # 15: Sector
        OutputFile.write(struct.pack("<ii64s", 0, 5, "Sector"))
        OutputFile.close()        
    def ProduceProbabilityTables(self, Charge, OutputFileName):
        FeatureCount = len(self.FeatureLines[0])
        MaxValues = [0]*FeatureCount
        for Line in self.FeatureLines:
            for Index in range(FeatureCount):
                MaxValues[Index] = max(MaxValues[Index], Line[Index])
        ValueCounts = []
        for Index in range(FeatureCount):
            ValueCounts.append(MaxValues[Index] + 1)
        OutputFile = open(OutputFileName, "wb")
        # Write out the network connectivity.  Hard-coded magicks!
        OutputFile.write(struct.pack("<i", 16)) # network size
        # Node record:
        #  Flags, ValueCount, Name, parent array, parent block sizes, and prob-table.
        # If Flags = 0, then no parent-array is written.
        #NoParentBlock = struct.pack("<iiiiiiii", -1, -1, -1, -1, 0, 0, 0, 0)
        # 0: Dynamic range
        OutputFile.write(struct.pack("<ii64s", 0, ValueCounts[0], "Dynamic range"))
        # 1: y fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[1], "y ion"))
        ParentBlock = struct.pack("<iiiiiiii", 14, 15, -1, -1, ValueCounts[15] * ValueCounts[1], ValueCounts[1], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 1, [14,15], "y")
        # 2: b fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[2], "b ion"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, 13, -1, ValueCounts[15] * ValueCounts[13] * ValueCounts[2], ValueCounts[13] * ValueCounts[2], ValueCounts[2], 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable3(OutputFile, Charge, ValueCounts, 2, [1,15,13], "b")
        # 3: y isotopic fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[3], "y+1 isotopic peak"))
        ParentBlock = struct.pack("<iiiiiiii", 1, -1, -1, -1, ValueCounts[3], 0, 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable1(OutputFile, Charge, ValueCounts, 3, [1,], "y")
        # 4: b isotopic fragment
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[4], "b+1 isotopic peak"))
        ParentBlock = struct.pack("<iiiiiiii", 2, -1, -1, -1, ValueCounts[4], 0, 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable1(OutputFile, Charge, ValueCounts, 4, [2,], "b")
        # 5: doubly-charged y
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[5], "doubly-charged y"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, -1, -1, ValueCounts[5]*ValueCounts[15], ValueCounts[5], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 5, [1,15], "y")
        # 6: doubly-charged b
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[6], "doubly-charged b"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[6]*ValueCounts[15], ValueCounts[6], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 6, [2,15], "b")
        # 7: y-H2O
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[7], "y-H2O"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, -1, -1, ValueCounts[7]*ValueCounts[15], ValueCounts[7], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 7, [1,15], "y")
        # 8: a
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[8], "a ion"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[8]*ValueCounts[15], ValueCounts[8], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 8, [2,15], "b")
        # 9: b-h2o
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[9], "b-H2O"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[9]*ValueCounts[15], ValueCounts[9], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 9, [2,15], "b")
        # 10: y-NH3
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[10], "y-NH3"))
        ParentBlock = struct.pack("<iiiiiiii", 1, 15, -1, -1, ValueCounts[10]*ValueCounts[15], ValueCounts[10], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 10, [1,15], "y")
        # 11: b-h2o
        OutputFile.write(struct.pack("<ii64s", 3, ValueCounts[11], "b-NH3"))
        ParentBlock = struct.pack("<iiiiiiii", 2, 15, -1, -1, ValueCounts[11]*ValueCounts[15], ValueCounts[11], 0, 0)
        OutputFile.write(ParentBlock)
        self.ProduceProbabilityTable2(OutputFile, Charge, ValueCounts, 9, [2,15], "b")
        # 12: charge
        OutputFile.write(struct.pack("<ii64s", 0, 5, "Charge"))
        # 13: FlankB
        OutputFile.write(struct.pack("<ii64s", 0, 6, "FlankB"))
        # 14: FlankY
        OutputFile.write(struct.pack("<ii64s", 0, 4, "FlankY"))
        # 15: Sector (0..4)
        OutputFile.write(struct.pack("<ii64s", 0, 5, "Sector"))
        OutputFile.close()
    def ProduceProbabilityTable1(self, OutputFile, Charge, ValueCounts, Feature, ParentFeatures, Side):
        "Produce a probability table for one parent"
        ProbTableSize = ValueCounts[ParentFeatures[0]] * ValueCounts[Feature]
        OutputFile.write(struct.pack("<i", ProbTableSize));
        Dict = {}
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for FVal in range(ValueCounts[Feature]):
                Dict[(Val1, FVal)] = 5 # add some 'padding probability', so that nothing has probability zero
                Dict[(Val1,)] = Dict.get((Val1,), 0) + 5
        for Line in self.FeatureLines:
            if Line[12]!=Charge:
                continue
            if Side=="y":
                if Line[0] in (0, 1):
                    continue
            elif Side == "b":
                if Line[0] in (0, 2):
                    continue
            ParentA = Line[ParentFeatures[0]]
            Dict[(ParentA, Line[Feature])] += 1
            if MORE_INTENSE_IS_ALWAYS_BETTER:
                for FeatureX in range(Line[Feature] - 1, 0, -1):
                    Dict[(ParentA, FeatureX)] += 1
            Dict[(ParentA, )] += 1
        print "\nProbability table for %s:"%self.Headers[Feature]
        print "%s\t%s\t"%(self.Headers[ParentFeatures[0]], self.Headers[Feature])
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for FVal in range(ValueCounts[Feature]):
                Probability = Dict[(Val1, FVal)] / float(Dict[(Val1, )])
                OutputFile.write(struct.pack("<f", math.log(Probability)))
                print "%d\t%d\t%d/%d (%.2f%%)\t%.4f"%(Val1, FVal, Dict[(Val1, FVal)], Dict[(Val1,)], 100*Probability, math.log(Probability))
            print        
    def ProduceProbabilityTable2(self, OutputFile, Charge, ValueCounts, Feature, ParentFeatures, Side):
        "Produce a probability table for two parents"
        ProbTableSize = ValueCounts[ParentFeatures[0]] * ValueCounts[ParentFeatures[1]] * ValueCounts[Feature]
        OutputFile.write(struct.pack("<i", ProbTableSize));
        Dict = {}
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for Val2 in range(ValueCounts[ParentFeatures[1]]):
                for FVal in range(ValueCounts[Feature]):
                    Dict[(Val1, Val2, FVal)] = 5 # add some 'padding probability', so that nothing has probability zero
                    Dict[(Val1, Val2)] = Dict.get((Val1, Val2), 0) + 5
        for Line in self.FeatureLines:
            if Line[12]!=Charge:
                continue
            if Side=="y":
                if Line[0] in (0, 1):
                    continue
            elif Side == "b":
                if Line[0] in (0, 2):
                    continue
            ParentA = Line[ParentFeatures[0]]
            ParentB = Line[ParentFeatures[1]]
            Dict[(ParentA,ParentB,Line[Feature])] += 1
            if MORE_INTENSE_IS_ALWAYS_BETTER:
                for FeatureX in range(Line[Feature] - 1, 0, -1):
                    Dict[(ParentA, ParentB, FeatureX)] += 1
            Dict[(ParentA, ParentB)] += 1
        print "\nProbability table for %s:"%self.Headers[Feature]
        print "%s\t%s\t%s\t"%(self.Headers[ParentFeatures[0]], self.Headers[ParentFeatures[1]], self.Headers[Feature])
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for Val2 in range(ValueCounts[ParentFeatures[1]]):
                for FVal in range(ValueCounts[Feature]):
                    Probability = Dict[(Val1, Val2, FVal)] / float(Dict[(Val1, Val2)])
                    OutputFile.write(struct.pack("<f", math.log(Probability)))
                    print "%d\t%d\t%d\t%d/%d (%.2f%%)\t%.4f"%(Val1, Val2, FVal, Dict[(Val1, Val2, FVal)], Dict[(Val1, Val2)], 100*Probability, math.log(Probability))
                print
    def ProduceProbabilityTable3(self, OutputFile, Charge, ValueCounts, Feature, ParentFeatures, Side):
        "Produce a probability table for three parents"
        ProbTableSize = ValueCounts[ParentFeatures[0]] * ValueCounts[ParentFeatures[1]] * ValueCounts[ParentFeatures[2]] * ValueCounts[Feature]
        OutputFile.write(struct.pack("<i", ProbTableSize));
        # PROTECTION FROM UNDERFITTING:
        # Because combinations of three parent features may be rare, first generate a probability table for
        # all combinations of the first two parent features.  We'll use that probability table to fill
        # the 'padding'.  In other words: If you haven't seen parent features (x, y, z) before, go by what you
        # know about conditioning on (x, y).  (Note: We always list the parents in order of informative-ness!)
        SubDict = {}
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for Val2 in range(ValueCounts[ParentFeatures[1]]):
                for FVal in range(ValueCounts[Feature]):
                    SubDict[(Val1, Val2, FVal)] = 1
                    SubDict[(Val1, Val2)] = SubDict.get((Val1, Val2), 0) + 1
        for Line in self.FeatureLines:
            if Line[12]!=Charge:
                continue
            if Side=="y":
                if Line[0] in (0, 1):
                    continue
            elif Side == "b":
                if Line[0] in (0, 2):
                    continue
            ParentA = Line[ParentFeatures[0]]
            ParentB = Line[ParentFeatures[1]]
            SubDict[(ParentA, ParentB, Line[Feature])] += 1
            if MORE_INTENSE_IS_ALWAYS_BETTER:
                for FeatureX in range(Line[Feature]-1, 0, -1):
                    SubDict[(ParentA, ParentB, FeatureX)] += 1
            
            SubDict[(ParentA, ParentB)] += 1
        print "Probability sub-table for %s with first 2 parents:"%self.Headers[Feature]
        print "%s\t%s\t%s\tProb"%(self.Headers[ParentFeatures[0]], self.Headers[ParentFeatures[1]], self.Headers[Feature])
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for Val2 in range(ValueCounts[ParentFeatures[1]]):
                for FVal in range(ValueCounts[Feature]):
                    SubDict[(Val1, Val2, FVal)] /= float(SubDict[(Val1, Val2)])
                    print "%s\t%s\t%s\t%.2f\t"%(Val1, Val2, FVal, 100*SubDict[(Val1, Val2, FVal)])
                print
                
        # Now proceed to the main probability table:
        Dict = {}
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for Val2 in range(ValueCounts[ParentFeatures[1]]):
                for Val3 in range(ValueCounts[ParentFeatures[2]]):
                    for FVal in range(ValueCounts[Feature]):
                        Dict[(Val1, Val2, Val3, FVal)] = 10 * SubDict[(Val1, Val2, FVal)] # add some 'padding probability', so that nothing has probability zero
                    Dict[(Val1, Val2, Val3)] = 10 #Dict.get((Val1, Val2, Val3), 0) + 5
        for Line in self.FeatureLines:
            if Line[12]!=Charge:
                continue
            if Side=="y":
                if Line[0] in (0, 1):
                    continue
            else:
                if Line[0] in (0, 2):
                    continue
            ParentA = Line[ParentFeatures[0]]
            ParentB = Line[ParentFeatures[1]]
            ParentC = Line[ParentFeatures[2]]
            Dict[(ParentA, ParentB, ParentC, Line[Feature])] += 1
            if MORE_INTENSE_IS_ALWAYS_BETTER:
                for FeatureX in range(Lines[Feature]-1, 0, -1):
                    Dict[(ParentA, ParentB, ParentC, FeatureX)] += 1
            Dict[(ParentA, ParentB, ParentC)] += 1
        print "\nProbability table for %s:"%self.Headers[Feature]
        print "%s\t%s\t%s\t%s\t"%(self.Headers[ParentFeatures[0]], self.Headers[ParentFeatures[1]],
                              self.Headers[ParentFeatures[2]], self.Headers[Feature])
        for Val1 in range(ValueCounts[ParentFeatures[0]]):
            for Val2 in range(ValueCounts[ParentFeatures[1]]):
                for Val3 in range(ValueCounts[ParentFeatures[2]]):
                    for FVal in range(ValueCounts[Feature]):
                        Probability = Dict[(Val1, Val2, Val3, FVal)] / float(Dict[(Val1, Val2, Val3)])
                        OutputFile.write(struct.pack("<f", math.log(Probability)))
                        print "%d\t%d\t%d\t%d\t%d/%d (%.2f%%)\t%.4f"%(Val1, Val2, Val3, FVal, Dict[(Val1, Val2, Val3, FVal)], Dict[(Val1, Val2, Val3)], 100*Probability, math.log(Probability))
                    print                
if __name__ == "__main__":
    # Who predicts y well?
    #Trainer = ScorpionTrainer()
    #InputFile = "TrainingFiles\\PEPPRM2.txt"
    #Trainer.ReadFeatures(InputFile)
    #Trainer.FindInformativeFlanks(2)
    #Trainer.FindGoodParent3(7, [1,13,14,15,18,19,20,21,22,23,24,25])
    #Trainer.FindGoodParent2(2, [13,14,15,18,19,20,21])
    ################################
    # Main, for PRM scoring:
    ModelType = sys.argv[1].lower()
    for Charge in (2, 3):
        print "\n\nCharge %s"%Charge
        Trainer = ScorpionTrainer()
        if ModelType == "prm":
            InputFile = "TrainingFiles\\PRM%d.txt"%Charge
            Trainer.ReadFeatures(InputFile)
            OutputFile = "PRM%d.dat"%Charge
            Trainer.ProducePRMProbabilityTables(Charge, OutputFile)
        elif ModelType == "pep":
            InputFile = "TrainingFiles\\PEPPRM%d.txt"%Charge
            Trainer.ReadFeatures(InputFile)
            OutputFile = "Ch%dBNPEP.dat"%Charge
            Trainer.ProduceProbabilityTables(Charge, OutputFile)
        else:
            print "wtf?", ModelType

    