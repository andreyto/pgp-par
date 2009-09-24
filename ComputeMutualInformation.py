"""
Simple script for Bayesian Network topology generation:
Parse a table of features.  Report the mutual information for each pair of features.
"""
import os
import getopt
import sys
import math

PaddingCount = 1

class MutualInformationFinder:
    def __init__(self):
        self.ParseFileName = None
        self.OutputFileName = "MutualInformation.txt"
        self.QuickParseFlag = 0
    def ScoreFeatureRelation(self, FeatureX, FeatureY):
        """
        Compute a score reflecting how much the probability tables p(y) differ from the
        probability tables p(y|x), for each value of x.
        """
        FeatureNameX = self.FeatureNames[FeatureX]
        FeatureNameY = self.FeatureNames[FeatureY]
        # Determine how good feature A is for predicting the value of feature B.
        # Sum, over all values x and y, of p(x) * |ln(p(y|x)) - ln(p(y))|
        XValues = self.OccurrenceTables[FeatureX].keys()
        XValues.sort()
        YValues = self.OccurrenceTables[FeatureY].keys()
        YValues.sort()
        # Compute OddsX:
        OddsX = {}
        OverallCount = 0
        for (Key, Count) in self.OccurrenceTables[FeatureX].items():
            OverallCount += Count #+= Count + PaddingCount
        for (Key, Count) in self.OccurrenceTables[FeatureX].items():
            OddsX[Key] = max(0.0001, Count / float(OverallCount))
        # Compute OddsY:
        OddsY = {}
        OverallCount = 0
        for (Key, Count) in self.OccurrenceTables[FeatureY].items():
            OverallCount += Count #+= Count + PaddingCount
        for (Key, Count) in self.OccurrenceTables[FeatureY].items():
            OddsY[Key] = max(0.0001, Count / float(OverallCount))
        # Iterate over all pairs of feature values:
        FeatureKey = (min(FeatureX, FeatureY), max(FeatureX, FeatureY))
        TotalScore = 0
        for ValueX in XValues:
            OverallCount = 0
            for ValueY in YValues:
                if (FeatureX < FeatureY):
                    ValueKey = (ValueX, ValueY)
                else:
                    ValueKey = (ValueY, ValueX)
                OverallCount += self.OccurrenceTables[FeatureKey].get(ValueKey, 0)
            # OverallCount is the denominator for computing p(y|x)
            for ValueY in YValues:
                if (FeatureX < FeatureY):
                    ValueKey = (ValueX, ValueY)
                else:
                    ValueKey = (ValueY, ValueX)
                Count = self.OccurrenceTables[FeatureKey].get(ValueKey, 0)
                ConditionalOdds = max(0.0001, Count / float(OverallCount))
                Score = abs(math.log(OddsY[ValueY]) - math.log(ConditionalOdds))
                TotalScore += OddsX[ValueX] * Score
        return TotalScore
    def Main(self):
        self.OutputFile = open(self.OutputFileName, "wb")
        Header = "FeatureIndexA\tFeatureIndexB\tFeatureNameA\tFeatureNameB\tEntropyA\tEntropyB\tJointEntropy\tMutualInformation\tSymmetricUncertainty\tEntropyReductionA->B\tEntropyReductionB->A\tParentScore\t"
        self.ParseFeatures()
        self.ComputeEntropies()
        self.OutputFile.write(Header + "\n")
        self.ReportMutualInformation()
        self.ReportVerboseProbabilityTables()
    def ParseFeatures(self):
        # MaxValue[Column] -> Maximum value observed for that column
        # OccurrenceTables[(ColumnA, ColumnB)] -> dictionary of counts
        # OccurrenceTables[ColumnA] -> Dictionary of counts
        self.MaxValue = {}
        self.OccurrenceTables = {}
        if not self.ParseFileName:
            print UsageInfo
            return
        File = open(self.ParseFileName, "rb")
        # The first file line gives the NAMES of the features:
        FileLine = File.readline()
        Bits = FileLine.strip().split("\t")
        self.FeatureCount = len(Bits)
        self.FeatureNames = Bits
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "Line %s..."%LineNumber
            if self.QuickParseFlag and LineNumber > 1000:
                break
            Bits = FileLine.strip().split("\t")
            Values = []
            ValidLine = 1
            for Index in range(self.FeatureCount):
                try:
                    Value = int(Bits[Index])
                except:
                    print Bits
                    ValidLine = 0
                    break
                Values.append(Value)
                self.MaxValue[Index] = max(self.MaxValue.get(Index, 0), Value)
            if not ValidLine:
                continue
            for IndexA in range(len(Bits)):
                ValueA = Values[IndexA]
                if not self.OccurrenceTables.has_key(IndexA):
                    self.OccurrenceTables[IndexA] = {}
                self.OccurrenceTables[IndexA][ValueA] = self.OccurrenceTables[IndexA].get(ValueA, 0) + 1
                for IndexB in range(IndexA + 1, self.FeatureCount):
                    Key = (IndexA, IndexB)
                    if not self.OccurrenceTables.has_key(Key):
                        self.OccurrenceTables[Key] = {}
                    ValueB = Values[IndexB]
                    SubKey = (ValueA, ValueB)
                    self.OccurrenceTables[Key][SubKey] = self.OccurrenceTables[Key].get(SubKey, 0) + 1
        File.close()
    def ComputeEntropies(self):
        # Compute the entropy of each feature:
        self.Entropies = []
        
        for Index in range(self.FeatureCount):
            EntropySum = 0
            OverallCount = 0
            Str = "%s\t%s\t"%(Index, self.FeatureNames[Index])
            for (Key, Count) in self.OccurrenceTables[Index].items():
                OverallCount += Count #+= Count + PaddingCount
            for (Key, Count) in self.OccurrenceTables[Index].items():
                if Count == 0:
                    Odds = 0.0001
                else:
                    Odds = Count / float(OverallCount)
                Str += "%s\t"%Odds
                #Odds = (Count + PaddingCount) / float(OverallCount)
                EntropySum -= Odds * math.log(Odds)
            self.Entropies.append(EntropySum)
            self.OutputFile.write(Str + "\n")
        self.OutputFile.write("\n")
    def ReportMutualInformation(self):
        # For each pair of features, compute the mutual information:
        for IndexA in range(self.FeatureCount):
            for IndexB in range(self.FeatureCount):
                if IndexB == IndexA:
                    continue
                EntropySum = 0
                OverallCount = 0
                if IndexA < IndexB:
                    Key = (IndexA, IndexB)
                else:
                    Key = (IndexB, IndexA)
                for (SubKey, Count) in self.OccurrenceTables[Key].items():
                    OverallCount += Count  # += (Count + PaddingCount)
                for (SubKey, Count) in self.OccurrenceTables[Key].items():
                    #Odds = (Count + PaddingCount) / float(OverallCount)
                    if Count == 0:
                        Odds = 0.0001
                    else:
                        Odds = Count / float(OverallCount)
                    EntropySum -= Odds * math.log(Odds)
                MutualInformation = self.Entropies[IndexA] + self.Entropies[IndexB] - EntropySum
                SymmetricUncertainty = 2 * MutualInformation / (self.Entropies[IndexA] + self.Entropies[IndexB])
                Str = "%s\t%s\t%s\t%s\t"%(IndexA, IndexB, self.FeatureNames[IndexA], self.FeatureNames[IndexB])
                Str += "%s\t%s\t"%(self.Entropies[IndexA], self.Entropies[IndexB])
                Str += "%s\t%s\t%s\t"%(EntropySum, MutualInformation, SymmetricUncertainty)
                CBA = MutualInformation / self.Entropies[IndexA]
                CAB = MutualInformation / self.Entropies[IndexB]
                Str += "%s\t%s\t"%(CAB, CBA)
                Score = self.ScoreFeatureRelation(IndexB, IndexA)
                Str += "%s\t"%Score
                self.OutputFile.write(Str + "\n")
    def ReportVerboseProbabilityTable(self, FeatureIndex):
        FeatureName = self.FeatureNames[FeatureIndex]
        FeatureValues = self.OccurrenceTables[FeatureIndex].keys()
        FeatureValues.sort()
        File = open("ProbabilityTables.%s.txt"%FeatureName, "wb")
        # Write a header:
        Str = "Parent\tValue\tParentOdds"
        for Value in FeatureValues:
            Str += "%s\t"%Value
        Str += "\t"
        for Value in FeatureValues:
            Str += "Diff%s\t"%Value
        Str += "\t"
        for Value in FeatureValues:
            Str += "AbsDiff%s\t"%Value
        File.write(Str + "\n")
        # Compute BaseOdds:
        BaseOdds = {}
        OverallCount = 0
        for (Key, Count) in self.OccurrenceTables[FeatureIndex].items():
            OverallCount += Count #+= Count + PaddingCount
        for (Key, Count) in self.OccurrenceTables[FeatureIndex].items():
            BaseOdds[Key] = max(0.0001, Count / float(OverallCount))
        Str = "Base\t\t\t"
        for Value in FeatureValues:
            Str += "%s\t"%BaseOdds[Value]
        File.write(Str + "\n")
        for OtherFeature in range(self.FeatureCount):
            if OtherFeature == FeatureIndex:
                continue
            # Compute odds for the other feature:
            OtherFeatureOdds = {}
            OverallCount = 0
            for (Key, Count) in self.OccurrenceTables[OtherFeature].items():
                OverallCount += Count
            for (Key, Count) in self.OccurrenceTables[OtherFeature].items():
                OtherFeatureOdds[Key] = max(0.0001, Count / float(OverallCount))
            #########################################
            FeatureKey = (min(FeatureIndex, OtherFeature), max(FeatureIndex, OtherFeature))
            OtherFeatureValues = self.OccurrenceTables[OtherFeature].keys()
            OtherFeatureValues.sort()
            for OtherValue in OtherFeatureValues:
                Str = "%s\t%s\t%s\t"%(self.FeatureNames[OtherFeature], OtherValue, OtherFeatureOdds[OtherValue])
                # Get the overall count:
                OverallCount = 0
                for Value in FeatureValues:
                    if (FeatureIndex < OtherFeature):
                        ValueKey = (Value, OtherValue)
                    else:
                        ValueKey = (OtherValue, Value)
                    OverallCount += self.OccurrenceTables[FeatureKey].get(ValueKey, 0)
                # Report the row:
                Diffs = ""
                AbsDiffs = ""
                for Value in FeatureValues:
                    if (FeatureIndex < OtherFeature):
                        ValueKey = (Value, OtherValue)
                    else:
                        ValueKey = (OtherValue, Value)
                    Count = self.OccurrenceTables[FeatureKey].get(ValueKey, 0)
                    Odds = Count / float(OverallCount)
                    Str += "%s\t"%Odds
                    Diff = Odds - BaseOdds[Value]
                    Diffs += "%s\t"%(Diff)
                    AbsDiffs += "%s\t"%(abs(Diff))
                Str += "\t"
                Str += Diffs
                Str += "\t"
                Str += AbsDiffs
                
                File.write(Str + "\n")
                
    def ReportVerboseProbabilityTables(self):
        for FeatureIndex in range(self.FeatureCount):
            FeatureName = self.FeatureNames[FeatureIndex]
            if FeatureName in ("Y", "B", "Y-H2O", "Y-NH3", "B-H2O", "B-NH3", "A"):
                self.ReportVerboseProbabilityTable(FeatureIndex)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:Q")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.ParseFileName = Value
            elif Option == "-w":
                self.OutputFileName = Value
            elif Option == "-Q":
                self.QuickParseFlag = 1
        

UsageInfo = """
ComputeMutualInformation options:
 -r [FileName]: Parse table (written out by DriveScorpion.py)
 -w [OutputFile]: Write output to the specified file
 -Q: Parse only the first few lines (for debugging only)
"""
        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(no psyco)"
    EntropyReporter = MutualInformationFinder()
    EntropyReporter.ParseCommandLine(sys.argv[1:])
    EntropyReporter.Main()