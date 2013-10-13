#!/usr/bin/env python
"""
This script, based on PeptideProphet, computes the estimated probability that
a match is correct.  This probability is derived upon an F-Score.  The F-Score
for a match is a weighted sum of the length-corrected MQScore and the delta score.

We fit the distribution of F-scores as a mixture of two distributions:
A GAMMA DISTRIBUTION for false matches (with lower mean)
A NORMAL DISTRIBUTION for true matches (with higher mean)
Therefore, the probability that a match with a given F-Score is correct depends
upon the overall distribution of F-Scores for the rest of the run.

============================================================================
P-values for searches with shuffled-database:

Given a score cutoff, let the number of valid-protein hits above the
cutoff be V, and let the number of invalid-protein hits above the
cutoff be I.

# PVALUE WITH REMOVAL:
Throw out the I hits from invalid proteins.
Let TDB and FDB be the true and false database fractions.  (Note: For
a 1:1 database, TDB and FDB are both 0.5, so TDB/FDB equals 1.0)  Even
after throwing out hits from invalid proteins, there are still some
chance hits to true proteins.  The estimated number of hits from V
that are actually false is equal to: I*(TDB/FDB)

This, the odds that a match above this score cutoff is correct: 
(V - I*(TDB/FDB)) / V

This formulation of pvalue is the one normally used (e.g. in an
unmodified search).  Normally, we have no reason to keep any matches
to shuffled proteins; we only generate them in the first place so that
we can count them.

# PVALUE WITHOUT REMOVAL (-H command-line option):
Retain all hits.  As above, the number of
hits from V that are false is I*(TDB/FDB).  Thus, the odds that a
match above the score cutoff is correct:
(V - I*(TDB/FDB)) / (I+V)

"""
import os
import sys
import math
import getopt
import traceback
import InspectResults
import SelectProteins
import Learning

from Utils import *

Initialize()
try:
    from PIL import Image
    from PIL import ImageDraw
    from PIL import ImageFont
    # Fonts don't seem to work on Linux.  (Tried pdf, pcf, and pil formats...but no luck)
    # So, we'll content ourselves with a default font if we must:
    try:
        TheFont = ImageFont.truetype("Times.ttf", 12)
    except:
        TheFont = ImageFont.load_default()
except:
    print "(PIL not installed - image generation not available)"
    Image = None

class Colors:
    White = (255, 255, 255)
    Grey = (155, 155, 155)
    Background = (255, 255, 255)
    Black = (0, 0, 0)
    Green = (0, 155, 0)
    Red = (155, 0, 0)
    Blue = (0, 0, 155)

class Defaults:
    "Default F-score distribution; a starting point for E/M model fitting."
    MeanTrue = 4.48
    VarianceTrue = 1.50
    MeanFalse = 0.19
    VarianceFalse = 0.17
    PriorProbabilityTrue = 0.25
    GammaOffset = 0.3
    MQScoreWeight = 0.3
    DeltaScoreWeight = 1.5
    ###########################
    BlindMeanTrue = 5.0
    BlindVarianceTrue = 11.8
    BlindMeanFalse = -0.8
    BlindVarianceFalse = 0.7
    BlindPriorProbabilityTrue = 0.18
    BlindGammaOffset = 6.0
    BlindMQScoreWeight = 0.3
    BlindDeltaScoreWeight = 1.5

BLIND_MOD_PENALTY = 1.0
MIN_MQSCORE = -10.0

# Parse the scores from at most this many output files.  
MAX_RESULTS_FILES_TO_PARSE = 2000

BIN_MULTIPLIER = 10.0
SQRT2PI = math.sqrt(2 * math.pi)

Cof = [76.18009172947146, -86.50532032941677,
    24.01409824083091, -1.231739572450155, 
    0.1208650973866179e-2, -0.5395239384952e-5]

def Gamma(Z):
    X = Z
    Y = Z
    Temp = X + 5.5
    Temp -= (X + 0.5) * math.log(Temp)
    Ser = 1.000000000190015
    for J in range(6):
        Y += 1
        Ser += Cof[J] / Y
    Z = -Temp + math.log(2.5066282746310005 * Ser / X)
    return math.exp(Z)

class Bag:
    pass

class PValueParser():
    def __init__(self):
        self.RetainBadMatches = 0
        self.LoadDistributionPath = None
        self.ScoreHistogram2 = {}
        self.ScoreHistogram3 = {}
        self.ShuffledScoreHistogram2 = {}
        self.ShuffledScoreHistogram3 = {}
        self.MinimumPeptideLength = 7
        self.VerboseFlag = 0
        self.GenerateImageFlag = 0
        self.MQScoreWeight = Defaults.MQScoreWeight
        self.DeltaScoreWeight = Defaults.DeltaScoreWeight
        self.GammaOffset = Defaults.GammaOffset
        self.BlindFlag = 0
        self.PValueCutoff = 0.1 # default
        # aminos -> location list
        self.PeptideDict = {}
        self.MaxDeltaScoreGap = -3.5
        self.DBPath = []
        self.PerformProteinSelection = 0
        self.ProteinPicker = None
        self.WriteTopMatchOnly = 0
        self.ShuffledDatabaseFraction = None
        self.RemoveShuffledMatches = 1
        # Overwrite existing files in -w target:
        self.OverwriteNewScoresFlag = 1
        self.ClusterInfoPath = None

    def ReadDeltaScoreDistribution(self, FilePath):
        """
        Read delta-scores from a file, to compute the average delta-score.
        If passed a directory, iterate over all results files in the directory.
        """
        #
        self.AllSpectrumCount2 = 0
        self.AllSpectrumCount3 = 0
        self.MeanDeltaScore2 = 0
        self.MeanDeltaScore3 = 0
        self.ReadDeltaScoreDistributionFromFile(FilePath)
        self.MeanDeltaScore2 /= max(1, self.AllSpectrumCount2)
        self.MeanDeltaScore3 /= max(1, self.AllSpectrumCount3)
        if self.VerboseFlag:
            print "Mean delta score ch1..2: %s over %s spectra"%(self.MeanDeltaScore2, self.AllSpectrumCount2)
            print "Mean delta score ch3: %s over %s spectra"%(self.MeanDeltaScore3, self.AllSpectrumCount3)
        if not self.MeanDeltaScore2:
            self.MeanDeltaScore2 = 0.001
        if not self.MeanDeltaScore3:
            self.MeanDeltaScore3 = 0.001
    def ReadDeltaScoreDistributionFromFile(self, FilePath):
        "Read delta-scores from a single file, to compute the average delta-score."
        print "Read delta-score distribution from %s..."%FilePath
        try:
            FileParser = InspectResults.Parser( FilePath, MaxFilesToParse=MAX_RESULTS_FILES_TO_PARSE )
        except:
            traceback.print_exc()
            return
        OldSpectrum = None
        for row in FileParser:
            try:
                Charge = row.Charge
                MQScore = row.MQScore
                DeltaScore = row.DeltaScore
                Peptide = GetPeptideFromModdedName(row.Annotation[2:-2])
                Spectrum = (row.SpectrumFile, row.ScanNumber)
            except:
                traceback.print_exc()
                print row
                continue # header line
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            Length = len(Peptide.Aminos)
            if Length < self.MinimumPeptideLength:
                continue
            if DeltaScore < 0:
                print "## Warning: DeltaScore < 0!", Spectrum, FilePath
                print row
                continue
            if Charge < 3:
                self.AllSpectrumCount2 += 1
                self.MeanDeltaScore2 += DeltaScore
                
            else:
                self.AllSpectrumCount3 += 1
                self.MeanDeltaScore3 += DeltaScore

    def ReadScoreDistributionFromFile(self, FilePath):
        """
        Read F-scores from a single file, to compute the score histogram.
        """
        print "Read score distribution from %s..."%FilePath
        try:
            resultsParser = InspectResults.Parser(FilePath, MaxFilesToParse=MAX_RESULTS_FILES_TO_PARSE)
        except:
            traceback.print_exc()
            return
        OldSpectrum = None
        for row in resultsParser:
            try:
                Charge = row.Charge
                MQScore = row.MQScore
                DeltaScore = row.DeltaScore
                Peptide = GetPeptideFromModdedName(row.Annotation[2:-2])
                Protein = row.ProteinName
                Spectrum = (row.SpectrumFile, row.ScanNumber)
            except:
                continue # header line
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            Length = len(Peptide.Aminos)
            if Length < self.MinimumPeptideLength:
                continue
            if (Charge < 3):
                MeanDeltaScore = self.MeanDeltaScore2
            else:
                MeanDeltaScore = self.MeanDeltaScore3
            WeightedScore = self.MQScoreWeight * MQScore + self.DeltaScoreWeight * (DeltaScore / MeanDeltaScore)
            ScoreBin = int(round(WeightedScore * BIN_MULTIPLIER))
            Hit = 1
            if self.ClusterInfoPath:
                # Get this cluster's size:
                ClusterFileName = row.SpectrumFile.replace("/","\\").split("\\")[-1]
                ScanNumber = row.ScanNumber
                ClusterSize = self.ClusterSizes.get((ClusterFileName, ScanNumber), None)
                if not ClusterSize:
                    print "* Warning: ClusterSize not known for %s, %s"%(ClusterFileName, ScanNumber)
                else:
                    Hit = ClusterSize
            if Charge < 3:
                self.ScoreHistogram2[ScoreBin] = self.ScoreHistogram2.get(ScoreBin, 0) + Hit
            else:
                self.ScoreHistogram3[ScoreBin] = self.ScoreHistogram3.get(ScoreBin, 0) + Hit
            if self.ShuffledDatabaseFraction:
                if Protein[:3] == "XXX":
                    if Charge < 3:
                        self.ShuffledScoreHistogram2[ScoreBin] = self.ShuffledScoreHistogram2.get(ScoreBin, 0) + Hit
                    else:
                        self.ShuffledScoreHistogram3[ScoreBin] = self.ShuffledScoreHistogram3.get(ScoreBin, 0) + Hit

    def ProduceScoreDistributionImage(self, ImagePath, Charge3Flag = 0):
        """
        Write out, to the specified path, an image with f-score on the X-axis
        and p-value on the Y-axis.  If we fit a mixture model, plot the true
        and false (model) distributions; if we fit using shuffled proteins,
        plot the empirical distributions.
        """
        if Image == None:
            return
        
        if Charge3Flag:
            ScoreHistogram = self.ScoreHistogram3
            ShuffledScoreHistogram = self.ShuffledScoreHistogram3
        else:
            ScoreHistogram = self.ScoreHistogram2
            ShuffledScoreHistogram = self.ShuffledScoreHistogram2
        # Image size:
        self.Width = 900
        self.Height = 500
        self.LeftPadding = 50
        self.RightPadding = 80
        self.BottomPadding = 40
        self.TopPadding = 10
        self.PlotWidth = self.Width - (self.LeftPadding + self.RightPadding)
        self.PlotHeight = self.Height - (self.TopPadding + self.BottomPadding)
        self.PlotImage = Image.new("RGB", (self.Width, self.Height), Colors.Background)
        self.Draw = ImageDraw.Draw(self.PlotImage)
        # Find largest and smallest bins and entries:
        self.MaxScoreHistogramEntry = 10
        self.MaxScoreHistogramEntryValid = 10
        self.MaxScoreHistogramEntryInvalid = 10
        self.MaxScoreBin = -9999
        self.MinScoreBin = 9999
        self.TotalHistogramEntries = 0
        for (Bin, Entry) in ScoreHistogram.items():
            self.MaxScoreHistogramEntry = max(Entry, self.MaxScoreHistogramEntry)
            InvalidCount = ShuffledScoreHistogram.get(Bin, 0)
            ValidCount = max(0, Entry - InvalidCount)
            self.MaxScoreHistogramEntryValid = max(self.MaxScoreHistogramEntryValid, ValidCount)
            self.MaxScoreHistogramEntryInvalid = max(self.MaxScoreHistogramEntryInvalid, InvalidCount)
            self.MaxScoreBin = max(self.MaxScoreBin, Bin)
            self.MinScoreBin = min(self.MinScoreBin, Bin)
            self.TotalHistogramEntries += Entry
            #print "Bin %s: Valid %s, invalid %s"%(Bin, ValidCount, InvalidCount)
        self.BinCount = self.MaxScoreBin - self.MinScoreBin + 1
        # Draw the Y axis:
        self.Draw.line((self.LeftPadding, self.TopPadding, self.LeftPadding, self.TopPadding + self.PlotHeight), Colors.Black)
        self.Draw.line((self.Width - self.RightPadding, self.TopPadding, self.Width - self.RightPadding, self.TopPadding + self.PlotHeight), Colors.Black)
        Fraction = 0
        while Fraction <= 1.0:
            Y = self.TopPadding + self.PlotHeight * (1.0 - Fraction)
            Label = str(int(round(Fraction * self.MaxScoreHistogramEntry)))
            self.Draw.text((self.LeftPadding - 5 - len(Label)*5, Y - 6), Label, Colors.Black)
            self.Draw.line((self.LeftPadding - 5, Y, self.LeftPadding, Y), Colors.Black)
            Label = str(Fraction)
            self.Draw.text((self.Width - self.RightPadding + 10, Y - 6), Label, Colors.Black)
            self.Draw.line((self.Width - self.RightPadding, Y, self.Width - self.RightPadding + 5, Y), Colors.Black)
            Fraction += 0.1
        # Draw the X axis:
        self.Draw.line((self.LeftPadding, self.Height - self.BottomPadding, self.Width - self.RightPadding, self.Height - self.BottomPadding), Colors.Black)
        Bin = self.MinScoreBin
        while Bin % 10 != 0:
            Bin += 1
        while Bin < self.MaxScoreBin:
            BinNumber = Bin - self.MinScoreBin
            X = self.LeftPadding + BinNumber * self.PlotWidth / float(self.BinCount)            
            self.Draw.line((X, self.Height - self.BottomPadding - 2, X, self.Height - self.BottomPadding + 2), Colors.Black)
            Label = "%.1f"%(Bin / BIN_MULTIPLIER)
            self.Draw.text((X - len(Label) * 2.5, self.Height - self.BottomPadding + 2), Label, Colors.Black)
            Bin += 10
        if self.ShuffledDatabaseFraction != None:
            self.ProduceImageShuffledDB(Charge3Flag)
        else:
            self.ProduceImageMixtureModel(Charge3Flag)
        self.PlotImage.save(ImagePath)
        # Free:
        self.PlotImage = None
        self.Draw = None
    def ProduceImageShuffledDB(self, Charge3Flag = 0):
        if not Image:
            return

        if Charge3Flag:
            ScoreHistogram = self.ScoreHistogram3
            ShuffledScoreHistogram = self.ShuffledScoreHistogram3
            OddsTrue = self.OddsTrue3
        else:
            ScoreHistogram = self.ScoreHistogram2
            ShuffledScoreHistogram = self.ShuffledScoreHistogram2
            OddsTrue = self.OddsTrue2
        # Draw the legend:
        Y = self.Height - self.BottomPadding + 20
        self.Draw.line((105, Y, 125, Y), Colors.Black)
        self.Draw.rectangle((113, Y-2, 118, Y+2), Colors.Black)
        self.Draw.text((130, Y - 5), "p-value", Colors.Black)
        Y = self.Height - self.BottomPadding + 30
        self.Draw.line((105, Y, 125, Y), Colors.Blue)
        self.Draw.rectangle((113, Y-2, 118, Y+2), Colors.Grey)
        self.Draw.text((130, Y - 5), "All hits", Colors.Grey)
        Y = self.Height - self.BottomPadding + 20
        self.Draw.line((305, Y, 325, Y), Colors.Red)
        self.Draw.rectangle((313, Y-2, 318, Y+2), Colors.Green)
        self.Draw.text((330, Y - 5), "Valid proteins", Colors.Green)
        Y = self.Height - self.BottomPadding + 30
        self.Draw.line((305, Y, 325, Y), Colors.Red)
        self.Draw.rectangle((313, Y-2, 318, Y+2), Colors.Red)
        self.Draw.text((330, Y - 5), "Invalid proteins", Colors.Red)
        # Loop over bins, plotting distributions:
        PrevYOdds = None
        PrevYAll = None
        PrevYTrue = None
        PrevYFalse = None
        PrevX = None
        for Bin in range(self.MinScoreBin, self.MaxScoreBin + 1):
            BinNumber = Bin - self.MinScoreBin
            XX = self.LeftPadding + (BinNumber * self.PlotWidth / float(self.BinCount))
            # p-value:
            PValue = 1.0 - OddsTrue[Bin]
            YOdds = self.Height - self.BottomPadding - self.PlotHeight * PValue
            self.Draw.rectangle((XX - 2, YOdds - 2, XX + 2, YOdds + 2), Colors.Black)
            if PrevYOdds != None:
                self.Draw.line((PrevX, PrevYOdds, XX, YOdds), Colors.Black)
            # Overall:
            Count = ScoreHistogram.get(Bin, 0)
            YAll = self.Height - self.BottomPadding - self.PlotHeight * Count / float(self.MaxScoreHistogramEntry)
            self.Draw.rectangle((XX - 2, YAll - 2, XX + 2, YAll + 2), Colors.Grey)
            if (PrevYAll):
                self.Draw.line((PrevX, PrevYAll, XX, YAll), Colors.Grey)
            # Invalid:
            CountInvalid = ShuffledScoreHistogram.get(Bin, 0)
            YFalse = self.Height - self.BottomPadding - self.PlotHeight * CountInvalid / float(self.MaxScoreHistogramEntryInvalid)
            self.Draw.rectangle((XX - 2, YFalse - 2, XX + 2, YFalse + 2), Colors.Red)
            if (PrevYFalse):
                self.Draw.line((PrevX, PrevYFalse, XX, YFalse), Colors.Red)
            # Valid:
            CountValid = Count - CountInvalid
            YTrue = self.Height - self.BottomPadding - self.PlotHeight * CountValid / float(self.MaxScoreHistogramEntryValid)
            self.Draw.rectangle((XX - 2, YTrue - 2, XX + 2, YTrue + 2), Colors.Green)
            #print "Bin %s: Valid %s/%s invalid %s/%s"%(Bin, CountValid, self.MaxScoreHistogramEntryValid, CountInvalid, self.MaxScoreHistogramEntryInvalid)
            if (PrevYTrue):
                self.Draw.line((PrevX, PrevYTrue, XX, YTrue), Colors.Green)
            # Remember these values, for linking to the next in the series:
            PrevX = XX
            PrevYOdds = YOdds
            PrevYAll = YAll
            PrevYFalse = YFalse
            PrevYTrue = YTrue
    def ProduceImageMixtureModel(self, Charge3Flag = 0):
        """
        Helper for ProduceScoreDistributionImage, if we're using a mixture
        model (not a shuffled database)
        """
        if Image == None:
            return
        if Charge3Flag:
            ScoreHistogram = self.ScoreHistogram3
            ShuffledScoreHistogram = self.ShuffledScoreHistogram3
            MixtureModel = self.MixtureModel3
            OddsTrue = self.OddsTrue3
        else:
            ScoreHistogram = self.ScoreHistogram2
            ShuffledScoreHistogram = self.ShuffledScoreHistogram2
            MixtureModel = self.MixtureModel2
            OddsTrue = self.OddsTrue2
        # Draw the legend:
        Y = self.Height - self.BottomPadding + 20
        self.Draw.line((105, Y, 125, Y), Colors.Black)
        self.Draw.rectangle((113, Y-2, 118, Y+2), Colors.Black)
        self.Draw.text((130, Y - 5), "Empirical score distribution", Colors.Black)
        Y = self.Height - self.BottomPadding + 30
        self.Draw.line((105, Y, 125, Y), Colors.Blue)
        self.Draw.rectangle((113, Y-2, 118, Y+2), Colors.Blue)
        self.Draw.text((130, Y - 5), "Probability true (1-pvalue)", Colors.Blue)
        Y = self.Height - self.BottomPadding + 20
        self.Draw.line((305, Y, 325, Y), Colors.Red)
        self.Draw.rectangle((313, Y-2, 318, Y+2), Colors.Red)
        self.Draw.text((330, Y - 5), "Gamma dist. (fit to false matches)", Colors.Red)
        Y = self.Height - self.BottomPadding + 30
        self.Draw.line((305, Y, 325, Y), Colors.Green)
        self.Draw.rectangle((313, Y-2, 318, Y+2), Colors.Green)
        self.Draw.text((330, Y - 5), "Normal dist. (fit to true matches)", Colors.Green)
        Y = self.Height - self.BottomPadding + 20
        self.Draw.line((555, Y, 575, Y), Colors.Grey)
        self.Draw.rectangle((563, Y-2, 568, Y+2), Colors.Grey)
        self.Draw.text((580, Y - 5), "Fitted mixture model", Colors.Grey)
        # Draw the plot of OBSERVED SCORES:
        PrevX = None
        PrevY = None
        for Bin in range(self.MinScoreBin, self.MaxScoreBin + 1):
            BinNumber = Bin - self.MinScoreBin
            X = self.LeftPadding + BinNumber * self.PlotWidth / float(self.BinCount)
            Count = ScoreHistogram.get(Bin, 0)
            Y = self.Height - self.BottomPadding - self.PlotHeight * Count / float(self.MaxScoreHistogramEntry)
            self.Draw.rectangle((X - 2, Y - 2, X + 2, Y + 2), Colors.Black)
            if PrevX != None:
                self.Draw.line((PrevX, PrevY, X, Y), Colors.Black)
            PrevX = X
            PrevY = Y
        #######################################################
        # Find the scaling factor for the MERGED distribution:
        ComboDistTotal = 0
        for Bin in range(self.MinScoreBin, self.MaxScoreBin + 1):
            TrueScore = Bin / BIN_MULTIPLIER
            Pow = - ((TrueScore - MixtureModel.MeanTrue)**2) / (2 * MixtureModel.VarianceTrue)
            TrueNormal = math.exp(Pow) / (MixtureModel.StdDevTrue * SQRT2PI)
            GX = max(0.01, TrueScore + MixtureModel.GammaOffset)
            FalseGamma = math.pow(GX, MixtureModel.KFalse - 1) * math.exp(-GX / MixtureModel.ThetaFalse) / MixtureModel.GammaDemonFalse
            ComboDist = TrueNormal * MixtureModel.PriorProbabilityTrue + (1.0 - MixtureModel.PriorProbabilityTrue) * FalseGamma
            ComboDistTotal += ComboDist
        YFittedScalingFactor = self.TotalHistogramEntries / ComboDistTotal
        #######################################################
        # Draw the plot of the FALSE HIT GAMMA and TRUE HIT NORMAL and MERGED distributions:
        PrevX = None
        PrevYNormal = None
        PrevYGamma = None
        PrevYOdds = None
        PrevYFitted = None
        for Bin in range(self.MinScoreBin, self.MaxScoreBin + 1):
            BinNumber = Bin - self.MinScoreBin
            XX = self.LeftPadding + (BinNumber * self.PlotWidth / float(self.BinCount))
            TrueScore = Bin / BIN_MULTIPLIER
            Pow = - ((TrueScore - MixtureModel.MeanTrue)**2) / (2 * MixtureModel.VarianceTrue)
            TrueNormal = math.exp(Pow) / (MixtureModel.StdDevTrue * SQRT2PI)
            GX = max(0.01, TrueScore + MixtureModel.GammaOffset)
            FalseGamma = math.pow(GX, MixtureModel.KFalse - 1) * math.exp(-GX / MixtureModel.ThetaFalse) / MixtureModel.GammaDemonFalse
            YNormal = self.Height - self.BottomPadding - self.PlotHeight * TrueNormal
            # Normal distribution:
            self.Draw.rectangle((XX - 2, YNormal - 2, XX + 2, YNormal + 2), Colors.Green)
            if PrevX != None:
                self.Draw.line((PrevX, PrevYNormal, XX, YNormal), Colors.Green)
            # Gamma distribution:
            YGamma = self.Height - self.BottomPadding - self.PlotHeight * FalseGamma
            self.Draw.rectangle((XX - 2, YGamma - 2, XX + 2, YGamma + 2), Colors.Red)
            if PrevX != None:
                self.Draw.line((PrevX, PrevYGamma, XX, YGamma), Colors.Red)
            # Fitted curve:
            ComboDist = TrueNormal * MixtureModel.PriorProbabilityTrue + (1.0 - MixtureModel.PriorProbabilityTrue) * FalseGamma
            YFitted = ComboDist * YFittedScalingFactor / self.MaxScoreHistogramEntry
            YFitted = self.Height - self.BottomPadding - YFitted * self.PlotHeight
            #print TrueNormal, FalseGamma, self.AllSpectrumCount, ComboDist, YFitted
            self.Draw.rectangle((XX - 2, YFitted - 2, XX + 2, YFitted + 2), Colors.Grey)
            if PrevX != None:
                self.Draw.line((PrevX, PrevYFitted, XX, YFitted), Colors.Grey)
            # P-Value:
            PValue = 1.0 - OddsTrue.get(Bin, 0)
            YOdds = self.Height - self.BottomPadding - self.PlotHeight * PValue
            self.Draw.rectangle((XX - 2, YOdds - 2, XX + 2, YOdds + 2), Colors.Blue)
            if PrevX != None:
                self.Draw.line((PrevX, PrevYOdds, XX, YOdds), Colors.Blue)
            # Remember these points' coords for drawing lines next time:
            PrevX = XX
            PrevYNormal = YNormal
            PrevYGamma = YGamma
            PrevYOdds = YOdds
            PrevYFitted = YFitted
    def FitMixtureModel(self):
        self.MixtureModel2 = Learning.MixtureModelClass()
        self.MixtureModel2.Model(None, self.ScoreHistogram2)
        self.OddsTrue2 = self.MixtureModel2.OddsTrue
        self.MixtureModel3 = Learning.MixtureModelClass()
        self.MixtureModel3.Model(None, self.ScoreHistogram3)
        self.OddsTrue3 = self.MixtureModel3.OddsTrue        
        return 1
    def SavePValueDistribution(self, Charge3Flag = 0):
        """
        Write out a p-value distribution derived from forward+shuffled database
        """
        if Charge3Flag:
            OddsTrue = self.OddsTrue3
            MeanDeltaScore = self.MeanDeltaScore3
            ScoreHistogram = self.ScoreHistogram3
            ShuffledScoreHistogram = self.ShuffledScoreHistogram3
        else:
            OddsTrue = self.OddsTrue2
            MeanDeltaScore = self.MeanDeltaScore2
            ScoreHistogram = self.ScoreHistogram2
            ShuffledScoreHistogram = self.ShuffledScoreHistogram2
        Keys = OddsTrue.keys()
        if not Keys:
            return
        MinBin = min(Keys)
        MaxBin = max(Keys)
        self.OutputDistributionFile.write("#MeanDeltaScore\t%s\n"%MeanDeltaScore)
        self.OutputDistributionFile.write("#BlindFlag\t%s\n"%self.BlindFlag)
        if self.ShuffledDatabaseFraction != None:
            Header = "#Bin\tFDR\tTotalHits\tHitsValid\tHitsInvalid\tPeptideFDR\tPeptidesValid\tPeptidesInvalid\tProteinFDR\tProteinsValid\tProteinsInvalid\t\n"
        else:
            Header = "#Bin\tFDR\tTotalHits\t\n"
        self.OutputDistributionFile.write(Header)
        if self.ShuffledDatabaseFraction != None:
            # Count the total number of true hits, false hits, true peptides, false peptides...
            CumulativeTrueHits = 0
            CumulativeFalseHits = 0
            for Bin in range(MinBin, MaxBin + 1):
                AllHits = ScoreHistogram.get(Bin, 0)
                FalseHits = ShuffledScoreHistogram.get(Bin, 0)
                CumulativeFalseHits += FalseHits
                CumulativeTrueHits += (AllHits - FalseHits)
            if self.ProteinPicker:
                ######################################################
                # Peptides:
                ValidPeptides = {}
                InvalidPeptides = {}
                CumulativeTruePeptides = 0
                CumulativeFalsePeptides = 0
                BestScoreByProtein = {}
                for (Peptide, Score) in self.ProteinPicker.BestScoresByPeptide.items():
                    Bin = int(round(Score / BIN_MULTIPLIER))
                    ProteinID = self.ProteinPicker.PeptideProteins.get(Peptide, None)
                    if not ProteinID:
                        print "*** Warning: Peptide '%s' was never assigned to a protein!"%Peptide
                        LocationList = self.ProteinPicker.PeptideDict[Peptide]
                        print LocationList
                        for (ProteinID, Pos) in LocationList:
                            print ProteinID, self.ProteinPicker.ProteinNames[ProteinID], self.ProteinPicker.ProteinPeptideCounts[ProteinID]
                        continue # shouldn't occur!
                    ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                    OldScore = BestScoreByProtein.get(ProteinID, -9999)
                    BestScoreByProtein[ProteinID] = max(OldScore, Score)
                    if ProteinName[:3] == "XXX":
                        InvalidPeptides[Bin] = InvalidPeptides.get(Bin, 0) + 1
                        CumulativeFalsePeptides += 1
                    else:
                        ValidPeptides[Bin] = ValidPeptides.get(Bin, 0) + 1
                        CumulativeTruePeptides += 1
                ######################################################
                # Proteins:
                ValidProteins = {}
                InvalidProteins = {}
                CumulativeTrueProteins = 0
                CumulativeFalseProteins = 0
                for (ProteinID, Score) in BestScoreByProtein.items():
                    ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                    if ProteinName[:3] == "XXX":
                        InvalidProteins[Bin] = InvalidProteins.get(Bin, 0) + 1
                        CumulativeFalseProteins += 1
                    else:
                        ValidProteins[Bin] = ValidProteins.get(Bin, 0) + 1
                        CumulativeTrueProteins += 1
        for Bin in range(MinBin, MaxBin + 1):
            FDR = 1.0 - OddsTrue[Bin]
            AllHits = ScoreHistogram.get(Bin, 0)
            self.OutputDistributionFile.write("%s\t%s\t%s\t"%(Bin, FDR, AllHits))
            if self.ShuffledDatabaseFraction != None:
                FalseHits = ShuffledScoreHistogram.get(Bin, 0)
                TrueHits = AllHits - FalseHits
                CumulativeTrueHits -= TrueHits
                CumulativeFalseHits -= FalseHits                
                self.OutputDistributionFile.write("%s\t%s\t"%(CumulativeTrueHits, CumulativeFalseHits))
                if self.ProteinPicker:
                    # Peptide FDR:
                    FalseWithinTrue = min(CumulativeTruePeptides, CumulativeFalsePeptides * self.ShuffledScalingFactor)
                    PeptideFDR = FalseWithinTrue / float(max(1, CumulativeTruePeptides))
                    self.OutputDistributionFile.write("%.4f\t%s\t%s\t"%(PeptideFDR, CumulativeTruePeptides, CumulativeFalsePeptides))
                    CumulativeTruePeptides -= ValidPeptides.get(Bin, 0)
                    CumulativeFalsePeptides -= InvalidPeptides.get(Bin, 0)
                    # Protein FDR:
                    FalseWithinTrue = min(CumulativeTrueProteins, CumulativeFalseProteins * self.ShuffledScalingFactor)
                    ProteinFDR = FalseWithinTrue / float(max(1, CumulativeTrueProteins))
                    self.OutputDistributionFile.write("%.4f\t%s\t%s\t"%(ProteinFDR, CumulativeTrueProteins, CumulativeFalseProteins))
                    CumulativeTrueProteins -= ValidProteins.get(Bin, 0)
                    CumulativeFalseProteins -= InvalidProteins.get(Bin, 0)
            self.OutputDistributionFile.write("\n")
    def LoadPValueDistribution(self, FileName, Charge3Flag = 0):
        OddsTrue = {}
        if Charge3Flag:
            self.OddsTrue3 = OddsTrue
        else:
            self.OddsTrue2 = OddsTrue
        File = open(FileName, "rb")
        for FileLine in File.xreadlines():
            Bits = list(FileLine.strip().split("\t"))
            if len(Bits) < 2:
                continue            
            if FileLine[0] == "#":
                # Header line.  Parse special lines:
                Name = Bits[0][1:]
                if Name == "BlindFlag":
                    self.BlindFlag = int(Bits[1])
                elif Name == "MeanDeltaScore":
                    if Charge3Flag:
                        self.MeanDeltaScore3 = float(Bits[1])
                    else:
                        self.MeanDeltaScore2 = float(Bits[1])
                else:
                    print "(Skipping comment '%s', not understood)"%Bits[0]
                continue
            Bin = int(Bits[0])
            OddsTrue[Bin] = float(Bits[1])
        File.close()
        if self.BlindFlag:
            self.MQScoreWeight = Defaults.BlindMQScoreWeight
            self.DeltaScoreWeight = Defaults.BlindDeltaScoreWeight
            self.GammaOffset = Defaults.BlindGammaOffset
        else:
            self.MQScoreWeight = Defaults.MQScoreWeight
            self.DeltaScoreWeight = Defaults.DeltaScoreWeight
            self.GammaOffset = Defaults.GammaOffset
    def WriteMatchesForSpectrum(self, MatchesForSpectrum, OutHandle):
        if self.WriteTopMatchOnly:
            MatchesForSpectrum = MatchesForSpectrum[0:1]
        for Match in MatchesForSpectrum:
            # If we have a shuffled database (-S option), then by default we get to throw shuffled-protein
            # matches away for free.  We don't get to keep runners-up to them, though!!
            if Match.ProteinName[:3] == "XXX" and self.ShuffledDatabaseFraction != None and self.RemoveShuffledMatches:
                break
            # Skip matches with poor delta-score:
            if Match.DeltaScore < self.MaxDeltaScoreGap and not self.RetainBadMatches:
                continue
            # Skip short matches:
            Length = len(Match.Peptide.Aminos)
            if Length < self.MinimumPeptideLength:
                continue
            if Match.Charge < 3:
                MeanDeltaScore = self.MeanDeltaScore2
            else:
                MeanDeltaScore = self.MeanDeltaScore3
            WeightedScore = self.MQScoreWeight * Match.MQScore + self.DeltaScoreWeight * (Match.DeltaScore / MeanDeltaScore)
            ScoreBin = int(round(WeightedScore * BIN_MULTIPLIER))
            if Match.Charge < 3:
                TrueOdds = self.OddsTrue2.get(ScoreBin, None)
            else:
                TrueOdds = self.OddsTrue3.get(ScoreBin, None)
            if TrueOdds == None:
                if ScoreBin < 0:
                    TrueOdds = 0.00001
                else:
                    TrueOdds = 0.99999
            else:
                TrueOdds = max(0.00001, min(TrueOdds, 0.99999))
            Match.PValue = (1.0 - TrueOdds)
            Match.Bits.FScore = "%s"%WeightedScore
            Match.Bits.PValue = "%s"%Match.PValue
            if self.ProteinPicker:
                # Replace the original protein with the "correct" one:
                ProteinID = self.ProteinPicker.PeptideProteins.get(Match.Peptide.Aminos, None)
                if ProteinID != None:
                    Match.Bits.ProteinID = str(ProteinID)
                    Match.Bits.ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
            if (not self.RetainBadMatches):
                if (Match.PValue > self.PValueCutoff):
                    continue
                # Sometimes things with a horrible MQScore get a good pvalue.
                # We want to exclude these.
                if Match.MQScore < MIN_MQSCORE:
                    continue
            self.LinesAcceptedCount += 1
            OutHandle.write(str(Match.Bits))

    def WriteFixedScores(self, OutputPath):
        self.TotalLinesAcceptedCount = 0
        self.TotalLinesSecondPass = 0
        self.WriteScoresPath = OutputPath
        # Make the output directory, if it doesn't exist already.
        # Assume: OutputPath is a directory if ReadScoresPath is a directory,
        # and OutputPath is a file if ReadScoresPath is a file.
        if os.path.isdir(self.ReadScoresPath):
            DirName = OutputPath
        else:
            DirName = os.path.split(OutputPath)[0]
        try:
            os.makedirs(DirName)
        except:
            pass
        self.WriteFixedScoresFile(self.ReadScoresPath)
        print "Total accepted lines: %s of %s"%(self.TotalLinesAcceptedCount, self.TotalLinesSecondPass)

    def WriteFixedScoresFile(self, Path):
        try:
            inspectParser = InspectResults.Parser(Path,
                    inputMirrorTo=self.WriteScoresPath)

            LineCount = 0
            self.LinesAcceptedCount = 0
            OldSpectrum = None
            MatchesForSpectrum = []
            for row in inspectParser:
                Match = Bag()
                try:
                    Match.Bits = row
                    Match.Charge = row.Charge
                    Match.MQScore = row.MQScore
                    Match.DeltaScore = row.DeltaScore
                    Match.Peptide = GetPeptideFromModdedName(row.Annotation[2:-2])
                    Match.ProteinName = row.ProteinName
                except:
                    continue
                LineCount += 1
                Spectrum = (row.SpectrumFile, row.ScanNumber)
                if Spectrum != OldSpectrum:
                    self.WriteMatchesForSpectrum( MatchesForSpectrum,
                            inspectParser.mirrorOutHandle
                            )
                    MatchesForSpectrum = []
                OldSpectrum = Spectrum
                MatchesForSpectrum.append(Match)
            # Finish the last spectrum:
            self.WriteMatchesForSpectrum( MatchesForSpectrum,
                            inspectParser.mirrorOutHandle
                            )
            print "%s\t%s\t%s\t"%(Path, LineCount, self.LinesAcceptedCount)
            self.TotalLinesAcceptedCount += self.LinesAcceptedCount
            self.TotalLinesSecondPass += LineCount
        except:
            traceback.print_exc()
            print "* Error filtering annotations from '%s' to '%s'"%(Path, OutputPath)
    def ComputePValuesWithShuffled(self, Charge3Flag = 0):
        """
        Set self.OddsTrue using results from a partially-shuffled database.
        Given a score cutoff we assume that, above the score cutoff, there are
        T hits from valid proteins and F hits from invalid proteins.
        
        # PVALUE WITH REMOVAL:
        Let TDB and FDB be the true and false database fractions (FDB = self.ShuffledDatabaseFraction)
        After filtering out all F hits from invalid proteins, there are still some
        chance hits to true proteins.  The estimated number of hits from T that are actually
        false is equal to F*(TDB/FDB).  This, the odds true for this cutoff is:
        1.0 - (F*(TDB/FDB) / T)
        
        # PVALUE WITHOUT REMOVAL:
        The odds true for the cutoff is simply T/(T+F).  
        """
        OddsTrue = {}
        if Charge3Flag:
            self.OddsTrue3 = OddsTrue
            ScoreHistogram = self.ScoreHistogram3
            ShuffledScoreHistogram = self.ShuffledScoreHistogram3
        else:
            self.OddsTrue2 = OddsTrue
            ScoreHistogram = self.ScoreHistogram2
            ShuffledScoreHistogram = self.ShuffledScoreHistogram2
        CumulativeHits = 0
        CumulativeHitsTrue = 0
        CumulativeHitsFalse = 0
        Keys = ScoreHistogram.keys()
        Keys.sort()
        if not Keys:
            # There are NO HITS for this charge state.
            return
        MinKey = Keys[0]
        MaxKey = Keys[-1]
        TrueFraction = 1.0 - self.ShuffledDatabaseFraction
        self.ShuffledScalingFactor = TrueFraction / self.ShuffledDatabaseFraction
        for Key in range(MaxKey, MinKey - 1, -1):
            AllHits = ScoreHistogram.get(Key, 0)
            FalseHits = ShuffledScoreHistogram.get(Key, 0)
            ValidHits = AllHits - FalseHits
            CumulativeHitsTrue += ValidHits
            CumulativeHitsFalse += FalseHits
            FalseWithinTrue = min(CumulativeHitsTrue, CumulativeHitsFalse * self.ShuffledScalingFactor)
            if self.RemoveShuffledMatches:
                # OddsTrue = (V - I*(TDB/FDB)) / V
                BinOddsTrue = max(0, CumulativeHitsTrue - FalseWithinTrue) / float(max(1, CumulativeHitsTrue))
            else:
                # OddsTrue = (V - I*(TDB/FDB)) / (I+V)
                BinOddsTrue = max(0, CumulativeHitsTrue - FalseWithinTrue) / float(max(1, CumulativeHitsTrue + CumulativeHitsFalse))
            if self.VerboseFlag:
                # Bin, true, false, cumtrue, cumfalse
                Str = "%s\t%s\t%s\t%s\t%s\t"%(Key, ValidHits, FalseHits, CumulativeHitsTrue, CumulativeHitsFalse)
                Str += "%.5f\t%.5f\t"%(BinOddsTrue, 1.0 - BinOddsTrue)
                print Str
            OddsTrue[Key] = BinOddsTrue
        if self.VerboseFlag:
            print "\n\n"
    def SelectProteins(self, PValueCutoff, ReadScoresPath):
        """
        Using SelectProteins, assign each peptide the the most reasonable "owner" protein.
        """
        # Select the F-score cutoff:
        FScoreCutoff2 = 9999
        FScoreCutoff3 = 9999
        for FScoreBin in self.OddsTrue2.keys():
            OddsTrue = self.OddsTrue2[FScoreBin]
            if (1.0 - OddsTrue) <= PValueCutoff:
                if (FScoreBin / BIN_MULTIPLIER) < FScoreCutoff2:
                    FScoreCutoff2 = FScoreBin / BIN_MULTIPLIER
        for FScoreBin in self.OddsTrue3.keys():
            OddsTrue = self.OddsTrue3[FScoreBin]
            if (1.0 - OddsTrue) <= PValueCutoff:
                if (FScoreBin / BIN_MULTIPLIER) < FScoreCutoff3:
                    FScoreCutoff3 = FScoreBin / BIN_MULTIPLIER
        self.ProteinPicker.FScoreCutoff2 = FScoreCutoff2
        self.ProteinPicker.FScoreCutoff3 = FScoreCutoff3
        self.ProteinPicker.MeanDeltaScore2 = self.MeanDeltaScore2
        self.ProteinPicker.MeanDeltaScore3 = self.MeanDeltaScore3
        self.ProteinPicker.ParseAnnotations( ReadScoresPath )
        # We've COUNTED the protein hits.  Now ask the picker to decide which
        # protein 'owns' each peptide:
        self.ProteinPicker.ChooseProteins()
    def SetOutputDistributionPath(self, Path):
        self.OutputDistributionPath = Path
        self.OutputDistributionFile = open(Path, "wb")
    def ParseClusterInfo(self):
        """
        Parse cluster-sizes from an info file.  
        """
        self.ClusterSizes = {}
        File = open(self.ClusterInfoPath, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split()
            try:
                ScanNumber = int(Bits[1])
                ClusterSize = int(Bits[2])
            except:
                print "* Skipping this line:", FileLine
            self.ClusterSizes[(Bits[0], ScanNumber)] = ClusterSize
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "l:s:r:w:m:bp:vixzd:a1S:HX:")
        OptionsSeen = {}
        self.SaveDistributionPath = "PValues.txt" # default
        self.ReadScoresPath = None
        self.WriteScoresPath = None
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-l":
                if not os.path.exists(Value):
                    print "** Error: can't read p-value distribution from file '%s'"%Value
                    return 0
                self.LoadDistributionPath = Value
            elif Option == "-p":
                self.PValueCutoff = float(Value)
            elif Option == "-s":
                self.SaveDistributionPath = Value
            elif Option == "-x":
                self.RetainBadMatches = 1
            elif Option == "-b":
                self.BlindFlag = 1
                self.MQScoreWeight = Defaults.BlindMQScoreWeight
                self.DeltaScoreWeight = Defaults.BlindDeltaScoreWeight
                self.GammaOffset = Defaults.BlindGammaOffset
            elif Option == "-r":
                self.ReadScoresPath  = Value
            elif Option == "-w":
                self.WriteScoresPath = Value
            elif Option == "-m":
                global MAX_RESULTS_FILES_TO_PARSE
                MAX_RESULTS_FILES_TO_PARSE = int(Value)
            elif Option == "-v":
                self.VerboseFlag = 1
            elif Option == "-i":
                self.GenerateImageFlag = 1
            elif Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath.append(Value)
            elif Option == "-a":
                self.PerformProteinSelection = 1
            elif Option == "-1":
                self.WriteTopMatchOnly = 1
            elif Option == "-S":
                self.ShuffledDatabaseFraction = float(Value)
                if self.ShuffledDatabaseFraction <= 0 or self.ShuffledDatabaseFraction >= 1:
                    print "* Invalid value for -S: %s"%Value
                    return 0
            elif Option == "-H":
                self.RemoveShuffledMatches = 0
            elif Option == "-X":
                # Undocumented option for CLUSTER searches:
                self.ClusterInfoPath = Value
                self.ParseClusterInfo()
            else:
                print "** Unknown option:", Option, Value
        # Check validity of options:
        if self.PerformProteinSelection and not self.DBPath:
            print "* Error: -a option requires -d option!"
            return 0
        # No major problems - return TRUE for success.
        return 1

UsageInfo = """
PValue.py - Compute probability that each match from a tandem MS
peptide database search is correct.  Write out an updated results file containing
only the high-quality results.

Parameters:
 -r [FILENAME] Read results from filename (and fit the probability mixture
    model to these results).  If the option value is a directory, we'll read
    all the results-files from the directory.
 -w [FILENAME] Write re-scored results to a file.
 -l [FILENAME] Load p-value distribution from a file (written out earlier
    with -s option)

Protein selection can be performed, replacing the protein identification
with a parsimonious set of protein IDs (using a simple iterative
approach).  The following options are required for protein selection:
 -a: Replace protein identifications with a "parsimonious" set of protein IDs.
     Requires -d option!
 -d [FILENAME] Database (.trie file) searched
 -S [FRACTION]: (see below)

Other options:
 -S [FRACTION]: The fraction of the database consisting of shuffled
    proteins.  For instance, if you use a 1:1 mix of valid and invalid
    proteins, use -S 0.5.  If this option is set, p-values will be set using
    the number of matches to shuffled proteins, whose names begin with XXX
 -s [FILENAME] Save p-value distribution to a file.
 -i Write a .png image of the distribution graph (requires PIL)
 -p [PVALUE] Cutoff for saving results; by default, 0.1
 -b Blind search (use different score/deltascore weighting)
 -x If the -x flag is passed, even "bad" matches are written out (no p-value
    filtering is performed)
 -1 Write only the top hit for each spectrum, even if "good" runners-up exist

Internal use only:
 -v Verbose output (for debugging)
 -H Retain matches to shuffled proteins.  Used for further processing ONLY.
    
Example:
  PValue.py -r ShewanellaResults -s ShewPValue.txt -w ShewanellaFiltered
     -p 0.05 -d database\Shew.trie -a 
"""

def Main(Parser = None):
    global MAX_RESULTS_FILES_TO_PARSE
    if not Parser:
        Parser = PValueParser()
        Result = Parser.ParseCommandLine(sys.argv[1:])
        if not Result:
            print UsageInfo
            return
    if Parser.DBPath and Parser.PerformProteinSelection:
        Parser.ProteinPicker = SelectProteins.ProteinSelector()
        Parser.ProteinPicker.LoadMultipleDB(Parser.DBPath)
    if Parser.LoadDistributionPath:
        print "Load p-value distribution from %s..."%Parser.LoadDistributionPath
        Parser.LoadPValueDistribution(Parser.LoadDistributionPath)
    elif Parser.ReadScoresPath:
        print "Read scores from search results at %s..."%Parser.ReadScoresPath
        Parser.ReadDeltaScoreDistribution(Parser.ReadScoresPath)
        Parser.SetOutputDistributionPath(Parser.SaveDistributionPath)
        ##############################
        # Loop for F-score methods
        Parser.ReadScoreDistributionFromFile( Parser.ReadScoresPath )
        if Parser.ShuffledDatabaseFraction != None:
            print "Compute PValues with shuffled..."
            Parser.ComputePValuesWithShuffled(0)
            Parser.ComputePValuesWithShuffled(1)
        else:
            Result = Parser.FitMixtureModel()
            if not Result:
                sys.exit(1)
        if Parser.PerformProteinSelection:
            Parser.SelectProteins(Parser.PValueCutoff, Parser.ReadScoresPath)
        print "Write p-value distribution to %s..."%Parser.SaveDistributionPath
        (Stub, Extension) = os.path.splitext(Parser.SaveDistributionPath)
        Parser.SavePValueDistribution(0)
        Parser.SavePValueDistribution(1)
        ##############################
        Parser.OutputDistributionFile.close()
    else:
        print "** Please specify either a distribution file or results file."
        print UsageInfo
        sys.exit(1)
    if Parser.GenerateImageFlag and Image:
        ImagePath = os.path.splitext(Parser.SaveDistributionPath)[0] + ".2.png"
        Parser.ProduceScoreDistributionImage(ImagePath, 0)
        ImagePath = os.path.splitext(Parser.SaveDistributionPath)[0] + ".3.png"
        Parser.ProduceScoreDistributionImage(ImagePath, 1)
    if Parser.WriteScoresPath:
        Parser.WriteFixedScores(Parser.WriteScoresPath)
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "psyco not found - running without optimization"
    #TestMain()
    Main()
