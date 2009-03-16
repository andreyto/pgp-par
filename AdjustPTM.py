"""
Merge and reconcile peptide species, after running ComputePTMFeatures
and TrainPTMFeatures.  Iterate over peptide species from "best" to
"worst".  For each species, consider whether there's another species
which (1) is better, and either (2a) is the same after mod-shifting
(including charge-state), or (2b) is compatible after mod-shifting.
Case (2a): Try merging both species into one big cluster, determine
whether the MQScore / PValue improves.  If so, keep the merge.
Case (2b): Try shifting the inferior species to match the superior.
If the MQScore / PValue improves (or at least stays CLOSE), then
keep the shift.

A note on directories:
Consensus spectra and clusters can take quite a bit of disk space.  Running
AdjustPTM changes spectra and clusters, but we'd like the liberty to
re-run AdjustPTM.  Therefore, AdjustPTM uses a set of "adjusted" directories.
It wipes these when it starts a merge/reconcile run; it writes to them; it
reads clusters and spectra from these directories first, if possible.
"""
import os
import struct
import shutil
import math
import sys
import MSSpectrum
import string
import traceback
import getopt
import cPickle
import BuildConsensusSpectrum
import Learning
import PyInspect
import ResultsParser
import SpectralSimilarity
from Utils import *
Initialize()
import TrainPTMFeatures
from TrainPTMFeatures import FormatBits
from TrainPTMFeatures import FeatureBits

PROFILING_RUN = 0

class AnticipatedChemistry:
    """
    Represents a chemical adduct which we expect to see often, with relatively low
    site specificity.  Examples: M+16, *.Q-17.  We want to flag the adducts
    to highlight the remaining NON-adducts!
    """
    def __init__(self):
        # For the "allowed" members, None means "no restriction".
        self.AllowedResidues = None
        self.AllowedPrefix = None
        self.AllowedSuffix = None
        self.Terminus = None
        self.Mass = 0
        self.Name = ""

class SiteClass:
    """
    Wrapper for one or more Species instances which represent the same modification-mass
    at the same database-position.
    """
    def __init__(self):
        self.ModDBPos = None
        self.ModMass = None
        self.SpeciesList = []
    def __str__(self):
        return "%+d on dbpos %s"%(self.ModMass, self.ModDBPos)
    
class PeptideSpeciesClass:
    """
    Represents one (modified) peptide species; different charge states are different
    peptide species.
    """
    def __init__(self):
        self.MergedFlag = 0
        self.ConsensusSpectrum = None
        self.ConsensusMSSpectrum = None
        self.BestSpectrum = None
        self.BestMSSpectrum = None
        self.Bits = None
        self.ConsensusModlessMSSpectrum = None
    def __str__(self):
        return "<peptide species: %s>"%self.Bits
    def FreeCachedSpectra(self):
        """
        Discard our PySpectrum and MSSpectrum objects, because we can't hold ALL
        such in memory at once:
        """
        self.BestSpectrum = None
        self.BestMSSpectrum = None
        self.ConsensusSpectrum = None
        self.ConsensusMSSpectrum = None
        self.ConsensusModlessMSSpectrum = None
    def GetBestSpectrum(self, Master):
        if self.BestSpectrum:
            return self.BestSpectrum
        FilePath = self.Bits[FormatBits.BestSpectrumPath]
        ColonBits = FilePath.split(":")
        try:
            FilePos = int(ColonBits[-1])
            FilePath = string.join(ColonBits[:-1], ":")
        except:
            FilePos = 0
        FilePath = Master.FixSpectrumPath(FilePath)
        self.BestSpectrum = PyInspect.Spectrum(FilePath, FilePos)
        return self.BestSpectrum
    def GetMemberListStr(self, Master):
        Path = os.path.join(Master.ClusterScanListDirAdjusted, self.Annotation[2], "%s.%s.txt"%(self.Annotation, self.Charge))
        if not os.path.exists(Path):
            Path = os.path.join(Master.ClusterScanListDir, self.Annotation[2], "%s.%s.txt"%(self.Annotation, self.Charge))
        if not os.path.exists(Path):
            # Punt!
            return ""
        File = open(Path, "rb")
        Text = File.read()
        File.close()
        return Text
    def GetConsensusSpectrumPath(self, Master):
        Path = os.path.join(Master.ConsensusSpectraDirAdjusted, self.Annotation[2], "%s.%s.dta"%(self.Annotation, self.Charge))
        if not os.path.exists(Path):
            Path = os.path.join(Master.ConsensusSpectraDir, self.Annotation[2], "%s.%s.dta"%(self.Annotation, self.Charge))
        return Path
    def GetConsensusSpectrum(self, Master):
        if self.ConsensusSpectrum:
            return self.ConsensusSpectrum
        Path = self.GetConsensusSpectrumPath(Master)
        self.ConsensusSpectrum = PyInspect.Spectrum(Path, 0)
        return self.ConsensusSpectrum
    def GetConsensusMSSpectrum(self, Master):
        if self.ConsensusMSSpectrum:
            return self.ConsensusMSSpectrum
        Path = self.GetConsensusSpectrumPath(Master)
        self.ConsensusMSSpectrum = MSSpectrum.SpectrumClass()
        self.ConsensusMSSpectrum.ReadPeaks(Path)
        self.ConsensusMSSpectrum.FilterPeaks()
        self.ConsensusMSSpectrum.RankPeaksByIntensity()
        return self.ConsensusMSSpectrum
    def GetConsensusModlessMSSpectrum(self, Master):
        if self.ConsensusModlessMSSpectrum:
            return self.ConsensusModlessMSSpectrum
        Path = os.path.join(Master.ConsensusSpectraDirAdjusted, self.ModlessAnnotation[2], "%s.%s.dta"%(self.ModlessAnnotation, self.Charge))
        if not os.path.exists(Path):
            Path = os.path.join(Master.ConsensusSpectraDir, self.ModlessAnnotation[2], "%s.%s.dta"%(self.ModlessAnnotation, self.Charge))
        self.ConsensusModlessMSSpectrum = MSSpectrum.SpectrumClass()
        self.ConsensusModlessMSSpectrum.ReadPeaks(Path)
        self.ConsensusModlessMSSpectrum.FilterPeaks()
        self.ConsensusModlessMSSpectrum.RankPeaksByIntensity()
        return self.ConsensusModlessMSSpectrum
    def ParseBits(self, Bits):
        self.Bits = Bits
        self.Annotation = Bits[FormatBits.Peptide]
        self.Peptide = GetPeptideFromModdedName(self.Annotation)
        self.ModlessAnnotation = "%s.%s.%s"%(self.Peptide.Prefix, self.Peptide.Aminos, self.Peptide.Suffix)
        self.Charge = int(Bits[FormatBits.Charge])
        self.ModDBPos = int(Bits[FormatBits.DBPos])
        ModIndex = self.Peptide.Modifications.keys()[0]
        self.DBPos = self.ModDBPos - ModIndex
        self.ModMass = self.Peptide.Modifications[ModIndex][0].Mass
        self.ModAA = self.Peptide.Aminos[ModIndex]
        self.DBEnd = self.DBPos + len(self.Peptide.Aminos)
        self.ConsensusMQScore = float(Bits[FormatBits.ConsensusMQScore])
        try:
            self.ModelScore = float(Bits[FormatBits.ModelScore])
            self.PValue = float(Bits[FormatBits.ModelPValue])
        except:
            self.ModelScore = None
            self.PValue = None
        # Parse old features:
        self.Features = []
        for FeatureIndex in range(FormatBits.FirstFeature, FormatBits.LastFeature + 1):
            try:
                self.Features.append(float(Bits[FeatureIndex]))
            except:
                self.Features.append(0)
        self.ComputePrefixes()
    def ComputePrefixes(self):
        # self.Prefixes[DBPos] is the mass that this species accumulates
        # *before* the specified residue.  Examples:
        # ~ Species.Prefixes[Species.DBPos] = 0 always,
        # ~ Species.Prefix[Species.DBPos + 1] is equal to the mass (with modification, if any)
        #    of the first residue
        self.Prefixes = {}
        self.Suffixes = {}
        ParentMass = self.Peptide.Masses[-1] + 19
        AccumulatedMass = 0
        for Pos in range(len(self.Peptide.Aminos)):
            self.Prefixes[self.DBPos + Pos] = AccumulatedMass
            self.Suffixes[self.DBPos + Pos] = ParentMass - AccumulatedMass
            AccumulatedMass += GetMass(self.Peptide.Aminos[Pos])
            for Mod in self.Peptide.Modifications.get(Pos, []):
                AccumulatedMass += Mod.Mass

class PTMAdjuster(ResultsParser.SpectrumOracleMixin):
    def __init__(self):
        self.HeaderLines = []
        self.CompatibilityTolerance = 3
        self.CachedClusterPath = None
        self.ConsensusClusterDir = "PTMScore\\Lens-99-10\\Cluster" # default
        self.ConsensusSpectraDir = "PTMScore\\Lens-99-10\\Spectra" # default
        self.SortByModel = 1
        self.PeptideDict = {} # keys: (Annotation, Charge)
        self.KnownChemistryFileName = None
        self.OutputModelFileName2 = None
        self.OutputModelFileName3 = None
        self.DBStart = None
        self.DBEnd = None
        self.SpectrumRoot = None
        self.CheckDirectoriesFlag = 0
        self.MergeBlockRunsFlag = 0
        self.KnownPTMVerboseOutputFileName = None
        self.MaxPeptideWindowWidth = 2500
        ResultsParser.SpectrumOracleMixin.__init__(self)
    def PerformMergeReconcileOnWindow(self, PerformMergeFlag):
        ###############################################################
        # Consider merging/reconciling these peptides:
        SortedList = []
        self.PeptideDict = {}
        for Species in self.WindowPeptides:
            if self.SortByModel:
                SortedList.append((Species.ModelScore, Species))
            else:
                SortedList.append((Species.ConsensusMQScore, Species))
            Key = (Species.Annotation, Species.Charge)
            self.PeptideDict[Key] = Species
        # SortedList lists species from BEST to WORST.
        SortedList.sort()
        SortedList.reverse()
        # Dual iteration over the peptides from the window: Species A has the
        # lower score, species B has the higher score.
        # Consider reconciling species A to species B:
        for IndexA in range(len(SortedList)):
            (ScoreA, SpeciesA) = SortedList[IndexA]
            Str = "(%s/%s) %s %s"%(IndexA, len(SortedList), SpeciesA.Charge, SpeciesA.Annotation)
            if PerformMergeFlag:
                print "M", Str
            else:
                print "C", Str
            if SpeciesA.MergedFlag:
                # A has already been merged into another species.
                continue
            for IndexB in range(IndexA):
                (ScoreB, SpeciesB) = SortedList[IndexB]
                if SpeciesB.MergedFlag:
                    # B has already been merged into another species.
                    continue
                # Compatibility checks.
                # Charge must be the same in order to MERGE (but not to RECONCILE):
                if SpeciesA.Charge != SpeciesB.Charge and PerformMergeFlag:
                    continue
                # Peptides must overlap:
                if SpeciesA.DBEnd <= SpeciesB.DBPos or SpeciesB.DBEnd <= SpeciesA.DBPos:
                    continue
                # To reconcile, Peptide A must cover the DBposition which is modified in B:
                if SpeciesB.ModDBPos >= SpeciesA.DBEnd or SpeciesB.ModDBPos < SpeciesA.DBPos:
                    if not PerformMergeFlag:
                        continue
                # First, look for a MERGE:
                # Prefix and suffix must be the same at some point:
                SamePrefixSuffix = 0
                for DBPos in SpeciesA.Prefixes.keys():
                    PMassA = SpeciesA.Prefixes[DBPos]
                    PMassB = SpeciesB.Prefixes.get(DBPos, -9999)
                    if abs(PMassA - PMassB) >= self.CompatibilityTolerance:
                        continue
                    SMassA = SpeciesA.Suffixes[DBPos]
                    SMassB = SpeciesB.Suffixes.get(DBPos, -9999)
                    if abs(SMassA - SMassB) >= self.CompatibilityTolerance:
                        continue
                    SamePrefixSuffix = 1
                    break
                if SamePrefixSuffix and SpeciesA.Charge == SpeciesB.Charge:
                    # Merge is possible.  If this is first-cycle, then do a merge;
                    # if not, then continue.
                    if PerformMergeFlag:
                        MergeFlag = self.AttemptMerge(SpeciesA, SpeciesB) # A into B
                        if MergeFlag:
                            SpeciesB.FreeCachedSpectra()
                            break
                        # We didn't merge A into B.  But perhaps we can merge B into A!
                        MergeFlag = self.AttemptMerge(SpeciesB, SpeciesA, 1) # B into A
                        if MergeFlag:
                            SpeciesB.FreeCachedSpectra()
                            break
                    continue
                # Merge is impossible.  If this is first-cycle, bail out:
                if PerformMergeFlag:
                    continue
                # If species A and B are already compatible, then there's nothing to do:
                if SpeciesA.ModDBPos == SpeciesB.ModDBPos and SpeciesA.ModMass == SpeciesB.ModMass:
                    print "(Already reconciled to %s)"%Species.Annotation
                    continue
                # Perhaps A could be conformed to B
                # if the modification-masses are similar (possibly after
                # an endpoint shift):
                if SpeciesA.DBPos < SpeciesB.DBPos:
                    ExtraPrefixA = GetMass(self.DB[SpeciesA.DBPos:SpeciesB.DBPos])
                else:
                    ExtraPrefixA = 0
                if SpeciesB.DBPos < SpeciesA.DBPos:
                    ExtraPrefixB = GetMass(self.DB[SpeciesB.DBPos:SpeciesA.DBPos])
                else:
                    ExtraPrefixB = 0
                if SpeciesA.DBEnd > SpeciesB.DBEnd:
                    ExtraSuffixA = GetMass(self.DB[SpeciesB.DBEnd:SpeciesA.DBEnd])
                else:
                    ExtraSuffixA = 0
                if SpeciesB.DBEnd > SpeciesA.DBEnd:
                    ExtraSuffixB = GetMass(self.DB[SpeciesA.DBEnd:SpeciesB.DBEnd])
                else:
                    ExtraSuffixB = 0
                # VERBOSE:
                for DBPos in SpeciesA.Prefixes.keys():
                    PMassA = SpeciesA.Prefixes[DBPos] + ExtraPrefixA
                    PMassB = SpeciesB.Prefixes.get(DBPos, -9999) - ExtraPrefixB
                    SMassA = SpeciesA.Suffixes[DBPos] + ExtraSuffixA
                    SMassB = SpeciesB.Suffixes.get(DBPos, -9999) - ExtraSuffixB
                    #print "DBPos %s: PreA %s PreB %s (%s)"%(DBPos, PMassA, PMassB, abs(PMassA-PMassB))
                    #print "   PostA %s PostB %s (%s)"%(SMassA, SMassB, abs(SMassA-SMassB))
                    if abs(PMassA - PMassB) >= self.CompatibilityTolerance:
                        continue
                    if abs(SMassA - SMassB) >= self.CompatibilityTolerance:
                        continue
                    SamePrefixSuffix = 1
                    break
                if not SamePrefixSuffix:
                    # irreconcilable, move on to try the next species:
                    continue
                self.ReconcileDetailOutput.write("%s\t%s\t%s\t%s\t\n"%(SpeciesA.ModDBPos, SpeciesB.ModDBPos, SpeciesA.ModMass, SpeciesB.ModMass))
                self.ReconcileDetailOutput.write("%s\t%s\t%s\t\t%s\t%s\t%s\t\n"%(SpeciesA.Annotation, SpeciesA.DBPos, SpeciesA.DBEnd, SpeciesB.Annotation, SpeciesB.DBPos, SpeciesB.DBEnd))                    
                # Species A *could* be reconciled with B.
                Result = self.AttemptReconcile(SpeciesA, SpeciesB)
                SpeciesB.FreeCachedSpectra()
                if Result:
                    # We reconciled to B.  Stop now, don't re-reconcile to
                    # another (INFERIOR) species:
                    break                    
            SpeciesA.FreeCachedSpectra()        
    def PerformAllMerges(self, PerformMergeFlag = 1):
        """
        Workhorse for the merge and reconcile procedure.  Double-loop over sites,
        from high to low scoring.  Consider either MERGING or RECONCILING the
        low-scoring site to match the high-scoring site.
        """
        self.HeaderLines = []
        self.ModTypeSpectrumCount = {}
        self.ModTypeSiteCount = {}
        # A list of peptides within our 'window'.  At each iteration,
        # we read peptides until we hit eof, hit a new protein-name, or hit a peptide 2500 characters
        # past the first one.  Then we attempt reconciliation/merging for peptides
        # which are "covered" by the window.  Then we advance the window.
        self.WindowPeptides = []
        EOFFlag = 0
        WroteHeaderFlag = 0
        NextProteinFirstPeptide = None
        while 1:
            if NextProteinFirstPeptide:
                self.WindowPeptides = [NextProteinFirstPeptide]
                NextProteinFirstPeptide = None
            if not len(self.WindowPeptides):
                WindowStart = 0
                WindowEnd = 0
                CurrentProteinName = None
                if EOFFlag:
                    break
            else:
                WindowStart = self.WindowPeptides[0].DBPos
                WindowEnd = self.WindowPeptides[-1].DBEnd
                CurrentProteinName = self.WindowPeptides[0].ProteinName
            ###############################################################
            # Parse some more peptides:
            while (not EOFFlag) and (WindowEnd < WindowStart + self.MaxPeptideWindowWidth):
                FileLine = self.InputFile.readline()
                if not FileLine:
                    EOFFlag = 1
                    break
                if FileLine[0] == "#":
                    self.HeaderLines.append(FileLine)
                    continue # skip comment line
                FileLine = FileLine.replace("\r", "").replace("\n", "")
                if not FileLine.strip():
                    continue # skip blank line
                Bits = FileLine.split("\t")
                try:
                    Species = PeptideSpeciesClass()
                    Species.ParseBits(Bits)
                    Species.ProteinName = Bits[FormatBits.ProteinName]
                except:
                    traceback.print_exc()
                    print Bits
                    continue
                # SKIP the species if it falls outside our block:
                if self.DBStart != None and Species.DBPos < self.DBStart:
                    continue
                if self.DBEnd != None and Species.DBPos >= self.DBEnd:
                    continue
                # If the species comes from a NEW protein, finish the window and save the new species
                # for next iteration:
                #print "CurrentProtein '%s', new species protein '%s'"%(CurrentProteinName, Species.ProteinName)
                if CurrentProteinName == None:
                    # We have no current-protein.  Start the list:
                    CurrentProteinName = Species.ProteinName
                    WindowStart = Species.DBPos
                else:
                    # Check whether this species matches the current protein:
                    if Species.ProteinName != CurrentProteinName:
                        NextProteinFirstPeptide = Species
                        break
                self.WindowPeptides.append(Species)
                WindowStart = min(WindowStart, Species.DBPos)
                WindowEnd = max(WindowEnd, Species.DBEnd)
            ###############################################################
            print "->Handling %s peptides in range %s...%s\n  %s"%(len(self.WindowPeptides), WindowStart, WindowEnd, CurrentProteinName)
            self.PerformMergeReconcileOnWindow(PerformMergeFlag)
            ###############################################################
            # Re-sort the window peptides, so that sites fall together:
            SortedList = []
            for Peptide in self.WindowPeptides:
                SortedList.append((Peptide.ModDBPos, Peptide.ModMass, Peptide))
            SortedList.sort()
            self.WindowPeptides = []
            for Tuple in SortedList:
                self.WindowPeptides.append(Tuple[-1])
            ###############################################################
            # Write file header now, if we haven't:
            if not WroteHeaderFlag:
                for HeaderLine in self.HeaderLines:
                    self.OutputFile.write(HeaderLine)
                WroteHeaderFlag = 1
            ###############################################################                
            # Write out, and free, peptides in (the early portion of) the window:
            PeptideIndex = 0
            while PeptideIndex < len(self.WindowPeptides):
                Species = self.WindowPeptides[PeptideIndex]
                ###print "  Species %s of %s: %s-%s (window %s-%s)"%(PeptideIndex, len(self.WindowPeptides), Species.DBPos, Species.DBEnd, WindowStart, WindowEnd)
                if Species.DBEnd > WindowEnd - 100 and (not EOFFlag and not NextProteinFirstPeptide):
                    ###print "  ->Leave it in the window"
                    PeptideIndex += 1
                    continue
                # This peptide can be dropped from the list.
                if Species.MergedFlag:
                    pass
                else:
                    # Update Species.Bits to reflect changes to Species.Features:
                    for Index in range(FormatBits.FirstFeature, FormatBits.LastFeature + 1):
                        Species.Bits[Index] = str(Species.Features[Index - FormatBits.FirstFeature])
                    String = string.join(Species.Bits, "\t")
                    self.OutputFile.write(String + "\n")
                    ModTypeKey = (Species.ModAA, Species.ModMass)
                    # Note NOW the number of sites and spectra for each modification-type.
                    self.ModTypeSiteCount[ModTypeKey] = self.ModTypeSiteCount.get(ModTypeKey, 0) + 1
                    self.ModTypeSpectrumCount[ModTypeKey] = self.ModTypeSpectrumCount.get(ModTypeKey, 0) + Species.Features[FeatureBits.SpectrumCount]
                ###print "  ->Drop it from the window"
                del self.WindowPeptides[PeptideIndex]
    def GetSiteScore(self, Site):
        Site.PValue = 1.0
        SortedSpeciesScores = []
        for Species in Site.SpeciesList:
            SortedSpeciesScores.append(Species.ModelScore)
            if Species.Charge > 2:
                Species.PValue = self.Model3.GetPValue(Species.ModelScore)
            else:
                Species.PValue = self.Model2.GetPValue(Species.ModelScore)
            Site.PValue *= Species.PValue
        SortedSpeciesScores.sort()
        SortedSpeciesScores.reverse()
        SiteScore = [-math.log(Site.PValue)]
        SiteScore.extend(SortedSpeciesScores)
        return SiteScore
    def LoadKnownModifications(self):
        """
        Parse the KnownPTM file, to get a list of chemical events which we
        already have a (hypothetical) annotation for.
        """
        if not self.KnownChemistryFileName:
            return
        self.KnownPTMs = []
        self.KnownPTMByMass = {}
        File = open(self.KnownChemistryFileName, "rb")
        # Load one PTM from each line of the file.
        # Line format:
        # Mass, name, AllowedResidues, Terminus, AllowedPrefix, AllowedSuffix
        # Example: -17, pyroglutamate formation, Q, N, "", ""
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = list(FileLine.strip().split("\t"))
            if len(Bits) == 1:
                continue
            if FileLine[0] == "#":
                continue
            while len(Bits)<6:
                Bits.append("")
            try:
                PTM = AnticipatedChemistry()
                PTM.Mass = int(Bits[0])
            except:
                print "** Skipping invalid line %s of %s"%(LineNumber, self.KnownChemistryFileName)
                print Bits
                continue
            PTM.Name = Bits[1]
            if len(Bits[2]) == 0 or Bits[2][0] == "*":
                PTM.AllowedResidues = None
            else:
                PTM.AllowedResidues = Bits[2]
            if len(Bits[3].strip()):
                PTM.Terminus = Bits[3].strip()
            if len(Bits[4].strip()):
                PTM.AllowedPrefix = Bits[4].strip()
            if len(Bits[5].strip()):
                PTM.AllowedSuffix = Bits[5].strip()
            self.KnownPTMs.append(PTM)
            if not self.KnownPTMByMass.has_key(PTM.Mass):
                self.KnownPTMByMass[PTM.Mass] = []
            self.KnownPTMByMass[PTM.Mass].append(PTM)
        File.close()
        print "Loaded %s known PTMs from '%s'."%(len(self.KnownPTMs), self.KnownChemistryFileName)
    def AttemptKnownPTM(self, Site):
        InitialScore = self.GetSiteScore(Site)
        # Initialize known-ptm information for each species:
        for Species in Site.SpeciesList:
            Species.BestKnownPTMAnnotation = ""
            Species.BestKnownPTMAnnotationName = ""
            Species.BestKnownPTMAnnotationPValue = ""
            Species.BestKnownPTMAnnotationScore = ""
            Species.BestKnownPTMAnnotationSitePValue = ""
        # Loop over species to find the allowed database residues.
        ResidueCounts = {}
        ResidueInitialCounts = {}
        for Species in Site.SpeciesList:
            for Pos in range(Species.DBPos, Species.DBPos + len(Species.Peptide.Aminos)):
                ResidueCounts[Pos] = ResidueCounts.get(Pos, 0) + 1
            ResidueInitialCounts[Species.DBPos] = ResidueInitialCounts.get(Species.DBPos, 0) + 1
        BestEditedSiteScore = None
        self.KPTMVerbose.write("\n===============================\n")
        self.KPTMVerbose.write("Site %s initial score (%.3f, %.3f)\n"%(Site, InitialScore[0], InitialScore[1]))
        Residues = ResidueCounts.items()
        FirstCoveredResidue = min(ResidueCounts.keys())
        LastCoveredResidue = max(ResidueCounts.keys())
        ###############################################################################
        # Consider shifting the modification to any (legal) residue, with any legal mass:
        ModMass = Site.SpeciesList[0].ModMass
        # Decide which endpoint-shifts we'll considered.
        # If the peptides don't all share the same N-terminus, then shifting the
        # N-terminus isn't allowed.  (That would perforce SPLIT this site!)
        Shifts = [None]
        if ResidueCounts[FirstCoveredResidue] == len(Site.SpeciesList):
            Shifts.extend(["N+1", "N+2"])
            if FirstCoveredResidue > 0:
                Shifts.append("N-1")
            if FirstCoveredResidue > 1:
                Shifts.append("N-2")
        if ResidueCounts[LastCoveredResidue] == len(Site.SpeciesList):
            Shifts.extend(["C-2", "C-1"])
            if LastCoveredResidue < len(self.DB) - 1:
                Shifts.append("C+1")
            if LastCoveredResidue < len(self.DB) - 2:
                Shifts.append("C+2")
        for Shift in Shifts:
            # CoreMass is equal to the modified mass, possibly modified due to
            # endpoint shifts:
            CoreMass = ModMass
            # AllowedDBList is a list of the database positions where the ptm could
            # be attached.  We've already selected the range the peptide will cover,
            # but the PTM could fall on various residues:
            AllowedDBList = list(range(FirstCoveredResidue, LastCoveredResidue + 1))
            if Shift in ("N-1", "N-2"):
                AllowedDBList.append(FirstCoveredResidue - 1)
                DropMass = Global.AminoMass.get(self.DB[FirstCoveredResidue - 1], None)
                if not DropMass:
                    continue
                CoreMass -= DropMass
            if Shift == "N-2":
                AllowedDBList.append(FirstCoveredResidue - 2)
                DropMass = Global.AminoMass.get(self.DB[FirstCoveredResidue - 2], None)
                if not DropMass:
                    continue
                CoreMass -= DropMass
            if Shift in ("N+1", "N+2"):
                AllowedDBList.remove(FirstCoveredResidue)
                CoreMass += Global.AminoMass.get(self.DB[FirstCoveredResidue], None)
            if Shift == "N+2":
                AllowedDBList.remove(FirstCoveredResidue + 1)
                CoreMass += Global.AminoMass.get(self.DB[FirstCoveredResidue + 1], None)
            if Shift in ("C+1", "C+2"):
                AllowedDBList.append(LastCoveredResidue + 1)
                DropMass = Global.AminoMass.get(self.DB[LastCoveredResidue + 1], None)
                if not DropMass:
                    continue
                CoreMass -= DropMass
            if Shift == "C+2":
                AllowedDBList.append(LastCoveredResidue + 2)
                DropMass = Global.AminoMass.get(self.DB[LastCoveredResidue + 2], None)
                if not DropMass:
                    continue
                CoreMass -= DropMass
            if Shift in ("C-1", "C-2"):
                AllowedDBList.remove(LastCoveredResidue)
                CoreMass += Global.AminoMass.get(self.DB[LastCoveredResidue], None)
            if Shift == "C-2":
                AllowedDBList.remove(LastCoveredResidue - 1)
                CoreMass += Global.AminoMass.get(self.DB[LastCoveredResidue - 1], None)
            if CoreMass < -250 or CoreMass > 250:
                continue
            CoreMass = int(round(CoreMass))
            ShiftablePeptides = []
            for DBPos in AllowedDBList:
                TryMassList = (CoreMass - 2, CoreMass - 1, CoreMass, CoreMass + 1, CoreMass + 2)
                for NearMass in TryMassList:
                    KnownPTMList = self.KnownPTMByMass.get(NearMass, [])
                    for KnownPTM in KnownPTMList:
                        # Determine whether this is a legal PTM placement.
                        # The amino acid type must be valid:
                        if KnownPTM.AllowedResidues != None and self.DB[DBPos] not in KnownPTM.AllowedResidues:
                            continue
                        # The terminus must be valid (for at least one peptide species):
                        if KnownPTM.Terminus in ("N", "^"):
                            if DBPos > FirstCoveredResidue:
                                continue
                        # The prefix and suffix residues must be valid:
                        if DBPos:
                            PrefixAA = self.DB[DBPos - 1]
                        else:
                            PrefixAA = "-"
                        if DBPos < len(self.DB) - 1:
                            SuffixAA = self.DB[DBPos + 1]
                        else:
                            SuffixAA = "-"
                        if KnownPTM.AllowedPrefix != None and PrefixAA not in KnownPTM.AllowedPrefix:
                            continue
                        if KnownPTM.AllowedSuffix != None and SuffixAA not in KnownPTM.AllowedSuffix:
                            continue
                        ############################################################
                        # Okay, this is a LEGAL placement.  Determine its score:
                        Score = self.TryShiftedSite(Site, NearMass, DBPos, KnownPTM, Shift)
                        if Score > BestEditedSiteScore:
                            BestEditedSiteScore = Score
                            self.RememberOptimalKnownPTM(Site)
            # Now, for this shift, let's try no modification at all...if
            # our mass is not too large.  Many spectra, especially from LTQ
            # data-sets, have spurious +6 PTMs:
            if abs(CoreMass) < 10.0:
                Score = self.TryShiftedSite(Site, 0, DBPos, None, Shift)
                if Score > BestEditedSiteScore:
                    BestEditedSiteScore = Score
                    self.RememberOptimalKnownPTM(Site)
        ###############################################
        # Loop over modification SHIFTS, POSITIONS and MASSES is now complete.
        # Clean up memory usage:
        for Species in Site.SpeciesList:
            Species.FreeCachedSpectra()
        # Edit the species bits, for output:
        for Species in Site.SpeciesList:
            while len(Species.Bits) <= FormatBits.KnownPTMSitePValue:
                Species.Bits.append("")
            Species.Bits[FormatBits.KnownPTMName] = Species.BestKnownPTMAnnotationName
            Species.Bits[FormatBits.KnownPTMAnnotation] = Species.BestKnownPTMAnnotation
            Species.Bits[FormatBits.KnownPTMScore] = str(Species.BestKnownPTMAnnotationScore)
            Species.Bits[FormatBits.KnownPTMSitePValue] = str(Species.BestKnownPTMAnnotationSitePValue)
        # Verbose output:
        if BestEditedSiteScore:
            # Let stdout know what we're up to:
            Species = Site.SpeciesList[0]
            self.KPTMVerbose.write("Result: PValue %.3f (versus %.3f)\n"%(Species.BestKnownPTMAnnotationSitePValue, Site.PValue))
            ScoreDiff = -math.log(Species.BestKnownPTMAnnotationSitePValue) + math.log(Site.PValue)
            self.KPTMVerbose.write("==>Score change: %s\n"%ScoreDiff)
            for Species in Site.SpeciesList:
                self.KPTMVerbose.write("  %s (original)\n"%Species.Annotation)
                self.KPTMVerbose.write("  %s (%s)\n"%(Species.BestKnownPTMAnnotation, Species.BestKnownPTMAnnotationName))
                self.KPTMVerbose.write("Score %s (vs %s)\n"%(Species.BestKnownPTMAnnotationScore, Species.ModelScore))
        ###############################################
        # And now, we can output the site:
        for Species in Site.SpeciesList:
            Str = string.join(Species.Bits, "\t")
            self.OutputFile.write(Str + "\n")
    def TryShiftedSite(self, Site, ModMass, ModDBPos, KnownPTM, Shift):
        """
        Try shifting this modification site to the specified database position
        and mass.  Consider the effects on the peptide-score of each peptide.
        Return the resulting site-score.
        The value of Shift can be None, N-1, N-2, N+1, N+2, C-1, C-2, C+1, C+2.
        """
        SitePValue = 1.0
        SortedSpeciesScores = []
        #print "try shift to %+d on %s%s shift %s (%s)"%(ModMass, self.DB[ModDBPos], ModDBPos, Shift, KnownPTM.Name)
        for Species in Site.SpeciesList:
            if Species.Charge > 2:
                Model = self.Model3
            else:
                Model = self.Model2
            # Default "null" values:
            Species.KnownPTMAnnotation = ""
            Species.KnownPTMAnnotationScore = ""
            Species.KnownPTMAnnotationPValue = ""
            Species.KnownPTMAnnotationName = ""
            DBStart = Species.DBPos
            DBEnd = Species.DBEnd
            if Shift == "N-1":
                DBStart -= 1
            elif Shift == "N-2":
                DBStart -= 2
            elif Shift == "N+1":
                DBStart += 1
            elif Shift == "N+2":
                DBStart += 2
            elif Shift == "C-1":
                DBEnd -= 1
            elif Shift == "C-2":
                DBEnd -= 2
            elif Shift == "C+1":
                DBEnd += 1
            elif Shift == "C+2":
                DBEnd += 2
            if DBEnd <= ModDBPos or DBStart > ModDBPos:
                continue
            if KnownPTM and KnownPTM.Terminus in ("N", "^") and ModDBPos != DBStart:
                continue
            EditedFeatures = Species.Features[:]
            BestSpectrum = Species.GetBestSpectrum(self)
            NewAminos = self.DB[DBStart:DBEnd]
            NewPrefix = self.DB[DBStart - 1]
            NewSuffix = self.DB[DBEnd]
            ModPos = ModDBPos - DBStart
            if ModMass == 0:
                NewAnnotation = NewAminos
            else:
                NewAnnotation = "%s%+d%s"%(NewAminos[:ModPos + 1], ModMass, NewAminos[ModPos + 1:])
            NewAnnotation = "%s.%s.%s"%(NewPrefix, NewAnnotation, NewSuffix)
            if NewAnnotation == Species.Annotation:
                # Shortcut - the annotation hasn't changed, so neither will the score!
                if Species.Charge > 2:
                    PValue = self.Model3.GetPValue(Species.ModelScore)
                else:
                    PValue = self.Model2.GetPValue(Species.ModelScore)
                SortedSpeciesScores.append(Species.ModelScore)
                SitePValue *= PValue
                # Store these temp values; if the site is GOOD then we'll edit the species:
                Species.KnownPTMAnnotation = NewAnnotation
                if KnownPTM:
                    Species.KnownPTMAnnotationName = KnownPTM.Name
                else:
                    Species.KnownPTMAnnotationName = "unmodified"
                Species.KnownPTMAnnotationScore = Species.ModelScore
                Species.KnownPTMAnnotationPValue = PValue
                continue
            # Best spectrum score:
            Tuple = BestSpectrum.ScorePeptideDetailed(NewAnnotation)
            EditedFeatures[FeatureBits.BestMQScore] = Tuple[0]
            #EditedFeatures[FeatureBits.BestMQScore] = BestSpectrum.ScorePeptide(NewAnnotation)
            # Delta-score:
            ScoreDiff = EditedFeatures[FeatureBits.BestMQScore] - Species.Features[FeatureBits.BestMQScore]
            EditedFeatures[FeatureBits.BestDeltaScore] += ScoreDiff
            # Consensus spectrum score:
            ConsensusSpectrum = Species.GetConsensusSpectrum(self)
            ScoreInfo = ConsensusSpectrum.ScorePeptideDetailed(NewAnnotation)
            EditedFeatures[FeatureBits.ConsensusMQScore] = ScoreInfo[0]
            ScoreDiff = EditedFeatures[FeatureBits.ConsensusMQScore] - Species.Features[FeatureBits.ConsensusMQScore]
            EditedFeatures[FeatureBits.DeltaVsBigDB] += ScoreDiff
            EditedFeatures[FeatureBits.PeptideLength] = ScoreInfo[1]
            EditedFeatures[FeatureBits.TotalCutScore] = ScoreInfo[2]
            EditedFeatures[FeatureBits.MedianCutScore] = ScoreInfo[3]
            EditedFeatures[FeatureBits.YPresent] = ScoreInfo[4]
            EditedFeatures[FeatureBits.BPresent] = ScoreInfo[5]
            EditedFeatures[FeatureBits.BYIntensity] = ScoreInfo[6]
            EditedFeatures[FeatureBits.NTT] = ScoreInfo[7]
##            EditedFeatures[FeatureBits.PRMScore] = ScoreInfo[1]
##            EditedFeatures[FeatureBits.BYPresence] = ScoreInfo[2]
##            EditedFeatures[FeatureBits.TopPeakExplanation] = ScoreInfo[3]
##            EditedFeatures[FeatureBits.NTT] = ScoreInfo[4]
            ModTypeKey = (self.DB[ModDBPos], ModMass)
            EditedFeatures[FeatureBits.SpectraThisModType] = self.ModTypeSpectrumCount.get(ModTypeKey, 1)
            EditedFeatures[FeatureBits.SitesThisModType] = self.ModTypeSiteCount.get(ModTypeKey, 1)
            EditedFeatures[FeatureBits.LogSpecThisType] = math.log(EditedFeatures[FeatureBits.SpectraThisModType])
            EditedFeatures[FeatureBits.LogSitesThisType] = math.log(EditedFeatures[FeatureBits.SitesThisModType])
            # Spectral similarity:
            # Spectral similarity:
            try:
                SisterAnnotationFlag = int(Species.Bits[FormatBits.SisterAnnotationFlag])
            except:
                SisterAnnotationFlag = 0
            if SisterAnnotationFlag:
                try:
                    ConsensusMSSpectrum = Species.GetConsensusMSSpectrum(self)
                    ModlessMSSpectrum = Species.GetConsensusModlessMSSpectrum(self)
                    Comparator = SpectralSimilarity.SpectralSimilarity(ConsensusMSSpectrum,
                        ModlessMSSpectrum, NewAnnotation, Species.ModlessAnnotation)
                    # COPIED from ComputePTMFeatures:
                    Comparator.LabelPeaks(0.5)
                    Similarity = Comparator.DotProduct(0.5)
                    EditedFeatures[FeatureBits.Dot] = Similarity
                    Similarity = Comparator.GetSharedPeakCount(0, 1)
                    EditedFeatures[FeatureBits.Shared01] = Similarity
                    Similarity = Comparator.GetSharedPeakCount(1, 1)
                    EditedFeatures[FeatureBits.Shared11] = Similarity
                    CorrelationCoefficient = Comparator.ComputeCorrelationCoefficient(1.0)
                    EditedFeatures[FeatureBits.Correlation] = Similarity
                except:
                    traceback.print_exc()
                    print "*** Unable to generate spectral-similarity features; continuing."
                    print Site, ModMass, ModDBPos, KnownPTM, Shift
                    print Species.Annotation, NewAnnotation
                    print "BITS:", Species.Bits
                    print "SisterAnnotationFlag:", SisterAnnotationFlag
            # Features are complete - score the altered peptide!
            ModelScore = self.ScoreInstance(Model, EditedFeatures)
##            # TEMP: Verbose output
##            self.KPTMVerbose.write("%s -> %s\n"%(Species.Annotation, NewAnnotation))
##            for Index in range(len(EditedFeatures)):
##                self.KPTMVerbose.write("%s: %.2f %.2f (%.4f)\n"%(Index, Species.Features[Index], EditedFeatures[Index], EditedFeatures[Index] - Species.Features[Index]))
##            self.KPTMVerbose.write("Score: %s versus old %s\n"%(ModelScore, Species.ModelScore))
            PValue = Model.GetPValue(ModelScore)
            SortedSpeciesScores.append(ModelScore)
            SitePValue *= PValue
            # Store these temp values; if the site is GOOD then we'll edit the species:
            Species.KnownPTMAnnotation = NewAnnotation
            if KnownPTM:
                Species.KnownPTMAnnotationName = KnownPTM.Name
            else:
                Species.KnownPTMAnnotationName = "Unmodified"
            Species.KnownPTMAnnotationScore = ModelScore
            Species.KnownPTMAnnotationPValue = PValue
        SortedSpeciesScores.sort()
        SortedSpeciesScores.reverse()
        for Species in Site.SpeciesList:
            Species.KnownPTMAnnotationSitePValue = SitePValue
        SiteScore = [-math.log(SitePValue)]
        SiteScore.extend(SortedSpeciesScores)
        if len(SortedSpeciesScores):
            self.KPTMVerbose.write("%s%+d on %s (%s): score (%.3f, %.3f)\n"%(self.DB[ModDBPos], ModMass, ModDBPos, Shift, SiteScore[0], SiteScore[1]))
        return SiteScore
    def ScoreInstance(self, Model, Features):
        NiceFeatures = []
        #for FeatureIndex in TrainPTMFeatures.ValidFeatureIndices:
        #    NiceFeatures.append(Features[FeatureIndex])
        return Model.ScoreInstance(Features)
    def RememberOptimalKnownPTM(self, Site):
        for Species in Site.SpeciesList:
            Species.BestKnownPTMAnnotation = Species.KnownPTMAnnotation
            Species.BestKnownPTMAnnotationScore = Species.KnownPTMAnnotationScore
            Species.BestKnownPTMAnnotationPValue = Species.KnownPTMAnnotationPValue
            Species.BestKnownPTMAnnotationSitePValue = Species.KnownPTMAnnotationSitePValue
            Species.BestKnownPTMAnnotationName = Species.KnownPTMAnnotationName
    def LoadCoverageLevels(self):
        "Load peptide coverage levels, written by ComputePTMFeatures"
        CoveragePath = os.path.join(self.TempFileDir, "Coverage.dat")
        try:
            CoverageFile = open(CoveragePath, "rb")
        except:
            traceback.print_exc()
            print "** WARNING: Coverage levels not found at '%s'"%CoveragePath
            self.Coverage = [1] * len(self.DB)
            self.ModCoverage = [1] * len(self.DB)
            return
        self.Coverage = []
        self.ModCoverage = []
        BlokSize = struct.calcsize("<II")
        for DBPos in range(len(self.DB)):
            Blok = CoverageFile.read(BlokSize)
            Tuple = struct.unpack("<II", Blok)
            self.Coverage.append(Tuple[0])
            self.ModCoverage.append(Tuple[1])
        CoverageFile.close()
    def SaveCoverageLevels(self):
        "Save peptide coverage levels, which may have changed during merge+reconcile"
        Dir = self.TempFileDir #os.path.split(self.OutputFileName)[0]
        if self.DBEnd != None:
            CoveragePath = os.path.join(Dir, "AdjustedCoverage.%s.%s.dat"%(self.DBStart, self.DBEnd))
        else:
            CoveragePath = os.path.join(Dir, "AdjustedCoverage.dat")
        CoverageFile = open(CoveragePath, "wb")
        for DBPos in range(len(self.DB)):
            if self.Coverage[DBPos] < 0 or self.Coverage[DBPos] >= 65535:
                print "* Coverage of %s is %s"%(DBPos, self.Coverage[DBPos])
            if self.ModCoverage[DBPos] < 0 or self.ModCoverage[DBPos] >= 65535:
                print "* ModCoverage of %s is %s"%(DBPos, self.ModCoverage[DBPos])
            Str = struct.pack("<II", self.Coverage[DBPos], self.ModCoverage[DBPos])
            CoverageFile.write(Str)
        CoverageFile.close()
    def MergeAndReconcile(self):
        """
        Iterate over peptide species, from best to worst.  Two iterations: In the first, we
        consider MERGING a peptide with a superior one.  In the second iteration, we consider
        RECONCILING each species to a superior one.
        """
        self.LoadCoverageLevels()
        #self.ParseOriginalSpectraForModType(self.InputFileName)
        if self.DBEnd != None:
            MergeFileName = "MergeDetails.%s.%s.txt"%(self.DBStart, self.DBEnd)
            ReconcileFileName = "ReconcileDetails.%s.%s.txt"%(self.DBStart, self.DBEnd)
        else:
            MergeFileName = "MergeDetails.txt"
            ReconcileFileName = "ReconcileDetails.txt"
        OutputDir = os.path.split(self.OutputFileName)[0]
        self.MergeDetailOutput = open(os.path.join(OutputDir, MergeFileName), "wb")
        Header = "Flag\tCharge\tMaster\tServant\tMasterMQScore\tServantMQScore\tServantRescore\tScoreChange\tNewConsScore\tOldModelScore\tNewModelScore\t"
        self.MergeDetailOutput.write(Header + "\n")
        self.ReconcileDetailOutput = open(os.path.join(OutputDir, ReconcileFileName), "wb")
        # Flow for file lines:
        # Initial input file -> MergeTemp -> ReconcileTemp -> output file
        ########################################################
        # FIRST cycle through points: Consider merging.
        self.InputFile = open(self.InputFileName, "rb")
        MergeTempPath = "%s.merge"%self.OutputFileName
        self.OutputFile = open(MergeTempPath, "wb")
        print ">>>PerformAllMerges 1: Read from %s, write to %s"%(self.InputFileName, MergeTempPath)
        self.PerformAllMerges(1)
        self.InputFile.close()
        self.OutputFile.close()
        ########################################################
        # SECOND cycle through points: Consider conforming.
        print "\n\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
        print "Reconcile:"
        self.InputFile = open(MergeTempPath, "rb")
        ConformTempPath = "%s.conform"%self.OutputFileName
        self.OutputFile = open(ConformTempPath, "wb")
        print ">>>PerformAllMerges 0: Read from %s, write to %s"%(MergeTempPath, ConformTempPath)
        self.PerformAllMerges(0)
        self.InputFile.close()
        self.OutputFile.close()
        self.SaveCoverageLevels()
        if self.DBStart != None:
            # We're handling just one block.  Therefore, we shouldn't re-score the peptides:
            return 
        ########################################################
        # At this point, we know the number of sites and spectra for each
        # modification type.  We need that information for when we
        # consider changing to known PTMs.  Let's pickle it.
        self.SaveSitesByType()
        ########################################################
        # THIRD cycle through points: Update sites-per-mod and spectra-per-mod,
        # and accumulate feature vectors for the model.
        #self.InputFile = open(ConformTempPath, "rb")
        self.ParseFeatureVectors(ConformTempPath)
        #self.InputFile.close()
        ########################################################
        # FOURTH cycle through points: Write the revised score!
        self.OutputFile = open(self.OutputFileName, "wb")
        for HeaderLine in self.HeaderLines:
            self.OutputFile.write(HeaderLine)
        self.ProcessSites(ConformTempPath, "rescore")
    def MergeBlockRuns(self):
        """
        Combine the output of several AdjustPTM runs for sub-blocks of the database.  
        """
        Directory = os.path.split(self.OutputFileName)[0]
        self.CombineBlockCoverage(self.TempFileDir)
        self.SaveCoverageLevels()
        # Populate self.ModTypeSpectrumCount and self.ModTypeSiteCount:
        self.LoadModSitesByTypeBlocks(Directory)
        self.SaveSitesByType()
        # Concatenate block files into one large file:
        ConcatenatedFileName = os.path.join(Directory, "ConcatenatedFeatures.txt")
        self.ConcatenateBlockOutputFiles(Directory, ConcatenatedFileName)
        # Parse feature-vectors, train model, output model:
        self.ParseFeatureVectors(ConcatenatedFileName)
        # Rescore:
        self.OutputFile = open(self.OutputFileName, "wb")
        for HeaderLine in self.HeaderLines:
            self.OutputFile.write(HeaderLine)
        self.ProcessSites(ConcatenatedFileName, "rescore")
        self.OutputFile.close()
    def ConcatenateBlockOutputFiles(self, Directory, OutputFileName):
        """
        Concatenate the merge-and-reconcile output files from various blocks of the database.
        """
        OutputFile = open(OutputFileName, "wb")
        FirstFileFlag = 1
        for FileName in os.listdir(Directory):
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension != ".conform":
                continue
            FilePath = os.path.join(Directory, FileName)
            File = open(FilePath, "rb")
            for FileLine in File.xreadlines():
                if FileLine[0] == "#":
                    # Header line.  Write it out iff this is the first file
                    if FirstFileFlag:
                        OutputFile.write(FileLine)
                else:
                    OutputFile.write(FileLine)
            File.close()
            FirstFileFlag = 0
        print "Concatenated results from '%s' into '%s'"%(Directory, OutputFileName)
    def LoadModSitesByTypeBlocks(self, Directory):
        """
        Iterate over block output files from this directory.  Populate
        self.ModTypeSiteCount and self.ModTypeSpectrumCount based on the contents.
        """
        self.ModTypeSpectrumCount = {}
        self.ModTypeSiteCount = {}        
        for FileName in os.listdir(Directory):
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension != ".conform":
                continue
            Path = os.path.join(Directory, FileName)
            File = open(Path, "rb")
            print "Read SitesByType from %s..."%Path
            for FileLine in File.xreadlines():
                Bits = FileLine.split("\t")
                if FileLine[0] == "#":
                    continue
                AA = Bits[FormatBits.ModifiedAA]
                Mass = int(Bits[FormatBits.ModificationMass])
                Spectra = int(float(Bits[FormatBits.SpectrumCount]))
                Key = (AA, Mass)
                self.ModTypeSiteCount[Key] = self.ModTypeSiteCount.get(Key, 0) + 1
                self.ModTypeSpectrumCount[Key] = self.ModTypeSpectrumCount.get(Key, 0) + Spectra
    def CombineBlockCoverage(self, Directory):
        # Load the original coverage, just for reference:
        self.LoadCoverageLevels()
        # Iterate over coverage output files from the individual blocks:
        for FileName in os.listdir(Directory):
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension != ".dat":
                continue
            Bits = FileName.split(".")
            # Names have the form AdjustedCoverage.START.END.dat
            if len(Bits) < 4:
                continue
            DBStart = int(Bits[1])
            DBEnd = int(Bits[2])
            Path = os.path.join(Directory, FileName)
            print "Read block coverage from %s..."%Path
            CoverageFile = open(Path, "rb")
            BlokSize = struct.calcsize("<II")
            for DBPos in range(len(self.DB)):
                Blok = CoverageFile.read(BlokSize)
                Tuple = struct.unpack("<II", Blok)
                if DBPos < DBStart or DBPos >= DBEnd:
                    continue
                self.Coverage[DBPos] = Tuple[0]
                self.ModCoverage[DBPos] = Tuple[1]
            CoverageFile.close()
    def ProcessSites(self, InputFileName, Command):
        """
        Parse modification-sites from the input file.  Once all the peptides
        for a site have been read, execute the command.
        """
        if Command == "knownptm":
            if self.KnownPTMVerboseOutputFileName:
                self.KPTMVerbose = open(self.KnownPTMVerboseOutputFileName, "wb")
            else:
                self.KPTMVerbose = sys.stdout
        CurrentSite = None
        InputFile = open(InputFileName, "rb")
        for FileLine in InputFile.xreadlines():
            FileLine = FileLine.replace("\r", "").replace("\n", "")
            if not FileLine:
                continue  # skip blank lines
            Bits = FileLine.split("\t")
            if FileLine[0] == "#" or len(Bits) < 2:
                continue
            try:
                Species = PeptideSpeciesClass()
                Species.ParseBits(Bits)
            except:
                traceback.print_exc()
                print Bits
                continue
            if CurrentSite == None or Species.ModDBPos != CurrentSite.ModDBPos or Species.ModMass != CurrentSite.ModMass:
                # Finish the previous (if any) site:
                if CurrentSite:
                    if Command == "rescore":
                        self.RescoreAndWriteSite(CurrentSite)
                    elif Command == "knownptm":
                        self.AttemptKnownPTM(CurrentSite)
                CurrentSite = SiteClass()
                CurrentSite.ModDBPos = Species.ModDBPos
                CurrentSite.ModMass = Species.ModMass
            # Add a species to the current site:
            CurrentSite.SpeciesList.append(Species)
        InputFile.close()
        # Finish the last site:
        if CurrentSite:
            if Command == "rescore":
                self.RescoreAndWriteSite(CurrentSite)
            elif Command == "knownptm":
                self.AttemptKnownPTM(CurrentSite)
    def ParseFeatureVectors(self, FileName):
        """
        Called after merge and reconcile.  Read feature-vectors, updating
        the spectrum/site counts for modification type, and obtain scores!
        """
        FeatureSet2 = Learning.FeatureSetClass()
        FeatureSet3 = Learning.FeatureSetClass()
        File = open(FileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if FileLine[0] == "#":
                continue
            if len(Bits) < 2:
                continue
            Charge = int(Bits[FormatBits.Charge])
            Vector = Learning.FeatureVector()
            Vector.Features = [] 
            for BitIndex in range(FormatBits.FirstFeature, FormatBits.LastFeature + 1):
                try:
                    Vector.Features.append(float(Bits[FeatureIndex]))
                except:
                    Vector.Features.append(0)
            # Tweak spectra-by-type and sites-by-type:
            ModTypeKey = (Bits[FormatBits.ModifiedAA], int(Bits[FormatBits.ModificationMass]))
            TotalSpectra = self.ModTypeSpectrumCount.get(ModTypeKey, 0)
            TotalSites = self.ModTypeSiteCount.get(ModTypeKey, 0)
            Vector.Features[FeatureBits.SpectraThisModType] = TotalSpectra
            Vector.Features[FeatureBits.SitesThisModType] = TotalSites
            Vector.Features[FeatureBits.LogSpecThisType] = math.log(TotalSpectra)
            Vector.Features[FeatureBits.LogSitesThisType] = math.log(TotalSites)
            if Charge > 2:
                FeatureSet = FeatureSet3
            else:
                FeatureSet = FeatureSet2
            if int(Bits[FormatBits.TrueProteinFlag]):
                Vector.TrueFlag = 1
                FeatureSet.TrueVectors.append(Vector)
            else:
                Vector.TrueFlag = 0
                FeatureSet.FalseVectors.append(Vector)
            FeatureSet.AllVectors.append(Vector)
        File.close()
        FeatureSet2.SetCounts()
        FeatureSet3.SetCounts()
        self.Model2.Test(FeatureSet2)
        self.Model3.Test(FeatureSet3)
        if self.OutputModelFileName2:
            self.Model2.SaveModel(self.OutputModelFileName2)
            self.Model3.SaveModel(self.OutputModelFileName3)
    def LoadCluster(self, Path):
        Builder = BuildConsensusSpectrum.ConsensusBuilder()
        Builder.UnpickleCluster(Path)
        return Builder
    def AttemptMerge(self, SpeciesA, SpeciesB, BWeak = 0):
        """
        Attermpt a merge of SpeciesA into SpeciesB.  If the merge is valid,
        perform the merge, set SpeciesA.MergedFlag, and RETURN TRUE.
        """
        print "AttemptMerge chg%s %s into %s"%(SpeciesA.Charge, SpeciesA.Annotation, SpeciesB.Annotation)
        ################################################
        # Condition 1: Consensus A seems to match annotation B reasonably well
        if BWeak:
            # Species A has the stronger score, so don't screw it up!
            ScoreLossLimit = 0.1
        else:
            ScoreLossLimit = 3
        SpectrumA = SpeciesA.GetConsensusSpectrum(self)
        Score = SpectrumA.ScorePeptide(SpeciesB.Annotation)
        DetailStr = "%s\t%s\t%s\t%s\t%s\t"%(SpeciesA.Charge, SpeciesB.Annotation, SpeciesA.Annotation, SpeciesB.Features[FeatureBits.ConsensusMQScore], SpeciesA.Features[FeatureBits.ConsensusMQScore])
        DetailStr += "%s\t%s\t"%(Score, SpeciesA.Features[FeatureBits.ConsensusMQScore] - Score)
        if (SpeciesA.Features[FeatureBits.ConsensusMQScore] - Score) > ScoreLossLimit:
            DetailStr = "FailAScore\t"+DetailStr
            self.MergeDetailOutput.write(DetailStr + "\n")
            return 0
        ################################################
        # Condition 2: The merged consensus spectrum is not significantly WORSE.
        # Load in ClusterA and ClusterB (we cache the species-A cluster)
        ClusterPathA = os.path.join(self.ConsensusClusterDir, SpeciesA.Annotation[2], "%s.%s.cls"%(SpeciesA.Annotation, SpeciesA.Charge))
        if ClusterPathA == self.CachedClusterPath:
            ClusterA = self.CachedCluster
        else:
            self.CachedCluster = self.LoadCluster(ClusterPathA)
            self.CachedClusterPath = ClusterPathA
            ClusterA = self.CachedCluster
        ClusterPathB = os.path.join(self.ConsensusClusterDir, SpeciesB.Annotation[2], "%s.%s.cls"%(SpeciesB.Annotation, SpeciesB.Charge))
        ClusterB = self.LoadCluster(ClusterPathB)
        # If we combine these two clusters into a single consensus
        # spectrum, what sort of MQScore do we end up with?
        TempConsensusPath = os.path.join(self.ConsensusSpectraDir, "Consensus.dta")
        ClusterB.AssimilateCluster(ClusterA)
        NewConsensusSpectrum = ClusterB.ProduceConsensusSpectrum()
        NewConsensusSpectrum.WritePeaks(TempConsensusPath)
        # Set the file members of the spectrum, since Label.py reads them:
        NewConsensusSpectrum.FilePath = TempConsensusPath
        NewConsensusSpectrum.FilePos = 0
        NewConsensusSpectrum.RankPeaksByIntensity()
        PySpectrum = PyInspect.Spectrum(TempConsensusPath, 0)
        ScoreInfo = PySpectrum.ScorePeptideDetailed(SpeciesB.Annotation)
        DetailStr += "%s\t"%ScoreInfo[0]
        if SpeciesB.Features[FeatureBits.ConsensusMQScore] - ScoreInfo[0] > 2:
            DetailStr = "FailBConsensus\t"+DetailStr
            self.MergeDetailOutput.write(DetailStr + "\n")
            return 0
        NewFeatures = SpeciesB.Features[:]
        # Update feature values for the merged guy:
        # Spectra:
        NewFeatures[FeatureBits.SpectrumCount] = SpeciesA.Features[FeatureBits.SpectrumCount] + SpeciesB.Features[FeatureBits.SpectrumCount]
        NewFeatures[FeatureBits.LogSpectrumCount] = math.log(NewFeatures[FeatureBits.SpectrumCount])
        NewFeatures[FeatureBits.ModlessSpectrumCount] = SpeciesB.Features[FeatureBits.ModlessSpectrumCount] # unchanged
        # BestMQ, BestDelta, PeptideCount:
        BestSpectrumA = SpeciesA.GetBestSpectrum(self)
        NewBestMQA = BestSpectrumA.ScorePeptide(SpeciesB.Annotation)
        NewBestDeltaA = SpeciesA.Features[3] + (NewBestMQA - SpeciesA.Features[2])
        print "A best MQ %.4f (was %.4f) delta %.4f (was %.4f)"%(NewBestMQA, SpeciesA.Features[2],
            NewBestDeltaA, SpeciesA.Features[3])
        # Best MQScore:
        NewFeatures[FeatureBits.BestMQScore] = max(NewBestMQA, SpeciesB.Features[FeatureBits.BestMQScore])
        # BestDelta:
        NewFeatures[FeatureBits.BestDeltaScore] = max(NewBestDeltaA, SpeciesB.Features[FeatureBits.BestDeltaScore])
        # Peptide length:
        NewFeatures[FeatureBits.PeptideLength] = SpeciesB.Features[FeatureBits.PeptideLength]
        # Peptide count:
        NewFeatures[FeatureBits.PeptideCount] = SpeciesB.Features[FeatureBits.PeptideCount]
        # Consensus scoring (Score, and score components):
        NewFeatures[FeatureBits.ConsensusMQScore] = ScoreInfo[0]
        NewFeatures[FeatureBits.PeptideLength] = ScoreInfo[1]
        NewFeatures[FeatureBits.TotalCutScore] = ScoreInfo[2]
        NewFeatures[FeatureBits.MedianCutScore] = ScoreInfo[3]
        NewFeatures[FeatureBits.YPresent] = ScoreInfo[4]
        NewFeatures[FeatureBits.BPresent] = ScoreInfo[5]
        NewFeatures[FeatureBits.BYIntensity] = ScoreInfo[6]
        NewFeatures[FeatureBits.NTT] = ScoreInfo[7]
        # Adjust delta-score by the difference in consensus-mq-score:
        NewFeatures[FeatureBits.DeltaVsBigDB] = SpeciesB.Features[FeatureBits.DeltaVsBigDB] + (ScoreInfo[0] - SpeciesB.Features[FeatureBits.ConsensusMQScore])
        # Similarity scores for the new consensus spectrum:
        if SpeciesB.Bits[FormatBits.SisterAnnotationFlag]:
            ModlessSpectrum = SpeciesB.GetConsensusModlessMSSpectrum(self)
            #print "Build comparator:"
            Comparator = SpectralSimilarity.SpectralSimilarity(NewConsensusSpectrum,
                ModlessSpectrum, SpeciesB.Annotation, SpeciesB.ModlessAnnotation)
            #print "Label peaks:"
            # COPIED from ComputePTMFeatures:
            Comparator.LabelPeaks(0.5)
            #print "Compute:"
            Similarity = Comparator.DotProduct(0.5)
            NewFeatures[FeatureBits.Dot] = Similarity
            Similarity = Comparator.GetSharedPeakCount(0, 1)
            NewFeatures[FeatureBits.Shared01] = Similarity
            Similarity = Comparator.GetSharedPeakCount(1, 1)
            NewFeatures[FeatureBits.Shared11] = Similarity
            CorrelationCoefficient = Comparator.ComputeCorrelationCoefficient(1.0)
            NewFeatures[FeatureBits.Correlation] = Similarity
        # Ask the trained model what it thinks of this new feature-vector:
        if SpeciesA.Charge > 2:
            Model = self.Model3
        else:
            Model = self.Model2
        NewModelScore = self.ScoreInstance(Model, NewFeatures)
        PValue = Model.GetPValue(NewModelScore)
        OldPValue = Model.GetPValue(SpeciesB.ModelScore)
        print "Score of the NEW FEATURES: %.3f (%.1f%%) versus %.3f (%.1f%%) old)"%(NewModelScore, 100 * PValue, SpeciesB.ModelScore, 100 * OldPValue)
        DetailStr += "%s\t%s\t"%(SpeciesB.ModelScore, NewModelScore)
        ################################################
        # Condition 3: Model score should not be dramatically worse!        
        if NewModelScore < SpeciesB.ModelScore - 0.5:
            DetailStr = "FailModelScore\t%s"%DetailStr
            self.MergeDetailOutput.write(DetailStr + "\n")
            return 0
        DetailStr = "MERGE\t" + DetailStr
        print "MERGE %s and %s"%(SpeciesA.Annotation, SpeciesB.Annotation)
        self.MergeDetailOutput.write(DetailStr + "\n")
        #######################################################
        # ALL CONDITIONS PASSED, NOW LET'S MERGE:
        SpeciesA.MergedFlag = 1 # this species won't be written out.
        # Remember the new "best spectrum", if it belonged to A:
        if (NewBestMQA > SpeciesB.Features[2]):
            SpeciesB.ConsensusMQScore = NewBestMQA
            SpeciesB.Bits[FormatBits.BestSpectrumPath] = SpeciesA.Bits[FormatBits.BestSpectrumPath]
        SpeciesB.Features = NewFeatures
        ############################################
        # Update our COVERAGE and MODDED-FRACTION:
        for DBPos in range(SpeciesA.DBPos, SpeciesA.DBEnd):
            if SpeciesA.Peptide.Modifications.keys():
                self.ModCoverage[DBPos] -= int(SpeciesA.Features[FeatureBits.SpectrumCount])
            else:
                self.Coverage[DBPos] -= int(SpeciesA.Features[FeatureBits.SpectrumCount])
        for DBPos in range(SpeciesB.DBPos, SpeciesB.DBEnd):
            if SpeciesB.Peptide.Modifications.keys():
                self.ModCoverage[DBPos] += int(SpeciesB.Features[FeatureBits.SpectrumCount])
            else:
                self.Coverage[DBPos] += int(SpeciesB.Features[FeatureBits.SpectrumCount])
        SpeciesB.Features[FeatureBits.ModdedFraction] = self.ModCoverage[SpeciesB.ModDBPos] / float(self.ModCoverage[SpeciesB.ModDBPos] + self.Coverage[SpeciesB.ModDBPos])
        ############################################    
        # Write the adjusted consensus cluster:
        ClusterPathB = os.path.join(self.ConsensusClusterDirAdjusted, SpeciesB.Annotation[2], "%s.%s.cls"%(SpeciesB.Annotation, SpeciesB.Charge))
        ClusterB.PickleCluster(ClusterPathB)
        SpeciesB.ConsensusMSSpectrum = NewConsensusSpectrum
        SpeciesB.ConsensusSpectrum = PySpectrum
        SpeciesB.ModelScore = NewModelScore
        # Write the adjusted consensus spectrum:
        ConsensusPath = os.path.join(self.ConsensusSpectraDirAdjusted, SpeciesB.Annotation[2], "%s.%s.dta"%(SpeciesB.Annotation, SpeciesB.Charge))
        NewConsensusSpectrum.WritePeaks(ConsensusPath)
        # Write the merged list of members:
        MemberListStr = ""
        try:
            MemberListStr += SpeciesA.GetMemberListStr(self)
            MemberListStr += SpeciesB.GetMemberListStr(self)
        except:
            print "* ERROR!"
            print SpeciesA
            print SpeciesB
            raise
        Path = os.path.join(self.ClusterScanListDirAdjusted, SpeciesB.Annotation[2], "%s.%s.txt"%(SpeciesB.Annotation, SpeciesB.Charge))
        ClusterMemberFile = open(Path, "wb")
        ClusterMemberFile.write(MemberListStr)
        ClusterMemberFile.close()
        # Remove peptide A from self.PeptideDict:
        Key = (SpeciesA.Annotation, SpeciesA.Charge)
        try:
            del self.PeptideDict[Key]
        except:
            pass
        return 1
    def AttemptReconcileFixEndpoints(self, SpeciesA, SpeciesB, OldDBPos, OldDBEnd):
        # Adjust endpoints, if necessary for reconciliation:
        if SpeciesA.ModMass > SpeciesB.ModMass:
            # Species A has a larger modification.  Maybe we can ADD 1-2 residues
            # and make the modification mass equal?
            # Try shifting N-terminus:
            ExtraMass = Global.AminoMass.get(self.DB[SpeciesA.DBPos - 1], 999999)
            FullMass = int(round(SpeciesA.ModMass - ExtraMass))
            if FullMass == SpeciesB.ModMass:
                return (OldDBPos - 1, OldDBEnd)
            if FullMass < SpeciesB.ModMass:
                ExtraMass += Global.AminoMass.get(self.DB[SpeciesA.DBPos - 2], 999999)
                FullMass = int(round(SpeciesA.ModMass - ExtraMass))
                if FullMass == SpeciesB.ModMass:
                    return (OldDBPos - 2, OldDBEnd)
            # Try shifting C-terminus:
            ExtraMass = Global.AminoMass.get(self.DB[SpeciesA.DBEnd], 999999)
            FullMass = int(round(SpeciesA.ModMass - ExtraMass))
            if abs(FullMass - SpeciesB.ModMass) < 2:
                return (OldDBPos, OldDBEnd + 1)
            if FullMass < SpeciesB.ModMass:
                ExtraMass += Global.AminoMass.get(self.DB[SpeciesA.DBEnd + 1], 999999)
                FullMass = int(round(SpeciesA.ModMass - ExtraMass))
                if abs(FullMass - SpeciesB.ModMass) < 2:
                    return (OldDBPos, OldDBEnd + 2)
        if SpeciesA.ModMass < SpeciesB.ModMass:
            # Species A has a smaller modification.  Maybe we can REMOVE 1-2 residues
            # and make the modification mass equal?
            # Try shifting N-terminus:
            ExtraMass = Global.AminoMass[self.DB[SpeciesA.DBPos]]
            FullMass = int(round(SpeciesA.ModMass + ExtraMass))
            if abs(FullMass - SpeciesB.ModMass) < 2:
                return (OldDBPos + 1, OldDBEnd)
            if FullMass > SpeciesB.ModMass:
                ExtraMass += Global.AminoMass[self.DB[SpeciesA.DBPos + 1]]
                FullMass = int(round(SpeciesA.ModMass + ExtraMass))
                if abs(FullMass - SpeciesB.ModMass) < 2:
                    return (OldDBPos + 2, OldDBEnd)
            # Try shifting C-terminus:
            ExtraMass = Global.AminoMass[self.DB[SpeciesA.DBEnd - 1]]
            FullMass = int(round(SpeciesA.ModMass + ExtraMass))
            if abs(FullMass - SpeciesB.ModMass) < 2:
                return (OldDBPos, OldDBEnd - 1)
            if FullMass < SpeciesB.ModMass:
                ExtraMass += Global.AminoMass[self.DB[SpeciesA.DBEnd - 2]]
                FullMass = int(round(SpeciesA.ModMass + ExtraMass))
                if abs(FullMass - SpeciesB.ModMass) < 2:
                    return (OldDBPos, OldDBEnd - 2)
        return (OldDBPos, OldDBEnd)
    def AttemptReconcile(self, SpeciesA, SpeciesB):
        """
        Attempt to reconcile species A with species B.  In other words,
        edit annotation A so that it carries the same modification as B, and
        on the same database position.  If the effects on match quality score
        (and/or model score) are an IMPROVEMENT (or at least, not a big
        disappointment), then perform the reconciliation, and return TRUE.
        """
        OldDBPos = SpeciesA.DBPos
        OldDBEnd = SpeciesA.DBEnd
        (NewDBPos, NewDBEnd) = self.AttemptReconcileFixEndpoints(SpeciesA, SpeciesB, SpeciesA.DBPos, SpeciesA.DBEnd)
        NewPrefix = self.DB[NewDBPos - 1]
        NewSuffix = self.DB[NewDBEnd]
        ConsensusSpectrum = SpeciesA.GetConsensusSpectrum(self)
        ModIndex = SpeciesB.ModDBPos - NewDBPos
        ModdedAnnotation = "%s%+d%s"%(self.DB[NewDBPos:SpeciesB.ModDBPos + 1], SpeciesB.ModMass, self.DB[SpeciesB.ModDBPos + 1:NewDBEnd])
        NewAnnotation = "%s.%s.%s"%(NewPrefix, ModdedAnnotation, NewSuffix)
        NewAnnotation = NewAnnotation.replace("*", "-")
        NewConsensusScore = ConsensusSpectrum.ScorePeptide(NewAnnotation)
        ScoreDiff = NewConsensusScore - SpeciesA.ConsensusMQScore
        self.ReconcileDetailOutput.write("%s\t%s\t%s\t%s\t%s\t\n"%(SpeciesA.Annotation, NewAnnotation, SpeciesA.ConsensusMQScore, NewConsensusScore, ScoreDiff))
        if ScoreDiff < -0.5:
            return 0
        OldAnnotation = SpeciesA.Annotation
        NewPeptide = GetPeptideFromModdedName(NewAnnotation)
        NewModlessAnnotation = "%s.%s.%s"%(NewPeptide.Prefix, self.DB[NewDBPos:NewDBEnd], NewPeptide.Suffix)
        # Compute new features of the 'reconciled peptide':
        NewFeatures = SpeciesA.Features[:]
        # Best spectrum MQScore and Delta-score:
        PySpectrum = SpeciesA.GetBestSpectrum(self)
        NewBestMQ = PySpectrum.ScorePeptide(NewAnnotation)
        NewFeatures[FeatureBits.BestDeltaScore] += (NewBestMQ - NewFeatures[FeatureBits.BestMQScore])
        NewFeatures[FeatureBits.BestMQScore] = NewBestMQ
        PeptideLength = NewDBEnd - NewDBPos
        NewFeatures[FeatureBits.PeptideLength] = PeptideLength
        NewFeatures[FeatureBits.LogPeptideLength] = math.log(PeptideLength)
        # Consensus score:
        PySpectrum = SpeciesA.GetConsensusSpectrum(self)
        ScoreInfo = PySpectrum.ScorePeptideDetailed(NewAnnotation)
        NewFeatures[FeatureBits.ConsensusMQScore] = ScoreInfo[0]
        NewFeatures[FeatureBits.PeptideLength] = ScoreInfo[1]
        NewFeatures[FeatureBits.TotalCutScore] = ScoreInfo[2]
        NewFeatures[FeatureBits.MedianCutScore] = ScoreInfo[3]
        NewFeatures[FeatureBits.YPresent] = ScoreInfo[4]
        NewFeatures[FeatureBits.BPresent] = ScoreInfo[5]
        NewFeatures[FeatureBits.BYIntensity] = ScoreInfo[6]
        NewFeatures[FeatureBits.NTT] = ScoreInfo[7]
        # Adjust delta-score by the difference in consensus-mq-score:
        NewFeatures[FeatureBits.DeltaVsBigDB] = SpeciesA.Features[FeatureBits.DeltaVsBigDB] + (ScoreInfo[0] - SpeciesA.Features[FeatureBits.ConsensusMQScore])
        # Adjust spectra, sites for this modification type:
        NewFeatures[FeatureBits.SpectraThisModType] = SpeciesB.Features[FeatureBits.SpectraThisModType]
        NewFeatures[FeatureBits.SitesThisModType] = SpeciesB.Features[FeatureBits.SitesThisModType]
        NewFeatures[FeatureBits.LogSpecThisType] = SpeciesB.Features[FeatureBits.LogSpecThisType]
        NewFeatures[FeatureBits.LogSitesThisType] = SpeciesB.Features[FeatureBits.LogSitesThisType]
        # Modless spectra:
        # - If our endpoint didn't change, then we keep our old modless spectrum
        # - If our endpoints changed and we're assimilating a target (ExistingSpecies),
        #   then we inherit *its* modless spectra
        # - Otherwise, we LOSE our modless spectra!
        Key = (NewAnnotation, SpeciesA.Charge)
        ExistingSpecies = self.PeptideDict.get(Key, None)
        # Initialize:
        ModlessSpectrumFlag = 0
        BestModlessSpectrumPath = ""
        BestModlessMQScore = ""
        if OldDBPos == NewDBPos and OldDBEnd == NewDBEnd:
            try:
                ModlessSpectrumFlag = int(SpeciesA.Bits[FormatBits.SisterAnnotationFlag])
                BestModlessSpectrumPath = SpeciesA.Bits[FormatBits.BestModlessSpectrumPath]
                BestModlessMQScore = float(SpeciesA.Bits[FormatBits.BestModlessMQScore])
                ModlessMSSpectrum = SpeciesA.GetConsensusModlessMSSpectrum(self)
            except:
                pass # the modless-bits weren't set; that's fine
        elif ExistingSpecies:
            try:
                ModlessSpectrumFlag = int(ExistingSpecies.Bits[FormatBits.SisterAnnotationFlag])
                BestModlessSpectrumPath = ExistingSpecies.Bits[FormatBits.BestModlessSpectrumPath]
                BestModlessMQScore = float(ExistingSpecies.Bits[FormatBits.BestModlessMQScore])
                ModlessMSSpectrum = ExistingSpecies.GetConsensusModlessMSSpectrum(self)
            except:
                pass # the modless-bits weren't set; that's fine
        else:
            ModlessSpectrumFlag = ""
            BestModlessSpectrumPath = ""
            BestModlessMQScore = ""
        if ModlessSpectrumFlag:
            MSSpectrum = SpeciesA.GetConsensusMSSpectrum(self)
            Comparator = SpectralSimilarity.SpectralSimilarity(MSSpectrum,
                ModlessMSSpectrum, NewAnnotation, NewModlessAnnotation)
            # COPIED from ComputePTMFeatures:
            Comparator.LabelPeaks(0.5)
            Similarity = Comparator.DotProduct(0.5)
            NewFeatures[FeatureBits.Dot] = Similarity
            Similarity = Comparator.GetSharedPeakCount(0, 1)
            NewFeatures[FeatureBits.Shared01] = Similarity
            Similarity = Comparator.GetSharedPeakCount(1, 1)
            NewFeatures[FeatureBits.Shared11] = Similarity
            CorrelationCoefficient = Comparator.ComputeCorrelationCoefficient(1.0)
            NewFeatures[FeatureBits.Correlation] = Similarity
        if SpeciesA.Charge > 2:
            Model = self.Model3
        else:
            Model = self.Model2
        NewModelScore = self.ScoreInstance(Model, NewFeatures)
        # Finalize:
        # If self.PeptideDict already has an entry, then this peptide
        # is the same as another after reconciliation.  But, we already
        # skipped the opportunity to merge with that other.
        if ExistingSpecies:
            self.ReconcileDetailOutput.write("# Existing species has score %s vs model %s\n"%(ExistingSpecies.ModelScore, SpeciesA.ModelScore))
            if ExistingSpecies.ModelScore < SpeciesA.ModelScore:
                ExistingSpecies.MergedFlag = 1
            else:
                # We want to reconcile to the master...but that would make us the same as another
                # peptide species, which we refused to merge with!
                self.ReconcileDetailOutput.write("# *-> We'd like to reconcile %s to %s, but...\n"%(SpeciesA.Annotation, NewAnnotation))
                self.ReconcileDetailOutput.write("# ...there's ALREADY a superior peptide at %s\n"%str(Key))
                return 0
        ################################################################################
        # All tests passed.  RECONCILE!
        self.ReconcileDetailOutput.write("> Reconcile %s to %s\n"%(SpeciesA.Annotation, NewAnnotation))
        # Copy over our consensus spectrum.  
        NewSpectrumPath = os.path.join(self.ConsensusSpectraDirAdjusted, NewAnnotation[2], "%s.%s.dta"%(NewAnnotation, SpeciesA.Charge))
        OldSpectrumPath = os.path.join(self.ConsensusSpectraDirAdjusted, SpeciesA.Annotation[2], "%s.%s.dta"%(SpeciesA.Annotation, SpeciesA.Charge))
        if os.path.exists(OldSpectrumPath):
            # If we've already adjusted once, move the *old* adjusted to the *new* adjusted:
            Command = "move \"%s\" \"%s\""%(OldSpectrumPath, NewSpectrumPath)
        else:
            OldSpectrumPath = os.path.join(self.ConsensusSpectraDir, SpeciesA.Annotation[2], "%s.%s.dta"%(SpeciesA.Annotation, SpeciesA.Charge))
            Command = "copy \"%s\" \"%s\""%(OldSpectrumPath, NewSpectrumPath)
        print Command
        os.system(Command)
        # Copy over the cluster member list:
        MemberListStr = SpeciesA.GetMemberListStr(self)
        Path = os.path.join(self.ClusterScanListDirAdjusted, NewAnnotation[2], "%s.%s.txt"%(NewAnnotation, SpeciesA.Charge))
        ClusterMemberFile = open(Path, "wb")
        ClusterMemberFile.write(MemberListStr)
        ClusterMemberFile.close()
        # Update features and such:
        SpeciesA.Features = NewFeatures
        SpeciesA.ConsensusMQScore = NewConsensusScore
        SpeciesA.ModelScore = NewModelScore
        SpeciesA.Peptide = NewPeptide
        SpeciesA.ComputePrefixes()
        SpeciesA.Annotation = NewAnnotation
        SpeciesA.ModDBPos = SpeciesB.ModDBPos
        SpeciesA.ModMass = SpeciesB.ModMass
        SpeciesA.ModAA = SpeciesB.ModAA
        SpeciesA.ModlessAnnotation = NewModlessAnnotation
        # Revise file bits:
        SpeciesA.Bits[FormatBits.DBPos] = str(SpeciesA.ModDBPos)
        SpeciesA.Bits[FormatBits.ModificationMass] = str(SpeciesA.ModMass)
        SpeciesA.Bits[FormatBits.ModifiedAA] = SpeciesB.Bits[FormatBits.ModifiedAA]
        SpeciesA.Bits[FormatBits.ModifiedResidueNumber] = SpeciesB.Bits[FormatBits.ModifiedResidueNumber]
        SpeciesA.Bits[FormatBits.Peptide] = NewAnnotation
        for FeatureIndex in range(len(SpeciesA.Features)):
            SpeciesA.Bits[FormatBits.FirstFeature + FeatureIndex] = str(SpeciesA.Features[FeatureIndex])
        SpeciesA.Bits[FormatBits.ModelScore] = str(NewModelScore)
        SpeciesA.Bits[FormatBits.ConsensusMQScore] = str(NewConsensusScore)
        # Bits for modless spectra:
        SpeciesA.Bits[FormatBits.SisterAnnotationFlag] = str(ModlessSpectrumFlag)
        SpeciesA.Bits[FormatBits.BestModlessSpectrumPath] = BestModlessSpectrumPath
        SpeciesA.Bits[FormatBits.BestModlessMQScore] = str(BestModlessMQScore)
        # Remove our old PeptideDict entry:
        try:
            del self.PeptideDict[(OldAnnotation, SpeciesA.Charge)]
        except:
            pass
        # Add a new PeptideDict entry:
        self.PeptideDict[Key] = SpeciesA
        ############################################
        # Update our COVERAGE and MODDED-FRACTION:
        for DBPos in range(OldDBPos, OldDBEnd):
            self.ModCoverage[DBPos] -= int(SpeciesA.Features[FeatureBits.SpectrumCount])
        for DBPos in range(SpeciesA.DBPos, SpeciesA.DBEnd):
            self.ModCoverage[DBPos] += int(SpeciesA.Features[FeatureBits.SpectrumCount])
        SpeciesA.Features[FeatureBits.ModdedFraction] = self.ModCoverage[SpeciesA.ModDBPos] / float(self.ModCoverage[SpeciesA.ModDBPos] + self.Coverage[SpeciesA.ModDBPos])
        ############################################    
        return 1
    def LoadModel(self):
        print "load %s"%self.SavedModelFileName2
        self.Model2 = Learning.LoadGeneralModel(self.SavedModelFileName2)
        print "load %s"%self.SavedModelFileName3
        self.Model3 = Learning.LoadGeneralModel(self.SavedModelFileName3)
        #print "Model.MixtureModel", self.Model.MixtureModel
    def GroupPeptidesBySite(self):
        self.Sites = {} # (ModDBPos, ModMass) -> site instance
        for Species in self.SpeciesList:
            if Species.MergedFlag:
                continue
            Key = (Species.ModDBPos, Species.ModMass)
            Site = self.Sites.get(Key, None)
            if not Site:
                Site = SiteClass()
                Site.ModDBPos = Species.ModDBPos
                Site.ModMass = Species.ModMass
                Site.ModAA = Species.ModAA
                self.Sites[Key] = Site
            Site.SpeciesList.append(Species)
    def ScorePTMSites(self):
        """
        Group peptide species by site, and compute the p-value (odds FALSE) for
        each site. 
        """
        self.GroupPeptidesBySite()
        for Site in self.Sites.values():
            Site.PValue = 1.0
            for Species in Site.SpeciesList:
                while len(Species.Bits) <= FormatBits.SitePValue:
                    Species.Bits.append("")
                if Species.Charge > 2:
                    Model = self.Model3
                else:
                    Model = self.Model2
                PeptidePValue = Model.GetPValue(Species.ModelScore)
                Species.Bits[FormatBits.ModelPValue] = str(PeptidePValue)
                Site.PValue *= PeptidePValue
            for Species in Site.SpeciesList:
                Species.Bits[FormatBits.SitePValue] = str(Site.PValue)
    def RescoreAndWriteSite(self, Site):
        Site.PValue = 1.0
        for Species in Site.SpeciesList:
            ModTypeKey = (Species.Bits[FormatBits.ModifiedAA], Species.ModMass)
            TotalSpectra = self.ModTypeSpectrumCount.get(ModTypeKey, 0)
            TotalSites = self.ModTypeSiteCount.get(ModTypeKey, 0)
            Species.Features[FeatureBits.SpectraThisModType] = TotalSpectra
            Species.Features[FeatureBits.SitesThisModType] = TotalSites
            Species.Features[FeatureBits.LogSpecThisType] = math.log(TotalSpectra)
            Species.Features[FeatureBits.LogSitesThisType] = math.log(TotalSites)
            Species.Bits[FormatBits.SpectraWithThisModType] = str(TotalSpectra)
            Species.Bits[FormatBits.SitesWithThisModType] = str(TotalSites)
            Species.Bits[FormatBits.LogSpectraThisModType] = str(math.log(TotalSpectra))
            Species.Bits[FormatBits.LogSitesThisModType] = str(math.log(TotalSites))
            DBPosition = int(Species.Bits[FormatBits.DBPos])
            # Update modded%:
            ModdedSpectra = self.ModCoverage[DBPosition]
            ModlessSpectra = self.Coverage[DBPosition]
            TotalSpectra = ModdedSpectra + ModlessSpectra
            if TotalSpectra <= 0:
                print "*** Warning: Site %s has no coverage at DB position %s"%(Species.Annotation, DBPosition)
            Species.Bits[FormatBits.ModdedFraction] = str(ModdedSpectra / float(max(1, TotalSpectra)))
            # Pad with empty bits if necessary:
            while len(Species.Bits) <= FormatBits.SitePValue:
                Species.Bits.append("")
            if Species.Charge > 2:
                Model = self.Model3
            else:
                Model = self.Model2
            PeptidePValue = Model.GetPValue(Species.ModelScore)
            Species.Bits[FormatBits.ModelPValue] = str(PeptidePValue)
            Site.PValue *= PeptidePValue
        for Species in Site.SpeciesList:
            Species.Bits[FormatBits.SitePValue] = str(Site.PValue)
            Str = string.join(Species.Bits, "\t")
            self.OutputFile.write(Str + "\n")
    def OutputPTMs(self):
        File = open(self.OutputFileName, "wb")
        for Line in self.HeaderLines:
            File.write(Line)
        # Sort the sites:
        SortedSites = []
        for Site in self.Sites.values():
            SortedSites.append((Site.PValue, Site))
        SortedSites.sort()
        CumulativeSiteCount = 0
        CumulativeTrueSiteCount = 0
        CumulativeSpeciesCount = 0
        # Report peptides, grouped by site, from best site to worst:
        for (PValue, Site) in SortedSites:
            BestSpecies = None
            for Species in Site.SpeciesList:
                if Species.MergedFlag:
                    continue
                CumulativeSpeciesCount += 1
                if (BestSpecies == None) or (Species.PValue < BestSpecies.PValue):
                    BestSpecies = Species
            if not BestSpecies:
                continue
            CumulativeSiteCount += 1
            if int(Species.Bits[FormatBits.TrueProteinFlag]):
                CumulativeTrueSiteCount += 1
            # The LENS way:
            FalseProteinCount = CumulativeSiteCount - CumulativeTrueSiteCount
            #FalseWithinTrue = FalseProteinCount * 0.01
            FalseWithinTrue = FalseProteinCount
            TrueCount = max(0, CumulativeTrueSiteCount - FalseWithinTrue)
            SiteCount = CumulativeTrueSiteCount
            # FDR:
            if CumulativeTrueSiteCount <= 0:
                FDR = 1.0
            else:
                # False discovery rate:
                # The number of spurious sites which come from valid proteins
                # divided by the number of sites that come from valid proteins
                FDR = FalseWithinTrue / float(CumulativeTrueSiteCount)
                FDR = min(1.0, FDR)
            print "pvalue %.6f sites%d species%d T%d F%d FWT %.3f SC %.3f FDR %.3f"%(\
                Site.PValue, CumulativeSiteCount, CumulativeSpeciesCount, 
                CumulativeTrueSiteCount, FalseProteinCount, FalseWithinTrue, SiteCount, FDR)
            for Species in Site.SpeciesList:
                if Species.MergedFlag:
                    continue
                try:
                    Str = string.join(Species.Bits, "\t")
                    Str += "\t%s\t"%FDR
                except:
                    traceback.print_exc()
                    print Species.Bits
                    print map(type, Species.Bits)
                File.write(Str + "\n")
        File.close()
    def ParseDB(self):
        DBFile = open(self.DBFileName, "rb")
        self.DB = DBFile.read()
        DBFile.close()
    def CheckSpectrumDirectories(self):
        """
        Create our adjusted-spectrum and adjusted-cluster directories, wiping out
        old ones beforehand if we must.  We do this for merge-and-reconcile;
        we DON'T do it for biochem tweaking
        """
        print "Prepare spectrum directories... (-z option, set for "
        print "single-block runs and the first block of multi-block runs...)"
        try:
            shutil.rmtree(self.ConsensusClusterDirAdjusted)
        except:
            pass
        try:
            shutil.rmtree(self.ConsensusSpectraDirAdjusted)
        except:
            pass
        try:
            shutil.rmtree(self.ClusterScanListDirAdjusted)
        except:
            pass
        MakeDirectory(self.ConsensusClusterDirAdjusted)
        MakeDirectory(self.ConsensusSpectraDirAdjusted)
        MakeDirectory(self.ClusterScanListDirAdjusted)
        Aminos = "ACDEFGHIKLMNOPQRSTUVWY"
        for Amino in Aminos:
            Dir = os.path.join(self.ConsensusClusterDirAdjusted, Amino)
            MakeDirectory(Dir)
            Dir = os.path.join(self.ConsensusSpectraDirAdjusted, Amino)
            MakeDirectory(Dir)
            Dir = os.path.join(self.ClusterScanListDirAdjusted, Amino)
            MakeDirectory(Dir)
            
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:d:c:m:k:M:x:y:X:zev:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFileName = Value
            elif Option == "-m":
                self.SavedModelFileName2 = "%s.2"%Value
                self.SavedModelFileName3 = "%s.3"%Value
            elif Option == "-M":
                self.OutputModelFileName2 = "%s.2"%Value
                self.OutputModelFileName3 = "%s.3"%Value
            elif Option == "-w":
                self.OutputFileName = Value
            elif Option == "-d":
                self.DBFileName = Value
            elif Option == "-k":
                self.KnownChemistryFileName = Value
            elif Option == "-c":
                self.TempFileDir = Value
                self.ConsensusClusterDir = os.path.join(Value, "Clusters")
                self.ConsensusSpectraDir = os.path.join(Value, "Spectra")
                self.ClusterScanListDir = os.path.join(Value, "ClusterMembers")
                self.ConsensusClusterDirAdjusted = os.path.join(Value, "ClustersAdjusted")
                self.ConsensusSpectraDirAdjusted = os.path.join(Value, "SpectraAdjusted")
                self.ClusterScanListDirAdjusted = os.path.join(Value, "ClusterMembersAdjusted")
            elif Option == "-x":
                self.DBStart = int(Value)
            elif Option == "-y":
                self.DBEnd = int(Value)
            elif Option == "-X":
                self.SpectrumRoot = Value
            elif Option == "-z":
                self.CheckDirectoriesFlag = 1
            elif Option == "-e":
                # mErge block runs:
                self.MergeBlockRunsFlag = 1
            elif Option == "-v":
                self.KnownPTMVerboseOutputFileName = Value
    def SaveSitesByType(self):
        "Save sites/spectra by modification type"
        PicklePath = os.path.join(self.TempFileDir, "SitesByModType.dat")
        SitesFile = open(PicklePath, "wb")
        cPickle.dump(self.ModTypeSiteCount, SitesFile)
        cPickle.dump(self.ModTypeSpectrumCount, SitesFile)
        SitesFile.close()
    def LoadSitesByType(self):
        "Load sites/spectra by modification type"
        PicklePath = os.path.join(self.TempFileDir, "SitesByModType.dat")
        SitesFile = open(PicklePath, "rb")
        self.ModTypeSiteCount = cPickle.load(SitesFile)
        self.ModTypeSpectrumCount = cPickle.load(SitesFile)
        SitesFile.close()
    def LoadHeaderLines(self, FileName):
        File = open(FileName, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber > 10:
                break
            if FileLine[0] == "#":
                self.HeaderLines.append(FileLine)
        File.close()
    def Main(self):
        if self.SpectrumRoot:
            self.PopulateSpectrumOracle(self.SpectrumRoot)
        print "Load model..."
        self.LoadModel()
        print "Parse database..."
        self.ParseDB()
        if self.MergeBlockRunsFlag:
            self.MergeBlockRuns()
        elif self.KnownChemistryFileName:
            print "Tweak sites to match KNOWN CHEMISTRY..."
            self.LoadKnownModifications()
            self.LoadSitesByType()
            self.LoadHeaderLines(self.InputFileName)
            self.OutputFile = open(self.OutputFileName, "wb")
            for HeaderLine in self.HeaderLines:
                self.OutputFile.write(HeaderLine)
            self.ProcessSites(self.InputFileName, "knownptm")
        else:
            print "MERGE and RECONCILE..."
            if self.CheckDirectoriesFlag:
                self.CheckSpectrumDirectories()
            self.MergeAndReconcile()
            # Re-compute spectra for mod type, since some peptides
            # have now been re-annotated:
            #self.ComputeTotalSpectraForModType()
            #self.RescorePeptides()
            #self.ScorePTMSites()
        
        #print "Write output to %s"%(self.OutputFileName)
        #self.OutputPTMs()        

UsageInfo = """
AdjustPTM: Merge, reconcile, and tweak PTM annotations.

Arguments:
-r [FILENAME]: Feature file (from TrainPTMFeatures) to read in
-w [FILENAME]: Output modded-peptide filename
-d [DBFILE]: Database searched
-c [DIR]: Cluster directory
-k [FILENAME]: Known chemistry filename.  If specified, consider altering sites
   to match known chemical adducts; report the best site-score attainable by using
   known chemical adducts.
-m [FILENAME]: Peptide scoring model INPUT filename
-M [FILENAME]: Peptide scoring model OUTPUT filename
-x [POS]: Database start position
-y [POS]: Database end position
"""
   
if __name__ == "__main__":
    if PROFILING_RUN:
        import profile
        profile.run("Main()")
    else:
        try:
            import psyco
            psyco.full()
        except:
            print "(Psyco not found - no optimization)"
        Adjutant = PTMAdjuster()
        Adjutant.ParseCommandLine(sys.argv[1:])
        Adjutant.Main()
