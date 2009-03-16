"""
Parse all the peptides, including modified peptides, from an unrestrictive search.
Store them in a text file.  We want this list so that we can report the residue-level
coverage of the database.
"""
import os
import getopt
import sys
import struct
import traceback
import ResultsParser
import StripPTM
from Utils import *
Initialize()

class PeptideGrabber(ResultsParser.ResultsParser):
    def __init__(self):
        # Hard-coded for now:
        self.ResultsDir = r"f:\ftproot\briggs\hek293\resultsblind1fixed"
        self.OutputFileDir = "."
        self.PValueCutoff = 0.05
        self.PeptideDict = {}
        self.GroupNames = ["All"]
        self.DBPath = None
        # Cache group names:
        self.GroupNamesForSpectrum = {}
        ResultsParser.ResultsParser.__init__(self)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:p:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.ResultsDir = Value
            elif Option == "-w":
                self.OutputFileDir = Value
            elif Option == "-p":
                self.PValueCutoff = float(Value)
            elif Option == "-d":
                self.DBPath = Value
    def ParseFile(self, FileName):
        """
        Callback to parse a single results-file.  Read all the top-scoring
        annotations, and count the number of occurrences of each annotation.
        """
        try:
            File = open(FileName, "rb")
            OldSpectrum = None
            LineNumber = 0
            for FileLine in File.xreadlines():
                LineNumber += 1
                if LineNumber % 1000 == 0:
                    print "Line %s..."%LineNumber
                Bits = FileLine.strip().split("\t")
                if FileLine[0] == "#":
                    continue
                if len(Bits) < self.Columns.PValue:
                    continue
                try:
                    Annotation = Bits[self.Columns.Annotation]
                    #Peptide = GetPeptideFromModdedName(Bits[self.Columns.Annotation])
                    PValue = float(Bits[self.Columns.PValue])
                    Spectrum = (Bits[self.Columns.SpectrumFile], Bits[self.Columns.ScanNumber])
                    SpectrumFileName = Bits[self.Columns.SpectrumFile]
                except:
                    traceback.print_exc()
                    continue
                if PValue > self.PValueCutoff:
                    continue
                if Spectrum == OldSpectrum:
                    continue
                OldSpectrum = Spectrum
                PeptideInfo = self.PeptideDict.get(Annotation, None)
                if not PeptideInfo:
                    PeptideInfo = Bag()
                    PeptideInfo.Annotation = Annotation
                    PeptideInfo.SpectrumCount = {}
                    self.PeptideDict[Annotation] = PeptideInfo
                # Update the spectrum count(s) for this peptide:
                Groups = self.GroupNamesForSpectrum.get(SpectrumFileName, None)
                if not Groups:
                    Groups = self.GetGroups(SpectrumFileName)
                    self.GroupNamesForSpectrum[SpectrumFileName] = Groups
                for Group in Groups:
                    PeptideInfo.SpectrumCount[Group] = PeptideInfo.SpectrumCount.get(Group, 0) + 1
                #self.PeptideDict[Annotation] = self.PeptideDict.get(Annotation, 0) + 1
            File.close()
        except:
            traceback.print_exc()
            print "* Unable to parse results from:", FileName
        print "So far: %s annotations"%len(self.PeptideDict.keys())
    def GetGroups(self, SpectrumFileName):
        """
        Return a list of the sub-samples this spectrum is a part of.
        Default behavior: Every spectrum is part of "all"; otherwise, you're
        part of a group if your spectrum file-name includes the group name
        """
        GroupNames = ["All"]
        for GroupName in self.GroupNames:
            if GroupName == "All":
                continue
            Pos = SpectrumFileName.find(GroupName)
            if Pos != -1:
                GroupNames.append(GroupName)
        return GroupNames
    def FixupPeptides(self):
        """
        Call StripNeedlessModifications on each of our peptides.  Also, set the
        Peptide member of all PeptideInfo objects.
        """
        Keys = self.PeptideDict.keys()
        Keys.sort()
        for Annotation in Keys:
            OldPeptideInfo = self.PeptideDict[Annotation]
            (DBPos, FixedAnnotation) = StripPTM.StripNeedlessModifications(self.DB, Annotation)
            OldPeptideInfo.DBPos = DBPos
            OldPeptideInfo.Peptide = GetPeptideFromModdedName(FixedAnnotation)
            if FixedAnnotation == Annotation:
                continue
            if self.PeptideDict.has_key(FixedAnnotation):
                # We're MERGING this peptide with a peptide that doesn't have
                # the superfluous PTM.  Update the spectrum count(s) for NewPeptideInfo, 
                # the assimilator.
                NewPeptideInfo = self.PeptideDict[FixedAnnotation]
                for (Key, Count) in OldPeptideInfo.SpectrumCount.items():
                    NewPeptideInfo.SpectrumCount[Key] = NewPeptideInfo.SpectrumCount.get(Key, 0) + Count
            else:
                self.PeptideDict[FixedAnnotation] = OldPeptideInfo
            del self.PeptideDict[Annotation]
    def OutputPeptides(self):
        """
        Output a table with one line for each peptide observed on our sample.
        """
        HeaderString = "#Annotation\tDBPos\tSpectra\tModifiedFlag\t\n"
        for Group in self.GroupNames:
            Path = os.path.join(self.OutputFileDir, "Peptides.%s.txt"%Group)
            OutputFile = open(Path, "wb")
            OutputFile.write(HeaderString)
            for (Annotation, PeptideInfo) in self.PeptideDict.items():
                Peptide = PeptideInfo.Peptide
                SpectrumCount = PeptideInfo.SpectrumCount.get(Group, 0)
                if SpectrumCount:
                    Str = "%s\t%s\t%s\t"%(Annotation, PeptideInfo.DBPos, SpectrumCount)
                    if len(Peptide.Modifications.keys()):
                        Str += "1\t"
                    else:
                        Str += "0\t"
                    OutputFile.write(Str + "\n")
            OutputFile.close()
    def Main(self):
        # Rudimentary checking:
        if os.path.exists(self.OutputFileDir):
            if not os.path.exists(self.OutputFileDir):
                print "** Error: Output file should be a directory!"
                print UsageInfo
                sys.exit(-1)
        else:
            MakeDirectory(self.OutputFileDir)
        # Parse the database:
        DBFile = open(self.DBPath, "rb")
        self.DB = DBFile.read()
        DBFile.close()
        # Read the peptides from p-valued results:
        self.ProcessResultsFiles(self.ResultsDir, self.ParseFile)
        # Fix unnecessary modifications:
        self.FixupPeptides()
        # Write the peptides:
        print "Output peptides to %s..."%self.OutputFileDir
        self.OutputPeptides()
        # Compute and report coverage:
        print "Compute and output coverage levels to %s..."%self.OutputFileDir
        self.ReportResidueCoverage()
    def LoadPeptides(self, PeptideFileName):
        """
        Load a 'pickled' collection of peptides from a text file.
        """
        File = open(PeptideFileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                Count = int(Bits[1])
                Annotation = Bits[0]
            except:
                continue
            PeptideInfo = Bag()
            PeptideInfo.SpectrumCount = Count
            PeptideInfo.Annotation = Annotation
            PeptideInfo.Peptide = GetPeptideFromModdedName(Annotation)
            self.PeptideDict[Annotation] = PeptideInfo
        File.close()
    def MeasureCoverage(self, DatabaseString):
        self.Coverage = {}
        self.ModdedCoverage = {}
        self.DB = DatabaseString
        for Annotation in self.PeptideDict.keys():
            Peptide = GetPeptideFromModdedName(Annotation)
            if len(Peptide.Modifications.keys()):
                ModdedFlag = 1
            else:
                ModdedFlag = 0
            LastPos = -1
            # Iterate over database locations for the peptide, and
            # count up the coverage:
            while 1:
                Pos = DatabaseString.find(Peptide.Aminos, LastPos + 1)
                if Pos == -1:
                    break
                for DBPos in range(Pos, Pos + len(Peptide.Aminos)):
                    if ModdedFlag:
                        self.ModdedCoverage[DBPos] = self.Coverage.get(DBPos, 0) + 1
                    else:
                        self.Coverage[DBPos] = self.Coverage.get(DBPos, 0) + 1
                LastPos = Pos
    def ReportCoverage(self):
        DBLength = len(self.DB)
        DBLength -= self.DB.count("*")
        CoverFlags = len(self.Coverage.keys())
        ModdedCoverFlags = len(self.ModdedCoverage.keys())
        print "Residue coverage: %s of %s (%.2f%%)"%(CoverFlags, DBLength, 100 * CoverFlags / float(DBLength))
        print "Modded coverage: %s of %s (%.2f%%)"%(ModdedCoverFlags, DBLength, 100 * ModdedCoverFlags / float(DBLength))
    def ReportResidueCoverage(self):
        for Group in self.GroupNames:
            ####################################################
            # Compute coverage:
            print "Coverage for group '%s'..."%Group
            Coverage = {}
            ModdedCoverage = {}
            Keys = self.PeptideDict.keys()
            for KeyIndex in range(len(Keys)):
                Annotation = Keys[KeyIndex]
                PeptideInfo = self.PeptideDict[Annotation]
                if KeyIndex % 1000 == 0:
                    print "Peptide %s/%s..."%(KeyIndex, len(Keys))
                SpectrumCount = PeptideInfo.SpectrumCount.get(Group, 0)
                if not SpectrumCount:
                    continue
                if len(PeptideInfo.Peptide.Modifications.keys()):
                    ModifiedFlag = 1
                else:
                    ModifiedFlag = 0
                for Pos in range(PeptideInfo.DBPos, PeptideInfo.DBPos + len(PeptideInfo.Peptide.Aminos)):
                    if ModifiedFlag:
                        ModdedCoverage[Pos] = ModdedCoverage.get(Pos, 0) + SpectrumCount
                    else:
                        Coverage[Pos] = Coverage.get(Pos, 0) + SpectrumCount
            ####################################################
            # Output coverage (binary format):
            CoveragePath = os.path.join(self.OutputFileDir, "Coverage.%s.dat"%Group)
            CoverageFile = open(CoveragePath, "wb")
            for DBPos in range(len(self.DB)):
                Str = struct.pack("<II", Coverage.get(DBPos, 0), ModdedCoverage.get(DBPos, 0))
                CoverageFile.write(Str)
            CoverageFile.close()
##            ####################################################
##            # Output coverage (text format):
##            CoveragePath = os.path.join(self.OutputFileDir, "Coverage.%s.txt"%Group)
##            CoverageFile = open(CoveragePath, "wb")
##            for DBPos in range(len(self.DB)):
##                Spectra = Coverage.get(DBPos, 0)
##                ModdedSpectra = ModdedCoverage.get(DBPos, 0)
##                if Spectra + ModdedSpectra == 0:
##                    continue
##                Str = "%s\t%s\t%s\t"%(DBPos, Spectra, ModdedSpectra)
##                CoverageFile.write(Str + "\n")
##            CoverageFile.close()
            pass
                
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(no psyco - no optimization performed)"
    Grabber = PeptideGrabber()
    Grabber.ParseCommandLine(sys.argv[1:])
    Grabber.Main()