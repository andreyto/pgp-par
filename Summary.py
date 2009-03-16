"""
The new Summary script, modified for lower memory usage.  This version
assumes that the top annotation is the correct one, and does not consider
re-assigning a spectrum to a close runner-up or to a homologous protein.

- Iterate over annotations.  Filter any that aren't a top-scoring annotation with
good p-value.
- Remember the best modified and the best unmodified annotation for each database interval
- Generate the protein report: one sub-report for each protein, sorted by coverage
- Sub-report includes a table of covered residues, followed by a list of peptide rows.
- Each peptide row gives a spectrum count and score for the best spectrum
- Peptide rows are sorted by protein position, then mass
"""
import os
import time
import sys
import struct
import traceback
import shutil
import getopt
import Label
import MakeImage
import MSSpectrum
import SelectProteins
import ResultsParser
from Utils import *
Initialize()

UsageInfo = """Summary.py - Parse search results, and generate either a webpage
  summarizing the results, or a filtered database for unrestrictive PTM search.

Required options:
 -r [FileName] - The name of the results file to parse.  If a directory is
    specified, then all .txt files within the directory will be combined into
    one report
 -d [FileName] - The name of the database file (.trie format) searched.
         (allows more than one database file; use multiple -d options)
 
Additional options:
 -b [FileName] - Second-pass database filename.  If specified, the proteins
    selected will be written out to a database (.fasta, .trie and .index files)
    suitable for unrestrictive search.
 -w [FileName] - Webpage directory.  If specified, a webpage will be written
    to the specified directory, summarizing the proteins identified, and the
    degree of coverage.

 (Note: Either -b or -w, or both, must be provided)
    
 -p [Value] - Cutoff p-value.  Annotations with inferior p-value are ignored.
    Defaults to 0.05.
 -e [Count] - Minimum number of peptides that a protein must annotate in order
    to add it to the report or the filtered database.  Defaults to 1.
 -m [Count] - Minimum number of spectra that a protein must annotate in order
    to add it to the report or the filtered database.  By default, this count
    is set to (SpectrumCount / ProteinsInDatabase) * 2.  If the protein
    database has already been filtered, set this parameter to 1.
 -v [Count] - Verbose spectrum output count.  If set, report [Count] spectra
    for each distinct peptide identified.  This option is slower and
    consumes more memory, but can be more informative.
 -i [SpectraPath] - For use if verbose spectrum output (-v) is enabled.
    Images will be generated for each annotation, if the Python Imaging
    Library (PIL) is installed.  This option generates many files on disk,
    so it's recommended that you set the summary file (-w option) in its own
    directory.  "SpectraPath" is the path to the folder with MS2 spectra
    
    
Examples:
    Summary.py -r Frac1Ouput.txt -d Database%ssprot.trie -p 0.1
    Summary.py -r Frac1Ouput.txt -d Database%ssprot.trie -w F1Summary\index.html -v
"""%(os.sep, os.sep)

class SummarizerClass(ResultsParser.ResultsParser):
    def __init__(self):
        # SharedProteinSpectra[ProteinA][ProteinB] = number of spectra for
        # which the annotation is found in both protein A and protein B.  When
        # we accept one of the two proteins, the other one loses some annotations!
        self.SharedProteinSpectra = {}
        # SharedProteinPeptides is similar to SharedProteinSpectra, but tracks
        # distinct peptide records
        self.SharedProteinPeptides = {} 
        self.AnnotationCounts = {} # ProteinID -> count
        # BestRepresentatives[(DBStart, DBEnd)] is a list of peptide
        # instances for the best spectra for that positon.
        self.BestRepresentatives = {}
        self.BestModRepresentatives = {}
        self.ResultsFileName = None
        self.DatabasePath = []
        self.SecondPassDatabasePath = None
        self.MinimumProteinHits = None
        self.SummaryPagePath = None
        self.MZXMLPath = None
        self.PValueCutoff = 0.05
        # Keys are peptides (after I->L substitution), values are lists of protein
        # record numbers.  We keep this dictionary so that we needn't repeatedly map
        # the same protein to the database.
        self.PeptideDict = {}
        self.VerboseProteinReportCount = 0
        self.GenerateSpectrumImagesFlag = 0
        # Very short peptides are uninformative...skip them.
        self.MinimumPeptideLength = 7
        self.SpectrumCount = 0
        self.IntervalHitCount = {}
        self.MinimumPeptides = 1
        ResultsParser.ResultsParser.__init__(self)
    def GetSpectrumPath(self, Path, SpectraPath):
        """
        This requires a bit of trickery, because sometimes the results are
        generated on a unix machine (creating a unix path), and this script
        is run on a windows machine (can't split a unix path). So I have going
        to hack things out on my own.
        """
        FileName = None
        if Path.find("/") >= 0:
            "results files made on a unix machine"
            LastSlash = Path.rfind("/")
            FileName = Path[LastSlash+1:]
        else:
            "results files made on a windows machine, hopefully.  Any other users, go home"
            LastBackSlash = Path.rfind("\\")
            FileName = Path[LastBackSlash+1:]
        if not FileName:
            print "unable to create a path to the spectrum file %s"%Path
            return Path
        return os.path.join(SpectraPath,FileName)
    def WriteSecondPassDatabase(self):
        """
        Write out the "present" proteins to our second-pass database.
        self.ProteinSelector is responsible for deciding which peptides
        belong to which proteins.
        """
        Bits = os.path.split(self.SecondPassDatabasePath)
        DBPathStub = os.path.splitext(self.SecondPassDatabasePath)[0]
        if len(Bits[0]) == 0:
            DBPath = os.path.join("Database", "%s.trie"%DBPathStub)
            IndexPath = os.path.join("Database", "%s.index"%DBPathStub)
            FastaPath = os.path.join("Database", "%s.fasta"%DBPathStub)
        else:
            DBPath = DBPathStub + ".trie"
            IndexPath = DBPathStub + ".index"
            FastaPath = DBPathStub + ".fasta"
        print "Writing second-pass database to %s..."%DBPath
        DBFile = open(DBPath, "wb")
        IndexFile = open(IndexPath, "wb")
        FASTAFile = open(FastaPath, "wb")
        DBFilePos = 0
        for (ProteinID, ScoreTuple) in self.ProteinSelector.SelectedProteins.items():
            (PeptideCount, SpectrumCount) = ScoreTuple
            if SpectrumCount < self.MinimumProteinHits:
                continue
            if PeptideCount < self.MinimumPeptides:
                continue
            # Let's write out the protein.  Write to the INDEX file, the
            # TRIE file, and a FASTA file.  (The FASTA file is just for
            # humans to read)
            ProteinName = self.ProteinSelector.ProteinNames[ProteinID]
            ProteinSequence = self.ProteinSelector.ProteinSequences[ProteinID]
            Str = struct.pack("<qi80s", 0, DBFilePos, ProteinName[:80])
            IndexFile.write(Str)
            DBFile.write("%s*"%ProteinSequence)
            DBFilePos += len(ProteinSequence) + 1
            Pos = 0
            FASTAFile.write(">%s\n"%ProteinName)
            while Pos < len(ProteinSequence):
                Chunk = ProteinSequence[Pos:Pos+70]
                FASTAFile.write(Chunk)
                FASTAFile.write("\n")
                Pos += 70
        IndexFile.close()
        FASTAFile.close()
        DBFile.close()
    def GetProteinHREF(self, ProteinID, ProteinName):
        # By default, just print the name with no hyperlinking.
        # Subclass can override to hyperlink to IPI, swiss-prot, etc
        return ProteinName
    def WriteSummaryPage(self):
        """
        Produce the protein report.  The index file contains the protein coverage
        information.  If verbose output is requested, also contains one row per peptide.
        """
        Dir = os.path.split(self.SummaryPagePath)[0]
        try:
            os.makedirs(Dir)
        except:
            pass
        # Populate PeptidesForProtein.  Keys are protein IDs.  Values are lists
        # of peptide annotations; the annotations are, in turn, keys for
        # self.ProteinSelector.BestRepresentatives
        self.PeptidesForProtein = {}
        for (Annotation, RepList) in self.ProteinSelector.BestRepresentatives.items():
            if len(RepList) < 1:
                continue
            Peptide = RepList[0][1]
            ProteinID = self.ProteinSelector.PeptideProteins.get(Peptide.Aminos, None)
            if ProteinID == None:
                continue
            PeptideList = self.PeptidesForProtein.get(ProteinID, [])
            PeptideList.append(Annotation)
            self.PeptidesForProtein[ProteinID] = PeptideList
        # Sort proteins from "best" to "worst".  For now, just sort by the
        # number of distinct peptides.
        SortedProteins = []
        for (ProteinID, AnnotationList) in self.PeptidesForProtein.items():
            AnnotationCount = len(AnnotationList)
            if AnnotationCount >= self.MinimumPeptides:
                SortedProteins.append((AnnotationCount, ProteinID))
        SortedProteins.sort()
        SortedProteins.reverse()
        # Start the index file:
        self.SummaryPageDir = os.path.split(self.SummaryPagePath)[0]
        # Ensure the directory exists:
        try:
            os.makedirs(self.SummaryPageDir)
        except:
            pass
        self.IndexFile = open(self.SummaryPagePath, "w")
        self.IndexFile.write("<html><title>Protein Report</title>\n")
        # Iterate over proteins, writing one record for each one:
        for (PeptideCount, ProteinID) in SortedProteins:
            ProteinName = self.ProteinSelector.ProteinNames[ProteinID]
            if ProteinName[:3]== "XXX":
                #print "Found a fake protein %s"%ProteinName
                continue
            #print "Write protein %s (%s), %s peptides"%(ProteinID, self.ProteinSelector.ProteinNames[ProteinID], PeptideCount)
            self.WriteProteinSummary(ProteinID)
        self.IndexFile.close()
    def WriteProteinSummary(self, ProteinID):
        """
        Write summary page section for a single protein.
        """
        ##########################################################
        # Determine coverage, and sort peptides by position within the protein:
        ProteinName = self.ProteinSelector.ProteinNames[ProteinID]
        ProteinSequence = self.ProteinSelector.ProteinSequences[ProteinID]
        Coverage = [0] * len(ProteinSequence)
        PeptideCount = len(self.PeptidesForProtein.get(ProteinID, []))
        for Annotation in self.PeptidesForProtein.get(ProteinID, []):
            Peptide = self.ProteinSelector.BestRepresentatives[Annotation][0][1]
            MatchPos = ProteinSequence.find(Peptide.Aminos)
            if MatchPos == -1:
                print "** Error: Peptide '%s' assigned to incompatible protein %s '%s'"%(Aminos, ProteinID, ProteinName)
            for Pos in range(MatchPos, MatchPos + len(Peptide.Aminos)):
                Coverage[Pos] += 1
        CoverFlags = 0
        for CoverageCount in Coverage:
            if CoverageCount:
                CoverFlags += 1
        CoverageRate = CoverFlags / float(len(ProteinSequence))
        # Write header:
        SpectrumCount = self.ProteinSelector.ProteinSpectrumCounts[ProteinID]
        HREF = self.GetProteinHREF(ProteinID, ProteinName)
        self.IndexFile.write("<h3>%s</h3>\n<b>%s peptides, %s spectra, %.1f%% coverage</b><br>\n"%(HREF, PeptideCount, SpectrumCount, CoverageRate*100))
        # Write protein sequence:
        ColorUncovered = "#aaaaaa"
        ColorCovered = "#000000"
        OldColor = ColorUncovered
        OldBoldFlag = 0
        BoldFlag = 0
        self.IndexFile.write("<tt>")
        for Pos in range(len(ProteinSequence)):
            ResidueNumber = Pos + 1
            if ResidueNumber%50 == 1:
                if BoldFlag:
                    self.IndexFile.write("</b>")
                    BoldFlag = 0                  
                self.IndexFile.write("</font><br>\n<font color=#000000>")
                OldColor = ColorCovered
                OldBoldFlag = 0
                for Foo in range(4 - len(str(ResidueNumber))):
                    self.IndexFile.write("&nbsp;")
                self.IndexFile.write("%d "%ResidueNumber)
            if ResidueNumber % 10 == 1:
                self.IndexFile.write(" ")
            if Coverage[Pos]:
                Color = ColorCovered
            else:
                Color = ColorUncovered
                BoldFlag = 0
            if Color != OldColor:
                self.IndexFile.write("</font><font color=%s>"%Color)
                OldColor = Color
            if BoldFlag != OldBoldFlag:
                if BoldFlag:
                    self.IndexFile.write("<b>")
                else:
                    self.IndexFile.write("</b>")
                OldBoldFlag = BoldFlag
            self.IndexFile.write("%s"%ProteinSequence[Pos])
        self.IndexFile.write("<br><br></font><font color=#000000></tt>\n\n")
        ###############################################
        # Write individual peptides, if requested:
        if self.VerboseProteinReportCount:
            # Write out peptides:
            self.WritePeptideHeader(ProteinID, self.IndexFile)
            SortedAnnotations = []
            for Annotation in self.PeptidesForProtein.get(ProteinID, []):
                RepresentativeList = self.ProteinSelector.BestRepresentatives[Annotation]
                Peptide = RepresentativeList[0][1]
                Pos = ProteinSequence.find(Peptide.Aminos)
                SortedAnnotations.append((Pos, Pos + len(Peptide.Aminos) - 1, Annotation))
            SortedAnnotations.sort()
            for (StartPos, EndPos, Annotation) in SortedAnnotations:
                IntervalString = "%s-%s"%(StartPos + 1, EndPos + 1)
                TotalHitCount = self.ProteinSelector.AnnotationSpectrumCounts[Annotation]
                RepresentativeList = self.ProteinSelector.BestRepresentatives[Annotation]
                RepresentativeList.reverse() # they're sorted from worst-to-best; fix that.
                for Index in range(len(RepresentativeList)):
                    Peptide = RepresentativeList[Index][-1]
                    self.WritePeptideLine(self.IndexFile, IntervalString, ProteinID, Index, Peptide, TotalHitCount)
            self.WritePeptideFooter(self.IndexFile)
        self.IndexFile.write("<hr>")
    def WritePeptideFooter(self, IndexFile):
        IndexFile.write("</table>\n")
    def WritePeptideHeader(self, ProteinID, IndexFile):
        IndexFile.write("<table><tr><td><b>Residues</b></td><td><b>Total Spectra</b></td><td><b>Peptide</b></td><td><b>p-value</b></td><td><b>MQScore</b></td><td><b>File</b></td><td><b>Scan</b></td></tr>")
    def WritePeptideLine(self, File, IntervalStr, ProteinID, SpectrumIndex, Peptide, TotalHitCount):
        Dir = os.path.split(self.SummaryPagePath)[0]
        Annotation = Peptide.GetModdedName()
        SpecFileName = Peptide.SpectrumFilePath.replace("/","\\").split("\\")[-1]
        if self.GenerateSpectrumImagesFlag:
            ImageFileName = "%s.%s.png"%(Annotation, SpectrumIndex)
            ImageFilePath = os.path.join(Dir, ImageFileName)
            Maker = MakeImage.MSImageMaker()
            MSpectrum = MSSpectrum.SpectrumClass()
            Path = self.GetSpectrumPath(Peptide.SpectrumFilePath, self.MZXMLPath)
            FileName = "%s:%s"%(Path, Peptide.SpectrumFilePos)
            try:
                #SpectrumFile = Label.OpenAndSeekFile(FileName)
                #print FileName
                #MSpectrum.ReadPeaksFromFile(SpectrumFile, FileName)
                #MSpectrum.RankPeaksByIntensity()
                #SpectrumFile.close()
                #Label.LabelSpectrum(MSpectrum, Peptide)
                #Maker.ConvertSpectrumToImage(MSpectrum, ImageFilePath, Peptide, Width = 500, Height = 380)
                Args = " -r %s -b %d -a %s -w %s -p"%(Path, int(Peptide.SpectrumFilePos), Annotation, ImageFilePath)
                ArgsList = Args.split()
                #print "Parsing Results for %s, scan %s"%(FileName, Scan)
                Dymo = Label.LabelClass()
                Dymo.ParseCommandLine(ArgsList)
                #Dymo.LoadModel(0, Dymo.PeptideHasPhosphorylation)
                Dymo.Main()
                
            except:
                traceback.print_exc()
            File.write("<tr><td>%s</td><td>%s</td><td><a href=\"%s\">%s</td><td>%s</td><td>%s</td>"%(IntervalStr, TotalHitCount, ImageFileName, Peptide.GetFullModdedName(), Peptide.PValue, Peptide.MQScore))
            File.write("<td>%s</td><td>%s</td></tr>\n"%(SpecFileName, Peptide.ScanNumber))
        else:
            File.write("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"%(IntervalStr, TotalHitCount, Peptide.GetFullModdedName(), Peptide.PValue, Peptide.MQScore))
            File.write("<td>%s</td><td>%s</td></tr>\n"%(SpecFileName, Peptide.ScanNumber))
    def SetMinimumProteinHits(self):
        ProteinCount = len(self.ProteinNames)
        self.MinimumProteinHits = (2 * self.SpectrumCount) / ProteinCount
        self.MinimumProteinHits = max(self.MinimumProteinHits, 2)
        print "%s spectra with a valid annotation, %s proteins in all"%(self.SpectrumCount, ProteinCount)
        print "Minimum hits required to accept an additional protein: ", self.MinimumProteinHits
    def Main(self):
        self.ProteinSelector = SelectProteins.ProteinSelector()
        self.ProteinSelector.PValueCutoff = self.PValueCutoff
        self.ProteinSelector.LoadMultipleDB(self.DatabasePath)
        # If we're expected to write out a summary page, then keep track
        # of the top N representatives for each annotation:
        if self.SummaryPagePath:
            if self.VerboseProteinReportCount:
                self.ProteinSelector.RetainRepresentativeCount = self.VerboseProteinReportCount
            else:
                self.ProteinSelector.RetainRepresentativeCount = 1
        self.ProcessResultsFiles(self.ResultsFileName, self.ProteinSelector.ParseAnnotations)
        self.ProteinSelector.ChooseProteins()
        if self.SecondPassDatabasePath:
            self.WriteSecondPassDatabase()
            print "Second-pass database written to:", self.SecondPassDatabasePath
        if self.SummaryPagePath:
            self.WriteSummaryPage()
            print "Summary page written to:", self.SummaryPagePath
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:b:p:w:m:v:i:e:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ResultsFileName = Value
            elif Option == "-d":
                # -d database
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DatabasePath.append(Value)
            elif Option == "-b":
                # -b Second-pass database
                self.SecondPassDatabasePath = Value
            elif Option == "-i":
                self.GenerateSpectrumImagesFlag = 1
                self.MZXMLPath = Value
            elif Option == "-m":
                # -m Minimum number of spectra for a new protein
                self.MinimumProteinHits = int(Value)
            elif Option == "-e":
                # -e Minimum number of peptides for a new protein
                self.MinimumPeptides = int(Value)
            elif Option == "-w":
                # -w Summary page filename
                self.SummaryPagePath = Value
            elif Option == "-p":
                # -p p-value cutoff
                self.PValueCutoff = float(Value)
            elif Option == "-v":
                # -v Verbose output flag
                self.VerboseProteinReportCount = int(Value)
        # Error out, if we didn't see required options:
        if not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-r"):
            print "** Please specify database (-d) and results file (-r)"
            print UsageInfo
            sys.exit(1)
        # If neither -b nor -w was specified, assume they want a summary:
        if not OptionsSeen.has_key("-b") and not OptionsSeen.has_key("-w"):
            self.SummaryPagePath = os.path.join("ProteinSummary", "index.html")
            print "** Summary page will be written to '%s'; use -w to override this"%Summarizer.SummaryPagePath
        print "Summary page path:", self.SummaryPagePath
        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Summarizer = SummarizerClass()
    Summarizer.ParseCommandLine(sys.argv[1:])
    StartTime = time.clock()
    Summarizer.Main()
    EndTime = time.clock()
    print "ELAPSED:", EndTime - StartTime
