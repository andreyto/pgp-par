"""
Summarize the best peptide candidates across a run.
Focus is on finding mutations and PTMs.
"""
from Utils import *
import urllib
import MSSpectrum
import MakeImage
import Label
import os
import struct
import sys
import shutil
import traceback
import math
Initialize()      
Global.FixedMods = {"C":57.0518} #%%% IKKB TEMP

#DBPath = r"c:\ftproot\swt\inspect\database\ipi.HUMAN.v3.08.trie"
DBPath = r"database\drosoph_fbgn_040630.trie"
DBFile = open(DBPath, "rb")
#DBIndexPath = r"c:\ftproot\swt\inspect\database\ipi.HUMAN.v3.08.index"
DBIndexPath = r"database\drosoph_fbgn_040630.index"
DBIndexFile = open(DBIndexPath, "rb")
#GIFile = open("RequestGIIDs.txt","w")

class ResultClass:
    def __init__(self, Candidate):
        self.Candidate = Candidate
        self.Better = None
        self.Worse = None
        self.ExplainedBY = 0
        self.Score = -99
        self.PValue = 1
        self.MQScore = -99
        self.DeltaCNOther = -99
        self.ExplainedIntensity = 0
        self.ExplainedPeaks = 0

def GetNCBIURL(ProteinName):
    Bits = ProteinName.split("|")
    if len(Bits)>1:
        Bit = Bits[1]
        if Bit.find(":")!=-1:
            Bit = Bit.split(":")[-1]
        return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=%s"%Bit
        
def GetSafeFileName(Str):
    return Str.replace("*", "-").replace(">","").replace("#","+")

KnownPhosphoPeptides = {}
##File = open("KnownPhosphoPeptides.txt", "r")
##for FileLine in File.xreadlines():
##    KnownPhosphoPeptides[FileLine.strip().replace("*","")] = 1
##File.close()

class PeptideCandidateClass:
    def __init__(self, Peptide):
        self.Sequence = Peptide.Aminos
        self.TopHits = []
        self.RunnersUp = []
        self.TopCount = 0
        self.RunnerUpCount = 0
        self.SummaryScore = 0
        self.King = None
        self.Serfs = []
        self.Peptide = Peptide
        self.Name = Peptide.GetModdedName()
        #self.Parse(Aminos, Mods) #%%% obsolete
        self.BestDeltaScore = 0
        self.ProteinName = None
        self.FixedProteinName = None
        self.FixedProteinURL = None
        self.ProteinPosition = None
        
    def GetAllSerfs(self, Level = 1):
        #print " "*Level, "GetAllSerfs %s '%s'"%(self.Score, self.Name)
        List = self.Serfs[:]
        for Serf in self.Serfs:
            for SubSerf in Serf.GetAllSerfs(Level+2):
                if SubSerf not in List:
                    List.append(SubSerf)
        return List
    def Parse(self, Aminos, ModString): #%%% obsolete
        self.Name = ModString
        if ModString[1] == ".":
            ModString = ModString[2:-2]
        try:
            self.Peptide = GetPeptideFromModdedName(ModString)
        except:
            pass
        return
        Mods = {}
        ModLists = {}
        ModBits = ModString.split(",")
        ModCount = 0
        if Aminos[1]!=".":
            Aminos = "-." + Aminos + ".-"
        for Bit in ModBits:
            Bit = Bit.strip()
            if not Bit:
                continue
            ModCount += 1
            CBits = Bit.split(":")
            Amino = int(CBits[0])
            Name = CBits[1][:4].lower()
            Mods[Amino] = Mods.get(Amino,"") + Name
        Str = Aminos[:2]
        for Amino in range(2, len(Aminos)-2):
            Str += Aminos[Amino] + Mods.get(Amino-2, "")
        Name = Str + Aminos[-2:]
        self.Name = Name
        self.Peptide = GetPeptideFromModdedName(self.Name[2:-2])


class SummarizerClass:
    LineLimit = None # 30000
    MaxAllowedScoreGap = 500 #0.1
    MaxAllowedMQScoreGap = 0.5 #0.1
    MinSummaryScore = 2.0
    MinKingSummaryScore = 1.0
    MaxKingBestPValue = 0.1
    def __init__(self, PhosphoSummary):
        self.PhosphoSummary = PhosphoSummary
        self.SkipNoPhosphateKings = 0
        self.WritePeptidePages = 1 # default yes yesssss my preciousss
        self.SpectrumRoot =r"C:\Program Files\Apache Group\Apache2\cgi-bin"
        self.UseSerfsFlag = 1
        self.FinalSort = "score"
    def FixFileName(self, Name, FilePos):
        Stub = Name.split("/")[-1]
        UnderBits = Stub.split("_")
        Dir = "%s_%s"%(UnderBits[0], UnderBits[1])
        return "e:\\ftproot\\haixu\\drosprot\\reg\\%s\\%s:%s"%(Dir, Stub, FilePos)
    def ProduceCandidateList(self, FileName, SpectrumDir):
        self.RunName = os.path.splitext(os.path.split(FileName)[1])[0]
        print "Summarize:", FileName
        File = open(FileName, "r")
        self.PepCandidates = {} # Keys are sequences, values are PeptideCandidate objects
        self.ProteinCandidates = {}
        AnnotationsSeenForSpectrum = {}
        LineNumber = 0
        CurrentResult = None
        CurrentBestResult = None
        CurrentFileName = None
        HitsForSpectrum = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber%500 == 0:
                print ".",
            if self.LineLimit and LineNumber > self.LineLimit:
                break
            Bits = FileLine.split("\t")
            FileName = self.FixFileName(Bits[0], Bits[20])
            #FileName = Bits[0][5:] # skip "data/" #%%%
            if len(Bits)<2:
                continue
            RecordNumber = int(Bits[18])
            # We'll remember only 2 hits per spectrum
            if FileName != CurrentFileName:
                CurrentResult = None
                CurrentBestResult = None
                CurrentFileName = FileName
                HitsForSpectrum = 0
                AnnotationsSeenForSpectrum = {}
            try:
                Score = float(Bits[5])
            except:
                continue
            MQScore = Score
            if (MQScore < 0.0):
                continue
            
            # Grab fields from the line:            
            Aminos = Bits[2][2:-2]
            if AnnotationsSeenForSpectrum.get(Aminos, None):
                continue
            AnnotationsSeenForSpectrum[Aminos] = 1
            RecordNumber = int(Bits[18])
            ProteinName = Bits[3]
            
            PValue = float(Bits[15])
            #PValue = float(Bits[6])
            #if (PValue > 0.08):
            #    continue
            DeltaCNOther = float(Bits[17])
            Key = Aminos
            # Get the candidate (or build a new one):
            Candidate = self.PepCandidates.get(Key, None)
            if not Candidate:
                Pep = GetPeptideFromModdedName(Aminos)
                Candidate = PeptideCandidateClass(Pep)
                Candidate.RecordNumber = RecordNumber
                self.PepCandidates[Key] = Candidate
                if not self.ProteinCandidates.has_key(RecordNumber):
                    self.ProteinCandidates[RecordNumber] = []
                self.ProteinCandidates[RecordNumber].append(Candidate)
            Result = ResultClass(Candidate)
            Result.FileName = FileName
            Result.ScanNumber = Bits[1]
            Candidate.ProteinName = ProteinName
            Result.Score = Score
            Result.MQScore = MQScore
            Result.PValue = PValue
            Result.DeltaCNOther = DeltaCNOther
            if (CurrentBestResult):
                Result.Better = CurrentBestResult
                Candidate.RunnersUp.append(Result)
            else:
                CurrentBestResult = Result
                Candidate.TopHits.append(Result)
            if (CurrentResult):
                CurrentResult.Worse = Result
            CurrentResult = Result
            #print Result.Score, ProteinName
    def FindKingCandidates(self):
        print "Score all candidates:"
        for Candidate in self.PepCandidates.values():
            Candidate.BestPValue = 1
            Candidate.Score = 0
            PValues = []
            for TopHit in Candidate.TopHits:
                Candidate.BestPValue = min(Candidate.BestPValue, TopHit.PValue)
                PValues.append(TopHit.PValue)
            PValues.sort()
            for Index in range(len(PValues)):
                Candidate.Score -= math.log(max(0.00001, Candidate.BestPValue)) / float(Index+1)
        if self.UseSerfsFlag:
            print "Attach kings to serfs..."
            for (RecordNumber, CandidateList) in self.ProteinCandidates.items():
                SortedList = []
                for Candidate in CandidateList:
                    if Candidate.Score > self.MinSummaryScore:
                        SortedList.append((Candidate.Score,Candidate))
                SortedList.sort() # From BEST to WORST
                for IndexA in range(len(SortedList)):
                    CandidateA = SortedList[IndexA][-1]
                    for IndexB in range(IndexA+1, len(SortedList)):
                        CandidateB = SortedList[IndexB][-1]
                        OverlapFlag = self.CheckPeptideOverlap(CandidateA, CandidateB)
                        if OverlapFlag:
                            CandidateA.King = CandidateB
                            CandidateB.Serfs.append(CandidateA)
                            break # a candidate can only be in one serfs list at a time!
    def CheckPeptideOverlap(self, CandidateA, CandidateB):
        AminosA = CandidateA.Peptide.Aminos
        AminosB = CandidateB.Peptide.Aminos
        HalfPoint = len(AminosA)/2
        Pos = AminosB.find(AminosA[HalfPoint:])
        if Pos!=-1:
            for X in range(20):
                if X>HalfPoint:
                    return 1
                if X>Pos:
                    return 1
                if AminosA[HalfPoint-X]!=AminosB[Pos-X]:
                    break
        Pos = AminosB.find(AminosA[:HalfPoint])
        if Pos!=-1:
            for X in range(20):
                if X+HalfPoint>=len(AminosA):
                    return 1
                if X+Pos+HalfPoint>=len(AminosB):
                    return 1
                if AminosA[X+HalfPoint]!=AminosB[X+Pos+HalfPoint]:
                    break
        return 0
    def Summarize(self, FileName, SpectrumDir):
        self.SpectrumRoot = SpectrumDir
        self.ProduceCandidateList(FileName, SpectrumDir)
        self.FindKingCandidates()
        self.ReportKingCandidates()
    def GetFixedProteinName(self, Candidate):
        try:
            IPIID = Candidate.ProteinName.split("|")[0][4:].split(".")[0] # barftacular.
            if IPIID[0] == ":":
                IPIID = IPIID[1:]
            Description = ""
            URL = "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-e+[IPI-acc:%s]+-vn+2"%IPIID
            TargetFile = "IPI\\%s.txt"%IPIID
            urllib.urlretrieve(URL, TargetFile)
            IPIFile = open(TargetFile, "r")
            for IPILine in IPIFile.xreadlines():
                if IPILine[:3] == "DE ":
                    Description += IPILine[3:].strip().replace("\r"," ").replace("\n"," ")
            IPIFile.close()
            Candidate.FixedProteinName = "%s - %s"%(IPIID, Description)
        except:
            Candidate.FixedProteinName = Candidate.ProteinName
    def GetProteinResidue(self, Candidate):
        # Get the residue number where a peptide begins within its protein.
        DBIndexFile.seek(Candidate.RecordNumber * 92 + 8)
        FilePosition = struct.unpack("<i", DBIndexFile.read(4))[0]
        print "RecordNumber %s is at fileposition %s"%(Candidate.RecordNumber, FilePosition)
        DBFile.seek(FilePosition)
        Sequence = DBFile.read(400000)
        Position = Sequence.find(Candidate.Peptide.Aminos)
        #print Sequence[:500]
        Candidate.ProteinPosition = Position + 1
        print Position
    def ReportKingCandidates(self):
        # Open the summary file:
        self.SummaryDir = r"C:\Program Files\Apache Group\Apache2\dev\www\Summary%s"%self.RunName
        try:
            os.makedirs(self.SummaryDir)
        except:
            pass #traceback.print_exc()
        if self.PhosphoSummary:
            HTMLFile = open(os.path.join(self.SummaryDir, "PhosphoSummary.html"), "w")
            ExcelFile = open(os.path.join(self.SummaryDir, "PhosphoSummary.xls"), "w")
        else:
            HTMLFile = open(os.path.join(self.SummaryDir, "Summary.html"), "w")
            ExcelFile = open(os.path.join(self.SummaryDir, "Summary.xls"), "w")
        self.WriteHTMLHeader(HTMLFile)
        self.WriteExcelHeader(ExcelFile)
        print "Produce sorted king list:"
        # Produce the list of (king) candidates:
        SortedList = []
        for Candidate in self.PepCandidates.values():
            if Candidate.Score < self.MinKingSummaryScore:
                continue
            if Candidate.BestPValue > self.MaxKingBestPValue:
                continue
            if Candidate.King:
                continue
            if not Candidate.TopHits:
                continue
            if self.SkipNoPhosphateKings:
                Phos = 0
                if Candidate.Name.find("phos")!=-1:
                    Phos = 1
                else:
                    for Serf in Candidate.GetAllSerfs():
                        if Serf.Name.find("phos")!=-1:
                            Phos = 1
                            break
                if not Phos:
                    continue
            if self.FinalSort == "score":
                SortedList.append((Candidate.Score, Candidate))
            else:
                SortedList.append((Candidate.ProteinName, Candidate))
        SortedList.sort()
        SortedList.reverse()
        self.LineNumber = 1
        # Print a list of candidates to stdout:
        for Index in range(len(SortedList)):
            (Dummy, Candidate) = SortedList[Index]
            print "%s of %s: %s %s (% hits, pvalue %.4f)"%(Index, len(SortedList), -Candidate.Score, Candidate.Name,
                                                           len(Candidate.TopHits), Candidate.BestPValue)
        # Report each one:
        for Index in range(len(SortedList)):
            (Dummy, Candidate) = SortedList[Index]
            print "%s of %s: %s %s (% hits, pvalue %.4f)"%(Index, len(SortedList), -Candidate.Score, Candidate.Name,
                                                           len(Candidate.TopHits), Candidate.BestPValue)
            #print Dummy, len(Candidate.TopHits), Candidate.TopHits[0].PValue, Candidate.Name
            self.ReportCandidate(Candidate, HTMLFile, ExcelFile)
    def GetCandidateMaxScores(self, Candidate):
        Candidate.BestScore = -99999
        Candidate.BestMQScore = -999
        Candidate.BestDeltaScore = -99999
        Candidate.BestExplainedIntensity = 0
        Candidate.BestExplainedPeaks = 0
        Candidate.BestBY = 0
        Candidate.BestPValue = 1.0
        
        for Result in Candidate.TopHits:            
            Candidate.BestBY = max(Result.ExplainedBY, Candidate.BestBY)
            Candidate.BestScore = max(Result.Score, Candidate.BestScore)
            Candidate.BestPValue = min(Result.PValue, Candidate.BestPValue)
            Candidate.BestMQScore = max(Result.MQScore, Candidate.BestMQScore)
            Candidate.BestDeltaScore = max(Result.DeltaCNOther, Candidate.BestDeltaScore)
            Candidate.BestExplainedIntensity = max(Result.ExplainedIntensity , Candidate.BestExplainedIntensity)
            Candidate.BestExplainedPeaks = max(Result.ExplainedPeaks  , Candidate.BestExplainedPeaks)
    def GetResultFileName(self, Result):
        FileName = Result.FileName
        #(Root,Edge) = os.path.split(FileName)
        #FileName = os.path.split(FileName)[1] #%%% temp for ikkb / lens
        #FileName = FileName.replace("/", "\\")
        return FileName        
    def GetResultDetails(self, Candidate, Result):
        # Look up the spectrum, and get summary statistics:
        #(Root, Name) = os.path.split(Result.FileName)
        #(Root, Dir) = os.path.split(Root)
        FilePath = os.path.join(self.SpectrumRoot, Result.FileName)
        #FilePath = os.path.join(self.SpectrumRoot, Dir, Name)
        print FilePath
        FileName = self.GetResultFileName(Result)
        FileNameAminos = Candidate.Peptide.Aminos.replace("*", "-").replace(">","").replace("#","+")
        Spectrum = MSSpectrum.SpectrumClass()
        #FilePath = FixFilePath(FileName)
        File = Label.OpenAndSeekFile(FileName)
        Spectrum.ReadPeaksDTA(File)
        File.close()
        Spectrum.FilterPeaks(50, 6)
        Label.LabelSpectrum(Spectrum, Candidate.Peptide)
        Result.ExplainedIntensity = Spectrum.GetExplainedIntensity()
        Result.ExplainedPeaks = Spectrum.GetExplainedPeaks(50)
        (Present, Count, IonPercent, CutsPresent, CutCount, CutPercent) = Spectrum.GetExplainedIons(Candidate.Peptide, 125, 2000)
        Result.IonString = "%s/%s (%.1f%%)"%(Present, Count, IonPercent*100)
        Result.ExplainedBY = IonPercent
        return Spectrum
    def ReportCandidateDetails(self, Candidate, HTMLFile):
        HTMLFile.write("""<title>Peptide results for %s</title>
        <h1>%s</h1><br>
        <h2>Peptide results: %s (%s times top hit, %s times runner up)</h2>
        This page presents an annotated spectrum graph (left) for each spectrum where this
        peptide was the top hit, as well as an annoteted graph (right) for the next-best candidate.
        """%(Candidate.Name, Candidate.FixedProteinName, Candidate.Name, len(Candidate.TopHits), len(Candidate.RunnersUp)))
        SortedList = []
        for Result in Candidate.TopHits:
            SortedList.append((-Result.MQScore, Result))
        SortedList.sort()
        for (Dummy, Result) in SortedList[:25]: # at most 25 per page
            self.ReportResultDetails(Candidate, Result, HTMLFile)
        self.GetCandidateMaxScores(Candidate)
        if len(Candidate.RunnersUp):
            HTMLFile.write("<hr>")
            HTMLFile.write("""<h1>Runner-up results</h1>
            For the following spectra, this peptide was among the top 10 hits.  The annotated spectrum
            for the top hit is shown on the left, the annotation for this peptide on the right.
            """)
            for Result in Candidate.RunnersUp[:25]: 
                self.ReportResultDetails(Candidate, Result, HTMLFile)
    def ReportResultDetails(self, Candidate, Result, HTMLFile):
        Spectrum = self.GetResultDetails(Candidate, Result)
        #return #%%%% shortcut!
        FileName = self.GetResultFileName(Result)
        # Write out an image for this guy:
        if FileName.find("\\"):
            SlashBits = FileName.split("\\") # c:\mzxml\tripure1\foo.dta
        else:
            SlashBits = FileName.split("/")
        if len(SlashBits)>1:
            FileStub = SlashBits[-2] + "-" + SlashBits[-1].split(".")[0]
        else:
            FileStub = FileName.split(".")[0]
        FileAminos = Candidate.Peptide.GetModdedName().replace("*", "-").replace(">","").replace("#","+")
        FileStub = os.path.split(FileName)[1]
        if FileStub.find(":")>2:
            FileStub = FileStub.split(":")[0]
        ImageFileName = "%s.%s.png"%(FileAminos, FileStub)
        ImageFilePath = os.path.join(self.SummaryDir, ImageFileName)
        Maker = MakeImage.MSImageMaker()
        Maker.ConvertSpectrumToImage(Spectrum, ImageFilePath, Candidate.Peptide, Width=500, Height=380) #%%%
        # Write out an image for the runner-up or best match:
        if not (Result.Better or Result.Worse):
            HTMLFile.write("<hr><table><tr>")
            HTMLFile.write("""<td><img src="%s"></td>"""%ImageFileName)
            HTMLFile.write("</tr><tr><td>%s</td></tr><tr><td>This hit</td></tr></table>\n"%Candidate.Name)
        else:
            if Result.Better:
                OtherResult = Result.Better
                OtherCandidate = Result.Better.Candidate
                OtherPeptide = OtherCandidate.Peptide
            else:
                OtherResult = Result.Worse
                OtherCandidate = Result.Worse.Candidate
                OtherPeptide = OtherCandidate.Peptide
            OtherFileNameAminos = OtherPeptide.GetModdedName().replace("*", "-").replace(">","").replace("#","+")
            OtherImageFileName = "%s.%s.png"%(OtherFileNameAminos, FileStub)
            OtherImageFilePath = os.path.join(self.SummaryDir, OtherImageFileName)
            Maker = MakeImage.MSImageMaker()
            Label.LabelSpectrum(Spectrum, OtherPeptide)
            Maker.ConvertSpectrumToImage(Spectrum, OtherImageFilePath, OtherPeptide, Width=500, Height=380) #%%%
            if Result.Better:
                HTMLFile.write("<hr><table><tr>")
                HTMLFile.write("""<td><img src="%s"></td>"""%OtherImageFileName)
                HTMLFile.write("""<td><img src="%s"></td>"""%ImageFileName)
                HTMLFile.write("""</tr><tr><td><center>%s</center></td><td><center>%s</center></td>"""%(OtherCandidate.Name, Candidate.Name))
                HTMLFile.write("""</tr><tr><td><center>Top hit</center></td><td><center>This peptide</center></td>""")
                HTMLFile.write("""</tr><tr><td><center>Score %.4f, pvalue %.4f</center></td><td><center>Score %.4f, pvalue %.4f</center></td>"""%\
                (OtherResult.Score, OtherResult.PValue, Result.Score, Result.PValue))
                HTMLFile.write("""</tr></table>\n""")
            else:
                HTMLFile.write("<hr><table><tr>")
                HTMLFile.write("""<td><img src="%s"></td>"""%ImageFileName)
                HTMLFile.write("""<td><img src="%s"></td>"""%OtherImageFileName)
                HTMLFile.write("""</tr><tr><td><center>%s</center></td><td><center>%s</center></td>"""%(Candidate.Name, OtherCandidate.Name))
                HTMLFile.write("""</tr><tr><td><center>This peptide</center></td><td><center>Runner-up</center></td>""")
                HTMLFile.write("""</tr><tr><td><center>Score %.4f, pvalue %.4f</center></td><td><center>Score %.4f, pvalue %.4f</center></td>"""%\
                (Result.Score, Result.PValue, OtherResult.Score, OtherResult.PValue))
                HTMLFile.write("""</tr></table>\n""")
        # And write some score stuff:
        HTMLFile.write("Spectrum %s scan %s:<br>\np-value %s, score %s, delta-score %s, mqscore %.2f, intensity %.2f%%, peaks %.2f%%, b/y %s<br>\n"%(Result.FileName, Result.ScanNumber,
            Result.PValue, Result.Score, Result.DeltaCNOther, Result.MQScore,
            100*Result.ExplainedIntensity, 100*Result.ExplainedPeaks, Result.IonString))
    def ReportCandidate(self, Candidate, HTMLFile, ExcelFile, IsSerf = 0):
        # Compute stuff:
        print Candidate.Name
        self.GetFixedProteinName(Candidate)
        self.GetProteinResidue(Candidate)
        
        Candidate.TrypticTermini = 2
        if (Candidate.Name[0] not in ("R","K","-","*")):
            Candidate.TrypticTermini -= 1
        if (Candidate.Peptide.Aminos[-1] not in ("R","K")):
            Candidate.TrypticTermini -= 1
        Candidate.KnownFlag = "" #default
        if Candidate.Name.find("pho")!=-1 and KnownPhosphoPeptides.has_key(Candidate.Peptide.Aminos):
            Candidate.KnownFlag = "Yes"
        # Individual results:
        FileName = os.path.join(self.SummaryDir, "%s.html"%GetSafeFileName(Candidate.Name))
        CandidateResultURL = "%s.html"%GetSafeFileName(Candidate.Name)
        if self.WritePeptidePages:
            CandidateHTMLFile = open(FileName, "w")
            self.ReportCandidateDetails(Candidate, CandidateHTMLFile)
            CandidateHTMLFile.close()
        else:
            for Index in range(len(Candidate.TopHits)):
                Result = Candidate.TopHits[Index]
                #print "Result %d of %d..."%(Index, len(Candidate.TopHits))
                self.GetResultDetails(Candidate, Result)
            #print "Get candidate max scores..."
            self.GetCandidateMaxScores(Candidate)
        #print "write to html:"
        # Report stuff:
        Bold = 0
        Color = "#FFFFFF"
        # Highlight phosphopeptide rows in color, if we're doing phoshpopeptide summary:
        if self.PhosphoSummary:
            if Candidate.Name.find("pho")!=-1:
                Bold = 1
                if Candidate.KnownFlag:
                    Color = "#66FF99"
                else:
                    Color = "#FF66CC"
            else:
                Bold = 0
        Color = "#FFFFFF" #%%%
        HTMLFile.write("<tr bgcolor=%s>\n"%Color)
        HTMLFile.write("<td>%s</td>\n"%self.LineNumber)
        if IsSerf:
            HTMLFile.write("""<td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;""")
        else:
            HTMLFile.write("""<td>""")
        URL = GetNCBIURL(Candidate.ProteinName)
        Name = Candidate.FixedProteinName
        #HTMLFile.write("%s\t"%Candidate.Peptide.GetModdedName())
        NameStr = """<a href="%s">%s</a></td>"""%(CandidateResultURL, Candidate.Peptide.GetModdedName())
        HTMLFile.write("%s"%NameStr)
        HTMLFile.write("""<td><a href="%s">%s</a></td>\n"""%(URL, Name))
        HTMLFile.write("<td>%s</td>\n"%len(Candidate.TopHits))
        HTMLFile.write("<td>%s</td>\n"%len(Candidate.RunnersUp))
        HTMLFile.write("<td>%s</td>\n"%Candidate.ProteinPosition)
        HTMLFile.write("<td>%s</td>\n"%Candidate.BestPValue)
        HTMLFile.write("<td>%s</td>\n"%Candidate.BestScore)
        HTMLFile.write("<td>%s</td>\n"%Candidate.BestDeltaScore)
        HTMLFile.write("<td>%s</td>\n"%Candidate.TrypticTermini)
        HTMLFile.write("<td>%.2f%%</td>\n"%(Candidate.BestExplainedIntensity*100))
        HTMLFile.write("<td>%.2f%%</td>\n"%(Candidate.BestExplainedPeaks*100))
        HTMLFile.write("<td>%.2f</td>\n"%(Candidate.BestBY*100))
        HTMLFile.write("<td>%s</td>\n"%(Candidate.TopHits[0].FileName))
        HTMLFile.write("</tr>\n")
        HTMLFile.flush()
        #print "Excel file write:"
        ExcelStr = "%s\t%s\t%s\t%s\t"%(self.LineNumber,Candidate.Name, Candidate.FixedProteinName, len(Candidate.TopHits))
        ExcelStr += "%s\t%s\t%s\t%s\t"%(len(Candidate.RunnersUp), Candidate.BestPValue, Candidate.BestScore, Candidate.BestDeltaScore)
        ExcelStr += "%s\t%s\t%s\t%s\t"%(Candidate.TrypticTermini, Candidate.BestExplainedIntensity*100, Candidate.BestExplainedPeaks*100, Candidate.BestBY*100)
        ExcelStr += "%s\t"%(Candidate.KnownFlag)
        ExcelFile.write(ExcelStr+"\n")
        ExcelFile.flush()
        self.LineNumber += 1
        # Serfs:
        if not Candidate.King:
            AllSerfs = Candidate.GetAllSerfs()
            SortedList = []
            for Serf in AllSerfs:
                SortedList.append((-Serf.Score, Serf))
            SortedList.sort()
            for (Dummy, Serf) in SortedList:
                if len(Serf.TopHits) and Serf.Score > self.MinSummaryScore:
                    self.ReportCandidate(Serf, HTMLFile, ExcelFile, 1)
        #print "Candidate complete."
    def WriteExcelHeader(self, ExcelFile):
        Header = "Line\tPeptide\tProtein\tTopHits\tRunner-up Hits\tBest p-value\tBest score\tBest delta-score\tTryptic Termini\tBest explained intensity\tBest explained peaks\tBest explained b/y\tFileName\t"
        ExcelFile.write(Header+"\n")
    def WriteHTMLHeader(self, HTMLFile):
        HTMLFile.write("""<title>%s results</title>
        <small><small>
        <table><tr>
        <td>Line</td>
        <td><a href="Help.html#FFPeptide">Peptide</a></td>
        <td><a href="Help.html#FFProtein">Protein</a></td>
        <td><a href="Help.html#FFTopHits">Top<br>Hits</a></td>
        <td><a href="Help.html#FFRunnerUpHits">Runner-up<br>Hits</a></td>
        <td><a href="Help.html#FFPosition">Protein<br>residue</a></td>
        <td><a href="Help.html#FFBestPValue">Best hit<br>p-value</a></td>
        <td><a href="Help.html#FFBestScore">Best<br>Score</a></td>
        <td><a href="Help.html#FFBestDeltaScore">Best<br>delta-score</a></td>
        <td><a href="Help.html#FFTrypticTermini">Tryptic<br>termini</a></td>
        <td><a href="Help.html#FFBestExplainedIntensity">Best<br>Explained<br>Intensity</a></td>
        <td><a href="Help.html#FFBestExplainedPeaks">Best<br>Explained<br>Peaks</a></td>
        <td><a href="Help.html#FFBestBY">Best<br>b/y</a></td>
        <td>FileName</td>
        </tr>
        """%self.RunName)
        

try:
    import psyco
    psyco.full()
except:
    print "<no psyco>"

# Summary of peptides in general:
#Summarizer = SummarizerClass(0)
#Summarizer.ProduceCandidateList(sys.argv[1], sys.argv[2])
#Summarizer.WriteFindings(0)

# Phosphopeptide summary:
# SummarizeResults.py AfCSOutput\ds030808_01_Exp236-D-C.txt AfCS\ds030808_01_Exp236-D-C
Summarizer = SummarizerClass(0)
Summarizer.SkipNoPhosphateKings = 0 # (special-purpose for AfCS)
Summarizer.MinKingSummaryScore = 9.0
Summarizer.MinSummaryScore = 1.0
Summarizer.WritePeptidePages = 1
Summarizer.UseSerfsFlag = 1 # (Set TRUE to group all peptides by locus)
Summarizer.FinalSort = "score" # can be PROTEIN or SCORE

# FileName, SpectrumDir
Summarizer.Summarize(sys.argv[1], sys.argv[2]) # FileName, SpectrumDir
#Summarizer.WriteFindings(1)

