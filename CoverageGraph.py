"Simple plot of how well-covered by search results each residue in a protein is"
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

import sys
import os
from Utils import *
Initialize()

class Colors:
    "Color scheme for web display.  Colorful."
    White = (255,255,255)
    Green = (0,255,0)
    Blue = (0,0,255)
    PaleBlue = (10,10,80)
    Red = (255,0,0)
    Grey = (155,155,155)
    ResidueNumber = (155, 155, 155)
    #Grey = (0,0,0)
    Background = (255, 255, 255)
    Coverage = (155, 155, 155)
    Mod = [(0,0,155),(0,155,155),(155,0,0),(155,0,155),(155,155,0),
           (0,0,55),(0,55,55),(55,0,0),(55,0,55),(55,55,0)]


class CoverageGraphMaker:
    Header = 20
    Footer = 10
    SideMargin = 10
    def __init__(self):
        pass
    def PlotCoverage(self, BlindAnnotationFile, FASTAFile, OutputDir):
        try:
            os.makedirs(OutputDir)
        except:
            pass
        self.BlindAnnotationFile = BlindAnnotationFile
        self.OutputDir = OutputDir
        File = open(FASTAFile, "r")
        self.Sequence = ""
        for FileLine in File.xreadlines():
            if FileLine[0] == ">":
                if self.Sequence:
                    self.GetCoverageForProtein(ProteinName, self.Sequence)
                    self.Sequence = ""
                ProteinName = FileLine[1:].strip()
            else:
                self.Sequence += FileLine.strip()
        if self.Sequence:
            self.GetCoverageForProtein(ProteinName, self.Sequence)
    def GetCoverageForProtein(self, ProteinName, Sequence):
        CoverageList = [] # Tuple of the form (StartPos, EndPos, Modifications)
        File = open(self.BlindAnnotationFile, "r")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            RawSequence = Bits[2]
            Pos = Sequence.find(RawSequence)
            if Pos==-1:
                continue
            # Convert to 1-based numbering:
            Coverage = [Pos + 1, Pos + 1 + len(RawSequence), {}]
            Peptide = GetPeptideFromModdedName(Bits[3])
            Coverage[2] = Peptide.Modifications
            Peptide.Modifications = None 
            CoverageList.append(Coverage)
        CoverageList.sort()
        self.PlotProtein(ProteinName, Sequence, CoverageList)
    def PlotProtein(self, ProteinName, Sequence, CoverageList):
        Width = len(Sequence) + self.SideMargin*2
        Height = len(CoverageList) + self.Header + self.Footer
        CoverageImage = Image.new("RGB", (Width, Height), Colors.Background)  # mode, size, [startcolor]
        Draw = ImageDraw.Draw(CoverageImage)
        Draw.line((self.SideMargin, self.Header, Width - self.SideMargin, self.Header), Colors.ResidueNumber)
        for X in range(10, len(self.Sequence), 20):
            Draw.text((X + self.SideMargin, 2), "%s"%X, Colors.ResidueNumber)
            Draw.line((X + self.SideMargin, self.Header, X + self.SideMargin, self.Header - 2), Colors.ResidueNumber)
        for Index in range(len(CoverageList)):
            Coverage = CoverageList[Index]
            Y = self.Header + Index
            Draw.line((self.SideMargin + Coverage[0], Y, self.SideMargin + Coverage[1], Y), Colors.Coverage)
            for (Key, List) in Coverage[2].items():
                X = self.SideMargin + Coverage[0] + Key
                Color = Colors.Mod[List[0].Mass % 10]
                Draw.line((X, Y, X+1, Y), Color)
        SafeProteinName = ProteinName.replace("|",".").replace("/",".").replace("\\",".").replace("[",".").replace("]",".").replace("\x01","#")
        OutputFileName = os.path.join(self.OutputDir, SafeProteinName) + ".png"
        CoverageImage.save(OutputFileName, "png")

if __name__ == "__main__":
    Bob = CoverageGraphMaker()
    BlindAnnotationFile = sys.argv[1]
    FASTAFile = os.path.join("Database", sys.argv[2])
    OutputDir = "CoveragePlots"
    Bob.PlotCoverage(BlindAnnotationFile, FASTAFile, OutputDir)

    