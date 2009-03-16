"""
MakeImage.py takes a spectrum labeled with ion types (built by Label.py) and produces a graph.
"""
from Utils import *
import traceback
import MSSpectrum
import math
import traceback

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
    traceback.print_exc()
    print "WARNING: Python Imaging Library (PIL) not installed.\n  Image creation is NOT available."
    Image = None


class WebColors:
    "Color scheme for web display.  Colorful."
    White = (255,255,255)
    Green = (0,255,0)
    Blue = (0,0,255)
    PaleBlue = (10,10,80)
    Red = (255,0,0)
    Grey = (155,155,155)
    #Grey = (0,0,0)
    Background = (255, 255, 255)
    Peak = (199, 199, 199)
    #Peak = (0,0,0)
    AnnotatedPeak = (0, 0, 0)
    Axis = (155, 155, 155)
    LabeledPeak = (0, 0, 0)
    PeakLabel = (200, 0, 0)
    BSeries = (55,55,200)
    BSeriesPale = (155,155,255)
    YSeries = (155,155,55)
    YSeriesPale = (200,200,100)

class PrintingColors:
    "Color scheme for printing.  Black-and-white, dark shades."
    White = (255,255,255)
    Green = (0,0,0)
    Blue = (0,0,0)
    PaleBlue = (80,80,80)
    Red = (0,0,0)
    Grey = (80,80,80)
    Background = (255, 255, 255)
    Peak = (200, 200, 200)
    AnnotatedPeak = (80, 80, 80)
    Axis = (80, 80, 80)
    LabeledPeak = (10, 10, 10)
    PeakLabel = (0, 0, 0)
    BSeries = (0,0,0)
    BSeriesPale = (180,180,180)
    YSeries = (0,0,0)
    YSeriesPale = (180,180,180)
    
#Colors = PrintingColors
Colors = WebColors

def SetColors(PrintFlag):
    global Colors
    if PrintFlag:
        Colors = PrintingColors
    else:
        Colors = WebColors
    print "COLORS SET!", Colors.YSeries
    
def RoundOff(Int):
    "Round an integer to the nearest power of ten"
    PowerOfTen = 10
    while (1):
        if PowerOfTen>Int:
            return Int
        Radix = Int % PowerOfTen
        if Radix == 0:
            PowerOfTen *= 10
            continue
        if Radix > PowerOfTen/2:
            return (Int-Radix) + PowerOfTen
        else:
            return (Int-Radix)
        
    
def GetTickWidth(MassWidth):
    "Look for a nice round number that divides the width into about 20 pieces."
    IdealSliceCount = 20
    # If 6-15, round to 10; if 26-35, round to 20...  If 151...
    #print "GetTickWidth:", MassWidth
    Width = int(round(MassWidth / IdealSliceCount))
    SliceCount = MassWidth / Width
    while (1):
        # Try rounding off another digit:
        RoundedWidth = RoundOff(Width)
        if RoundedWidth==Width:
            break # we can't round off any more!
        RoundedSliceCount = MassWidth / RoundedWidth
        if (RoundedSliceCount < IdealSliceCount / 2) or (RoundedSliceCount > IdealSliceCount * 2):
            break # Slices are too skinny or too fat
        Width = RoundedWidth
    return Width


class BaseImageMaker:
    """
    Graph generating class.  Not very MS-specific.  Subclassed by the spectrum plotter.
    """
    LeftPad = 30
    RightPad = 3
    UpperPad = 50
    LowerPad = 20
    def __init__(self, Width = 600, Height = 400):
        """
        Width is the total image width, in pixels.  The plot body ranges from
        self.LeftPad to self.Width.
        Height is the total image height, in pixels.  The plot body ranges from
        self.UpperPad to (self.Height-self.LowerPad)
        """
        self.MinX = 0
        self.MaxX = 100
        self.MinY = 0
        self.MaxY = 1
        self.Width = Width
        self.Height = Height
        self.BaseLine = self.Height - self.LowerPad
        self.YBreak = None
    def GetNiceTickMark(self, Width, GoodTickCount):
        """
        Given a width (axis size), compute a good interval for major tick marks.
        We want around 10 ticks, but we round the tick size up or down in order
        to get cleaner numbers.
        """
        if (Width <= 0):
            return 0.1
        Tick = Width / float(GoodTickCount)
        MinimumPower = int(math.log10(Tick)) - 1
        PossibleTicks = []
        for Power in range(MinimumPower, MinimumPower+5):
            PossibleTicks.append(10**Power)
            PossibleTicks.append(10**Power * 2)
            PossibleTicks.append(10**Power * 5)
        for Index in range(len(PossibleTicks) - 1):
            if PossibleTicks[Index+1] > Tick:
                return PossibleTicks[Index+1]
        return 0.1 # hacky fall-back case!
    def GetValueName(self, Value):
        """
        Take a number - potentially a very large or small one - and format it so as
        to use few characters, so as to make a usable axis label.
        1.2002 -> 1.2
        0.149 -> 0.15
        10000 -> 1e5
        -0.00005 -> -5e-5
        """
        if Value == 0:
            return "0"
        if Value < 0:
            Sign = "-"
            Value *= -1
        else:
            Sign = ""
        if Value > 1000:
            Exp = int(math.log10(Value))
            Abscissa = 10**(math.log10(Value) - int(math.log10(Value)))
            return "%s%.1fe%d"%(Sign,Abscissa, Exp)
        if Value < 0.01:
            Exp = math.floor(math.log10(Value))
            Abscissa = 10**(math.log10(Value) - Exp)
            return "%s%.1fe%d"%(Sign,Abscissa, Exp)
        if abs(Value) < 1:
            return "%s%.2f"%(Sign, Value)
        else:
            return "%s%.1f"%(Sign, Value)
    def DrawYAxis(self):
        if self.YBreak:
            self.DrawYAxisHelper(self.MinLowY, self.MaxLowY, 8)
            self.DrawYAxisHelper((self.MinHighY + self.MaxHighY) / 2.0, self.MaxHighY, 1)
            self.Draw.line((self.LeftPad, self.UpperPad, self.LeftPad, self.YBreak - 2), Colors.Axis)
            self.Draw.line((self.LeftPad - 4, self.YBreak + 4, self.LeftPad + 4, self.YBreak), Colors.Axis)
            self.Draw.line((self.LeftPad - 4, self.YBreak, self.LeftPad + 4, self.YBreak - 4), Colors.Axis)
            self.Draw.line((self.LeftPad, self.YBreak + 2, self.LeftPad, self.BaseLine), Colors.Axis)
        else:
            self.DrawYAxisHelper(self.MinY, self.MaxY, 10)
            self.Draw.line((self.LeftPad, self.UpperPad, self.LeftPad, self.BaseLine), Colors.Axis)
    def DrawYAxisHelper(self, MinVal, MaxVal, TickCount):
        "Draw the y-axis (including tick-marks and labels)"
        TickSize = self.GetNiceTickMark(MaxVal - MinVal, TickCount)
        ##print "Y axis: %s ticks from %s to %s; tick size %s"%(TickCount, MinVal, MaxVal, TickCount)
        MaxPowers = int(math.log10(MaxVal)) + 1
        LastY = None
        try:
            IntensityLevel = MinVal
            while IntensityLevel < MaxVal:
                Y = self.GetY(IntensityLevel)
                self.Draw.line((self.LeftPad-2,Y,self.LeftPad,Y),Colors.Axis)
                Str = self.GetValueName(IntensityLevel)
                self.Draw.text((0, Y-5), Str, Colors.Axis, font = TheFont)
                LastY = Y
                IntensityLevel += TickSize
        except:
            traceback.print_exc()
            print 0, self.MinY, self.MaxY, TickSize
            raise
    def BreakY(self, SmallestOfBig, BiggestOfSmall):
        self.MinLowY = self.MinY
        self.MaxLowY = BiggestOfSmall * 1.05
        self.MinHighY = max(self.MaxLowY * 1.001, SmallestOfBig * 0.95)
        self.MaxHighY = self.MaxY
        self.BrokenY = 1
        self.YBreak = int((self.Height - self.UpperPad - self.LowerPad) * 0.2)
        #print "Y break is:", self.YBreak
    def GetY(self, YValue):
        if self.YBreak:
            if YValue > self.MaxLowY:
                YPercent = (YValue - self.MinHighY) / max(1, self.MaxHighY - self.MinHighY)
                return self.YBreak - int((self.YBreak - self.UpperPad)*YPercent)
            else:
                YPercent = YValue / max(1, self.MaxLowY - self.MinLowY)
                return self.BaseLine - int((self.BaseLine - self.YBreak)*YPercent)
        YPercent = YValue / float(max(1, self.MaxY - self.MinY))
        return self.BaseLine - int((self.BaseLine - self.UpperPad) * YPercent)
    def GetX(self, XValue):
        XPercent = (XValue - self.MinX) / max(1, (self.MaxX - self.MinX))
        TotalWidth = self.Width - (self.LeftPad + self.RightPad)
        return self.LeftPad + int(XPercent * TotalWidth)
    def DrawTickMarks(self):
        TickPos = 0
        while TickPos < self.MaxX-1.0:
            if TickPos < self.MinX:
                TickPos += self.TickWidth
                continue
            X = self.GetX(TickPos)
            # Draw text, unless it would go over the edge:
            TextWidth = len(str(TickPos))*6
            TextX = X - TextWidth/2
            if TextX < self.Width - TextWidth:
                self.Draw.line((X, self.BaseLine, X, self.BaseLine+3), Colors.Axis)
                self.Draw.text((TextX, self.BaseLine+2), str(TickPos), Colors.Axis, font = TheFont)
            TickPos += self.TickWidth
    
class MSImageMaker(BaseImageMaker):
    def __init__(self, *args, **kw):
        # Some options:
        self.YBreakThreshold = 0.3333
        #self.Labels = {} Label -> (PeakX, PeakY, Label, PeakIntensity, IntensityRank)
        self.IntensityRank = {} # PeakIndex -> Rank
        BaseImageMaker.__init__(self, *args, **kw)
    def ConvertPeakAnnotationToImage(self, PeakAnnotationList, OutputFileName, Peptide = None, Width = 600, Height = 400):
        if not Image: # catch for no PIL
            return None
        self.Width = Width
        self.Height = Height
        self.BaseLine = self.Height - self.LowerPad
        self.PeakAnnotationList = PeakAnnotationList
        self.GetPeakDemographics(Peptide) #computes min, max, intensityrank
        self.PlotImage = Image.new("RGB", (Width, Height), Colors.Background)  # mode, size, [startcolor]
        self.Draw = ImageDraw.Draw(self.PlotImage)
        self.RoofLine = Height * 0.5
        MassWidth = self.Width - (self.LeftPad + self.RightPad)
        # Draw baseline
        self.Draw.line((self.LeftPad, self.BaseLine, Width - self.RightPad, self.BaseLine), Colors.Axis)
        # Draw x axis tickmarks (and labels):
        self.TickWidth = 200
        self.DrawTickMarks()
        # Draw y axis:
        self.DrawYAxis()
        # Draw peaks, with labels
        self.DrawPeaks()
        if Peptide:
            self.DrawBSeries(Peptide)
            self.DrawYSeries(Peptide)
        self.DrawPeakLabels()
        self.PlotImage.save(OutputFileName, "png")
    def GetPeakDemographics(self, Peptide = None):
        """
        Because this version uses a list of peak annotations given by
        PyInspect, and not peak objects, they don't come with an associated
        rank.  here I do a quick bubble sort and rank stuff (Ali, I wrote this bubble sort)
        """
        Intensities = []
        self.MinX = 1000
        self.MaxX = 0
        # Sort peak, from most to least intense:
        PeaksByIntensity = []
        for PeakIndex in range(len(self.PeakAnnotationList)):
            Tuple = self.PeakAnnotationList[PeakIndex]
            PeaksByIntensity.append((Tuple[1], Tuple[0], PeakIndex))
            ## set min and max masses
            if Tuple[0] > self.MaxX:
                self.MaxX = Tuple[0]
            if Tuple[0] < self.MinX:
                self.MinX = Tuple[0]
            self.MaxY = max(self.MaxY, Tuple[1])
        PeaksByIntensity.sort()
        PeaksByIntensity.reverse()
        ###
        self.IntensityRank = [None] * len(self.PeakAnnotationList)
        for IntensityRank in range(len(PeaksByIntensity)):
            self.IntensityRank[Tuple[-1]] = IntensityRank
        ## possibly reset x min and max
        if Peptide:
            self.MaxX = max(self.MaxX, Peptide.Masses[-1] + 10)
            self.MinX = 0
        FullMass = self.MaxX - self.MinX
        # Move the left and right edges a little bit, so that peaks don't
        # hit the edges of the image.
        self.MinX = max(0, self.MinX - FullMass * 0.05)
        self.MaxX = self.MaxX + FullMass * 0.05
        self.MaxY = max(10.0, self.MaxY)
        ##############
        ## do some Y breakage
        Intensity1 = PeaksByIntensity[0][0]
        if len(PeaksByIntensity) < 2:
            return # why we would have a spectrum with only one peak, I will never know.
        Intensity2 = PeaksByIntensity[1][0]
        if (Intensity2 / Intensity1 < self.YBreakThreshold):
            self.BreakY(Intensity1, Intensity2)
            return
        if len(PeaksByIntensity) < 3:
            return
        Intensity3 = PeaksByIntensity[2][0]
        if Intensity3 / Intensity2 < self.YBreakThreshold:
            self.BreakY(Intensity2, Intensity3)
            return
    def DrawPeaks(self):
        self.Labels = {}
        MaxIntensity = 0
        # First, no-ion-type peaks:
        for PeakTuple in self.PeakAnnotationList:
            (Mass, Intensity, Label, AminoIndex) = PeakTuple
            MaxIntensity = max(MaxIntensity, Intensity)
            if Label: # don't draw in grey anything which has been labeled
                continue
            PeakX = self.GetX(Mass)
            PeakY = self.GetY(Intensity)
            self.Draw.line((PeakX, PeakY, PeakX, self.BaseLine), Colors.Peak)
        MinLabelIntensity = MaxIntensity * 0.25 # any annotated peak above this threshold receives a text label, even if it's a neutral loss.
        # Next, peaks with assigned ions:
        for PeakIndex in range(len(self.PeakAnnotationList)):
            (Mass, Intensity, Label, AminoIndex) = self.PeakAnnotationList[PeakIndex]
            if not Label: #skip all the unlabeled peaks this round
                continue
            PeakX = self.GetX(Mass)
            PeakY = self.GetY(Intensity)
            self.Draw.line((PeakX, PeakY, PeakX, self.BaseLine), Colors.AnnotatedPeak)
            TextLabelPeakNames = ("B", "Y", "GlcNAc", "Y2", "B2")
            if Label in TextLabelPeakNames:  
                PeptideIndex = AminoIndex
                if Label in ("B", "Y"):
                    Label = Label.lower()
                if Label in ("B2", "Y2"):
                    Label = "%s2"%Label[0].lower() 
                if AminoIndex != None:
                    Label = "%s %s"%(Label, AminoIndex)
                OldLabelInfo = self.Labels.get(Label, None)
                if OldLabelInfo!=None and OldLabelInfo[-2] > Intensity:
                    continue
                self.Labels[Label] = (PeakX, PeakY, Label, Intensity, self.IntensityRank[PeakIndex])
            if Label[0] == "M":
                ##special case for phorphorylated spectra. the label
                ## Parent loss has been changed to M-p or M-p-h2o
                OldLabelInfo = self.Labels.get(Label, None)
                if OldLabelInfo != None and OldLabelInfo[-2] > Intensity:
                    continue
                self.Labels[Label] = (PeakX, PeakY, Label, Intensity, self.IntensityRank[PeakIndex])
    def CollideRectangles(self, X1, Y1, X2, Y2, Rectangles):
        for (CX1, CY1, CX2, CY2) in Rectangles:
            if CX1 <= X1 <= CX2 and CY1 <= Y1 <= CY2:
                return (CX1, CY1, CX2, CY2)
            if CX1 <= X2 <= CX2 and CY1 <= Y1 <= CY2:
                return (CX1, CY1, CX2, CY2)
            if CX1 <= X1 <= CX2 and CY1 <= Y2 <= CY2:
                return (CX1, CY1, CX2, CY2)
            if CX1 <= X2 <= CX2 and CY1 <= Y2 <= CY2:
                return (CX1, CY1, CX2, CY2)
        return None
    def DrawPeakLabels(self):
        # Sort labels by priority.  b and y take precedence over all else;
        # intense peaks take precedence over others.
        SortedLabels = []
        for (X, Y, Label, Intensity, Rank) in self.Labels.values():
            if Rank < 10:
                Priority = 0
            else:
                Priority = 1
            SortedLabels.append([Priority, Y, X, Label, Intensity, None, None])
        SortedLabels.sort()
        SortedLabels = SortedLabels[:25]
        DirtyRectangles = []
        for List in SortedLabels:
            (IsBY, Y, X, Label, Intensity, Dummy, Dummy) = List
            (Width, Height) = self.Draw.textsize(Label, font = TheFont)
            Height *= 2 # for superscript and subscript
            X1 = X - Width/2
            X2 = X + Width/2
            Y1 = Y - Height
            Y2 = Y
            Tuple = self.CollideRectangles(X1, Y1, X2, Y2, DirtyRectangles)
            if Tuple == None:
                List[5] = (X1, Y1, X2, Y2)
                DirtyRectangles.append((X1, Y1, X2, Y2))
                continue
            (CX1, CY1, CX2, CY2) = Tuple
            # Try moving this label off to the side:
            if (X1 + X2) / 2 < (CX1 + CX2) / 2:
                Move = (X2 - CX1) + 1
                X1 -= Move
                X2 -= Move
                Y1 -= 5
                Y2 -= 5
            else:
                Move = (CX2 - X1) + 1
                X1 += Move
                X2 += Move
                Y1 -= 5
                Y2 -= 5
            Tuple = self.CollideRectangles(X1, Y1, X2, Y2, DirtyRectangles)
            if Tuple == None:
                List[5] = (X1, Y1, X2, Y2)
                DirtyRectangles.append((X1, Y1, X2, Y2))
                List[6] = ((X1 + X2) / 2, Y2, X, Y)
                continue
        for Index in range(len(SortedLabels)-1, -1, -1):
            List = SortedLabels[Index]
            if List[5]:
                (X1, Y1, X2, Y2) = List[5]
                self.Draw.line((List[2], List[1], List[2], self.BaseLine), Colors.LabeledPeak) # color the peak
        for Index in range(len(SortedLabels)-1, -1, -1):
            List = SortedLabels[Index]
            if List[5]:
                (X1, Y1, X2, Y2) = List[5]
                # Most peaks are drawn using superscripts and subscripts:
                PeakName = List[3]
                if PeakName == "GlcNAc":
                    self.Draw.text((X1, Y1 + 5), PeakName, Colors.PeakLabel, font = TheFont)
                else:
                    self.Draw.text((X1, Y1 + 5), PeakName[0], Colors.PeakLabel, font = TheFont)
                    NumIndex = len(PeakName) - 1
                    while PeakName[NumIndex] in ("0123456789"):
                        NumIndex -= 1
                    if NumIndex > 0:
                        SuperScript = PeakName[1:NumIndex+1].strip()
                    else:
                        SuperScript = ""
                    if SuperScript:
                        self.Draw.text((X1+7, Y1+10), SuperScript, Colors.PeakLabel, font = TheFont)
                        self.Draw.text((X1+7, Y1), PeakName[NumIndex+1:], Colors.PeakLabel, font = TheFont)
                    else:
                        self.Draw.text((X1+7, Y1+5), PeakName[NumIndex+1:], Colors.PeakLabel, font = TheFont)
            # Draw the dotted line from the label to its peak:                
            if List[6]:
                self.Draw.line(List[6], Colors.LabeledPeak)
    def DrawDottedLine(self, X1, Y1, X2, Y2, Color):
        Distance = math.sqrt((X2-X1)**2 + (Y2-Y1)**2)
        if Distance == 0:
            return
        DX = (X2-X1)/Distance
        DY = (Y2-Y1)/Distance
        OldLineLength = 0
        Dot = 1
        while (1):
            LineLength = min(Distance, OldLineLength + 5)
            XA = int(X1 + DX*OldLineLength)
            XB = int(X1 + DX*LineLength)
            YA = int(Y1 + DY*OldLineLength)
            YB = int(Y1 + DY*LineLength)
            if Dot:
                self.Draw.line((XA, YA, XB, YB), Color)
            OldLineLength = LineLength
            Dot = not Dot
            if (LineLength == Distance):
                break
    def DrawBSeries(self, Peptide):
        BHeight = 17
        if getattr(Peptide, "Seed", None):
            SeedIndex = Peptide.Aminos.rfind(Peptide.Seed)
        else:
            SeedIndex = -999
        self.Draw.text((10, BHeight-7), "b", Colors.BSeries, font = TheFont)
        # First, draw tickmarks for the b peaks
        for MassIndex in range(len(Peptide.Masses)):
            Mass = Peptide.Masses[MassIndex]
            X = self.GetX(Mass)
            Label = ("b %s"%MassIndex)
            if self.Labels.has_key(Label):
                PeakIntensity = self.Labels[Label][-2]
                Y = self.GetY(PeakIntensity)
                self.DrawDottedLine(X, BHeight-2, X, Y, Colors.BSeriesPale)
            else:
                self.Draw.line((X, BHeight-2, X, BHeight+2), Colors.BSeriesPale)
        # Now draw horizontal lines, and amino labels:
        for AminoIndex in range(len(Peptide.Aminos)):
            LabelA = "b %s"%AminoIndex
            LabelB = "b %s"%(AminoIndex+1)
            XA = self.GetX(Peptide.Masses[AminoIndex])
            XB = self.GetX(Peptide.Masses[AminoIndex+1])
            HasLabelA = self.Labels.has_key(LabelA)
            if AminoIndex == 0:
                HasLabelA = 1
            HasLabelB = self.Labels.has_key(LabelB)
            if AminoIndex ==len(Peptide.Aminos)-1:
                HasLabelB = 1
            if HasLabelA and HasLabelB:
                self.Draw.line((XA, BHeight, XB, BHeight), Colors.BSeries)
            else:
                self.DrawDottedLine(XA, BHeight, XB, BHeight, Colors.BSeriesPale)
            if AminoIndex in (SeedIndex, SeedIndex+1, SeedIndex+2):
                self.Draw.line((XA, BHeight-1, XB, BHeight-1), Colors.BSeries)
                self.Draw.line((XA, BHeight, XB, BHeight), Colors.BSeries)
                self.Draw.line((XA, BHeight+1, XB, BHeight+1), Colors.BSeries)
                
            X = (XA+XB)/2 - 3
            Str = Peptide.Aminos[AminoIndex]
            if Peptide.Modifications.get(AminoIndex):
                self.Draw.text((X-4, BHeight-14), "%s*"%Str, Colors.BSeries, font = TheFont)
            else:
                self.Draw.text((X, BHeight-14), Str, Colors.BSeries, font = TheFont)
    def DrawYSeries(self, Peptide): # Copied and modded from DrawBSeries
        YHeight = 34
        if getattr(Peptide, "Seed", None):
            SeedIndex = Peptide.Aminos.rfind(Peptide.Seed)
        else:
            SeedIndex = -999
        
        self.Draw.text((10, YHeight-7), "y", Colors.YSeries, font = TheFont)
        # First, draw tickmarks for the y peaks
        PM = Peptide.Masses[-1] + 19 # parent mass
        for MassIndex in range(len(Peptide.Masses)):
            Mass = PM - Peptide.Masses[MassIndex]
            X = self.GetX(Mass)
            Label = "y %s"%(len(Peptide.Masses) - MassIndex - 1)
            #Label = "y %s"%(MassIndex)
            #print "Y series y %d, mass %f, label %s"%(MassIndex, Mass, Label)
            if self.Labels.has_key(Label):
                PeakIntensity = self.Labels[Label][-2]
                Y = self.GetY(PeakIntensity)
                self.DrawDottedLine(X, YHeight-2, X, Y, Colors.YSeriesPale)
            else:
                self.Draw.line((X, YHeight-2, X, YHeight+2), Colors.YSeriesPale)
        # Now draw horizontal lines, and amino labels:
        for AminoIndex in range(len(Peptide.Aminos)):
            LabelA = "y %s"%(len(Peptide.Aminos) - AminoIndex)
            LabelB = "y %s"%(len(Peptide.Aminos) - AminoIndex - 1)
            XA = self.GetX(PM - Peptide.Masses[AminoIndex])
            XB = self.GetX(PM - Peptide.Masses[AminoIndex+1])
            HasLabelA = self.Labels.has_key(LabelA)
            if AminoIndex == 0:
                HasLabelA = 1
            HasLabelB = self.Labels.has_key(LabelB)
            if AminoIndex ==len(Peptide.Aminos)-1:
                HasLabelB = 1
            if HasLabelA and HasLabelB:
                self.Draw.line((XA, YHeight, XB, YHeight), Colors.YSeries)
            else:
                self.DrawDottedLine(XA, YHeight, XB, YHeight, Colors.YSeriesPale)
            if AminoIndex in (SeedIndex, SeedIndex+1, SeedIndex+2):
                self.Draw.line((XA, YHeight-1, XB, YHeight-1), Colors.YSeries)
                self.Draw.line((XA, YHeight, XB, YHeight), Colors.YSeries)
                self.Draw.line((XA, YHeight+1, XB, YHeight+1), Colors.YSeries)
                
            X = (XA+XB)/2 - 3
            Str = Peptide.Aminos[AminoIndex]
            if Peptide.Modifications.get(AminoIndex):
                self.Draw.text((X-4, YHeight-14), "%s*"%Str, Colors.YSeries, font = TheFont)
            else:
                self.Draw.text((X, YHeight-14), Str, Colors.YSeries, font = TheFont)

UsageInfo = """
Usage:
   MakeImage.py <LabeledSpectrum> [<OutputFileName>]
"""

if __name__ == "__main__":
    if len(sys.argv)<2:
        print UsageInfo
        sys.exit(1)
    InputFileName = sys.argv[1]
    if len(sys.argv)>2:
        OutputFileName = sys.argv[2]
    else:
        OutputFileName = os.path.splitext(InputFileName)[0] + ".png"
    Maker = MSImageMaker()
    Maker.ConvertSpectrumFileToImage(InputFileName, OutputFileName)
