"""
Generic annotation-track plotter
"""
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
try:
    TheFont = ImageFont.truetype("Arial.ttf", 12)
except:
    TheFont = ImageFont.load_default()

class FeatureTypes:
    Plain = "plain"
    Dashed = "dashed"
    Fat = "fat"
    FatDashed = "fatdashed"

class Colors:
    "Bag of constants specifying a color scheme"
    Background = (255, 255, 255)
    Feature = (155, 155, 155)
    Exon = (155, 155, 155)
    SNP = (55, 55, 55)
    Mismatch = (255, 100, 100)
    Peptide = (200, 55, 55)
    Peptide2 = (100, 0, 0)
    PeptideSpliceGhost = (220, 155, 155)
    ESTWrong = (100, 255, 100)
    LinkCorrect = (200, 0, 0)
    LinkWrong = (125, 125, 255)
    GF = (0, 0, 100)
    Label = (0, 0, 0)
    Grey = (155, 155, 155)
    Black = (0, 0, 0)
    Border = (200, 200, 200)
    SpanningEdge = (0, 155, 155)

class AnnotationFeature:
    def __init__(self, Start, End):
        self.Start = Start
        self.End = End
        self.Color = None
        self.Label = None
        self.Blocks = []
        self.Type = FeatureTypes.Plain
    def __cmp__(self, OtherFeature):
        if self.Start < OtherFeature.Start:
            return -1
        if self.Start > OtherFeature.Start:
            return 1
        if self.End < OtherFeature.End:
            return -1
        if self.End > OtherFeature.End:
            return 1
        if id(self) < id(OtherFeature):
            return -1
        if id(self) > id(OtherFeature):
            return 1
        return 0
    
class AnnotationTrack:
    def __init__(self, Name):
        self.Name = Name
        self.Features = []
        self.Y = 100
        self.LabelLineFlag = 0
        self.Color = Colors.Feature
        self.Start = None
        self.End = None
        self.RequiredFlag = 0
    def AddFeature(self, Start, End):
        Feature = AnnotationFeature(Start, End)
        self.Features.append(Feature)
        if (self.Start == None or self.Start > Start):
            self.Start = Start
        if (self.End == None or self.End < End):
            self.End = End
        return Feature
    def SortFeatures(self):
        self.Features.sort()
    def GetStart(self):
        return self.Features[0].Start
    def GetEnd(self):
        End = 0
        for Feature in self.Features:
            End = max(Feature.End, End)
        return End
##    def MergeFeatures(self):
##        """
##        Examine our feature-list.  If any two features overlap, MERGE them into a single large feature.
##        """
##        pass
class PlotBlockClass:
    def __init__(self):
        self.Start = None
        self.End = None
        self.Left = None
        self.Right = None
        
class PlotterClass:
    def __init__(self):
        self.TitleLines = []
        self.Tracks = []
        self.Blocks = []
        self.BlockEdgePadding = 20
        self.GapSizeForNewBlock = 100 #%%% temp
        self.PixelsPerPos = 1.0 # 0.1
        self.LabelPositionsFlag = 1 # by default, DO indicate residue numbers along chromosome
    def AddHeaderLine(self, Str):
        self.TitleLines.append(Str)
    def AddTrack(self, Name):
        Track = AnnotationTrack(Name)
        self.Tracks.append(Track)
        return Track
    def CreateBlocks(self):
        """
        A -block- is a contiguous region of the genome.  We start by creating a block
        for the first feature, and we continue adding features to this block until
        the next feature comes too far after the end of the current block.  
        """
        RequiredFeatures = []
        AllFeatures = [] 
        for Track in self.Tracks:
            if Track.RequiredFlag:
                for Feature in Track.Features:
                    if Feature.Type not in (FeatureTypes.Dashed, FeatureTypes.FatDashed):
                        RequiredFeatures.append(Feature)
            AllFeatures.extend(Track.Features)
        if len(RequiredFeatures) == 0:
            # No required features?  Assume all features are required, then:
            RequiredFeatures = AllFeatures
        if not AllFeatures:
            print "* Warning: No features to plot!"
        AllFeatures.sort()
        RequiredFeatures.sort()
        Block = PlotBlockClass()
        Block.Start = RequiredFeatures[0].Start - self.BlockEdgePadding
        Block.End = RequiredFeatures[0].End + self.BlockEdgePadding
        Block.Left = 100
        self.Blocks = [Block]
        print "Create blocks spanning %s required features"%len(RequiredFeatures)
        for FeatureIndex in range(len(RequiredFeatures)):
            Feature = RequiredFeatures[FeatureIndex]
            print "req feature %s: %s-%s %s %s"%(FeatureIndex, Feature.Start, Feature.End, Feature.Color, Feature.Type)
            if Feature.Start < (Block.End + self.GapSizeForNewBlock):
                Feature.Blocks = [Block]
                Block.End = max(Block.End, Feature.End + self.BlockEdgePadding)
            else:
                # End the old block:
                BlockWidth = Block.End - Block.Start
                Block.Right = Block.Left + BlockWidth * self.PixelsPerPos
                # Start a new block!
                NewBlock = PlotBlockClass()
                NewBlock.Left = Block.Right + 1
                NewBlock.Start = Feature.Start - self.BlockEdgePadding
                NewBlock.End = Feature.End + self.BlockEdgePadding
                self.Blocks.append(NewBlock)
                Block = NewBlock
                Feature.Blocks = [Block]
        # Set the right-edge of the last block:
        BlockWidth = Block.End - Block.Start
        Block.Right = Block.Left + BlockWidth * self.PixelsPerPos
        print "Blocks (%s):"%len(self.Blocks)
        for Block in self.Blocks:
            print " ", Block.Start, Block.End
        # Assign all non-required features to blocks:
        for Feature in AllFeatures:
            for Block in self.Blocks:
                if Feature.End > Block.Start and Feature.Start < Block.End:
                    if Block not in Feature.Blocks:
                        Feature.Blocks.append(Block)
                
    def DrawTitleLines(self):
        Y = 2
        for TitleLine in self.TitleLines:
            self.Draw.text((5, Y), TitleLine, Colors.Label)
            Y += 10
    def DrawBlockTimeline(self):
        for BlockIndex in range(len(self.Blocks)):
            Block = self.Blocks[BlockIndex]
            self.Draw.line((Block.Left + 2, self.GenomeY, Block.Right - 2, self.GenomeY), Colors.Label)
            # Sigil: ----//----
            self.Draw.line((Block.Left, self.GenomeY + 2, Block.Left + 4, self.GenomeY - 2), Colors.Label)
            self.Draw.line((Block.Right - 4, self.GenomeY + 2, Block.Right, self.GenomeY - 2), Colors.Label)
            if self.LabelPositionsFlag:
                # Draw genome, and label the genome positions
                LabelGenomePos = ((Block.Start + 51) / 100) * 100
                LabelY = self.GenomeY - 15 - (BlockIndex % 2) * 10
                while LabelGenomePos < Block.End:
                    X = Block.Left + (LabelGenomePos - Block.Start) 
                    # Don't draw the label if we're close to the edge of the exon group.
                    # (We don't want the label to spill over *too* far into the adjoining group)
                    if (X > Block.Left + 5 and X < Block.Right - 5):
                        self.Draw.line((X, self.GenomeY, X, LabelY + 10), Colors.Border)
                        Str = str(LabelGenomePos)
                        self.Draw.text((X - 3*len(Str), LabelY), Str, Colors.Label)
                    LabelGenomePos += 100
            # Vertical line between blocks:
            #if (BlockIndex < len(self.Blocks) - 1):
            #    self.Draw.line((Block.Right, 0, Block.Right, self.Height), Colors.Border)
    def DrawTrack(self, Track):
        if Track.Name:
            self.Draw.text((1, Track.Y - 5), Track.Name, Track.Color)
        for Feature in Track.Features:
            for Block in Feature.Blocks:
                # Plot the feature inside this block, even if the feature extends farther:
                X1 = Block.Left + max(0, (Feature.Start - Block.Start) * self.PixelsPerPos)
                X2 = Block.Right - max(0, (Block.End - Feature.End) * self.PixelsPerPos)
                Color = Feature.Color
                if Color == None:
                    Color = Track.Color
                # DRAW the feature:
                if Feature.Type == FeatureTypes.Plain:
                    # Be sure that the feature is visible:
                    X2 = max(X2, X1 + 1)
                    self.Draw.line((X1, Track.Y, X2, Track.Y), Color)
                    self.Draw.line((X1, Track.Y + 1, X2, Track.Y + 1), Color)
                    self.Draw.line((X1, Track.Y + 2, X2, Track.Y + 2), Color)
                elif Feature.Type == FeatureTypes.Dashed:
                    X = X1 + 2
                    while (X < X2):
##                        self.Draw.line((X, Track.Y, X + 2, Track.Y), Color)
##                        self.Draw.line((X, Track.Y + 1, X + 2, Track.Y + 1), Color)
##                        X += 4
                        self.Draw.line((X, Track.Y, X + 1, Track.Y), Color)
                        self.Draw.line((X, Track.Y + 1, X + 1, Track.Y + 1), Color)
                        self.Draw.line((X, Track.Y + 2, X + 1, Track.Y + 2), Color)
                        X += 2
                elif Feature.Type == FeatureTypes.Fat:
                    for YY in range(-4, 5):
                        self.Draw.line((X1, Track.Y + YY, X2, Track.Y + YY), Color)
                elif Feature.Type == FeatureTypes.FatDashed:
                    X = X1 + 2
                    while (X < X2):
                        #for YY in range(-4, 5):
                        self.Draw.line((X, Track.Y - 4, X, Track.Y + 5), Color)
                        X += 2
                        
                # LABEL the feature:
                if Feature.Label:
                    if Track.LabelLineFlag:
                        LabelPos = (X1 - 2*len(Feature.Label), Track.Y + 3)
                        #LabelPos = (X1 - 30, Track.Y + 3)
                    else:
                        LabelPos = (X2 + 3, Track.Y - 5)
                    self.Draw.text(LabelPos, Feature.Label, Color)
    def Plot(self, OutputFileName):
        """
        Plot a region of the genome, using our tracks.
        """
        print
        print "="*78
        for Line in self.TitleLines:
            print Line
        self.Start = None
        self.End = None
        for Track in self.Tracks:
            if (self.Start == None or self.Start > Track.Start):
                self.Start = Track.Start
            if (self.End == None or self.End > Track.End):
                self.End = Track.End
        self.CreateBlocks()
        # Prune any tracks whose features "slipped through the cracks" between blocks:
        for Track in self.Tracks[:]:
            ValidFeatureCount = 0
            for Feature in Track.Features:
                if Feature.Blocks:
                    ValidFeatureCount += 1
            if not ValidFeatureCount and not Track.RequiredFlag:
                self.Tracks.remove(Track)
        # Decide upon vertical layout and image dimensions:
        Height = 2
        # Some vertical space for the title lines:
        for Line in self.TitleLines:
            Height += 10
        Height += 25 # position labels
        self.GenomeY = Height
        Height += 20 # Spacing between genome 'timeline' and first track
        for Track in self.Tracks:
            Track.Y = Height
            if Track.LabelLineFlag:
                Height += 15
            else:
                Height += 10
        self.Height = int(round(Height))
        self.Width = self.Blocks[-1].Right + 250
        self.Width = max(self.Width, 500)
        self.Width = int(round(self.Width))
        print "->-> Create image %s x %s"%(self.Width, self.Height)
        self.Image = Image.new("RGB", (self.Width, self.Height), Colors.Background)
        
        self.Draw = ImageDraw.Draw(self.Image)
        self.DrawTitleLines()
        self.DrawBlockTimeline()
        for Track in self.Tracks:
            self.DrawTrack(Track)
        self.Image.save(OutputFileName, "png")