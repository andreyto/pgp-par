import sys
import os
import getopt
import PyInspect
import MakeImage
from Utils import *
import MSSpectrum
import GetByteOffset
import PLSUtils
UsageInfo = """
Label.py - Generate a labeled spectrum, given a peptide.

Required Options
 -r [FileName] Spectrum file
 -b [Offset] The byte offset in the file for the spectrum, as reported in
     the Inspect output, it can be left blank for single-spectrum
 -a [Peptide] The annotation for the spectrum
 -c [Charge] The charge of the peptide

Additional Options
 -w [FileName] Output file name.  Default is to temp.png
 -v [FileName] Write verbose scoring details to the specified file
 -d [Width]: Image width
 -h [Height]: Image height
 -s [ScanNumber]: Scan number
 -x: Use black and white (for Printing)
Example:
    Label.py -r Sample346.mzXML -b 38289818 -a R.A+226LLAAFDFPFR.K
"""

class LabelClass:
    def __init__(self):
        self.SpectrumPath = None
        self.SpectrumFilePos = 0
        self.Peptide = None
        self.OutputFileName = "temp.png"
        self.VerboseFileName = None
        self.LabeledPeaks = None
        self.InspectFeatures = None
        self.InspectFeatureNames = ["MQScore", "Length", "Total Cut Score", "Median Cut Score", "Y present", "B present", "Intensity in BY", "NTT"]
        self.AutoPopUp = 1
        self.PeptideHasPhosphorylation = 0
        self.InstrumentType = "ESI-ION-TRAP" # or QTOF or FT-HYBRID
        self.Charge = 0 #guessed or set by user.
        self.ImageWidth = 600
        self.ImageHeight = 400
        self.ScanNumber = None
        self.DoPLS = 0 #don't do this by default
    def ParseCommandLineSimple(self, Arguments):
        self.SpectrumPath = Arguments[0]
        ColonPos = self.SpectrumPath.rfind(":")
        try:
            self.SpectrumFilePos = int(self.SpectrumPath[ColonPos + 1:])
            self.SpectrumPath = self.SpectrumPath[:ColonPos]
        except:
            self.SpectrumFilePos = 0
        self.Peptide = GetPeptideFromModdedName(Arguments[1])
        if Arguments[1].find("phos") > 0:
            self.PeptideHasPhosphorylation = 1
        if len(Arguments) > 2:
            self.OutputFileName = Arguments[2]
    def ParseCommandLine(self, Arguments):
        # Hack:
        if len(Arguments) > 0 and Arguments[0][1] == ":":
            return self.ParseCommandLineSimple(Arguments)
        (Options, Args) = getopt.getopt(Arguments, "r:b:s:a:w:v:pi:c:d:h:s:xP")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-a":
                Annotation = Value
                self.Peptide = GetPeptideFromModdedName(Annotation)
                if Annotation.find("phos") > 0:
                    self.PeptideHasPhosphorylation = 1
            elif Option == "-r":
                self.SpectrumPath = Value
            elif Option == "-b":
                self.SpectrumFilePos = int(Value)
            elif Option == "-w":
                self.OutputFileName = Value
            elif Option == "-v":
                self.VerboseFileName = Value
            elif Option == "-i":
                self.InstrumentType = Value
            elif Option == "-c":
                self.Charge = int(Value)
            elif Option == "-p":
                self.AutoPopUp = 0
                #secret option to supress the image from popping up to the screen
            elif Option == "-d":
                self.ImageWidth = int(Value)
            elif Option == "-h":
                self.ImageHeight = int(Value)
            elif Option == "-s":
                self.ScanNumber = int(Value)
            elif Option == "-x":
                MakeImage.SetColors(1)
            elif Option == "-P":
                self.DoPLS = 1
            else:
                raise ValueError, "* Unknown option %s"%Option
        # Filename and annotation are required.  (Byte position is optional,
        # since there are many single-scan .dta files out there)
        if not OptionsSeen.has_key("-a") or not OptionsSeen.has_key("-r"):
            print UsageInfo
            sys.exit(1)
    def Main(self):
        if self.ScanNumber != None:     # scan number is provided in input
            # get byte offset using scan number
            Abacus = GetByteOffset.Abacus()
            self.ScanOffset = Abacus.GetByteOffset(self.SpectrumPath)
            self.SpectrumFilePos = self.ScanOffset[self.ScanNumber]
            print "ByteOffset # = %s"%self.SpectrumFilePos
        self.LabelPeaks()
        #self.ConvertDoublyChargedPeakLabelsOLD()
        self.LabeledPeaks = self.ConvertDoublyChargedPeakLabels(self.LabeledPeaks, self.Peptide)
        #self.ConvertYPeakNumberingOLD()
        self.LabeledPeaks = self.ConvertYPeakNumbering(self.LabeledPeaks)
        self.ConvertParentLossLabels()
        self.VerboseOutput()
        Maker = MakeImage.MSImageMaker(Width = self.ImageWidth, Height = self.ImageHeight)
        Maker.ConvertPeakAnnotationToImage(self.LabeledPeaks, self.OutputFileName, self.Peptide,
            Width = self.ImageWidth, Height = self.ImageHeight)
        if self.AutoPopUp:
            os.startfile(self.OutputFileName)
    def VerboseOutput(self):
        """
        Extra output for the curious
        1. Inspect scoring features
        """
        if not self.VerboseFileName:
            return
        ##1. Inspect scoring features
        VerboseHandle = open(self.VerboseFileName, "wb")
        VerboseHandle.write("M/z %f\n"%self.MZ)
        VerboseHandle.write("Annotation %s\n"%self.Peptide.GetFullModdedName())
        VerboseHandle.write("ParentMass: Hypothetical, Observered, Error: %.2f, %.2f, %.2f\n"%(self.HypotheticalParentMass, self.ObservedParentMass, self.ObservedParentMassError))
        for Index in range(len(self.InspectFeatures)):
            VerboseHandle.write("%s\t%.3f\n"%(self.InspectFeatureNames[Index],self.InspectFeatures[Index]))
        ##If they want the Phosphate localization Score do it here.
        if self.DoPLS:
            PLS = self.CalculatePLS()
            if PLS:
                #here it is possible that we get a different peptide as winner
                VerboseHandle.write("Phosphate Localization Score: %.3f\n"%PLS[0])
                if len(PLS) >  1:
                    VerboseHandle.write("WARNING: Better annotation than input. %.4f, %s"%(PLS[1], PLS[2]))
            else:
                VerboseHandle.write("Phosphate Localization Score: N/A\n")

        ## 2. Write out the peaks found, not found
        String = self.GetFoundPeaksTable()
        VerboseHandle.write("\n\nPeaksFoundTable\n%s\n\n"%String)

        ## 3. Dump the peak list
        VerboseHandle.write("Mass\tIntensity\tLabel\tAminoIndex\n")
        for Tuple in self.LabeledPeaks:
            Label = Tuple[2] 
            if not Label:
                Label = "UnLabeled"
            Str = "%f\t%f\t%s\t%d\n"%(Tuple[0], Tuple[1], Label, Tuple[3])
            VerboseHandle.write(Str)
        VerboseHandle.write("\n\n")
            
        VerboseHandle.close()

    def GetFoundPeaksTable(self):
        ## Get masses from Peptide object
        ## mark the ones found with bold something
        IonMasses = {} #key = IonName, value\tvalue\t
        IonsFound = {} #key = IonNameAndIndex, value = 1 if found
        IonNamesSorted = ("b2", "b", "y", "y2")
        ##Mark Ions found
        for Tuple in self.LabeledPeaks:
            Label = Tuple[2]
            if not Label:
                continue
            Label = Label.lower()
            if Label in IonNamesSorted:
                Index = Tuple[3]
                Key = "%s:%s"%(Label,Index)
                IonsFound[Key] = 1
                #print "%s\t%s"%(Key, Tuple)
        ##get predicted values
        ReturnString = ""
        for IonName in IonNamesSorted:
            IonMasses[IonName] = "%s\t"%IonName
            Ion = Global.AllIonDict[IonName]
            for Index in range(1, len(self.Peptide.Masses)-1):
                Mass = self.Peptide.Masses[Index] # offset by 1, since mass 0 is in there.
                MassForIonType = Ion.GetPeakMass(Mass,self.Peptide.GetParentMass())
                if IonName[0] == "b":
                    Key = "%s:%s"%(IonName,Index)
                else:
                    NewIndex = len(self.Peptide.Aminos) - Index
                    Key = "%s:%s"%(IonName,NewIndex)
                if IonsFound.has_key(Key):
                    #print "I found %s at mass %s"%(Key,MassForIonType)
                    IonMasses[IonName]+="F%.3f\t"%MassForIonType
                else:
                    IonMasses[IonName]+="%.3f\t"%MassForIonType
                #print "%s, %s"%(Key, MassForIonType)
            ## Make String now that we are done with this IonName
            ReturnString += "%s\n"%IonMasses[IonName]
            if IonName == "b":
                ReturnString += " \t%s\n"%self.Peptide.GetModdedName()
        return ReturnString

    def CalculatePLS(self):
        """This function calculates the Phosphate Localization Score which is the ambuiguity score
        for phosphorylation placement to a specific amino acid in the sequence.  This
        is reported in Albuquerque et al.  Mol Cell Prot 2008
        """
        ##1. Get all the potential annotations for the peptide
        ##      If there are none, then the score is "N/A"
        ##2. Get the score of each alternate annotation
        ##3. determine winner and runner up
        ##4. calcualte PLS
        Abacus = PLSUtils.PLSClass()
        PotentialAnnotations = Abacus.GetAlternateAnnotations(self.Peptide)
        
        if len(PotentialAnnotations) == 0:
            return None
        ## 2. Try each individual annotation, keeping track of the top and runner up
        BestAlternateMQScore = -10
        RunnerUpAlternateMQScore = -10
        BestAlternatePeakList = None
        RunnerUpAlternatePeakList = None
        BestAlternatePeptide = None
        RunnerUpAlternatePeptide = None
        for Annotation in PotentialAnnotations:
            NewPeptide = GetPeptideFromModdedName(Annotation) # needed for peak label conversion
            NewPeptide.Prefix = self.Peptide.Prefix
            NewPeptide.Suffix = self.Peptide.Suffix
            ##have to load it each time, it was not getting the correct results if I reused the same object
            PySpectrum = PyInspect.Spectrum(self.SpectrumPath, self.SpectrumFilePos)
            PySpectrum.SetParentMass(self.HypotheticalParentMass, self.Charge)
            Features = PySpectrum.ScorePeptideDetailed(Annotation, self.Charge)
            PeakAnnotations = PySpectrum.LabelPeaks(Annotation, self.Charge)
            MQScore = Features[0]
            #print "The score is %s, %s"%(MQScore, Annotation)
            if MQScore > BestAlternateMQScore:
                ##swap with runner up
                RunnerUpAlternateMQScore = BestAlternateMQScore
                RunnerUpAlternatePeakList = BestAlternatePeakList
                RunnerUpAlternatePeptide = BestAlternatePeptide
                BestAlternateMQScore = MQScore
                BestAlternatePeakList = PeakAnnotations
                BestAlternatePeptide = NewPeptide
            elif MQScore > RunnerUpAlternateMQScore:
                RunnerUpAlternateMQScore = MQScore
                RunnerUpAlternatePeakList = PeakAnnotations
                RunnerUpAlternatePeptide = NewPeptide
        ## 3. Determine ther real winner and runner up.  we assume that the original inspect
        ## annotation is right, unless something beats it by say 0.3 units HARD CODED MAGIC!!!!!!!!!!!!!!!!
        ## this is hard coded magic, but I tested it and it performs better than range(0,1,0.1)
        ## then we swap out the top annotation.
        if BestAlternateMQScore > (self.MQScore + 0.3):
            TopPeakList = BestAlternatePeakList
            TopPeptide = BestAlternatePeptide
            TopMQScore = BestAlternateMQScore
            #now also consider the fate of the RU score
            if RunnerUpAlternateMQScore > (self.MQScore + 0.3):
                RunnerUpPeakList = RunnerUpAlternatePeakList
                RunnerUpPeptide = RunnerUpAlternatePeptide
                RunnerUpMQScore = RunnerUpAlternateMQScore
            else:
                RunnerUpPeakList = self.LabeledPeaks
                RunnerUpPeptide = self.Peptide
                RunnerUpMQScore = self.MQScore
        else:
            TopPeakList = self.LabeledPeaks
            TopPeptide = self.Peptide
            TopMQScore = self.MQScore
            RunnerUpPeakList = BestAlternatePeakList
            RunnerUpPeptide = BestAlternatePeptide
            RunnerUpMQScore = BestAlternateMQScore
        ## 4. Find the distinguishing peaks between the top 2 
        #print "Winner ", TopPeptide.GetFullModdedName()
        #print "runner up", RunnerUpPeptide.GetFullModdedName()
        TopPeakList = self.ConvertDoublyChargedPeakLabels(TopPeakList, TopPeptide)
        TopPeakList = self.ConvertYPeakNumbering(TopPeakList)
        RunnerUpPeakList = self.ConvertDoublyChargedPeakLabels(RunnerUpPeakList, RunnerUpPeptide)
        RunnerUpPeakList = self.ConvertYPeakNumbering(RunnerUpPeakList)
        DistinguishingPeakList = Abacus.GetDistinguishingPeaks(TopPeptide, RunnerUpPeptide)
        #print "finding peaks for %s"%TopPeptide.GetModdedName()
        nWinner = Abacus.GetSupportingPeaks(TopPeakList, DistinguishingPeakList)
        #print "finding peaks for %s"%RunnerUpPeptide.GetModdedName()
        nRunnerUp = Abacus.GetSupportingPeaks(RunnerUpPeakList, DistinguishingPeakList)
        ## 4.5 Here we take a slight detour.  If nWinner < nRunnerUp, then PLS predicts something
        ## different from Inspect.  This happens, scoring functions will have different opinions
        ## We simply swap the two and call it a day.
        #print "Getting the ambuiguity score with %s, %s, %s (top, ru, total)"%(nWinner, nRunnerUp, len(DistinguishingPeakList))
        AmbuigityScore = Abacus.ComputePLS(len(DistinguishingPeakList), nWinner, nRunnerUp)
        if AmbuigityScore < 0:
            ## means that nWinner < nRunnerUp
            AmbuigityScore *= -1
            ##now we shamelessly dump the top guy
            TopMQScore = RunnerUpMQScore
            TopPeptide = RunnerUpPeptide
        #print "Ascore is %s"%AmbuigityScore
        if not TopMQScore == self.MQScore:
            print "WARNING::Top score was %.2f for peptide %s"%(TopMQScore, TopPeptide.GetModdedName())
            print "\tInput was %.2f and %s"%(self.MQScore, self.Peptide.GetModdedName())
            return (AmbuigityScore, TopMQScore, TopPeptide.GetFullModdedName())
        return (AmbuigityScore,)

            
    def ConvertParentLossLabels(self):
        """
        Special case for phorphorylated spectra. Change the label
        'Parent loss' to M-p or M-p-h2o.
        """
        PhosLoss = 98.0 / self.Charge
        PhosWaterLoss = 116.0 / self.Charge
        PhosLabel = "M-p"
        PhosWaterLabel = "M-p-h2o"
        Error = 3.0
        for Index in range(len(self.LabeledPeaks)):
            Tuple = self.LabeledPeaks[Index]
            Label = Tuple[2]
            if not Label == "Parent loss":
                continue
            Mass = Tuple[0]
            Diff = abs(Mass - self.MZ)
            MaybePhosLoss = abs(Diff - PhosLoss)
            #print Tuple
            if MaybePhosLoss < Error:
                NewTuple = (Tuple[0], Tuple[1], PhosLabel, Tuple[3])
                self.LabeledPeaks[Index] = NewTuple
                #print self.LabeledPeaks[Index]
                continue
            MaybePhosWaterLoss = abs(Diff - PhosWaterLoss)
            if MaybePhosWaterLoss < Error:
                NewTuple = (Tuple[0], Tuple[1], PhosWaterLabel, Tuple[3])
                self.LabeledPeaks[Index] = NewTuple
                #print self.LabeledPeaks[Index]
    def ConvertDoublyChargedPeakLabelsOLD(self):
        """
        The inspect output does not distinguish between single, and doubly charged peaks.
        so in order for labeling to go well, we have to rewrite the labels as B2, Y2, etc
        """
        for Index in range(len(self.LabeledPeaks)):
            Tuple = self.LabeledPeaks[Index]
            Label = Tuple[2]
            TupleMass = Tuple[0]
            AminoIndex = Tuple[3]
            NewLabel = None
            PeptideMass = self.Peptide.Masses[AminoIndex]
            if Label == "B":
                PeptideMass += 1.0
                if abs(PeptideMass - TupleMass) > 5:
                    NewLabel = "B2"
            if Label == "Y":
                PeptideMass = self.Peptide.GetParentMass() - PeptideMass
                if abs(PeptideMass - TupleMass) > 5:
                    #a doublycharged peak, no isotope or error is this big
                    NewLabel = "Y2"
            if Label == "Y loss":  #hacky, but I can't think of a good way
                YPeptideMass = self.Peptide.GetParentMass() - PeptideMass
                Found = 0
                for CommonLoss in [17, 18, 98]:
                    YLossMass = YPeptideMass - CommonLoss
                    if abs(YLossMass - TupleMass) < 5:
                        Found = 1
                        break
                if not Found:
                    NewLabel = "Y2 Loss"
            if Label == "B loss":  #hacky, but I can't think of a good way
                BPeptideMass = PeptideMass + 1
                Found = 0
                for CommonLoss in [17, 18, 98]:
                    BLossMass = BPeptideMass - CommonLoss
                    if abs(BLossMass - TupleMass) < 5:
                        Found = 1
                        break
                if not Found:
                    NewLabel = "B2 Loss"

            if NewLabel:
                NewTuple = (Tuple[0], Tuple[1], NewLabel, Tuple[3])
                self.LabeledPeaks[Index] = NewTuple
    def ConvertDoublyChargedPeakLabels(self, Peaks, Peptide):
        """
        The inspect output does not distinguish between single, and doubly charged peaks.
        so in order for labeling to go well, we have to rewrite the labels as B2, Y2, etc
        """
        for Index in range(len(Peaks)):
            Tuple = Peaks[Index]
            Label = Tuple[2]
            TupleMass = Tuple[0]
            AminoIndex = Tuple[3]
            if abs(Tuple[0] - 402) < 1:
                Verbose = 1
            else:
                Verbose = 0
            NewLabel = None
            PeptideMass = Peptide.Masses[AminoIndex]
            if Label == "B":
                PeptideMass += 1.0
                if abs(PeptideMass - TupleMass) > 5:
                    NewLabel = "B2"
            if Label == "Y":
                PeptideMass = Peptide.GetParentMass() - PeptideMass
                if abs(PeptideMass - TupleMass) > 5:
                    #a doublycharged peak, no isotope or error is this big
                    NewLabel = "Y2"
            if Label == "Y loss":  #hacky, but I can't think of a good way
                YPeptideMass = Peptide.GetParentMass() - PeptideMass
                Found = 0
                for CommonLoss in [17, 18, 98]:
                    YLossMass = YPeptideMass - CommonLoss
                    if abs(YLossMass - TupleMass) < 5:
                        Found = 1
                        break
                if not Found:
                    NewLabel = "Y2 Loss"
            if Label == "B loss":  #hacky, but I can't think of a good way
                BPeptideMass = PeptideMass + 1
                Found = 0
                for CommonLoss in [17, 18, 98]:
                    BLossMass = BPeptideMass - CommonLoss
                    if abs(BLossMass - TupleMass) < 5:
                        Found = 1
                        break
                if not Found:
                    NewLabel = "B2 Loss"

            if NewLabel:
                NewTuple = (Tuple[0], Tuple[1], NewLabel, Tuple[3])
                Peaks[Index] = NewTuple
        return Peaks

    def ConvertYPeakNumberingOLD(self):
        """
        The amino indicies are numbered from the N- to C-terminus, but MakeImage numbers
        its Y peaks from y1 (nearest the C-terminus) upwards.  We re-number them here.
        """
        TempList = self.LabeledPeaks
        self.LabeledPeaks = [] #clean it out
        for Tuple in TempList:
            Label = Tuple[2]
            if not Label: #not a labeled peak.  proceede normally
                self.LabeledPeaks.append(Tuple)        
                continue
            if not Label[0] == "Y":
                self.LabeledPeaks.append(Tuple)
                continue
            ## should only have y derivates here.  switch indices
            AminoIndex = Tuple[-1]
            #print Tuple
            NewIndex = len(self.Peptide.Aminos) - AminoIndex
            NewTuple = (Tuple[0], Tuple[1], Tuple[2], NewIndex)
            self.LabeledPeaks.append(NewTuple)
    def ConvertYPeakNumbering(self, Peaks):
        """ SAME, just takes a parameter.  I know it's messy.
        The amino indicies are numbered from the N- to C-terminus, but MakeImage numbers
        its Y peaks from y1 (nearest the C-terminus) upwards.  We re-number them here.
        """
        TempList = Peaks
        Peaks = [] #clean it out
        for Tuple in TempList:
            Label = Tuple[2]
            if not Label: #not a labeled peak.  proceede normally
                Peaks.append(Tuple)        
                continue
            if not Label[0] == "Y":
                Peaks.append(Tuple)
                continue
            ## should only have y derivates here.  switch indices
            AminoIndex = Tuple[-1]
            #print Tuple
            NewIndex = len(self.Peptide.Aminos) - AminoIndex
            NewTuple = (Tuple[0], Tuple[1], Tuple[2], NewIndex)
            Peaks.append(NewTuple)
        return Peaks
            
    def LabelPeaks(self):
        """
        Uses PyInspect to label peaks in the spectrum according to Inspect's scoring
        PyInspect will always be current, so let's use it.
        """
        ## load a spectrum, set charge, parent mass, then label the peaks
        PySpectrum = PyInspect.Spectrum(self.SpectrumPath, self.SpectrumFilePos)
        self.MZ = PySpectrum.GetMZ()
        #print "m/z is %f"%self.MZ
        ParentMass = self.Peptide.GetParentMass()
        if not self.Charge:  ## Guess charge if not input
            BestDiff = 99999
            for Charge in range(1, 5):
                ParentMassFromCharge = self.MZ * Charge - (Charge - 1)*1.0078
                Diff = abs(ParentMass - ParentMassFromCharge)
                if Diff < BestDiff:
                    BestDiff = Diff
                    BestCharge = Charge
                    BestMass = ParentMassFromCharge
            self.Charge = BestCharge
            print "Appears to be charge %d with mass %.2f (oracle %.2f, error %.2f)"%(self.Charge, BestMass, ParentMass, BestDiff)
            if BestDiff > 5:
                print "\n** WARNING: Parent mass is off by %.2f!\n"%BestDiff
        else: #chage given, calculate observed mass
            BestMass = self.MZ * self.Charge - (self.Charge - 1)*1.0078
        self.HypotheticalParentMass = ParentMass
        self.ObservedParentMass = BestMass
        self.ObservedParentMassError = abs(self.HypotheticalParentMass - self.ObservedParentMass)
        PySpectrum.SetParentMass(ParentMass, self.Charge)
        Annotation = self.Peptide.GetModdedName() # lacks prefix/suffix
        ## self.LabeledPeaks is list of (Mass, intensity, ion, amino index)
        self.LabeledPeaks = PySpectrum.LabelPeaks(Annotation, self.Charge) 
        if self.VerboseFileName:
            self.InspectFeatures = PySpectrum.ScorePeptideDetailed(Annotation, self.Charge)
            print "The MQScore for %s is %f"%(Annotation, self.InspectFeatures[0])
            self.MQScore = self.InspectFeatures[0] 


def LabelSpectrum(Spectrum, Peptide, PeakTolerance):
    Labeler = LabelClass()
    Labeler.Peptide = Peptide
    Labeler.SpectrumPath = Spectrum.FilePath
    Labeler.SpectrumFilePos = Spectrum.FilePos
    #print "Label.LabelSpectrum(%s:%s)"%(Spectrum.FilePath, Spectrum.FilePos)
    Labeler.LabelPeaks()
    # Paired iteration through Spectrum.Peaks and Labeler.LabeledPeaks:
    IndexA = 0
    IndexB = 0
    while IndexA < len(Spectrum.Peaks) and IndexB < len(Labeler.LabeledPeaks):
        Diff = Spectrum.Peaks[IndexA].Mass - Labeler.LabeledPeaks[IndexB][0]
        if Diff > 0.01:
            # Mass A is too large; let B catch up
            IndexB += 1
            continue
        if Diff < 0.01:
            # Mass A is too small; iterate forward
            IndexA += 1
            continue
        Spectrum.Peaks[PeakIndex].IonType = Labeler.LabeledPeaks[IndexB][2]
        Spectrum.Peaks[PeakIndex].AminoIndex = Labeler.LabeledPeaks[IndexB][3]
    return Spectrum

if __name__ == "__main__":
    Dymo = LabelClass()
    Dymo.ParseCommandLine(sys.argv[1:])
    Dymo.Main()
