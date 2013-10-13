"""
Compare 'annotated' spectra to 'non-annotated' spectra.
For now, assume that anything with an MQScore above a cutoff is successfully annotated, and anything with
MQScore below the cutoff is not annotated.  (This is rather cheesy...)
"""
import sys
import os
import MSSpectrum
MQScoreCutoff = 0


def PrintSpectrumQuality(SpectrumFile, MQScore):
    # Fix up the file path, if needed:
    Bits = SpectrumFile.replace("/","\\").split("\\")
    #SpectrumFile = "c:\\ms\\ISB\\%s\\%s"%(Bits[-2], Bits[-1])
    #SpectrumFile = "c:\\ms\\LarryDavidOld\\%s\\%s"%(Bits[-2], Bits[-1])
    FileName = Bits[-1]
    UnderBits = FileName.split("_")
    if FileName[:6] == "040914":
        SpectrumFile = "c:\\ms\\Venom\\Venom\\0914_croat\\%s\\%s"%(UnderBits[2], FileName)
    elif FileName[:6] == "040903":
        SpectrumFile = "c:\\ms\\Venom\\Venom\\0914_croat\\%s\\%s"%(UnderBits[2], FileName)
        
    #File = open(SpectrumFile, "r")
    Spectrum = MSSpectrum.SpectrumClass()
    Spectrum.ReadPeaks(SpectrumFile)
    TotalIntensity  = 0
    PeakCount = 0
    MaxPeakIntensity = 0
    IntensityList = []
    for Peak in Spectrum.Peaks:
        TotalIntensity  += Peak.Intensity
        IntensityList.append(Peak.Intensity)
        MaxPeakIntensity = max(MaxPeakIntensity, Peak.Intensity)
        PeakCount += 1
    IntensityList.sort()
    GrassLevel = Spectrum.ApplyWindowFilter([], (100,), 6)
    FilteredPeakCount = len(Spectrum.Peaks)
    if GrassLevel < 0:
        GrassLevel = IntensityList[PeakCount / 3]
    BestPeakToNoise = IntensityList[-1] / GrassLevel
    HighIntensity = 0
    HighIntensityCount = min(10, len(Spectrum.Peaks))
    for Intensity in IntensityList[-HighIntensityCount:]:
        HighIntensity += Intensity
    SignalToNoise = HighIntensity / (HighIntensityCount*GrassLevel)
    if MQScore > MQScoreCutoff:
        Annotated = 1
    else:
        Annotated = 0
    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%d\t"%(Annotated, SpectrumFile, MQScore,
                          PeakCount, TotalIntensity, TotalIntensity / PeakCount, BestPeakToNoise,
                          SignalToNoise, FilteredPeakCount, FilteredPeakCount / float(PeakCount), MaxPeakIntensity,
                                                                  Spectrum.ParentMass, Spectrum.Charge)

AnnotationFile = open(sys.argv[1], "r")
OldSpectrumFile = None
for FileLine in AnnotationFile.xreadlines():
    Bits = FileLine.split("\t")
    SpectrumFile = Bits[0]
    if SpectrumFile != OldSpectrumFile:
        OldSpectrumFile = SpectrumFile
        MQScore = float(Bits[5])
        PrintSpectrumQuality(SpectrumFile, MQScore)