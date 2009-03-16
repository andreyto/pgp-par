import os
from Utils import *
Initialize()

PMCFile = open("PMCTraining.txt", "w")
PMCFile1 = open("PMCTraining1.txt", "w")
PMCFile2 = open("PMCTraining2.txt", "w")
PMCFile3 = open("PMCTraining3.txt", "w")


def GetTrainingData(SourceFile, SourceDir):
    """
    For training of the charge state SVM and PMC SVM.  File names should indicate the true charge state.
    Some of the file names may be wrong, so...correct them, based upon a file of annotations.
    """
    # Read peptide annotations from the "oracle"
    File = open(SourceFile, "r")
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        FileName = Bits[0]
        Score = float(Bits[2])
        if Score < 2:
            continue
        Slashy = FileName.replace("/","\\").split("\\")
        if Slashy[-2] == "ikkb_z2":
            FilePath = os.path.join(SourceDir, Slashy[-1])
        else:
            FilePath = os.path.join(SourceDir, Slashy[-2], Slashy[-1])
        if not os.path.exists(FilePath):
            print "Missing:", FilePath
            continue
        Peptide = GetPeptideFromModdedName(Bits[1])
        TruePM = Peptide.Masses[-1] + 19
        SpectrumFile = open(FilePath, "r")
        TargetPath = os.path.join("c:\\ms\\PMCTraining", Slashy[-1])
        HeaderLine = SpectrumFile.readline()
        #File.close()
        HBits = HeaderLine.split()
        PM = float(HBits[0]) - 1.0078
        HCharge = int(HBits[1])
        Diff = abs(PM - TruePM)
        if TargetPath[-5] not in ("123"):
            TargetPath = TargetPath[:-4] + ".%d.dta"%HCharge
        if Diff > 100:
            #print "Bad charge, probably:", FileName
            MZ = (PM + (1.0078 * (HCharge-1))) / HCharge
            PM1 = MZ
            PM2 = MZ*2 - 1.0078
            PM3 = MZ*3 - 2*1.0078
            Diff1 = abs(TruePM - PM1)
            Diff2 = abs(TruePM - PM2)
            Diff3 = abs(TruePM - PM3)
            MinDiff = min(Diff1, Diff2, Diff3)
            if Diff1 == MinDiff:
                TargetPath = TargetPath[:-5] + "1" + FilePath[-4:]
                Charge = 1
                PM = PM1
            elif Diff2 == MinDiff:
                TargetPath = TargetPath[:-5] + "2" + FilePath[-4:]
                Charge = 2
                PM = PM2
            elif Diff3 == MinDiff:
                TargetPath = TargetPath[:-5] + "3" + FilePath[-4:]
                Charge = 3
                PM = PM3
            print "Corrected charge %s to charge %s; new filepath %s"%(HCharge, Charge, TargetPath)
            HeaderLine = "%.2f %d\n"%(PM + 1, Charge)
        else:
            Charge = HCharge
        TargetFile = open(TargetPath, "w")
        TargetFile.write(HeaderLine)
        TargetFile.write(SpectrumFile.read())
        SpectrumFile.close()
        TargetFile.close()
        Str = "%s\t%s\t%s\t%s\t\n"%(TargetPath, TruePM, Bits[1], FilePath)
        PMCFile.write(Str)
        if Charge == 1:
            PMCFile1.write(Str)
        elif Charge == 2:
            PMCFile2.write(Str)
        elif Charge == 3:
            PMCFile3.write(Str)

##FixFileNames(r"C:\ms\ChargeStateTesting")
##
def CopyISBPeptides():
    File = open("SergeiPeptides.txt", "r")
    for FileLine in File.xreadlines():
        Bits = FileLine.split()
        if len(Bits)<4:
            continue
        TrueCharge = int(Bits[1])
        Path = Bits[0][2:]
        if (Path[-2] == "/"):
            Path = Path[:-3] + Bits[1]
        Path += ".dta"
        Stub = Path.split(".")[0]
        Path = os.path.join("C:\\source\\bafna\\AllSpectra\\%s\\%s"%(Stub, Path))
        Command = """copy %s C:\\ms\\ChargeStateTraining"""%Path
        print Command
        os.system(Command)

if __name__ == "__main__":
    GetTrainingData("AriAnnotations.ISB.txt", "c:\\source\\bafna\\AllSpectra")
    GetTrainingData("AriAnnotations.Lens.txt", "c:\\ms\\LarryDavidOld")
    GetTrainingData("AriAnnotations.IKKB.txt", "c:\\ms\\ikkb_z2")
    GetTrainingData("AriAnnotations.OMICS04.txt", "c:\\ms\\OMICS04")
    