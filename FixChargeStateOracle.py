from Utils import *
Initialize()

Global.FixedMods = {"C":57.0518}

def FixFileNames(Dir):
    """
    For training of the charge state SVM.  File names must indicate the true charge state.
    Some of the file names may be wrong, so...correct them, based upon a file of annotations.
    """
    # Read peptide annotations from the "oracle"
    TruePMDict = {}
    File = open("AriAnnotations.Lens.txt", "r")
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        FileName = Bits[0]
        FileName = FileName.replace("/","\\").split("\\")[-1]
        Peptide = GetPeptideFromModdedName(Bits[1])
        TruePMDict[FileName] = Peptide.Masses[-1] + 19
    File.close()
    # Iterate over files in the directory.  For each one with incorrect apparent charge,
    # fix it:
    for FileName in os.listdir(Dir):
        TruePM = TruePMDict.get(FileName, None)
        if not TruePM:
            continue
        FilePath = os.path.join(Dir, FileName)
        File = open(FilePath, "r")
        FileLine = File.readline()
        File.close()
        Bits = FileLine.split()
        PM = float(Bits[0])
        Charge = int(Bits[1])
        Diff = abs(PM - TruePM)
        if Diff > 100:
            #print "Bad charge, probably:", FileName
            MZ = PM / Charge
            PM1 = MZ
            PM2 = MZ*2
            PM3 = MZ*3
            Diff1 = abs(TruePM - PM1)
            Diff2 = abs(TruePM - PM2)
            Diff3 = abs(TruePM - PM3)
            MinDiff = min(Diff1, Diff2, Diff3)
            if Diff1 == MinDiff:
                NewName = FilePath[:-5] + "1" + FilePath[-4:]
                Charge = 1
            elif Diff2 == MinDiff:
                NewName = FilePath[:-5] + "2" + FilePath[-4:]
                Charge = 2
            elif Diff3 == MinDiff:
                NewName = FilePath[:-5] + "3" + FilePath[-4:]
                Charge = 3
            #Command = "move %s %s"%(FilePath, NewName)
            #print Command
            File = open(FilePath, "r")
            DiscardLine = File.readline()
            NewFile = open(NewName, "w")
            NewFile.write("%.2f %d\n"%(TruePM + 1, Charge))
            NewFile.write(File.read())
            File.close()
            NewFile.close()
            Command = "del %s"%FilePath
            print Command
            os.system(Command)
        else:
            #print "Charge ok:", FileName, PM, Diff
            pass


FixFileNames(r"C:\ms\ChargeStateTesting")         