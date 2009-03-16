from Utils import *
import Global
Initialize()

ErrorHistogram = {}
RelativeErrorHistogram = {}
IntegerMasses = {}
Aminos = "ACDEFGHIKLMNPQRSTVWY"
for AA in Aminos:
    IntegerMasses[AA] = int(round(Global.AminoMass[AA]))

IntegerMasses["C"] += 57
File = open("XGVKG1.txt", "rb")
LineNumber = 0
OldScan = None
for FileLine in File.xreadlines():
    LineNumber += 1
    if LineNumber > 50000:
        break
    Bits = FileLine.split("\t")
    try:
        Peptide = GetPeptideFromModdedName(Bits[2][2:-2])
    except:
        continue
    Scan = (Bits[0], Bits[1])
    if Scan == OldScan:
        continue
    OldScan = Scan
    TrueMass = Peptide.Masses[-1]
    if not TrueMass:
        continue
    IntMass = 0
    for AA in Peptide.Aminos:
        IntMass += IntegerMasses[AA]
    Diff = TrueMass - IntMass
    Bin = int(round(Diff * 10))
    ErrorHistogram[Bin] = ErrorHistogram.get(Bin, 0) + 1
    RelativeErrorHistogram[Bin] = RelativeErrorHistogram.get(Bin,0) + 1
    Ratio = IntMass / TrueMass
    print "%s\t%s\t%s\t%s\t"%(TrueMass, IntMass, Diff, Ratio)

##Keys =  ErrorHistogram.keys()
##Keys.sort()
##for Key in Keys:
##    print "%s\t%s\t"%(Key, ErrorHistogram[Key])

    