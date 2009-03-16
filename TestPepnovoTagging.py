"""
Test pepnovo's tagging accuracy.
"""
import os
import sys
from Utils import *
Global.FixedMods = {"C":57.0518} #%%% USUALLY, but not always...
Initialize()

if len(sys.argv)>1:
    OracleFile = sys.argv[1]
    SpectrumDir = sys.argv[2]
else:
    OracleFile = "c:\\ms\\TrainingSetL.txt"
    SpectrumDir = "c:\\ms\\TrainingSetL"

FinalOutputFile = open("PepNovoTagging.txt", "w")

Histogram = {}
SpectrumCount = 0
File = open(OracleFile, "r")
for FileLine in File.xreadlines():
    Bits = FileLine.strip().split("\t")
    FileName = Bits[0]
    SpectrumPath = os.path.join(SpectrumDir, FileName)
    Charge = int(Bits[1])
    TruePeptide = GetPeptideFromModdedName(Bits[3])
    if len(TruePeptide.Modifications.keys()):
        continue # skip modded guys
    TruePeptide.Aminos = TruePeptide.Aminos.replace("I","L").replace("Q","K") #so that I/L subs are allowed
    try:
        os.remove("PepTags.txt")
    except:
        pass
    Command = """pepnovo_tags -dta "%s" -model LTQ_tryp.txt -num_tags 100 > PepTags.txt"""%SpectrumPath
    print Command
    os.system(Command)
    if not os.path.exists("PepTags.txt"):
        print "*** No tags found!"
        continue
    TagFile = open("PepTags.txt")
    Found = None
    TagIndex = 0
    for TagFileLine in TagFile.xreadlines():
        TagBits = TagFileLine.strip().split()
        if len(TagBits)<3:
            continue
        try:
            PrefixMass = float(TagBits[0])
        except:
            continue
        Pep = TagBits[1].replace("I","L").replace("Q","K")
        for Index in range(len(TruePeptide.Aminos)):
            if TruePeptide.Aminos[Index:Index+3] == Pep and abs(PrefixMass - TruePeptide.Masses[Index]) < 2:
                Found = TagIndex
                Histogram[TagIndex] = Histogram.get(TagIndex, 0) + 1
                print "Tag %s is correct."%TagIndex
                break
        TagIndex += 1
        if Found!=None:
            break
    FinalOutputFile.write("%s\t%s\t%s\t%s\t%s\n"%(Bits[0], Bits[1], Bits[2], Bits[3], Found))
    FinalOutputFile.flush()
    Histogram[None] = Histogram.get(None, 0) + 1
    SpectrumCount += 1
    
Cumulative = 0
for Index in range(100):
    Cumulative += Histogram.get(Index, 0)
    print "%s\t%s\t%s\t%s\t"%(Index, Histogram.get(Index, 0), Cumulative, 100*Cumulative/SpectrumCount)