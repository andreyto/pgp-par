"""
Combine training data from multiple files.  Randomly select n rows from each one.
"""
import random
random.seed(1)

FileNames = [("SVMTrain\\Ch2SimMod1Blind.Correct.txt",800),
             ("SVMTrain\\SimMod0LargeDB.Correct.txt",800),
             ("SVMTrain\\SimMod0Output.Correct.txt",150),
             ("SVMTrain\\SimMod0WrongDBOutput.Correct.txt",150),
             ]
OutputFileA = open("SVMTrain.txt","w")
OutputFileB = open("SVMTest.txt","w")
for (FileName, KeepCount) in FileNames:
    GoodLines = []
    BadLines = []
    File = open(FileName, "r")
    for FileLine in File.xreadlines():
        if FileLine[:2] == "+1":
            GoodLines.append(FileLine)
        elif FileLine[:2] == "-1":
            BadLines.append(FileLine)
    File.close()
    random.shuffle(GoodLines)
    random.shuffle(BadLines)
    KeepLines = min(KeepCount, max(len(GoodLines), len(BadLines)))
    print "Keep %s lines from %s"%(KeepLines, FileName)
    for Line in GoodLines[:KeepLines/2]:
        OutputFileA.write(Line)
    for Line in GoodLines[KeepLines/2:KeepLines]:
        OutputFileB.write(Line)
    for Line in BadLines[:KeepLines/2]:
        OutputFileA.write(Line)
    for Line in BadLines[KeepLines/2:KeepLines]:
        OutputFileB.write(Line)
OutputFileA.close()
OutputFileB.close()