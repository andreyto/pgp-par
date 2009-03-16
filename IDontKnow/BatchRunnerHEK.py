import BatchRunner
import os
import sys
import traceback

#RequiredString = "H293CoCl2-total-try-40uL-std-a-200ug-2D34-LTQ1"
RequiredString = "H293CoCl2-total-try-4uL-std-a-200ug-2D34-LTQ1"


def Main():
    Runner = BatchRunner.BatchRunner()
    Runner.ResultsDirectory = "f:\\ftproot\\briggs\\HEK293\\ResultsBlind1"
    # Add a list of mzxml files:
    Dir = "f:\\ftproot\\briggs\\HEK293"
    for SubDirectory in os.listdir(Dir):
        SubDir = os.path.join(Dir, SubDirectory)
        if not os.path.isdir(SubDir):
            continue
        DirMZXMLCount = 0
        for FileName in os.listdir(SubDir):
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension.lower() != ".mzxml":
                continue
            # Only look for ONE group at a time:
            if Stub.find(RequiredString) == -1:
                continue
            Path = os.path.join(SubDir, FileName)
            Runner.MZXMLPaths.append(Path)
            DirMZXMLCount += 1
        print "%s\t%s\t"%(SubDirectory, DirMZXMLCount)
    Runner.MZXMLPaths.sort()
    Runner.SetGridNBCR()
    #####################################################
    # Now we have the list of mzxml files to be searched.
    Runner.InitFileInfo()
    
    Runner.LoadScanCounts()
    Runner.GetScanCounts() # ONCE ONLY
    #Runner.AssessProgress()
    Runner.Main()
    
if __name__ == "__main__":
    Main()
