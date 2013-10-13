"""
Merge a directory full of .dta files into a single, more managable, .mgf file.
"""
import os
import sys
import getopt
            
UsageInfo = """
MergeDTAFiles.py - Merge a directory of .dta files into a single .mgf file.

Options:
 -r [InputDir]: .dta files will be read from this directory
    (and all subdirectories!)
 -w [OutputDir]: .mgf files will be written into this directory
 -s: If set, then scans will be arbitrarily renumbered.
 -m: If set, then keep only one charge state if the input .dta files include
     redundant files with different charge states.
"""

def MergeDirFilesMS2(OutputFile, DirName):
    ScanCount = 0
    for FileName in os.listdir(DirName):
        Path = os.path.join(DirName, FileName)
        if os.path.isdir(Path):
            ScanCount += MergeDirFiles(OutputFile, Path)
        else:
            Dotty = FileName.split(".")
            # Only merge .dta files.  (TODO: allow .pkl or other abnormal file extensions)
            if Dotty[-1] != "dta":
                continue
            ScanCount += 1
            ScanNumber = Dotty[-3]
            OutputFile.write(":%s.%s.0\n"%(ScanNumber, ScanNumber))
            DTAFile = open(FileName, "r")
            OutputFile.write(DTAFile.read().strip())
            OutputFile.write("\n")
    return ScanCount


class MergeMaster:
    def __init__(self):
        self.RenumberScans = 0
        self.NextScanNumber = 1
        self.DTADirectory = None
        self.OutputDirectory = os.getcwd()
        self.OverallScanCount = 0
        self.MGFFileCount = 0
        self.MergeChargeStatesFlag = 0
    def MergeDirFilesMGF(self, DirName):
        SortedFileNames = []
        # First, get a list of all the scan numbers, so that we can write out
        # the DTA files in order by scan number.
        ScanNumbersSeen = {}
        FileNames = os.listdir(DirName)
        FileNames.sort()
        for FileName in FileNames:
            Path = os.path.join(DirName, FileName)
            if os.path.isdir(Path):
                self.MergeDirFilesMGF(Path)
            else:
                Dotty = FileName.split(".")
                if Dotty[-1] != "dta":
                    continue
                # ASSUMPTION: dta filenames have this structure:
                # Stuff.ScanNumber.Charge.dta
                ScanNumber = int(Dotty[-3])
                if self.MergeChargeStatesFlag and ScanNumbersSeen.has_key(ScanNumber):
                    continue
                ScanNumbersSeen[ScanNumber] = 1
                SortedFileNames.append((ScanNumber, FileName))
        print DirName, len(SortedFileNames)
        # If there were no .dta files in this directory, then return:
        if len(SortedFileNames) == 0:
            return
        # Output file name is taken from the directory name:
        Stub = os.path.split(DirName)[1]
        OutputPath = os.path.join(self.OutputDirectory, "%s.mgf"%Stub)
        print "Merge %s scans into %s..."%(len(SortedFileNames), OutputPath)
        OutputFile = open(OutputPath, "wb")
        self.MGFFileCount += 1
        self.NextScanNumber = 1
        SortedFileNames.sort()
        for (ScanNumber, FileName) in SortedFileNames:
            if self.RenumberScans:
                ScanNumber = self.NextScanNumber
                self.NextScanNumber += 1
            Path = os.path.join(DirName, FileName)
            # SCAN HEADER:
            OutputFile.write("#\n")
            OutputFile.write("BEGIN IONS\n")
            OutputFile.write("TITLE=%s\n"%FileName)
            OutputFile.write("SCAN=%s\n"%ScanNumber)
            DTAFile = open(Path, "r")
            HeaderLine = DTAFile.readline()
            try:
                (ParentMass, Charge) = HeaderLine.split()
                Charge = int(Charge)
                ParentMass = float(ParentMass)            
            except:
                print "*** Error parsing top line of:", Path
                print HeaderLine.strip()
                traceback.print_exc()
                continue
            OutputFile.write("CHARGE=%s\n"%Charge)
            if Charge == 0:
                MZ = ParentMass
            else:
                MZ = (ParentMass + (Charge - 1) * 1.0078) / Charge
            # IMPORTANT NOTE: The .mgf format uses the M/Z value for "PEPMASS",
            # not the parent mass!
            OutputFile.write("PEPMASS=%s\n"%MZ)
            # Scan body:
            OutputFile.write(DTAFile.read().strip())
            OutputFile.write("\n")
            # SCAN FOOTER:
            OutputFile.write("END IONS\n")
        OutputFile.close()
        self.OverallScanCount += len(SortedFileNames)
        return len(SortedFileNames)
    def MergeFiles(self):
        self.MergeDirFilesMGF(self.DTADirectory)
        print "Merged a total of %s scans into %s .mgf files"%(self.OverallScanCount, self.MGFFileCount)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:sw:m")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r directory of .dta/.pkl files
                if not os.path.exists(Value):
                    print "** Error: couldn't find directory '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DTADirectory = Value
            elif Option == "-m":
                self.MergeChargeStatesFlag = 1
            elif Option == "-w":
                # -w output directory for .mgf files
                if not os.path.exists(Value):
                    try:
                        os.makedirs(Value)
                    except:
                        pass
                self.OutputDirectory = Value
            elif Option == "-s":
                self.RenumberScans = 1
        if not self.DTADirectory:
            print "** Please specify the name of a .dta file directory."
            print UsageInfo
            sys.exit(-1)

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Merger = MergeMaster()
    Merger.ParseCommandLine(sys.argv[1:])
    Merger.MergeFiles()