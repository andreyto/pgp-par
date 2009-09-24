"""
SystemTest.py is the master test script for the Inspect toolkit.
It should be run after the inspect executable has been installed to
the current directory (and built, if necessary).

Run with no command-line arguments to perform a full test.
"""
import os
import sys
import traceback
try:
    from Utils import *
    Initialize()
except:
    print "** Error: Unable to load Utils!"
if hasattr(os, "sysconf"):
    IS_WINDOWS = 0
else:
    IS_WINDOWS = 1

if IS_WINDOWS:
    INSPECT_EXECUTABLE = "inspect"
else:
    INSPECT_EXECUTABLE = "./inspect"

SystemTestDir = "SystemTest"

class InspectRunner:
    def __init__(self):
        self.ErrorCount = 0
        self.TestsRun = 0
        self.TempOutputName = "SystemTestTemp.txt"
    def RunTestSearch(self, InputFileName, DesiredPeptide):
        "Run inspect, and verify that the desired peptide is the top match."
        Command = "%s -i %s -o %s"%(INSPECT_EXECUTABLE, InputFileName, self.TempOutputName)
        self.TestsRun += 1
        try:
            # Remove old output before running test:
            if os.path.exists(self.TempOutputName):
                os.remove(self.TempOutputName)
            # Run inspect:
            os.system(Command)
        except:
            traceback.print_exc()
            self.ErrorCount += 1
            return
        self.VerifyTestSearchResults(InputFileName, self.TempOutputName, DesiredPeptide)
    def VerifyTestSearchResults(self, InputFileName, OutputFileName, DesiredPeptide):
        if not os.path.exists(OutputFileName):
            print "** Error: No test output written for input '%s' to %s"%(InputFileName, OutputFileName)
            self.ErrorCount += 1
            return    
        File = open(OutputFileName, "rb")
        GoodHitPosition = None
        HitIndex = 0
        TopHit = None
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                Score = float(Bits[5])
            except:
                continue # header line
            if not TopHit:
                TopHit = Bits[2][2:-2]
            HitIndex += 1
            if Bits[2][2:-2] == DesiredPeptide:
                GoodHitPosition = HitIndex
                break
        if GoodHitPosition == 1:
            print "Test '%s' passed - top hit was '%s'"%(InputFileName, DesiredPeptide)
            return
        self.ErrorCount += 1
        print "** Error for test '%s':\n  Top hit was '%s'\n  Desired hit '%s' was seen at position %s"%(InputFileName, TopHit, DesiredPeptide, GoodHitPosition)
    def Summarize(self):
        print 
        print "-=- "*18
        print "System test summary: Ran %s tests, encountered %s errors."%(self.TestsRun, self.ErrorCount)
    def TestParentMassCorrection(self):
        # TestSpectra.pkl:
        # K.VLVLDTDYK.K, K.CLMEGAGDVAFVK.H, R.TPEVDDEALEK.F
        InputFileName = os.path.join(SystemTestDir, "TestPMC.txt")
        Command = "%s -i %s -o %s"%(INSPECT_EXECUTABLE, InputFileName, self.TempOutputName)
        self.TestsRun += 1
        os.system(Command)
        if not os.path.exists(self.TempOutputName):
            print "Error: TestPMC produced no output!"
            self.ErrorCount += 1
            return
        File = open(self.TempOutputName, "rb")
        Bits = File.readline().split("\t")
        File.close()
        if len(Bits) < 5:
            print "* Error: TestPMC produced invalid output!"
            self.ErrorCount += 1
            return
        Mass = float(Bits[3])
        Charge = int(Bits[4])
        DesiredCharge = 2
        DesiredMass = 1065.6
        if Charge != DesiredCharge or abs(Mass - DesiredMass) > 1.0:
            print "* Error: TestPMC produced invalid charge+mass (%s, %s), should be (%s, %s)"%(Charge, Mass, DesiredCharge, DesiredMass)
        else:
            print "TestPMC successful: Parent mass %s within tolerance"%Mass
        print "Parent mass correction complete."
            
    def TestTagging(self, InputFileName, Annotation, TagLength = None):
        """
        Run inspect in tag-generation mode.  Verify that one or more of the
        tags are correct for the target peptide (Annotation, a string).
        """
        DesiredPeptide = GetPeptideFromModdedName(Annotation)
        Command = "%s -i %s -o %s"%(INSPECT_EXECUTABLE, InputFileName, self.TempOutputName)
        self.TestsRun += 1
        try:
            # Remove old output before running test:
            if os.path.exists(self.TempOutputName):
                os.remove(self.TempOutputName)
            # Run inspect:
            os.system(Command)
        except:
            traceback.print_exc()
            self.ErrorCount += 1
            return
        if not os.path.exists(self.TempOutputName):
            print "** Error: No test output written for input '%s' to %s"%(InputFileName, self.TempOutputName)
            self.ErrorCount += 1
            return
        ValidTagCount = 0
        TagCount = 0
        File = open(self.TempOutputName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if FileLine[0] == "#" or len(Bits) < 7:
                continue
            #print Bits
            TagAminos = Bits[5]
            if TagLength != None and len(TagAminos) != TagLength:
                print "* Error in test '%s': Tag has length %s != %s"%(InputFileName, len(TagAminos), TagLength)
                self.ErrorCount += 1
            Tag = PeptideClass(Bits[5])
            Tag.PrefixMass = float(Bits[4])
            Tag.SuffixMass = float(Bits[6])
            if DesiredPeptide.IsValidTag(Tag):
                ValidTagCount += 1
            TagCount += 1
        if not ValidTagCount:
            print "* Test '%s' failed: No valid tags among %s attempts, for peptide %s"%(InputFileName, TagCount, DesiredPeptide.GetModdedName())
            self.ErrorCount += 1
        else:
            print "Tag test successful - found %s valid tags"%ValidTagCount
    def TestMS2DBConstruction(self):
        InputFileName = os.path.join(SystemTestDir, "BuildSimpleChromosome.txt")
        TempDBPath = "Temp.ms2db"
        Command = "%s -i %s -o %s"%(INSPECT_EXECUTABLE, InputFileName, TempDBPath)
        try:
            os.system(Command)
        except:
            print Command
            traceback.print_exc()
            self.ErrorCount += 1
        try:
            File = open(TempDBPath)
        except:
            print "** MS2DB test failed: No db constructed"
            self.ErrorCount += 1
            return
        MS2DB = File.read()
        File.close()
        Pos = MS2DB.find("AAPSASSTPTII")
        if Pos == -1:
            print "** MS2DB test failed: Expected peptide '%s' not present"%"SQVSFWLNRESVP"
            self.ErrorCount += 1
        # Now that the database has been constructed, let's search it:
        TempInputFileName = "TempMS2DB.in"
        TempScriptFile = open(TempInputFileName, "wb")
        TempScriptFile.write("db,%s\n"%os.path.abspath(TempDBPath))
        TempScriptFile.write("spectra,SystemTest/TestSpectrum.dta\n")
        TempScriptFile.write("protease,None\n")
        TempScriptFile.write("mod,+57,C,fix\n")
        TempScriptFile.close()
        Command = "%s -i %s -o %s"%(INSPECT_EXECUTABLE, TempInputFileName, self.TempOutputName)
        self.TestsRun += 1
        try:
            # Remove old output before running test:
            if os.path.exists(self.TempOutputName):
                os.remove(self.TempOutputName)
            # Run inspect:
            os.system(Command)
        except:
            traceback.print_exc()
            self.ErrorCount += 1
            return
        self.VerifyTestSearchResults(TempInputFileName, self.TempOutputName, "VKEAMAPK")
        try:
            os.remove(TempInputFileName)
        except:
            pass
        print "MS2DB search complete"
    def RunTests(self):
        self.TestMS2DBConstruction()
        self.TestTagging(os.path.join(SystemTestDir, "TestInputTag1.txt"), "VKEAMAPK", TagLength = 1)
        self.TestTagging(os.path.join(SystemTestDir, "TestInputTag3.txt"), "VKEAMAPK", TagLength = 3)
        self.RunTestSearch(os.path.join(SystemTestDir, "TestInput.txt"), "VKEAMGuserPK")
        self.RunTestSearch(os.path.join(SystemTestDir, "TestInputMod.txt"), "VKEAMG+14PK")
        self.RunTestSearch(os.path.join(SystemTestDir, "TestMS2.txt"), "AAEAATTDLTYR")
        self.TestParentMassCorrection()
        self.Summarize()
    
if __name__ == "__main__":
    Runner = InspectRunner()
    Runner.RunTests()