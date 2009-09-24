"""
PTMAnalysis.py
This wrapper script automates the steps involved in processing raw Insepct
results into finalized PTM site identifications. 
1. PValue.py - obtain a consistent pvalue across all runs
2. ComputePTMFeatures.py - group individual annotations into peptides and
    compute some features for each.
3. BuildMGF.py - builds a single .mgf file out of all the consensus spectra
    created in step 2 (in preparation for the Inspect run)
4. Inspect run.  Search clustered spectra (unmodified) against a large database
5. PTMSearchBigDB.py - Integrates the results of the inspect search against a
    large database, findig a delta-score.
6. TrainPTMFeatures.py - Computes the PValue of each site from a model
7. AdjustPTM.py - attempts to merge and reconcile sites
8. AdjustPTM.py (known chemistry) - attempts to find known explanations for
    the site

Depending on the size of your dataset, this program may take quite a while
(up to 1 day per million spectra).  It is reasonably easy to parallelize
the time-intensive steps (2) and (4) above; doing so is up to the user, since
compute clusters are so heterogenous.

To perform only a subset of the steps, you can designate start and stop steps.
This is a useful feature if perhaps the program crashes at some point, or if
you would like to run some steps on a grid, and some locally.
To start from a step other than one, the program assumes that all the previous
steps have been executed, and that their outputs are in the proper directories.
See the functions below for the expected directory and file names.
"""
import os
import getopt
import sys
import traceback
#import PValue
#import ComputePTMFeatures
#import TrainPTMFeatures
#import AdjustPTM
#import PTMSearchBigDB
#import BuildMGF
from Utils import *

UsageInfo = """
PTMAnalysis.py - produces a set of PTM sites.

Required options:
 -r [FileName] - The name of the results file to parse.  If a directory is
    specified, then all .txt files within the directory will be combined into
    one report
 -d [FileName] - The name of the database file (.trie format) searched.
 -w [FileName] - The final output filename
 -s [Dir] - Directory containing the MS/MS spectra
 -B [FilePath] - Path to the large database (for unmodified "decoy" search).
 -S [Value] - The fraction of sequences in the database that are shuffled.
 -k [FileName]: Known chemistry filename.  If specified, consider altering
   sites to match known chemical adducts; report the best site-score
   attainable by using known chemical adducts.

Additional options:
 -p [Value] - The pvalue cutoff for spectra annotations.  Default to 0.1
 -q [Value] - The pvalue cutoff for PTM site annotations. Default to 0.05
 -t [Path] - The path where all intermediate results will be written.
        Default to PTMTempFiles
 -i [Instrument] - The type of instrument that the Spectra were created on.
     Default is ESI-ION-TRAP.  Valid values are: QTOF, FT-Hybrid
     
Advanced options: to run only a subset of the steps
 -x [StartStep] - Start with step X (assume that previous steps are done)
 -y [StopStep] - Stop with step Y (inclusive)

Protein selection can be performed, replacing the protein identification
with a parsimonious set of protein IDs (using a simple iterative
approach).  The following options are required for protein selection:
 -a: Replace protein identifications with a "parsimonious" set of protein IDs.

"""

class Steps:
    First = 1
    PValue = 1
    ComputePTMFeatures = 2
    BuildMGF = 3
    RunInspect = 4
    SearchBigDB = 5
    TrainPTMFeatures = 6
    AdjustPeptides = 7
    AdjustToKnownChemistry = 8
    Last = 8

class WrapperClass:
    def __init__ (self):
        self.SpectrumPValue = 0.1 # for the pvalue.py program
        self.SitePValue = 0.05 # for the final output
        self.BasePath = "PTMTempFiles"
        self.PValuePath = None
        self.PValueOutput = None
        self.ComputePTMPath = None
        self.ComputePTMOutput = None
        self.TrainPTMPath = None
        self.TrainPTMOutput = None
        self.TrainPTMModelOutput = None
        self.AdjustDir = None
        self.AdjustOutputPath = None
        self.AdjustModelOutput = None
        self.SearchBigDBPath = None
        self.SearchBigDBOutput = None
        self.InputResults = None
        self.DatabaseFile = None
        self.PercentShuffled = None
        self.SpectraDir = None
        self.SelectProteins = 0 #default
        self.FinalOutputFile = None
        self.BigDB = None
        self.Instrument = "ESI-ION-TRAP"
        self.BuildMGFPath = None
        self.InspectOutDir = None
        self.MGFPath = None
        self.SpawnFlag = 0
        self.StartStep = Steps.First
        self.StopStep = Steps.Last
        self.KnownChemistryFileName = None
        self.TrainPTMModelType = "svm" # default
    def SetupDirectories(self):
        """
        Below the basepath there will be a group of directories, one for
        each major step
        """
        self.PValuePath = os.path.join(self.BasePath, "PValued")
        self.ComputePTMPath = os.path.join(self.BasePath, "ComputePTMFeatures")
        self.TrainPTMPath = os.path.join(self.BasePath, "TrainPTMFeatures")
        self.AdjustDir = os.path.join(self.BasePath, "AdjustPTM")
        self.SearchBigDBPath = os.path.join(self.BasePath, "SearchBigDB")
        self.BuildMGFPath = os.path.join(self.BasePath, "BuildMGF")
        self.InspectOutDir = os.path.join(self.BasePath, "InspectOut")
        print "Making temporary directories in %s for all intermediate output"%self.BasePath
        MakeDirectory(self.PValuePath)
        MakeDirectory(self.ComputePTMPath)
        MakeDirectory(self.TrainPTMPath)
        MakeDirectory(self.AdjustDir)
        MakeDirectory(self.SearchBigDBPath)
        MakeDirectory(self.BuildMGFPath)
        MakeDirectory(self.InspectOutDir)
    def RunPValue(self):
        """
        PValue.py
            -r InputResults
            -w OutputResults
            -S Percent of database shuffled (optional)
            -p pvalue cutoff
            -s Distribution file
            -i Distribution image
            -H write out results from shuffled protein
        """
        self.PValueOutput = self.PValuePath # default, a directory for directory input
        DistributionFile = os.path.join(self.PValuePath, "Distribution.txt")
        if not os.path.isdir(self.InputResults): # the InputResults is a single file
            FileName = os.path.split(self.InputResults)[1]
            self.PValueOutput = os.path.join(self.PValueOutput, FileName)
        PValueArgs = ""
        if self.PercentShuffled:
            PValueArgs = "-r %s -w %s -S %f -p %f -s %s -i -H" %(self.InputResults, self.PValueOutput, self.PercentShuffled, self.SpectrumPValue, DistributionFile)
        else:
            PValueArgs = "-r %s -w %s -p %f -s %s -i -H" %(self.InputResults, self.PValueOutput, self.SpectrumPValue, DistributionFile)
        if self.SelectProteins:
            PValueArgs += " -a -d %s"%self.DatabaseFile
        PValueArgs += " -b "
        if self.StartStep <= Steps.PValue and Steps.PValue <= self.StopStep:
            print "Step %s: PValue"%Steps.PValue
            print "Arguments: %s"%PValueArgs
            if self.SpawnFlag:
                Command = "python PValue.py %s"%PValueArgs
                print Command
                os.system(Command)
            else:
                import PValue
                ArgsList = PValueArgs.split()
                Parser = PValue.PValueParser()
                Parser.ParseCommandLine(ArgsList)
                PValue.Main(Parser)
                del PValue
        else:
            print "Skipping Step %s: PValue"%Steps.PValue
    def RunComputePTMFeatures(self):
        """
        ComputePTMFeatures.py
            -r InputResults
            -w OutputDir
            -d Database
            -s spectra
        """
        self.ComputePTMOutput = os.path.join(self.ComputePTMPath, "PTMFeatures.txt")
        Args = " -r %s -w %s -d %s -s %s"%(self.PValueOutput, self.ComputePTMPath, self.DatabaseFile, self.SpectraDir)
        if self.StartStep <= Steps.ComputePTMFeatures and Steps.ComputePTMFeatures <= self.StopStep:
            print "Step %s: ComputePTMFeatures"%Steps.ComputePTMFeatures
            print "Arguments: %s"%Args
            if self.SpawnFlag:
                Command = "ComputePTMFeatures.py %s"%Args
                print Command
                os.system(Command)
            else:
                import ComputePTMFeatures
                ArgsList = Args.split()
                Computer = ComputePTMFeatures.PTMFeatureComputer()
                Computer.ParseCommandLine(ArgsList)
                Computer.ComputeFeaturesMain()
                del ComputePTMFeatures
        else:
            print "Skipping Step %s: ComputePTMFeatures"%Steps.ComputePTMFeatures
    def RunBuildMGF(self):
        """
        BuildMGF.py
            -d PTM feature directory
            -m MGF file to make
        """
        self.MGFPath = os.path.join(self.BuildMGFPath, "spectra.mgf")
        Args = " -d %s -m %s"%(self.ComputePTMPath, self.MGFPath)
        
        if self.StartStep <= Steps.BuildMGF and Steps.BuildMGF <= self.StopStep:
            print "Step %s: BuildMGF"%Steps.BuildMGF
            print "Arguments:  %s"%Args
            if self.SpawnFlag:
                Command = "BuildMGF.py %s"%Args
                print Command
                os.system(Command)
            else:
                ArgsList = Args.split()
                import BuildMGF
                Builder = BuildMGF.MGFBuilder()
                Builder.ParseCommandLine(ArgsList)
                Builder.Main()
                del BuildMGF
        else:
            print "Skipping Step %s: BuildMGF"%Steps.BuildMGF
    def RunInspect(self):
        """
        Given that the mgf file was previously created, here we create an
        input file for Inspect and run it.
        """
        InspectExe = None
        if sys.platform == "win32":
            InspectExe = "Inspect.exe"
        else:
            InspectExe = "./inspect"
        InspectIn = "BigSearch.in"
        self.InspectOutFile = os.path.join(self.InspectOutDir, "Results.txt")
        Command = "%s -i %s -o %s"%(InspectExe, InspectIn, self.InspectOutFile)
        InFileCommands = "spectra,%s\n"%self.MGFPath
        InFileCommands += "instrument,%s\n"%self.Instrument
        InFileCommands += "protease,Trypsin\n"
        InFileCommands += "DB,%s\n"%self.BigDB
        InFileCommands += "mod,57,C,fix\n"
        InFileCommands += "tagcount,25\n"
        Handle = open(InspectIn, "wb")
        Handle.write(InFileCommands)
        Handle.close()
        if self.StartStep <= Steps.RunInspect and Steps.RunInspect <= self.StopStep:
            print "Step %s: Run Inspect, searching consensus spectra"%Steps.RunInspect
            print "Arguments: %s"%Command
            os.system(Command)
        else:
            print "Skipping Step %s: Run Inspect, searching consensus spectra"%Steps.RunInspect
    def RunPTMSearchBigDB(self):
        """
        PTMSearchBigDB.py
            -d PTM feature directory
            -w Outputfile to write
            -r Inspect Search Results (default to directory, not file)
        """
        self.SearchBigDBOutput = os.path.join(self.SearchBigDBPath, "Results.txt")
        # we use the InspectOutDir in case there are multiple files within the directory.
        Args = " -d %s -w %s -r %s"%(self.ComputePTMPath, self.SearchBigDBOutput, self.InspectOutDir)
        if self.StartStep <= Steps.SearchBigDB and Steps.SearchBigDB <= self.StopStep:
            print "Step %s: PTMSearchBigDB"%Steps.SearchBigDB
            print "Arguments:  %s"%Args
            if self.SpawnFlag:
                Command = "PTMSearchBigDB.py %s"%Args
                print Command
                os.system(Command)
            else:
                ArgsList = Args.split()
                import PTMSearchBigDB
                Searcher = PTMSearchBigDB.PTMSearcher()
                Searcher.ParseCommandLine(ArgsList)
                Searcher.Main()
                del PTMSearchBigDB
        else:
            print "Skipping Step %s: PTMSearchBigDB"%Steps.SearchBigDB
    def RunTrainPTMFeatures(self):
        """
        TrainPTMFeatures.py
            -u InputResults
            -v OutputResults
            -m ModelType
            -w OutputModel
        """
        self.TrainPTMOutput = os.path.join(self.TrainPTMPath, "Results.txt")
        self.TrainPTMModelOutput = os.path.join(self.TrainPTMPath, "model.%s.txt"%self.TrainPTMModelType)
        Args = "-u %s -v %s -m %s -w %s"%(self.SearchBigDBOutput, self.TrainPTMOutput, self.TrainPTMModelType, self.TrainPTMModelOutput)
        if self.StartStep <= Steps.TrainPTMFeatures and Steps.TrainPTMFeatures <= self.StopStep:
            print "Step %s: TrainPTMFeatures"%Steps.TrainPTMFeatures
            print "Arguments: %s"%Args
            if self.SpawnFlag:
                Command = "TrainPTMFeatures.py %s"%Args
                print Command
                os.system(Command)
            else:
                ArgsList = Args.split()
                import TrainPTMFeatures
                Trainer = TrainPTMFeatures.PTMFeatureTrainer()
                Trainer.ParseCommandLine(ArgsList)
                Trainer.TrainModel()
                del TrainPTMFeatures
        else:
            print "Skipping Step %s: TrainPTMFeatures"%Steps.TrainPTMFeatures
    def AdjustPeptides(self):
        """
        AdjustPTM.py
            -r InputResults
            -w OutputResults
            -d Database File
            -c Cluster directory from ComputePTMFeatures
            -m model INPUT filename
            -M model OUTPUT filename
        """
        self.AdjustOutputPath = os.path.join(self.AdjustDir, "Results.txt")
        self.AdjustModelOutput = os.path.join(self.AdjustDir, "model.%s.txt"%self.TrainPTMModelType)
        Args = "-r %s -w %s -d %s -c %s -m %s -M %s -z "%(self.TrainPTMOutput, self.AdjustOutputPath, self.DatabaseFile, self.ComputePTMPath, self.TrainPTMModelOutput, self.AdjustModelOutput)
        if self.StartStep <= Steps.AdjustPeptides and Steps.AdjustPeptides <= self.StopStep:
            print "Step %s: AdjustPTM"%Steps.AdjustPeptides
            print "Arguments: %s"%Args
            if self.SpawnFlag:
                Command = "AdjustPTM.py %s"%Args
                print Command
                os.system(Command)
            else:
                ArgsList = Args.split()
                import AdjustPTM
                Adjutant = AdjustPTM.PTMAdjuster()
                Adjutant.ParseCommandLine(ArgsList)
                Adjutant.Main()
                del AdjustPTM
        else:
            print "Skipping Step %s: AdjustPTM"%Steps.AdjustPeptides
    def AdjustToKnownChemistry(self):
        """
        AdjustPTM.py
            -r InputResults
            -w OutputResults
            -m input Model
            -d Database
            -c Clusters 
            -k Known PTM file
            -v verbose output file
        """
        if not self.KnownChemistryFileName:
            print "* Skipping AdjustToKnownChemistry: Requires a file (-k) of 'common' modifications."
            return
        KnownPTMVerboseOutputPath = os.path.join(self.BasePath, "KnownPTMOutput.txt")
        Args = "-r %s -w %s -m %s -d %s -c %s -k %s -v %s"%(self.AdjustOutputPath, self.FinalOutputFile, self.AdjustModelOutput, self.DatabaseFile, self.ComputePTMPath, self.KnownChemistryFileName, KnownPTMVerboseOutputPath)
        
        if self.StartStep <= Steps.AdjustToKnownChemistry and Steps.AdjustToKnownChemistry <= self.StopStep:
            print "Step %s: AdjustToKnownChemistry"%Steps.AdjustToKnownChemistry
            print "Arguments: %s"%Args
            if self.SpawnFlag == 1:
                Command = "AdjustPTM.py %s"%Args
                print Command
                os.system(Command)
            else:
                ArgsList = Args.split()
                import AdjustPTM
                Adjutant = AdjustPTM.PTMAdjuster()
                Adjutant.ParseCommandLine(ArgsList)
                Adjutant.Main()
                del AdjustPTM
        else:
            print "Skipping Step %s: AdjustToKnownChemistry"%Steps.AdjustToKnownChemistry
    def ParseCommandLine(self, Arguments):
        "Args is a list of arguments only (does not include sys.argv[0] == script name)"
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:s:p:q:t:S:B:i:x:y:k:bZm:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.InputResults = Value
            elif Option == "-m":
                self.TrainPTMModelType = Value # "svm" or "lda"
            elif Option == "-d":
                self.DatabaseFile = Value
            elif Option == "-w":
                self.FinalOutputFile = Value
            elif Option == "-s":
                self.SpectraDir = Value
            elif Option == "-p":
                self.SpectrumPValue = float(Value)
            elif Option == "-q":
                self.SitePValue = float(Value)
            elif Option == "-t":
                self.BasePath = Value
            elif Option == "-S":
                self.PercentShuffled = float (Value)
            elif Option == "-B":
                self.BigDB = Value
            elif Option == "-x":
                self.StartStep = int (Value)
            elif Option == "-y":
                self.StopStep = int (Value)
            elif Option == "-k":
                self.KnownChemistryFileName = Value
            elif Option == "-Z":
                self.SpawnFlag = 1
            else:
                print "** Unknown option:", Option, Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-s") or not OptionsSeen.has_key("-S"):
            print "Missing required options (r, d, w, s, S)"
            print UsageInfo
            sys.exit(-1)
        if not self.BigDB:
            print "Missing large DB (-B)"
            print UsageInfo
            sys.exit(-1)

def Main():
    "Main control box for the script"
    Wrap = WrapperClass()
    Wrap.ParseCommandLine(sys.argv[1:])
    Wrap.SetupDirectories()
    #now we run the parts of the scripts one after another
    print "\n*** Starting to run components ***"
    Wrap.RunPValue()
    Wrap.RunComputePTMFeatures()
    Wrap.RunBuildMGF()
    Wrap.RunInspect()
    Wrap.RunPTMSearchBigDB()
    Wrap.RunTrainPTMFeatures()
    Wrap.AdjustPeptides()
    Wrap.AdjustToKnownChemistry()
    

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "psyco not found - running without optimization"
    Main()    
        