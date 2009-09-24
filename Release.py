"""
RELEASE CHECKLIST:
- svn update
- Revise documentation
- Build inspect.exe (release version)
- Run SystemTest.py
- run Release.py
- unzip file to a different directory (ideally on different OS)
- make clean; make all; python SystemTest.py
- put release on peptide in /var/www/prod/www
- update links.  Run this on peptide:
  python ~/UpdateVersionNumber.py Inspect.XXXXXXXXX
- set bugs to RELEASED

This script helps prepare a release:
Copy files into the "zip" directory from the working directory.
To add a new file to the release, put a copy of it in the "zip" directory.
This script assumes the "zip" directory contains (some version of) each file desired.
"""
import os
import sys
import time
import shutil

class ReleaseTool:
    ReleaseDir = "zip"
    def __init__(self):
        pass
    def StampDocsWithDate(self):
        """
        Insert the version number into the documentation.
        Returns the number of errors encountered (0 if successful)
        """
        Path = os.path.join(self.ReleaseDir, "docs", "index.html")
        try:
            File = open(Path, "rb")
        except:
            print "* Can't open docs at '%s'"%Path
            self.ErrorCount += 1
            return
        HTML = File.read()
        File.close()
        if HTML.find("$VERSION_NUMBER")==-1:
            # this is ok:
            return
        HTML = HTML.replace("$VERSION_NUMBER", self.VersionNumber)
        File = open(Path, "wb")
        File.write(HTML)
        File.close()
        return
    def CopyFilesToDistributionDirectory(self):
        """
        The file "ReleaseFiles.txt" lists all the files which go into the external release.
        We don't release EVERYTHING from subversion, because that would be too much.
        """
        # Keep track of which files we've copied over into the release folder, so that
        # if we spy anything NOT on the list, we can complain loudly.  The keys are pathts
        # in the form of tuples, e.g. ("docs", "Searching.html")
        ReleasedFileNames = {}
        File = open("ReleaseFiles.txt", "rb")
        for FileLine in File.xreadlines():
            FileName = FileLine.strip()
            if not FileName:
                continue
            if FileName[0] == "#":
                continue # skip comments!
            if sys.platform == "win32":
                FileName = FileName.replace("/", "\\")
            else:
                FileName = FileName.replace("\\", "/")
            Key = tuple(FileName.replace("\\", "/").split("/"))
            ReleasedFileNames[Key] = 1
            if not os.path.exists(FileName):
                print "* Error: ReleaseFiles.txt told me to release '%s' but it doesn't exist"%FileName
                self.ErrorCount += 1
                continue
            TargetPath = os.path.join(self.ReleaseDir, FileName)
            # Check that the release directory exists:
            TargetDir = os.path.split(TargetPath)[0]
            if not os.path.exists(TargetDir):
                try:
                    os.makedirs(TargetDir)
                except:
                    print "* Error: Unable to build release directory '%s'"%TargetDir
                    self.ErrorCount += 1
                    continue
            # If the file is present and up-to-date, then don't bother recopying it:
            if os.path.exists(TargetPath):
                ReleaseModTime = os.stat(TargetPath).st_mtime
                DevModTime = os.stat(FileName).st_mtime
                if Key == ("docs", "index.html"):
                    ReleaseModTime = 0 # ALWAYS copy this over, then we'll edit the date.
                if ReleaseModTime > DevModTime:
                    print Key
                    print "* WARNING: Release file '%s' more recently modified than dev!  Copy manually to override."%TargetPath
                    self.ErrorCount += 1
                    continue
                if ReleaseModTime == DevModTime:
                    continue
            # Copy file to release directory.  Do *not* use shutil.copyfile, because
            # that changes the modification-time.  Use "copy" or "cp":
            if sys.platform == "win32":
                Command = "copy \"%s\" \"%s\""%(FileName, TargetPath)
            else:
                Command = "cp \"%s\" \"%s\""%(FileName, TargetPath)
            #print Command
            try:
                #shutil.copyfile(FileName, TargetPath)
                os.system(Command)
            except:
                print "* Error copying '%s' to '%s'"%(FileName, TargetPath)
                traceback.print_exc()
                self.ErrorCount += 1
                continue
        # Now, let's iterate over the files which are in the release directory.  If we see any files
        # that aren't on the list, then we'll complain.
        self.CheckForExtraReleaseFiles(ReleasedFileNames, self.ReleaseDir, [])
    def CheckForExtraReleaseFiles(self, ReleasedFileNames, Dir, PathBits):
        FileNames = os.listdir(Dir)
        for FileName in FileNames:
            Path = os.path.join(Dir, FileName)
            # Recurse into subdirectories:
            if os.path.isdir(Path):
                NewPathBits = PathBits[:]
                NewPathBits.append(FileName)
                self.CheckForExtraReleaseFiles(ReleasedFileNames, Path, NewPathBits)
                continue
            NewBits = PathBits[:]
            NewBits.append(FileName)
            if not ReleasedFileNames.has_key(tuple(NewBits)):
                print "* Error: Release directory contains unexpected file %s"%Path
                self.ErrorCount += 1
            else:
                pass
        
    def BuildRelease(self):
        Tuple = time.localtime()
        self.VersionNumber = "%s%02d%02d"%(Tuple[0], Tuple[1], Tuple[2])
        self.ErrorCount = 0
        self.CopyFilesToDistributionDirectory()
        #self.ErrorCount += self.PopulateDirectory(self.ReleaseDir, ".")
        #self.ErrorCount += self.PopulateDirectory(os.path.join(self.ReleaseDir, "docs"), "docs")
        self.StampDocsWithDate()
        if self.ErrorCount:
            print "%s Errors encountered in preparing release."%self.ErrorCount
            return 
        # Create a zip:
        ZipFileName = "Inspect.%s.zip"%(self.VersionNumber)
        os.chdir("zip")
        Command = "zip -r ..\\%s *"%ZipFileName
        print Command
        os.system(Command)
        os.chdir("..")
        print "Created %s"%ZipFileName
        print "pscp %s stanner@peptide.ucsd.edu:/var/www/prod/www"%ZipFileName
        
        
if __name__ == "__main__":
    Tool = ReleaseTool()
    Tool.BuildRelease()