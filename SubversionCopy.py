"""
SubversionCopy: Copy files from one directory tree to another.  Copy only the files
which are under version control.  Use case: Running the test suite on new code, *before*
checking the code in!
- Iterate over all the files under revision control in the current directory tree
- Copy each file to the target directory-tree, if it's not
  there (with current timestamp) already
"""
import sys
import os
import getopt
import shutil
import xml.dom.minidom
from Utils import *

UsageInfo = """
SubversionCopy: Copy revision-controlled files from one directory-tree
  to another.
Arguments:
  -r [DIR]: Source directory, defaults to current working directory
  -w [DIR]: Target directory
"""

class SubversionMirror:
    def __init__(self):
        pass
    def Mirror(self, SourceDir, TargetDir):
        EntryFilePath = os.path.join(SourceDir, ".svn", "entries")
        if not os.path.exists(EntryFilePath):
            return
        MakeDirectory(TargetDir)
        EntryFile = open(EntryFilePath, "rb")
        Text = EntryFile.read()
        EntryFile.close()
        DOM = xml.dom.minidom.parseString(Text)
        Entries = DOM.getElementsByTagName("entry")
        for Node in Entries:
            FileKind = str(Node.getAttribute("kind"))
            if FileKind != "file":
                continue
            # If the file is scheduled for deletion, then we needn't copy it over:
            Schedule = str(Node.getAttribute("schedule"))
            if Schedule == "delete":
                continue
            Baleeted = str(Node.getAttribute("deleted"))
            if Baleeted == "true":
                continue
            FileName = str(Node.getAttribute("name"))
            SourcePath = os.path.join(SourceDir, FileName)
            if not os.path.exists(SourcePath):
                print "* Warning: Revision-controlled file missing from '%s'"%SourcePath
                continue
            SourceModTime = os.stat(SourcePath).st_mtime
            TargetPath = os.path.join(TargetDir, FileName)
            if os.path.exists(TargetPath) and os.stat(TargetPath).st_mtime == SourceModTime:
                continue # it's up-to-date!
            #print "Copy %s to %s..."%(SourcePath, TargetPath)
            shutil.copyfile(SourcePath, TargetPath)
        # Handle sub-directories recursively:
        for FileName in os.listdir(SourceDir):
            SubDir = os.path.join(SourceDir, FileName)
            if os.path.isdir(SubDir):
                TargetSubDir = os.path.join(TargetDir, FileName)
                self.Mirror(SubDir, TargetSubDir)
    def Execute(self):
        self.Mirror(self.SourceDir, self.TargetDir)
    def ParseCommandLine(self, Arguments):
        self.TargetDir = None
        self.SourceDir = None
        (Options, Args) = getopt.getopt(sys.argv[1:], "r:w:")
        for (Option, Value) in Options:
            if Option == "-r":
                self.SourceDir = Value
            elif Option == "-w":
                self.TargetDir = Value
        if not self.SourceDir:
            self.SourceDir = os.getcwd()
        if not self.TargetDir:
            print UsageInfo
            sys.exit(-1)
            
if __name__ == "__main__":
    Mirror = SubversionMirror()
    Mirror.ParseCommandLine(sys.argv[1:])
    Mirror.Execute()