"""
This script finds the files and directories on a filesystem that are taking up lots of space.
"""
import os
import sys
import stat

REPORT_FILES = 100
REPORT_DIRECTORIES = 100

class MaidClass:
    def __init__(self):
        self.BigFiles = [(0, None)]
        self.BigDirectories = []
        self.BigTrees = []
        self.MP3Files = []
        self.FullSize = 0
        self.ImageFileList = open("c:\\ImageFileList.txt", "w")
    def Scan(self, RootDir):
        "Scan the files in (and below) this directory; return the tree size (in bytes)"
        SubDirs = []
        DirSize = 0
        TreeSize = 0
        if RootDir == "/proc":
            # Skip this directory, scanning it leads to infinite loopiness
            return 0
        # Add the size of each standard file, and queue up subdirectories for recursive processing:
        try:
            FileNames = os.listdir(RootDir)
        except:
            # There's no directory there!
            # Maybe some clever bastard deleted the directory out from under us.  Assume so, and carry on:
            return 0
        for FileName in FileNames:
            Path = os.path.join(RootDir, FileName)
            try:
                FileStat = os.stat(Path)
            except:
                # Assume it's a bogus file.  Skip it.
                continue
            Mode = FileStat[stat.ST_MODE]
            if stat.S_ISDIR(Mode):
                SubDirs.append(Path)
            elif stat.S_ISREG(Mode):
                if Path[-4:].lower() in (".mp3", ".ogg", ".xm", ".spc", ".m4u"):
                    self.MP3Files.append(Path)
                if Path[-4:].lower() in (".jpg", ".gif", ".bmp"):
                    self.ImageFileList.write(Path+"\n")
                FileSize = FileStat.st_size
                DirSize += FileSize
                # Add a big file to the "most wanted" list:
                if (len(self.BigFiles)<REPORT_FILES or FileSize > self.BigFiles[0][0]):
                    self.BigFiles.append((FileSize, Path))
                    self.BigFiles.sort()
                    self.BigFiles = self.BigFiles[-REPORT_FILES:]
        if (len(self.BigDirectories)<REPORT_DIRECTORIES or DirSize > self.BigDirectories[0][0]):
            self.BigDirectories.append((DirSize, RootDir))
            self.BigDirectories.sort()
            self.BigDirectories = self.BigDirectories[-REPORT_DIRECTORIES:]
        TreeSize = DirSize
        for Path in SubDirs:
            TreeSize += self.Scan(Path)
        if (len(self.BigTrees)<REPORT_DIRECTORIES or DirSize > self.BigTrees[0][0]):
            self.BigTrees.append((TreeSize, RootDir))
            self.BigTrees.sort()
            self.BigTrees = self.BigTrees[-REPORT_DIRECTORIES:]
        print RootDir
        return TreeSize
    def Report(self):
        print "-=- "*10
        print "Here are the biggest files:"
        self.BigFiles.reverse()
        for (Size, Path) in self.BigFiles:
            print Size, Path
        print "\n\n"
        print "-=- "*10
        print "Here are the biggest directories:"
        self.BigDirectories.reverse()
        for (Size, Path) in self.BigDirectories:
            print Size, Path
        print "\n\n"
        print "-=- "*10
        print "Here are the biggest directory subtrees:"
        self.BigTrees.reverse()
        for (Size, Path) in self.BigTrees:
            print Size, Path
        self.MP3Files.sort()
        DriveLetter = os.getcwd()[0]
        PlayList = open("%s:\\Everything.m3u"%DriveLetter, "wb")
        for MP3File in self.MP3Files:
            PlayList.write(MP3File + "\n")
        PlayList.close()
if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1] == "kill":
        # KILL!
        print "Killing files."
        for FileName in os.listdir("."):
            os.unlink(FileName)
        sys.exit(1)
    Betty = MaidClass()
    FullSize = Betty.Scan(os.path.abspath("."))
    print "Full size: %s (%.2fMb)"%(FullSize, FullSize / float(1024*1024))
    Betty.Report()
    