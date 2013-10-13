"""
Check some coding standards
"""
import os
import sys
import re

def CheckCodeInDir(Dir):
    for FileName in os.listdir(Dir):
        (Stub, Extension) = os.path.splitext(FileName)
        if Extension not in (".py", ".c", ".h"):
            continue
        Path = os.path.join(Dir, FileName)
        File = open(Path, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            # Let's ignore what happens in comments; they could be valid English text
            # or ASCII art.
            StripLine = FileLine.strip()
            if not StripLine:
                continue
            if StripLine[0] == "#":
                continue
            if StripLine[:2] == "//":
                continue
            # Comma should always be followed by whitespace:
            Match = re.search(",[^ \t\r\n\]]", FileLine)
            if Match:
                print ",%s [%s] %s"%(FileName, LineNumber, FileLine[:60].strip())                
                continue
            # Equal sign should be followed by whitespace, or part of ==
            Match = re.search("=[^= \t\r\n]", FileLine)
            if Match:
                print "=%s [%s] %s"%(FileName, LineNumber, FileLine[:60].strip())
                continue
            # Minus, when not unary, should be followed by whitespace
            Match = re.search("[^'\"]\-[^0123456789 \t\r\n=>]", FileLine)
            if Match:
                print "-%s [%s] %s"%(FileName, LineNumber, FileLine[:60].strip())
                continue
            # Plus, when not unary, should be followed by whitespace
            Match = re.search("\+[^0123456789 \t\r\n=d+]", FileLine)
            if Match:
                print "+%s [%s] %s"%(FileName, LineNumber, FileLine[:60].strip())
                continue
            # Division should be followed by whitespace
            Match = re.search("[^<]/[^ \t\r\n=/]", FileLine)
            if Match:
                print "/%s [%s] %s"%(FileName, LineNumber, FileLine[:60].strip())
                continue
            # Pointer types should be referenced as a lexical item
            Match = re.search("char \*", FileLine)
            if Match:
                print "*%s [%s] %s"%(FileName, LineNumber, FileLine[:60].strip())
                continue
            

CheckCodeInDir(".")
CheckCodeInDir("PyInspect")