import os
import sys
import traceback
import struct

Dir = "PATargets"
DBFile = open(os.path.join(Dir, "PATargets.dat"), "wb")
IndexFile = open(os.path.join(Dir, "PATargets.index"), "wb")
FilePos = 0
for FileName in os.listdir(Dir):
    (Stub, Extension) = os.path.splitext(FileName)
    if Stub.lower() == "patargets":
        continue
    if Extension.lower() != ".dat":
        continue
    IndexFile.write(struct.pack("<qi80s", 0, FilePos, Stub))
    File = open(os.path.join(Dir, FileName), "rb")
    Stuff = File.read()
    File.close()
    FilePos += len(Stuff)
    DBFile.write(Stuff)

    
    