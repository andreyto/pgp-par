import os

for FileName in os.listdir("."):
    (Stub, Extension) = os.path.splitext(FileName)
    if Extension in (".c", ".h"):
        File = open(FileName, "rb")
        Text = File.read()
        File.close()
        Text = Text.replace("Copyright 2005", "Copyright 2006")
        File = open(FileName, "wb")
        File.write(Text)
        File.close()
        