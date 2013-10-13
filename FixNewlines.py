"""
Change \r\n to \n, for stupid apache.  The 'pound-bang hack" (having #!/usr/bin/python on the first
line of a CGI script) only works if that CGI script has lines ending in \n, not \r\n.  (To test,
try running Main.py from the shell.  If it works, the line terminators are ok)
"""
import os
import sys
import traceback
from Utils import *
##
### If a filename was specified, then fix that file and exit:
##if len(sys.argv)>1:
##    FileName = sys.argv[1]
##    File = open(FileName, "r")
##    Text = File.read()
##    File.close()
##    FixedText = Text.replace("\r\n", "\n")
##    File = open(FileName, "w")
##    File.write(FixedText)
##    File.close()
##    sys.exit()
##
### By default: Fix all .py files in the current directory.    
##for FileName in os.listdir("."):
##    if FileName == "FixNewlines.py":
##        continue
##    Extension = os.path.splitext(FileName)[1].lower()
##    if Extension == ".py" and FileName[0]!=".":
##        File = open(FileName, "r")
##        Text = File.read()
##        if IS_WINDOWS:
##            #Text = Text.replace("#!/usr/bin/python", "#!python")
##            Text = Text.replace("#!/usr/bin/env python", "#!python")
##        else:
##            #Text = Text.replace("#!python", "#!/usr/bin/python")
##            Text = Text.replace("#!python", "#!/usr/bin/env python")
##        File.close()
##        FixedText = Text.replace("\r\n", "\n")
##        File = open(FileName, "w")
##        File.write(FixedText)
##        File.close()

# Also, add newlines at the end of all our .c and .h files,
# to reduce bogus build warnings:
for FileName in os.listdir("."):
    Extension = os.path.splitext(FileName)[1].lower()
    if Extension in (".c", ".h"):
        File = open(FileName, "rb")
        Text = File.read()
        if Text[-1] != '\n':
            Text += "\n"
        File.close()
        File = open(FileName, "wb")
        File.write(Text)
        File.close()
        
