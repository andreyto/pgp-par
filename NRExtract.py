"""
Extract peptpides with good species from NR.
Examples of name line:
>gi|4034483|emb|CAA10169.1| RAGE protein [Homo sapiens]
>gi|445136|prf||1908436B heat shock protein 16.9
    ^^^^^^ ^^^  ^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^
    ncbiid  |      |     everything after space is description
            |      |
         source    |
       database    source db ID

The ncbiID gives a usable link to the NCBI database.  The other ID may link to EMBL, to SWISS-PROT,
or somewhere else.  For now, let's not worry about it!
"""

WritingThis = 0
File = open("Database\\nr.fasta", "rb")
OutFile = open("Database\\MouseNR2.fasta", "wb")
for FileLine in File.xreadlines():
    if FileLine[0] == ">":
        WritingThis = 0
        Stuff = FileLine.lower()
        if Stuff.find("musculus")!=-1 or Stuff.find("mouse")!=-1:
            WritingThis = 1
        if FileLine.find("56205289")!=-1:
            print
            print FileLine, Stuff, WritingThis
    if WritingThis:
        OutFile.write(FileLine)
File.close()
OutFile.close()