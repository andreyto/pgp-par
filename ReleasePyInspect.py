"""
Script to build PyInspect
"""
import sys
import os

import distutils.core

PyInspectFileNames = [
    "PyInspect/PyInspect.c", "PyInspect/PySpectrum.c", "PyInspect/PyUtils.c",
    "base64.c", "BN.c", "BuildMS2DB.c", "ChargeState.c", "CMemLeak.c",
    "Errors.c", "ExonGraphAlign.c", "FreeMod.c", "IonScoring.c", "LDA.c",
    "Mods.c", "MS2DB.c", "ParentMass.c", "ParseInput.c", "ParseXML.c", "PValue.c",
    "Run.c", "Score.c", "Scorpion.c", "SNP.c",
    "Spectrum.c", "Spliced.c", "SpliceDB.c",
    "SpliceScan.c", "SVM.c", "Tagger.c", "Trie.c", "Utils.c"
    ]

def Main(Arguments):
    print "Prepping PyInspect..."
    if sys.platform == "win32":
        LibraryList = ["libexpat"]
    else:
        LibraryList = ["expat"]
        
    PyInspectExtension = distutils.core.Extension('PyInspect',
        sources = PyInspectFileNames,
        include_dirs = [".", "expat/lib"],
        library_dirs = ["expat/lib/release"], 
        libraries = LibraryList)

    distutils.core.setup(name = 'PyInspect', version = '1.0', ext_modules=[PyInspectExtension],
        script_args = Arguments)

if __name__ == "__main__":
    Main(sys.argv[1:])
