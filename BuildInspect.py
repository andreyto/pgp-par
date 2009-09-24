"""
Python script to build Inspect.  An alternative to makefiles.
"""
import sys
import distutils
import distutils.command.build
import distutils.ccompiler

def BuildInspect(BuildNow = 0):
    InspectSourceFiles = [
        "base64.c", "BN.c", "BuildMS2DB.c", "ChargeState.c", "CMemLeak.c",
        "Errors.c", "ExonGraphAlign.c", 
        "FreeMod.c", "IonScoring.c", "LDA.c", "main.c", "Mods.c",
        "MS2DB.c", "ParentMass.c", "ParseInput.c", 
        "ParseXML.c", "PValue.c", 
        "Run.c", "Score.c", "Scorpion.c", "SNP.c", "Spectrum.c", "Spliced.c", 
        "SpliceDB.c", "SpliceScan.c", "SVM.c", "Tagger.c", "Trie.c", "Utils.c"
        ]
    ExtraIncludeDirectories = ["expat\\lib",]
    class MyBuildClass(distutils.command.build.build):
        def build_opt(self):
            CC = distutils.ccompiler.new_compiler()
            #if sys.platform != 'win32':
            #    CC.add_library('m')
            #import os.path
            print dir(CC)
            CC.library_dirs.append("expat/lib/release")
            if sys.platform == "win32":
                CC.add_library("libexpat")
            else:
                CC.add_library("expat") # not "libexpat", that won't work on Linux.
                CC.add_library("m")
            CC.set_include_dirs(ExtraIncludeDirectories)
            opt_obj = CC.compile(InspectSourceFiles)
            CC.link_executable(opt_obj, "inspect")
        def run(self):
            self.build_opt()
            distutils.command.build.build.run(self)
    if BuildNow:
        Dist = distutils.dist.Distribution()
        Dist.parse_config_files()
        Dist.cmdclass["build"] = MyBuildClass
        Dist.commands = ["build"]
        Dist.run_commands()
    else:
        distutils.core.setup(cmdclass = {"build":MyBuildClass,})
    

if __name__ == "__main__":
    #sys.argv = ["", "build"]
    BuildInspect()
    