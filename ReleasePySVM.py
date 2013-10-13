"""
Script to build PySVM
"""
import sys
import distutils.core

def Main(Arguments):
    print "Prepping PySVM..."
    PySVMFileNames = ['PySVM/PySVM.c',"PySVM/svm-predict.c", "PySVM/svm.cpp",]
    PySVMExtension = distutils.core.Extension('PySVM', sources = PySVMFileNames)
    distutils.core.setup(name = 'PySVM', version = '1.0', ext_modules = [PySVMExtension],
          script_args = Arguments)

if __name__ == "__main__":
    Main(sys.argv[1:])
    
