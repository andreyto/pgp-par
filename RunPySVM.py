"""
Wrapper for PySVM
"""
import os
import sys
import traceback
try:
    import PySVM
except:
    print "(Warning: PySVM not imported - SVM training not available)"

def Predict(FeaturePath, ModelPath, OutputPath):
    PySVM.LoadModel(ModelPath)
    InputFile = open(FeaturePath, "rb")
    OutputFile = open(OutputPath, "wb")
    for FileLine in InputFile.xreadlines():
        Bits = FileLine.split()
        FeatureVector = []
        for Bit in Bits[1:]:
            ColonPos = Bit.find(":")
            if ColonPos == -1:
                continue
            FeatureIndex = int(Bit[:ColonPos]) - 1
            while len(FeatureVector) <= FeatureIndex:
                FeatureVector.append(0)
            FeatureVector[FeatureIndex] = float(Bit[ColonPos + 1:])
        Score = PySVM.Score(FeatureVector)
        OutputFile.write("%s\n"%Score)
    InputFile.close()
    OutputFile.close()
    

if __name__ == "__main__":
    Predict("TestFeatures.SVMScaled.txt", "SVM.model", "SVMPrediction.pytxt")