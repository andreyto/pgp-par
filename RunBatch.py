"""
Simple batch processing
"""
import os

def SpliceFindBatch():
    ProteinCount = 57366
    X = 41000
    while X < ProteinCount:
        Command = "inspect splicefind database\\IPIv314.trie estsplicedb\genome.full50.dat %s %s"%(X, X+1000)
        print Command
        os.system(Command)
        Command = "move SSDatabaseScan.txt SSDB.Full50.%s.txt"%X
        print Command
        os.system(Command)
        X += 1000

def RunATSpliceFindBatch():
    ProteinCount = 33862
    X = 0
    while X < ProteinCount:
        Command = "inspect splicefind database\\arabidopsis.trie atestsplicedb\genome.dat %s %s"%(X, X+1000)
        print Command
        os.system(Command)
        Command = "move SSDatabaseScan.txt SSDB.AT.%s.txt"%X
        print Command
        os.system(Command)
        X += 1000
        
def GeneParseBatch():
    for X in range(1, 49):
        Command = "ParseGeneID.py %s 0 > PGI.%s.0.txt"%(X,X)
        print Command
        os.system(Command)
        Command = "ParseGeneID.py %s 1 > PGI.%s.1.txt"%(X,X)
        print Command
        os.system(Command)    

def BuildSpliceDBBatch():
    for X in range(1, 49):
        Command = "inspect %s 0"%X
        print Command
        os.system(Command)
        Command = "inspect %s 1"%X
        print Command
        os.system(Command)    

def BuildReversedExonGraphBatch():
    for X in range(1, 49):
        Command = "BuildReversedExonGraph.py ESTSpliceDB\\%s+.dat ReverseSpliceDB\\%s+.dat"%(X, X)
        print Command
        os.system(Command)
        Command = "BuildReversedExonGraph.py ESTSpliceDB\\%s-.dat ReverseSpliceDB\\%s-.dat"%(X, X)
        print Command
        os.system(Command)    

def RunGeneMapperBatch():
    for X in range(1, 49):
        Command = "GeneMapper.py %s"%X
        print Command
        os.system(Command)

def RunESTCoverageBatch():
    for X in range(1, 49):
        Command = "ReportESTGFProteinCoverage.py %s"%X
        print Command
        os.system(Command)

def ESTRescueBatch():
    for X in range(1, 49):
        Command = "SpliceRescueEST.py %s 0"%X
        print Command
        os.system(Command)
        Command = "SpliceRescueEST.py %s 1"%X
        print Command
        os.system(Command)    

    
if __name__ == "__main__":
    #BuildSpliceDBBatch()
    #SpliceFindBatch()
    #GeneParseBatch()
    #BuildReversedExonGraphBatch()
    #RunGeneMapperBatch()
    #RunESTCoverageBatch()
    #RunATSpliceFindBatch()
    ESTRescueBatch()