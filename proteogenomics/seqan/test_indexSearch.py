from IndexSearch import IndexSearch
import sys, os


def testSearchs(db):
    pep = 'ERIERVEELL'

    #search = IndexSearch.IndexSearch(db)
    search = IndexSearch(db)
    for p in search.find(pep):
        print p

    os.system( "grep 'Vm' /proc/%d/status" % os.getpid())

    print 'Next pep'
    pep = 'ERIER'
    for p in search.find(pep):
        print p

db = sys.stdin.readline()
testSearchs(db)

os.system( "grep 'Vm' /proc/%d/status" % os.getpid())
