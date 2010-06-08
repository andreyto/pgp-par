import IndexSearch
import sys

db = sys.stdin.readline()

pep = 'ERIERVEELL'

search = IndexSearch.IndexSearch(db)
for p in search.find(pep):
    print p

print 'Next pep'
pep = 'ERIER'
for p in search.find(pep):
    print p
