import ExactSearch
import sys

db = sys.stdin.readline()

pep = 'ERIERVEELL'

search = ExactSearch.ExactSearch()
for p in search.find(db,pep):
    print p

print 'Next pep'
search = ExactSearch.ExactSearch()
pep = 'ERIER'
for p in search.find(db,pep):
    print p
