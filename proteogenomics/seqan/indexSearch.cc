#include <iostream>
#include <seqan/index.h>
#include <seqan/sequence.h>

using namespace seqan;

typedef String<AminoAcid> alphabet;

typedef Index< alphabet, Index_QGram< UngappedShape< 6 > > > MyIndex;

int main(int argc, char* argv[])
{
    ::std::string charTrie;
    ::std::cin >> charTrie;

    const alphabet trie(charTrie);
   // alphabet pep = "ERIERVEELL";

    MyIndex myIndex( trie );

    for(int i=1;i<argc;i++) {
        const alphabet fullPep( argv[i] );
        const int fullLen = length(fullPep);
        const alphabet pep = prefix(fullPep, 6);

        ::std::cout << pep << ::std::endl;
        ::std::vector<int> positions;

        Finder< MyIndex > finder(myIndex);
        Pattern<alphabet> pattern(pep);

        while( find(finder, pattern) ) {
            //positions.push_back( position(finder) );
            ::std::cout << position(finder) << ", ";
            ::std::cout << ::std::endl;
        }
/*
        ::std::vector<int>::const_iterator iter = positions.begin();
        while( iter != positions.end() ) {
            int place = *iter;
            iter++;
            const alphabet tpep = infixWithLength(trie,place,fullLen);
            if (tpep == fullPep) {
                ::std::cout << place << ", ";
            } else {
                ::std::cout << tpep << ", ";
            }
        }
        ::std::cout << ::std::endl;
*/
    }
}
