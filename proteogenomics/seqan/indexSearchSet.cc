#include <iostream>
#include <seqan/index.h>
#include <seqan/sequence.h>

using namespace seqan;

typedef String<AminoAcid> alphabet;
typedef StringSet<alphabet> MyStrSet;

typedef Index< alphabet, Index_QGram< UngappedShape< 6 > > > MyIndex;

int main(int argc, char* argv[])
{
    MyStrSet allInput;
    Index< MyStrSet > myIndex;
    String<char> defline;
    alphabet inputSeq;
    ::std::vector< String<char> > defs;
    ::std::fstream inFasta;

    inFasta.open(argv[1], ::std::ios_base::in | ::std::ios_base::binary);

    int i = 0;
    while( ! inFasta.eof()) {
        readMeta( inFasta, defline, Fasta());
        defs.push_back( defline );
        read( inFasta, inputSeq, Fasta());
        appendValue( allInput, inputSeq );
//        ::std::cout << allInput[i] << ::std::endl;
        i++;
    }
    ::std::cout << "Num sequences:" << length(allInput) << " defs " << length(defs)
        << ::std::endl;

    exit(0);
   // alphabet pep = "ERIERVEELL";
/*
    for(int i=1;i<argc;i++) {
        const alphabet fullPep( argv[i] );
        const int fullLen = length(fullPep);
        const alphabet pep = prefix(fullPep, 6);

        ::std::cout << pep << ::std::endl;
        ::std::vector<int> positions;

        Finder< MyIndex > finder(trie);
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
    }
*/
}
