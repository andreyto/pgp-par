#include <iostream>
#include <seqan/find.h>
#include <seqan/sequence.h>

//typedef seqan::String<seqan::AminoAcid, seqan::Packed<> > alphabet;
typedef seqan::String<char> alphabet;

int main(int argc, char* argv[])
{
    std::string charTrie;
    std::cin >> charTrie;

    alphabet trie(charTrie);
   // alphabet pep = "ERIERVEELL";

    seqan::Finder<alphabet> finder(trie);
    for(int i=1;i<argc;i++) {
        alphabet pep = argv[i];
        std::cout << pep << std::endl;
        seqan::Pattern<alphabet, seqan::Horspool> pattern(pep);
        while( seqan::find(finder, pattern) ) {
            std::cout << position(finder) << ", ";
        }
        std::cout << std::endl;
    }
}
