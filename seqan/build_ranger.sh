#!/bin/bash
#build IndexSearch python extension on TACC Ranger
SEQAN_HOME=$WORK/PGP/seqan/seqan-1.3.1/
#optimized build
SEQAN_CFLAGS="-O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0"
#debugging build
#SEQAN_CFLAGS="-g -O0 -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=0"
gcc -shared -fPIC $SEQAN_FLAGS -L$TACC_BOOST_LIB -lboost_python -I$TACC_BOOST_INC -I$TACC_PYTHON_DIR/include/python2.7 -L$TACC_PYTHON_LIB -I$SEQAN_HOME -o IndexSearch.so IndexSearch.cpp 

