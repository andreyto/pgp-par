#!/bin/bash

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   Copyright J. Craig Venter Institute 2012
#   See docs.PGP/COPYING file distributed along with the proteogenomics 
#   package for the copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

set -ex
#build IndexSearch python extension on TACC Ranger
[ -n "$PGP_SEQAN_HOME" ] || exit 1
#optimized build
SEQAN_CFLAGS="-O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0"
#debugging build
#SEQAN_CFLAGS="-g -O0 -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=0"
gcc -shared -fPIC $SEQAN_CFLAGS -L$TACC_BOOST_LIB -lboost_python -I$TACC_BOOST_INC -I$TACC_PYTHON_DIR/include/python2.7 -L$TACC_PYTHON_LIB -I$PGP_SEQAN_HOME -o IndexSearch.so IndexSearch.cpp
