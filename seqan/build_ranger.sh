SEQAN_HOME=$WORK/PGP/seqan/seqan-1.3.1/
gcc -shared -fPIC -L$TACC_BOOST_LIB -lboost_python -I$TACC_BOOST_INC -I$TACC_PYTHON_DIR/include/python2.7 -L$TACC_PYTHON_LIB -I$SEQAN_HOME -o IndexSearch.so IndexSearch.cpp 

