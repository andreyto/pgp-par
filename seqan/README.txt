This directory is about linking the seqan (http://www.seqan.de) libraries
to python. Initial this is done via the boost python wrappers:
http://www.boost.org/doc/libs/1_43_0/libs/python/doc/index.html

About the files
    The IndexSearch.cpp is the file linking seqan & python via boost.

    indexSearch.cc is a C++ only test program. I built it as such: 
    g++ -O3 -I$das/devel/seqan/projects/library -o indexSearch indexSearch.cc
    Where $das was /usr/local/devel/DAS/users/eventer. Any copy of seqan
    should be sufficient.

    ExactSearch.* and indexSearchSet.* are other explorations and or dead
    ends and may be deleted.

How to Build the module
    A copy of the seqan libraries and boost libraries must be available.
    The boost build tool bjam is also needed to build the module.
    By default it appears the user needs write access to the boost libraries for
    the bjam build.

    The user-config.jam file points to the python to use. This file can be put
    in the users home dir.

    The Jamroot and boost-build.jam files control the bjam build. This is where
    the seqan and boost library paths are specified.

    To use a newer g++ /usr/local/packages/gcc-4.4.3/bin should be in the users
    PATH before /usr/local/bin and /usr/bin.

    Building the module is then just: bjam --toolset=gcc

    This builds the modules into the bin/gcc-x.x.x/debug/* tree
    
Running the code 
    The users PYTHONPATH env var must include the path to the IndexSearch.so.
    We should be able to symlink or copy that from the build path above to the
    main proteogenomics dir.
    I built with a non default gcc, so the LD_LIBRARY_PATH env variable needed to
    include:
    /usr/local/packages/gcc-4.4.3/lib64:/usr/local/packages/boost-1.43.0/lib

    Then the users should be able to run the test_indexSearch.py program.
