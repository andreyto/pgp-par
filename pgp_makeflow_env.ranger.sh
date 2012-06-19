#!/bin/bash
#This is responsible for setting up the proper environment
#as expected by the PGP components
#E.g. PValue.py uses extensions compiled against a specific version of
#Boost Python library.
#This should be sourced from task scripts generated to be called from
#the makeflow makefile
export PGP_HOME=/work/01241/atovchig/PGP/proteogenomics
export PGP_PEPNOVO_HOME=$PGP_HOME/../vendor/productionPepNovo
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/boost-1.43.0/lib:/usr/local/packages/atlas/lib
export PYTHONPATH=$PGP_HOME:$PYTHONPATH
export PGP_JAVA=$JAVA_HOME/bin/java
export PGP_PYTHON=$PYTHONHOME/bin/python

