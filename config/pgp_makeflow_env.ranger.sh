#!/bin/bash
#This is responsible for setting up the proper environment
#as expected by the PGP components
#E.g. PValue.py uses extensions compiled against a specific version of
#Boost Python library.
#This should be sourced from task scripts generated to be called from
#the makeflow makefile
export PGP_ROOT=/work/01241/atovchig/PGP
export PGP_HOME=$PGP_ROOT/proteogenomics
export PGP_VENDOR_HOME=$PGP_ROOT/vendor
export PGP_VENDOR_BIN=$PGP_VENDOR_HOME/bin
export PATH=$PGP_VENDOR_BIN:$PATH
export PGP_PYTHON_PREFIX=$PGP_VENDOR_HOME
export PYTHONPATH=$PGP_HOME:$PYTHONPATH

[ -n "$JAVA_HOME" ] || echo "JAVA_HOME undefined"
export PGP_JAVA=$JAVA_HOME/bin/java

[ -n "$PYTHONHOME" ] || echo "PYTHONHOME undefined"
export PGP_PYTHON=$PYTHONHOME/bin/python

if [ -f "$PGP_PYTHON" ]; then

    export PGP_PY_VER=`$PGP_PYTHON -c 'from distutils.sysconfig import *; print get_python_version()'`

    export PGP_PYCOMMON=$PGP_PYTHON_PREFIX/lib/python${PGP_PY_VER}/site-packages
    export PGP_PYDIST=$PGP_PYTHON_PREFIX/lib/python${PGP_PY_VER}/dist-packages
    export PYTHONPATH=${PGP_PYDIST}:${PGP_PYCOMMON}:${PYTHONPATH}
fi

export PGP_PEPNOVO_HOME=$PGP_VENDOR_HOME/productionPepNovo

