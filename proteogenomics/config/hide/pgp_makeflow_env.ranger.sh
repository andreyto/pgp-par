#!/bin/bash

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   Copyright J. Craig Venter Institute 2012
#   See docs.PGP/COPYING file distributed along with the proteogenomics 
#   package for the copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#This is responsible for setting up the proper environment
#as expected by the PGP components
#E.g. PValue.py uses extensions compiled against a specific version of
#Boost Python library.
#This should be sourced from task scripts generated to be called from
#the makeflow makefile
export PGP_ROOT=@PGP_ROOT@
export PGP_HOME=@PGP_HOME@
export PATH=@CCTOOLS_INSTALL_DIR@/bin:$PATH
export PYTHONPATH=$PGP_HOME:@PY_INSTALL_PREFIX@:$PYTHONPATH

export PGP_JAVA=@JAVA@

export PGP_PYTHON=@PYTHON@

if [ -f "$PGP_PYTHON" ]; then

    export PGP_PY_VER=`$PGP_PYTHON -c 'from distutils.sysconfig import *; print get_python_version()'`

    export PGP_PYCOMMON=$PGP_PYTHON_PREFIX/lib/python${PGP_PY_VER}/site-packages
    export PGP_PYDIST=$PGP_PYTHON_PREFIX/lib/python${PGP_PY_VER}/dist-packages
    export PYTHONPATH=${PGP_PYDIST}:${PGP_PYCOMMON}:${PYTHONPATH}
fi

export PGP_PEPNOVO_HOME=@PEPNOVO_INSTALL_DIR@

