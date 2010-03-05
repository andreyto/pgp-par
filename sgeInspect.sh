#!/bin/bash -x

exe=$BASH_SOURCE
exe_name=${exe##*/} # basename
exe_path=${exe%/*}  # dirname

set -o errexit

inspectIn=jobs/$SGE_TASK_ID.in

$exe_path/inspect -i $inspectIn -o ResultsX/$SGE_TASK_ID.txt -r $exe_path

rundir=$PWD
inspectOut=$rundir/ResultsX/$SGE_TASK_ID.txt 
rescoreOut=${inspectOut/.txt/.res}
PepNovoDir=/usr/local/depot/projects/PGP/productionPepNovo

mzxml=`awk -F, '/spectra/ {print $2}' $inspectIn`
specCnt=`echo $spectra | wc -l`

if [ 1 != $specCnt ]
then
    echo "Error, more then one input spectra in $inspectIn"
    exit 1
fi

cd $PepNovoDir
ulimit -c 0 # PepNovo cores have filled up the disk before

if ./PepNovo_bin -model CID_IT_TRYP -PTMs 'C+57' -file $mzxml -rescore_inspect $inspectOut $rescoreOut
then
    echo "PepNovo success"
else
    echo "PepNovo failure"
    rescoreOut=''
fi
    
bzip2 $inspectOut $rescoreOut

touch $rundir/Done/$SGE_TASK_ID
