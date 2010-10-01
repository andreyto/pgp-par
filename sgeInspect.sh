#!/bin/bash -x

exe=$BASH_SOURCE
exe_name=${exe##*/} # basename
exe_path=${exe%/*}  # dirname

set -o errexit

inspectIn=jobs/$SGE_TASK_ID.in

$exe_path/inspect -i $inspectIn -o ResultsX/$SGE_TASK_ID.txt -r $exe_path

rundir=$PWD
results=$rundir/ResultsX
inspectOut=$results/$SGE_TASK_ID.txt
pepnovoOut=${inspectOut/.txt//pepnovo/.res}
pepOutDir=$results/pepnovo
PepNovoDir=/usr/local/depot/projects/PGP/productionPepNovo

mzxml=`awk -F, '/spectra/ {print $2}' $inspectIn`
specCnt=`echo $spectra | wc -l`

if [ 1 != $specCnt ]
then
    echo "Error, more then one input spectra in $inspectIn"
    exit 1
fi
if [ ! -d $pepOutDir ]
then 
    mkdir $pepOutDir
fi

cd $PepNovoDir
ulimit -c 0 # PepNovo cores have filled up the disk before

if ./PepNovo_bin -model CID_IT_TRYP -PTMs 'C+57' -file $mzxml -rescore_inspect $inspectOut $pepnovoOut
then
    echo "PepNovo success"
else
    echo "PepNovo failure, copying inspect output into $pepOutDir"
    cp -f $inspectOut $pepnovoOut
fi

bzip2 $inspectOut $pepnovoOut 
