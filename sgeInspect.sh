#!/bin/bash -x

exe=$BASH_SOURCE
exe_name=${exe##*/} # basename
exe_path=${exe%/*}  # dirname

set -o errexit

inspectIn=jobs/$SGE_TASK_ID.in

$exe_path/inspect -i $inspectIn -o ResultsX/$SGE_TASK_ID.txt -r $exe_path

rundir=$PWD
inspectOut=$rundir/ResultsX/$SGE_TASK_ID.txt 
rescoreOut=${inspectOut/.txt/.msgf}
PepNovoDir=/usr/local/depot/projects/PGP/productionPepNovo
MSGF=/usr/local/depot/projects/PGP/MSGF/MSGF.jar
java='/usr/local/java/1.6.0/bin/java -Xmx1000M'

mzxml=`awk -F, '/spectra/ {print $2}' $inspectIn`
specCnt=`echo $spectra | wc -l`

if [ 1 != $specCnt ]
then
    echo "Error, more then one input spectra in $inspectIn"
    exit 1
fi

if $java -jar $MSGF -inspect $inspectOut -d mzxml -x 0 > $rescoreOut
then
    echo "MSGF success"
else
    echo "MSGF failure"
    rescoreOut=''
fi
    
bzip2 $inspectOut $rescoreOut

touch $rundir/Done/$SGE_TASK_ID
