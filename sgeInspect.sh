#!/bin/bash -x

exe=$BASH_SOURCE
exe_name=${exe##*/} # basename
exe_path=${exe%/*}  # dirname

set -o errexit

inspectIn=jobs/$SGE_TASK_ID.in

$exe_path/inspect -i $inspectIn -o ResultsX/$SGE_TASK_ID.txt -r $exe_path

rundir=$PWD
results=$rundir/ResultsX
filter=$results/filter
inspectOut=$results/$SGE_TASK_ID.txt
pepnovoOut=${inspectOut/.txt/.res}
pvalueIn=$results/$SGE_TASK_ID.txt
pvalueOut=${inspectOut/.txt/.pval}
msgfOut=${inspectOut/.txt/.msgf}
filterOut=${msgfOut/ResultsX/ResultsX/filter}

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

cd $PepNovoDir
ulimit -c 0 # PepNovo cores have filled up the disk before

if ./PepNovo_bin -model CID_IT_TRYP -PTMs 'C+57' -file $mzxml -rescore_inspect $inspectOut $pepnovoOut
then
    echo "PepNovo success"
    pvalueIn=$pepnovoOut
else
    echo "PepNovo failure"
    pepnovoOut=''
fi

# PValue the inspect or pepnovo results
$exe_path/PValue.py -r $pvalueIn -w $pvalueOut -p 0.1 -S 0.5 -1 -H &> pvalue.log 

if $java -jar $MSGF -inspect $pvalueOut -d mzxml -x 0 > $msgfOut
then
    echo "MSGF success"
    if [ ! -d $filter ]
    then
        mkdir $filter
    fi
    $exe_path/FilterInspectMSGF.py -p 1e-6 -r $msgfOut -o $filterOut
else
    echo "MSGF failure"
    rescoreOut=''
    filterOut=''
fi
    
bzip2 $inspectOut $pepnovoOut $pvalueOut $msgfOut $filterOut

touch $rundir/Done/$SGE_TASK_ID
