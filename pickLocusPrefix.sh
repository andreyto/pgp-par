#!/bin/bash
# get the most frequent prefix on the locus_tag, and round the count up to the next 1000
for i in NC_*.original.tbl
do
  num_prefix=(`grep locus_tag $i | awk '{print $2}' |\
   perl -lne 'print $1 if /^(\D+)/' |\
    sort | uniq -c | sort -k1n,1 | tail -n 1 |\
     awk '{print int(($1+1000)/1000)*1000,$2}'`)

  nc=${i%%.*} # strip everything after 1st .
  num=${num_prefix[0]}
  prefix=${num_prefix[1]}
  gff=$nc*.starts.gff
  rejectFile=$nc.reject.list
  reject=''
  if [ -s $rejectFile ]
  then
      reject="-x $rejectFile"
  fi
  setTbl="setTblStartsFromGFF.py"
  logFile="$nc.tbl.log"
  echo $setTbl -g $gff -t $i -l $prefix -s $num $reject -o $nc.new.tbl
  # &> $logFile
  if [ -s $logFile ]
  then
      echo $i had warnings
  fi
done
