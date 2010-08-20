#!/bin/bash
# Reverse sort the tbl's so largest is 1st, so we get the highest locus number 1st
origInName=`\ls -s NC_[0-9]*.original.tbl 2>/dev/null | sort -nr | awk '{print $2}'`
noOrig=`\ls -s NC_[0-9]*[0-9].tbl 2>/dev/null | sort -nr | awk '{print $2}'`

# Keep track of previous, in case the locus prefix is shared with the plasmids
prev='prev'
prevNum=1000
for i in $origInName $noOrig
do
  if [ ! -s $i ]
  then
      continue
  fi
# get the most frequent prefix on the locus_tag, and round the count up to the next 1000
  num_prefix=(`grep locus_tag $i | awk '{print $2}' |\
   perl -lne 'print $1 if /^(\D+)/' |\
    sort | uniq -c | sort -k1n,1 | tail -n 1 |\
     awk '{print int(($1+1000)/1000)*1000,$2}'`)

  nc=${i%%.*} # strip everything after 1st .
  num=${num_prefix[0]}
  prefix=${num_prefix[1]}
  # seen this locus before, use a higher number to avoid duplicates
  if [ $prefix = $prev ]
  then
      num=$(( prevNum + 1000 ))
  fi
  gff=$nc*.starts.gff
  if [ ! -s $gff ]
  then
      echo No File $gff
      continue
  fi
  prev=$prefix
  prevNum=$num
  rejectFile=$nc.reject.list
  reject=''
  if [ -s $rejectFile ]
  then
      reject="-x $rejectFile"
  fi
  setTbl="setTblStartsFromGFF.py"
  logFile="$nc.tbl.log"
  echo -g $gff -t $i -l $prefix -s $num $reject -o $nc.new.tbl
  # &> $logFile
  if [ -s $logFile ]
  then
      echo $i had warnings
  fi
done
