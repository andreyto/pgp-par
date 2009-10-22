#!/bin/bash

exe=$BASH_SOURCE
exe_name=${exe##*/} # basename
exe_path=${exe%/*}  # dirname

$exe_path/inspect -i $SGE_TASK_ID.in -o $SGE_TASK_ID.txt -r $exe_path
