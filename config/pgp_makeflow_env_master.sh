#!/bin/bash
#source this file on the master (job submission) node before you
#the software or submit jobs
#get the dir of this file - works on BASH > 3.
config_dir=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
master_env=$config_dir/pgp_env_master.sh
#Ignore errors from sourcing the master host environment
#[ -f $master_env ] && (. $master_env || true)
. $config_dir/pgp_makeflow_env.sh

