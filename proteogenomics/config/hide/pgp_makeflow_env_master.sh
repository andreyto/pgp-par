#!/bin/bash

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   Copyright J. Craig Venter Institute 2012
#   See docs.PGP/COPYING file distributed along with the proteogenomics 
#   package for the copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#source this file on the master (job submission) node before you submit the jobs
#get the dir of this file - works on BASH > 3.
config_dir=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
master_env=$config_dir/pgp_login_env.sh
#source pgp_env_master.sh  and pgp_makeflow_env.sh
if [ -f $master_env ]; then
        source $master_env;
        if [ $? -eq 0 ]; then
                source $config_dir/pgp_makeflow_env.sh
        else
                echo "$master_env was not sourced successfully ..."
        fi
else
        echo "$master_env does'nt exist"
fi
