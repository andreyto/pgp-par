#!/bin/bash
# Proteogenomics pipeline (PGP) wrapper script

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   Copyright J. Craig Venter Institute 2012
#   See docs.PGP/COPYING file distributed along with the proteogenomics 
#   package for the copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#set -ex
set -e
#usage message
usage()
{
cat << EOF
Usage: $0 options

This script executes the first (installs the package) of the three seequential steps in the Proteogenomics pipeline (PGP) in TACC Ranger or other grid environments. The instructions on how to run the following 2 steps are given in the PGPdoc.pdf provided with this script. All options except for -h are mandatory.

REQUIREMENTS: 1). git utility for downloading proteogenomics package, 2). BASH Unix shell environment, 3). Installed Python, JAVA & Perl languages.

OPTIONS:
        -h      help message    
        -e      computing environment name (= ranger is the currently supported name)        
        -r      installation directory 
EOF
}

#parse options and their arguments
export ENVIRONMENT=
export ROOT=
while getopts "he:r:" OPTION
do
    case $OPTION in
        h)
           usage
           exit 1
           ;;
        e)
           ENVIRONMENT=$OPTARG
           ;;
        r)
           ROOT="$OPTARG"
           ;;
        ?)
          #printf "Invalid option: %s\n\n" -$OPTARG #uncomment this if the error mode is SILENT
           usage
           exit
           ;;
    esac
done

#exit if any of the options have missing arguments
if [[ -z $ENVIRONMENT ]] || [[ -z $ROOT ]]
then
    usage
    exit 1
fi

ROOT=$(cd $ROOT && pwd)


#log
export DATETIME=`date '+Date: %F Time: %T'`;

echo "-- runPGP started: ${DATETIME} --";
echo;

#Step1: Install the proteogenomics package & setup variables
##common to all environments
echo " Started installing proteogenomics package ..."
echo "Creating directory for the installation ..."
export PGP_ROOT="$ROOT"
mkdir -p $PGP_ROOT
cd $PGP_ROOT
#checkout the source with 'git'- it will create 'proteogenomics' directory
git clone git@bitbucket.org:andreyto/proteogenomics.git
mkdir -p build
pushd build
#get and unpack the dependencies
wget https://bitbucket.org/andreyto/proteogenomics/downloads/vendor.tar.gz
tar -zxf vendor.tar.gz
popd #back to PGP_ROOT
config_dir=$PGP_ROOT/proteogenomics/config #set config_dir variable for the location of shell scripts

#setup run-environment
cp $config_dir/pgp_makeflow.$ENVIRONMENT.ini $config_dir/pgp_makeflow.ini
perl -i  -pe 's/(pgp_root=|PGP_ROOT=).*/$1$ENV{PGP_ROOT}/g' $config_dir/pgp_makeflow.ini #change PGP_ROOT 
cp $config_dir/pgp_makeflow_env.$ENVIRONMENT.sh $config_dir/pgp_makeflow_env.sh
perl -i  -pe 's/(pgp_root=|PGP_ROOT=).*/$1$ENV{PGP_ROOT}/g' $config_dir/pgp_makeflow_env.sh #change PGP_ROOT 
cp $config_dir/pgp_login_env.$ENVIRONMENT.sh $config_dir/pgp_login_env.sh # no PGP_ROOT here
source $config_dir/pgp_makeflow_env_master.sh #setup ranger-specific run-environment
#Build & install dependencies: for the environment
(cd $PGP_HOME && make) || (echo "Failed to build Inspect binary!" && exit 1)
export PGP_VENDOR_BUILD=$PGP_ROOT/build/vendor
cd $PGP_VENDOR_BUILD
$config_dir/install_makeflow
#unpack seqan and build the extension
PGP_SEQAN_VER=1.3.1
unzip seqan-$PGP_SEQAN_VER.zip
export PGP_SEQAN_HOME=$PGP_VENDOR_BUILD/seqan-$PGP_SEQAN_VER
(cd $PGP_HOME/seqan && ./build_ranger.sh)
# copy IndexSearch.so into 'proteogenomics' directory
cp $PGP_HOME/seqan/IndexSearch.so $PGP_HOME/IndexSearch.so
#install pepnovo
$config_dir/install_pepnovo
#no build for MSGF, source is not available, just unpack the jar
(cd $PGP_VENDOR_HOME && tar -zxf $PGP_VENDOR_BUILD/MSGF.tar.gz)
#install python vendor (the python Bio module)
$config_dir/install_python_vendor

echo "BUILD IS COMPLETE!!!"
echo
echo "  Prepare data for running the analysis in Step 3 ..."
echo
#Instructions for Step 2 in the pipeline
echo "  Step 2: Generate the Makeflow input for the pipeline in a qsub script. A template (prepPGPdata.$ENVIRONMENT.qsub) is provided in $PGP_HOME/config directory ..." 
echo
echo "  REQUIRED steps for running 'prepPGPdata.$ENVIRONMENT.qsub':"
echo
echo "  1). source $config_dir/pgp_makeflow_env_master.sh ..."
echo
echo "  2). modify MANDATORY SGE (qsub) options as necessary (particularly, -A, -pe  & -l h_rt) in $config_dir/prepPGPdata.$ENVIRONMENT.qsub script ..." 
echo
echo "  3). provide values (paths) to 'DB' (input nucleotide/potein directory) and 'OUT' (output directory where data is saved) variables in $config_dir/prepPGPdata.$ENVIRONMENT.qsub script ..."
echo
echo "  4). cd to $PGP_ROOT (optional) and run $config_dir/prepPGPdata.$ENVIRONMENT.qsub ..."
echo
echo "  installPGP.sh finished: ${DATETIME} --"
echo
