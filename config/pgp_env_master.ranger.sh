#!/bin/bash
#This sets up the system environment needed to
#both build and run the PGP pipeline on TACC Ranger.
#Source (not execute) this file on your login node.
#This SHOULD NOT be sourced on compute nodes - you
#will get errors. Instead, your qsub script should use -V
#option to propagate the environment variables from login
#to compute nodes. It is best not to source this file
#from your Bash startup profiles - you might end up
#sourcing it on compute nodes with your MPI jobs not 
#running or even worse, getting stuck eating your allocation
#time.

if [ -z "PGP_ENV_MASTER_ENTERED" ]; then
export PGP_ENV_MASTER_ENTERED=1

#we need gcc compiler toolchain as the only one that has
#a working combination of boost and python builds
#associated with it
module swap pgi gcc
module unload openmpi
module load mvapich
#the latest boost/1.48 was compiled agains python with wrong unicode size
module load boost/1.46.1
module load python
#cctools make fail with default Globus 5
module add globus/4.0.8
#1.6 that PGP MSGF needs is default java here
module add java

fi

