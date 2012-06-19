#!/bin/bash
#testing locally w/o qsub
export MPI_QUEUE_PORT=10000
makeflow -T mpi-queue -J 1000 pgp_proc.mkf

