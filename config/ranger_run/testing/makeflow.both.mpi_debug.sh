#!/bin/bash
#testing locally w/o qsub
export MPI_QUEUE_PORT=10000
mpirun -np 3 -all-local -hostfile hostfile /work/01241/atovchig/packages/bin/mpi_queue_worker localhost $MPI_QUEUE_PORT &
makeflow -T mpi-queue -J 1000 pgp_proc.mkf
wait

