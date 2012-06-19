#!/bin/bash
#testing locally w/o qsub
export MPI_QUEUE_PORT=10000
mpirun -np 3 -all-local -hostfile hostfile /work/01241/atovchig/packages/bin/mpi_queue_worker login3 $MPI_QUEUE_PORT
