#!/bin/bash
# The name of the job.
#PBS -N debug-fourstep
# Asking for 4 cores.
#PBS -l nodes=2
# Mail user if job aborts (a) or ends (e) (currently commented out).
# PBS -m ae
# Set the job to run 24 hours 10 minutes and 5 seconds.
#PBS -l walltime=0:10:0
# Set output file name.
#PBS -o test.txt
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
ulimit -c unlimited
module load gcc/4.9.1
module load openmpi-gcc/1.8.3
module load boost/1.54.0-gcc
module unload gcc/4.8.2
module unload gcc/4.7.2 
export LD_LIBRARY_PATH=clp/build/lib/:$LD_LIBRARY_PATH 
mpiexec ./main markshare_5_0.mps.gz
