#!/bin/bash
# The name of the job.
#PBS -N csched2.bar
#PBS -l nodes=64
# Mail user if job aborts (a) or ends (e) (currently commented out).
#PBS -m ae
#PBS -l walltime=2:00:0
# Set output file name.
#PBS -o csched2.bar.txt
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
ulimit -c unlimited
ulimit -s unlimited
module load gcc/4.9.1
module load openmpi-gcc/1.8.3
module load boost/1.54.0-gcc
module unload gcc/4.8.2
module unload gcc/4.7.2 
export PBSCOREDUMP=1
export LD_LIBRARY_PATH=/home/bgbnbigben/project/clp/build/lib/:$LD_LIBRARY_PATH 
mpiexec ./main csched2.bar
