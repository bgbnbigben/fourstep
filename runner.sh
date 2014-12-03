filename="30_70_45_095_100.mps.gz a1c1s1.mps.gz aflow40b.mps.gz berlin_5_8_0.mps.gz core2536-691.mps.gz core4872-1529.mps.gz csched010.mps.gz neos-506422.mps.gz neos-847302.mps.gz"

filename="alan.bar batch.bar batchdes.bar csched1a.bar csched1.bar csched2a.bar csched2.bar elf.bar eniplac.bar enpro48.bar enpro48pb.bar enpro56.bar enpro56pb.bar ex1221.bar"

for file in $filename; do
    mkdir -p ${file}/rtrees;
    cp main ${file}
#    cp runjob.sh ${file%.gz}
    cp baron/$file ${file}
    cd ${file} && cat <<EOS | qsub -
#!/bin/bash
# The name of the job.
#PBS -N $file
#PBS -l nodes=64
# Mail user if job aborts (a) or ends (e) (currently commented out).
#PBS -m ae
#PBS -l walltime=2:00:0
# Set output file name.
#PBS -o $file.txt
cd \$PBS_O_WORKDIR
ulimit -c unlimited
ulimit -s unlimited
ulimit -f unlimited
module load gcc/4.9.1
module load openmpi-gcc/1.8.3
module load boost/1.54.0-gcc
module unload gcc/4.8.2
module unload gcc/4.7.2 
export PBSCOREDUMP=1
export LD_LIBRARY_PATH=/home/bgbnbigben/project/clp/build/lib/:$LD_LIBRARY_PATH 
mpiexec /home/bgbnbigben/project/${file}/main $file
EOS
cd /home/bgbnbigben/project/
done;
