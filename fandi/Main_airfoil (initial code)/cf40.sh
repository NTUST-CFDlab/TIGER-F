#!/bin/bash
#PBS -P MST108362
#PBS -N Dy_st1
#PBS -l select=1:ncpus=40:mpiprocs=5:ompthreads=8
#PBS -l walltime=96:00:00
#PBS -q cf40
#PBS -j oe
#PBS -M gsn9409@gmail.com 
#PBS -m be 



module purge
module load intel/2018_u1
module list
cd $PBS_O_WORKDIR


cat $PBS_NODEFILE
echo $PBS_O_WORKDIR
date

mpirun -PSM2 ./gogo


#qstat -xs
