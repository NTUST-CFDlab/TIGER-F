#!/bin/bash
#PBS -P GOV112050 
#PBS -N fd_00184
#PBS -l select=30:ncpus=32:mpiprocs=1:ompthreads=32
#PBS -l walltime=48:00:00	
#PBS -q cf1200
#PBS -o jobreport.out
#PBS -e jobreport.err
#PBS -M fandi.ds@gmail.com
#PBS -m be 


module purge
module load intel/2018_u1
module list
cd $PBS_O_WORKDIR


cat $PBS_NODEFILE
echo $PBS_O_WORKDIR
date

mpirun -PSM2 ./run


#qstat -xs

### select= nodes:CPUs/node:process/node:cores/process

###Partition type
### (cf : Memory 384 GB/node)
### serial -->    1 node ,        1 CPU, 'Max walltime(t): 96:00:00'
### cf40   -->    1 node ,     2-40 CPU, 'Max walltime(t): 96:00:00'
### cf160  -->  1-4 nodes,    2-160 CPU, 'Max walltime(t): 96:00:00'
### cf1200 --> 5-30 nodes, 161-1.2k CPU, 'Max walltime(t): 48:00:00'
### (ct : Memory 192 GB/node)
### ct160  -->  1-4 nodes,    2-160 CPU, 'Max walltime(t): 96:00:00'
### ct400  --> 5-10 nodes,  161-400 CPU, 'Max walltime(t): 96:00:00'
### ct800 --> 11-20 nodes,  401-800 CPU, 'Max walltime(t): 72:00:00'
### ct2k  --> 21-50 nodes,   801-2k CPU, 'Max walltime(t): 48:00:00'
### ct6k --> 51-150 nodes,  2001-6k CPU, 'Max walltime(t): 24:00:00'
