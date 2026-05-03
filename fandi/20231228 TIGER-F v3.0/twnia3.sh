#!/bin/bash
#SBATCH -A GOV112050       
#SBATCH -J fd_00188    
#SBATCH -p ct224               # Partition type [Take your needs by your case. See below as reference to set up.]
#SBATCH -n 4                  # Number of MPI tasks (i.e. processes)
#SBATCH -c 56                  # Number of cores per MPI task
#SBATCH -N 4                  # Maximum number of nodes to be allocated
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fandi.ds@gmail.com  # Where to send mail.  Set this to your email address
#SBATCH -t 96:00:00           # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log       # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err        # (-e) Path to the standard error ouput


module purge
module load compiler/intel/2021   IntelMPI/2021


### Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
### export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpiexec.hydra ./run

###Partition type
### ctest --> max nodes(N):1  , 'Max walltime(t): 00:30:00'
### ct56  --> max nodes(N):1  , 'Max walltime(t): 96:00:00'
### ct224 --> max nodes(N):4  , 'Max walltime(t): 96:00:00'
### ct560 --> max nodes(N):10 , 'Max walltime(t): 96:00:00'
### ct2k  --> max nodes(N):40 , 'Max walltime(t): 72:00:00'
### ct8k  --> max nodes(N):150, 'Max walltime(t): 72:00:00'
