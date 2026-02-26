#!/bin/bash
#SBATCH -A ACD114141       
#SBATCH -J fd_00363     
#SBATCH -p ct56               # Partition type [Take your needs by your case. See below as reference to set up.]
#SBATCH --nodes=2             # (-N) Total # of nodes to be allocated
#SBATCH --ntasks-per-node=2   # Number of MPI proceses (one process for each NUMA node is recommended!!!)
#SBATCH -c 28                 # Number of cores per MPI task
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fandi.ds@gmail.com  # Where to send mail.  Set this to your email address
#SBATCH -t 4-00:00:00         # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log       # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err        # (-e) Path to the standard error ouput


###module purge
###module load intel/2023 intelmpi


### Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
### export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export UCX_NET_DEVICES=mlx5_0:1
#export OMP_PROC_BIND=false  

# Launch MPI code
mpiexec.hydra -np $SLURM_NTASKS -genv OMP_PROC_BIND=close ./run


###Partition type
### ctest --> max nodes(N):1  , 'Max walltime(t):    2:00:00'
### ct56  --> max nodes(N):1  , 'Max walltime(t): 4-00:00:00'
### ct224 --> max nodes(N):4  , 'Max walltime(t): 4-00:00:00'
### ct560 --> max nodes(N):10 , 'Max walltime(t): 4-00:00:00'
### ct2k  --> max nodes(N):40 , 'Max walltime(t): 3-00:00:00'
### ct8k  --> max nodes(N):150, 'Max walltime(t): 3-00:00:00'
### use "sinfo" command for updates