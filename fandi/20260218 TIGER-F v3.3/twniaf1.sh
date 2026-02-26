#!/bin/bash
#SBATCH -A ACD114141          # (-A) iService Project ID       
#SBATCH -J fd_00363           # (-J) Job name    
#SBATCH -p ct112              # (-p) Slurm partition [Take your needs by your case. See below as reference to set up.]
#SBATCH --nodes=1             # (-N) Total # of nodes to be allocated
#SBATCH --ntasks-per-node=8   # Number of MPI proceses (one process for each NUMA node is recommended!!!)
#SBATCH -c 14                 # Number of cores per MPI task
#SBATCH --time=4-00:00:00     # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log       # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err        # (-e) Path to the standard error ouput
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fandi.ds@gmail.com  # Where to send mail.  Set this to your email address


###module purge
###module load intel/2023_2


### Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
### export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export UCX_NET_DEVICES=mlx5_0:1
#export OMP_PROC_BIND=false  

# Launch MPI code
mpiexec.hydra -np $SLURM_NTASKS -genv OMP_PROC_BIND=close ./run


###Partition type
### ct112 --> max nodes(N):1 , 'Max walltime(t): 4-00:00:00'
### ct448 --> max nodes(N):4 , 'Max walltime(t): 4-00:00:00'
### ct1k  --> max nodes(N):10, 'Max walltime(t): 2-16:00:00'
### ct2k  --> max nodes(N):20, 'Max walltime(t): 2-00:00:00'
### ct4k  --> max nodes(N):40, 'Max walltime(t): 2-00:00:00'
### ct8k  --> max nodes(N):80, 'Max walltime(t): 2-00:00:00'
### use "sinfo" command for updates