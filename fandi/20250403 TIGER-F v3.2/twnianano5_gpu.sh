#!/bin/bash
#SBATCH -A GOV113064          # (-A) iService Project ID
#SBATCH -J fd_00363           # (-J) Job name
#SBATCH -p normal             # (-p) Slurm partition [Take your needs by your case. See below as reference to set up.]
#SBATCH --nodes=1             # (-N) Total # of nodes to be allocated
#SBATCH --gpus-per-node=1     # Number of GPUs per node (must be the same as ntasks-per-node)
#SBATCH --ntasks-per-node=1   # Number of MPI proceses (must be the same as gpus-per-node)
#SBATCH --cpus-per-task=1     # (-c) Number of cores per MPI task (limited to 4 cores/GPU)
#SBATCH --time=2-00:00:00     # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log       # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err        # (-e) Path to the standard error ouput
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fandi.ds@gmail.com  # Where to send mail.  Set this to your email address




module list
nvidia-smi


### Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
### export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export UCX_NET_DEVICES=mlx5_0:1
export UCX_IB_GPU_DIRECT_RDMA=1
export OMP_PROC_BIND=false       

# Launch MPI code
mpirun --bind-to none -np $SLURM_NTASKS ./run


###PARTITION   TIMELIMIT
### gtest      30:00
### gp1d       1-00:00:00
### gp2d*      2-00:00:00
### gp4d       4-00:00:00

