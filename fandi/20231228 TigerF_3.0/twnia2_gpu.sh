#!/bin/bash
#SBATCH -A GOV112050          # Allocation name
#SBATCH -J fd_00330           # Job name
#SBATCH -p gp1d               # Partition type [Take your needs by your case. See below as reference to set up.]
#SBATCH --nodes=2             # Total # of nodes
#SBATCH --ntasks-per-node=1   # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=4     # Number of cores per MPI task
#SBATCH --gpus-per-node=1     # Number of GPUs per node (must be the same as ntasks-per-node)
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fandi.ds@gmail.com  # Where to send mail.  Set this to your email address
#SBATCH -t 1-00:00:00         # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log       # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err        # (-e) Path to the standard error ouput



module purge
module load nvhpc
module list
nvidia-smi


### Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
### export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export UCX_NET_DEVICES=mlx5_2:1
export OMP_PROC_BIND=false       

# Launch MPI code
mpirun --bind-to none -np $SLURM_NTASKS ./run


###PARTITION   TIMELIMIT
### gtest      30:00
### gp1d       1-00:00:00
### gp2d*      2-00:00:00
### gp4d       4-00:00:00

