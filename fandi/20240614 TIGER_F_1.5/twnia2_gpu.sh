#!/bin/bash
#SBATCH -A GOV112050          # Allocation name
#SBATCH -J fd_00185           # Job name
#SBATCH -p gp4d               # Partition type [Take your needs by your case. See below as reference to set up.]
#SBATCH -n 2                  # Number of MPI tasks (i.e. processes)
#SBATCH -c 56                 # Number of cores per MPI task
#SBATCH -N 2                  # Total # of nodes
#SBATCH --gpus-per-node=1     # Number of GPUs per node
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fandi.ds@gmail.com  # Where to send mail.  Set this to your email address
#SBATCH -t 96:00:00           # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log       # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err        # (-e) Path to the standard error ouput


module purge
module load nvhpc
module list


### Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
### export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Launch MPI code
mpirun -np $SLURM_NTASKS ./run


###PARTITION   TIMELIMIT
### gtest      30:00
### gp1d       1-00:00:00
### gp2d*      2-00:00:00
### gp4d       4-00:00:00
