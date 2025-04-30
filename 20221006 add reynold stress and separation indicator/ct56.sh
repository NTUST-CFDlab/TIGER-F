#!/bin/bash
#SBATCH -A GOV111056       
#SBATCH -J fd_00047    
#SBATCH -p ct56                     # Partiotion name [Take your needs by your case, take the line 7 as reference to set up.]
#SBATCH -n 8                        # Number of MPI tasks (i.e. processes)
#SBATCH -c 7                        # Number of cores per MPI task
#SBATCH -N 1                        # Maximum number of nodes to be allocated [ctest can setup the maxinmum nodes number:20 'Maximum walltime: 30 minutes' , ct56 can setup the maxinmum nodes number:1 'Maximum walltime: 96 hours', ct224 can setup the maxinmum nodes number:4 'Maximum walltime: 96 hours' ]
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fandi.ds@gmail.com  # Where to send mail.  Set this to your email address
#SBATCH -t 96:00:00                    # (-t) Wall time limit (days-hrs:min:sec)
#SBATCH --output=%j.log                # (-o) Path to the standard output and error files relative to the working directory
#SBATCH --error=%j.err                 # (-e) Path to the standard error ouput


# Set OMP_NUM_THREADS to the number of CPUs per task we asked for.
# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpiexec.hydra ./run
