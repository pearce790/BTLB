#!/bin/bash
#SBATCH --job-name PrimSim         # Set a name for your job.
#SBATCH --partition short          # Slurm partition to use
#SBATCH --ntasks 1                 # Number of tasks to run. By default, one CPU core per task
#SBATCH --time 0-01:00             # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=20000        # Memory limit for each tasks (in MB)
#SBATCH --array=1-180              # The number of iterations for this job
#SBATCH -o sim_%j.out      	   # File to which STDOUT will be written
#SBATCH -e sim_%j.err              # File to which STDERR will be written

module load R/4.2.2
Rscript PrimarySimulation.R $SLURM_ARRAY_TASK_ID