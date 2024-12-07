#!/bin/bash
#SBATCH --job-name ConfSim      # Set a name for your job.
#SBATCH --partition short       # Slurm partition to use
#SBATCH --ntasks 1              # Number of tasks to run. By default, one CPU core per task
#SBATCH --time 0-02:00          # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=20000     # Memory limit for each tasks (in MB)
#SBATCH --array=1-100           # The number of iterations for this job
#SBATCH -o conf_%j.out      	# File to which STDOUT will be written
#SBATCH -e conf_%j.err      	# File to which STDERR will be written

module load R/4.2.2
Rscript ConferenceSimulation.R $SLURM_ARRAY_TASK_ID