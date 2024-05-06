#!/bin/bash
#SBATCH --job-name sim             # Set a name for your job.
#SBATCH --partition short           # Slurm partition to use
#SBATCH --ntasks 1                  # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-08:00              # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=20000         # Memory limit for each tasks (in MB)
#SBATCH --array=1-100               # The number of iterations for this job
#SBATCH -o outputs/btlb_%j.out      # File to which STDOUT will be written
#SBATCH -e outputs/btlb_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL        # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mpp790@uw.edu # Email to which notifications will be sent

source /etc/profile.d/z00_lmod.sh
module load R
Rscript 20240429_ConferenceSimulation.R $SLURM_ARRAY_TASK_ID