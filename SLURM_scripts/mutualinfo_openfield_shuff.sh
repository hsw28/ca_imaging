#!/bin/bash
#SBATCH --account=p32072
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --job-name="sample_job_\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=MI_shuff.out ## use the jobid (A) and the specific job index (a) to name your log file
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=hsw@northwestern.edu  ## your email
#SBATCH --cpus-per-task=8

# Load MATLAB module
module load matlab/r2018b

# Run MATLAB script

matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('/home/hsw967/Programming/ca_imaging')); addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/')); mutualinfo_openfield_shuff_4SLURM; exit;"
