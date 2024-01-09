#!/bin/bash
#SBATCH --account=p32072
#SBATCH --partition=normal
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1 ## number of jobs to run "in parallel"
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --job-name="sample_job_\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=MI_trace_shuff.out ## use the jobid (A) and the specific job index (a) to name your log file
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=hsw@northwestern.edu  ## your email
#SBATCH --cpus-per-task=8

# Load MATLAB module
module load matlab/r2018b

# Run MATLAB script
matlab -nodisplay -nosplash -nodesktop -r   "addpath('/path/to/your/directory');
                                            SLURM_mutual_info_shuff; exit;"
