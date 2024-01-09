#!/bin/bash
#SBATCH --account=p32072 ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal
#SBATCH --nodes=1 ## Never need to change this
#SBATCH --ntasks-per-node=8 ## Never need to change this
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem=64G ## how much RAM do you need per computer (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name="shuff_\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=openfield_shuff.%A_%a.out ## use the jobid (A) and the specific job index (a) to name your log file
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=hsw@northwestern.edu  ## your email



module load matlab/r2018b

matlab -batch "addpath(genpath('/home/hsw967/Programming/ca_imaging')); addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/')); mutualinfo_openfield_shuff_4SLURM" -nodisplay -nosplash -nodesktop
