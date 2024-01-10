#!/bin/bash
#SBATCH --account=p32072 ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=short
#SBATCH --nodes=1 ## Never need to change this
#SBATCH --ntasks-per-node=52 ## Never need to change this
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem=16G ## how much RAM do you need per computer (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name="shuff_\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=openfield_shuff.%A_%a.out ## use the jobid (A) and the specific job index (a) to name your log file
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=hsw@northwestern.edu  ## your email

module load matlab/r2018b

matlab -batch mutualinfo_openfield_shuff_4SLURM" -nodisplay -nosplash -nodesktop
