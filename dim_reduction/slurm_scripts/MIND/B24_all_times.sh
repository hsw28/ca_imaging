#!/bin/bash
#SBATCH --account=p32072 ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal ### PARTITION (buyin, short, normal, w10001, etc)
#SBATCH --nodes=1 ## Never need to change this
#SBATCH --ntasks-per-node=32 ##
#SBATCH --time=15:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem=80G ## how much RAM you need per computer (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --output=matlab_job_result_B24_all.out ## standard out and standard error goes to this file


## job commands; simple is the MATLAB .m file, specified without the .m extension
module load matlab/r2018b
matlab -batch "B24_all_times"
