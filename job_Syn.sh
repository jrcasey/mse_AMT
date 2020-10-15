#!/bin/bash
#
#SBATCH -J jrcaseyjarray # A single job name for the array
#SBATCH -o jrcaseyjarray_%A_%a.out # Standard output
#SBATCH -e jrcaseyjarray_%A_%a.err # Standard error

cd /home/jrcasey/mse_AMT/

echo ${SLURM_ARRAY_TASK_ID}
export SLURM_ARRAY_TASK_ID
module load mit/matlab/2019a
matlab -nodesktop -nosplash -nojvm < toolbox/wrappers/Synthetic_Simulation_Wrapper.m
#matlab -nodesktop -nosplash -nojvm < AMT_Wrapper2(${SLURM_ARRAY_TASK_ID})

sleep 3
