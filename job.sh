#!/bin/bash
#SBATCH --job-name="mse_MESO_SCOPE"
#SBATCH --partition=
#SBATCH --account=
#SBATCH --nodes=
#SBATCH --tasks-per-node=
#SBATCH --time=XX:XX:XX
#SBATCH --mem=
#SBATCH --mail-user=jrcasey@mit.edu
#SBATCH --mail-type=ALL
#SBATCH --workdir="path/to/mse"
#SBATCH --output=mse_MESO_SCOPE-%j.out
#SBATCH --error==mse_MESO_SCOPE-%j.error
#SBATCH --array 1-19584


. /path/to/modules
module load matlab-2019a 
WORKDIR=/path/to/mse
MATLAB=/path/to/matlab/bin
OPTIONS="-nodesktop -nosplash"

/bin/echo "Executing matlab script"
cd $WORKDIR && $MATLAB $OPTIONS < ./toolbox/MESO_SCOPE_WRAPPER.m &> Solution.result