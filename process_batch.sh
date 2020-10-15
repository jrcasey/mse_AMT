#!/bin/bash

# Run post-processing
cd ~/mse_AMT/
module load mit/matlab/2019a
matlab -nodesktop -nosplash -nojvm < toolbox/core/compile_AMT_subjobs.m

# Clean up
cd ~/
find . -name "matlab_crash*" -exec rm {} \;
rm ~/mse_AMT/*.out
rm ~/mse_AMT/*.err
# rm -r /nobackup1/jrcasey/*