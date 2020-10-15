# Connecting to the MIT server
ssh -F $HOME/eofe-cluster/linux/config eofe7.mit.edu

# Submitting a batch job array using SLURM
## For a specified interval (e.g., iterations 15-29)
sbatch --array=15-29 -p sched_mit_darwin2 --time=12:00:00 job2.sh

## Looping through intervals
sh Job_loop.sh

#### Cancel a job
scancel -u jrcasey

#### compress results for download
zip -r All_Solutions.zip /nobackup1/jrcasey/

#Downloading data from MIT server locally
## Copying a single file to a local directory
scp -i /Users/jrcasey/eofe-cluster/linux/eofe-key jrcasey@eofe7.mit.edu:~/mse_AMT/All_Solutions.zip ~/Documents/MATLAB/CBIOMES/Data/Environmental_Data/Cruises/AMT13/

## Copying all results files to a local directory
scp -i /Users/jrcasey/eofe-cluster/linux/eofe-key -r jrcasey@eofe7.mit.edu:/nobackup1/jrcasey/. ~/Documents/MATLAB/CBIOMES/Data/Environmental_Data/Cruises/AMT13/mse_Results/

# Clean up scratch and home dir
rm -r /nobackup1/jrcasey/*

### Clean up mse (while logged in)
rm ~/mse_AMT/*.out
rm ~/mse_AMT/*.err
find . -name "matlab_crash*" -exec rm {} \;

##### List too long?
cd ~/
find . -name "*.out" -print0 | xargs -0 rm
find . -name "*.err" -print0 | xargs -0 rm
cd ~/mse_AMT/
find . -name "matlab_crash*" -exec rm {} \;

#### How many files?
ls -f . | wc -l

#### List top n files from sorted list
ls -ltu | head -10
