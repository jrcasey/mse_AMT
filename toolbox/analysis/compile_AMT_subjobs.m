%% Post processing on server
cd ~/mse_AMT/
addpath(genpath('~/mse_AMT'))
ResultsDirectory = '/nobackup1/jrcasey/';
load('data/output/Gridding.mat');
load('data/output/FileNames.mat');
load('data/output/CruiseData.mat');
load('data/output/PanGEM.mat');

[FullSolution] = get_AMT_Results_server(ResultsDirectory,PanGEM,FileNames,Gridding,CruiseData);

save('~/mse_AMT/data/output/FullSolution.mat','FullSolution');

%% Post processing on local machine for AMT cruise
cd /Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/
addpath(genpath('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/'))
%ResultsDirectory = '/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/Solution_20200408/';
ResultsDirectory = '/Users/jrcasey/Documents/MATLAB/CBIOMES/Data/Environmental_Data/Cruises/AMT13/Solution_20201015/';

load('data/output/Gridding.mat');
load('data/output/FileNames.mat');
load('data/output/CruiseData.mat');
load('data/GEM/PanGEM.mat');


[FullSolution] = get_AMT_Results_server(ResultsDirectory,PanGEM,FileNames,Gridding,CruiseData);

save('data/output/FullSolution.mat','FullSolution');

%% Post processing on local machine for Synthetic data
cd /Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/
addpath(genpath('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/'))
%ResultsDirectory = '/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/Solution_20200408/';
ResultsDirectory = '/Users/jrcasey/Documents/MATLAB/CBIOMES/Data/Environmental_Data/Cruises/AMT13/Solution_20200428/';

load('data/output/Synthetic/Gridding.mat');
load('data/output/Synthetic/FileNames.mat');
load('data/output/Synthetic/CruiseData.mat');
load('data/output/Synthetic/PanGEM.mat');


[FullSolution] = get_AMT_Results_server(ResultsDirectory,PanGEM,FileNames,Gridding,CruiseData);

save('data/output/Synthetic/FullSolution.mat','FullSolution');