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