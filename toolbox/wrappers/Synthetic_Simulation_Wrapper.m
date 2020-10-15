%% Synthetic Simulation Wrapper
% Runs and analyzes the MSE simulation for the AMT13 cruise.
% Set up currently to run on Centos7
%% set root dir (for server only)
%cd Documents/MATLAB/GitHub/mse_AMT/
%addpath '/nfs/cnhlab001/cnh/projects/jrcasey/mse/test/mosek/9.1/toolbox/r2015a'	
addpath(genpath('~/mse_AMT/'));	
cd ~/mse_AMT/	
%% Version
version = strcat({'_v'},datestr(date,'yyyymmdd'))

%% Version control and database paths (save a copy of this to the server)
FileNames = struct;
% Current PanGEM
FileNames.PanGEM_Path = 'data/GEM/PanGEM.mat';
% Organism database
FileNames.orgDB_Path = 'data/db/orgDatabase.csv';
% Current strain models
FileNames.StrMod_Path = 'data/GEM/targStrMod3.mat';
% List of strains to be analyzed
FileNames.strainList_Path = 'data/db/strainList.mat';
% Current OGTDat
FileNames.OGTDat_Path = 'data/db/OGTDat.csv';
% Cruise data
FileNames.CruiseDB_filename = 'data/envData/Synthetic_Cruise.csv';
% HyperPro profiles
FileNames.IrrDat_fileName = 'data/envData/IrrDat.mat';
% TpDat path
FileNames.TpDat_fileName = 'data/db/TpDat.csv';
% PhysOpt and PigOpt constraints path
FileNames.PhysOptPigOptConstraints_fileName = 'data/db/BOFConstraints.csv';
% PigDB path
FileNames.PigDB_fileName = 'data/db/AbsorptionDatabase.csv';
% Destination for solution
FileNames.destination_fileName = strcat('nobackups/jrcasey/Solution',version,'.mat');

%% Gridding (save a copy of this to the server)
Gridding = struct;
% Stations
Gridding.stationsVec =[{'PO4'},{'NH3'},{'NO2'},{'NO3'},{'Light'}];
Gridding.stationsVec2 = [1 2 3 4 5];
Gridding.nStations = numel(Gridding.stationsVec);
% Depth (m)
Gridding.minZ = 10;
Gridding.maxZ = 300;
Gridding.intervalZ = 10;
Gridding.depthVec = Gridding.minZ:Gridding.intervalZ:Gridding.maxZ;
Gridding.nZ = numel(Gridding.depthVec);
% Wavelength (nm)
Gridding.minLambda = 400; % nm
Gridding.maxLambda = 686; %nm
Gridding.bandwidth = 2; % nm
Gridding.lambdaVec = Gridding.minLambda:Gridding.bandwidth:Gridding.maxLambda; % nm wavelength

%% Options
Options = struct;
% Set options
Options.saveSolution = false;
Options.saveStrMod = false;
Options.maxIter_TpOpt = 1000;
Options.maxIter_physOpt = 1000;
%Options.saveCruiseData = true;

%% Get strains to analyze
load(FileNames.strainList_Path);
Gridding.strNameVec = strList;
Gridding.nStr = numel(Gridding.strNameVec);

%% Find index from SLURM
% preallocate a 3d matrix of dimensions nStations, nZ, nStr
idxMat = zeros(Gridding.nStations,Gridding.nZ,Gridding.nStr);
nIterations = size(idxMat,1).*size(idxMat,2).*size(idxMat,3);
% get job array index
% job_array_idx = str2num(getenv('SLURM_ARRAY_TASK_ID'));
% avoiding maxjobid on slurm (comment out for normal runs)
job_array_idx = str2num(getenv('SLURM_ARRAY_TASK_ID')) + 5760;
% % avoiding maxjobid on slurm (comment out for normal runs)
% job_array_idx_temp = str2num(getenv('SLURM_ARRAY_TASK_ID'));
% load('data/output/missingFileNo.mat');
% job_array_idx = missingFileNo(job_array_idx_temp);
% locate coordinates
[i,j,k] = ind2sub(size(idxMat),job_array_idx);

%% Run simulation
station = Gridding.stationsVec{i};
depth = Gridding.depthVec(j);
strName = Gridding.strNameVec{k};

%Solution = struct;
tic;
[Solution] = Synthetic_Simulation(strName, station, depth, Gridding, FileNames, Options);
dt = toc;
Solution.runtime = dt;
Solution.strName = strName;
Solution.z = depth;
Solution.station = station;

if ~Options.saveStrMod
    Solution.StrMod = [];
end

%% Save solution
%save(cell2str(FileNames.destination_fileName),'FullSolution');
save(strcat('/nobackup1/jrcasey/','Solution_',num2str(job_array_idx),'.mat'),'Solution');

% change path to nobackup/jrcasey/ for saving the full output.