%% MESO-SCOPE Simulation
% Runs and analyzes the MSE simulation for the MESO-SCOPE cruise.
% Set up currently to run on Centos7
%% set root dir (for server only)
%cd Documents/MATLAB/GitHub/mse/
%addpath '/nfs/cnhlab001/cnh/projects/jrcasey/mse/test/mosek/9.1/toolbox/r2015a'	
addpath(genpath('~/'));	
cd ~/mse/	
%% Version
version = strcat({'_v'},datestr(date,'yyyymmdd'))

%% Version control and database paths
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
FileNames.CruiseDB_filename = 'data/envData/MESO-SCOPE.csv';
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

%% Gridding
Gridding = struct;
% Stations
Gridding.stationsVec = 4:15;
Gridding.nStations = numel(Gridding.stationsVec);
% Depth (m)
Gridding.minZ = 15;
Gridding.maxZ = 181;
Gridding.intervalZ = 5;
Gridding.depthVec = Gridding.minZ:Gridding.intervalZ:Gridding.maxZ;
Gridding.nZ = numel(Gridding.depthVec);
% Wavelength (nm)
Gridding.minLambda = 300; % nm
Gridding.maxLambda = 850; %nm
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
job_array_idx = str2num(getenv('SLURM_ARRAY_TASK_ID'));
% locate coordinates
[i,j,k] = ind2sub(size(idxMat),job_array_idx);

%% Run simulation
station = Gridding.stationsVec(i);
depth = Gridding.depthVec(j);
strName = Gridding.strNameVec{k};

%Solution = struct;
tic;
[Solution] = AMT_Simulation(strName, station, depth, Gridding, FileNames, Options);
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