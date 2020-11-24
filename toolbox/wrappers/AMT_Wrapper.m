%% AMT Cruise Simulation
% Runs and analyzes the MSE simulation for the AMT13 cruise.
% Set up currently to run on Centos7
%% set root dir (for server only)
%cd Documents/MATLAB/GitHub/mse_AMT/
%addpath '/nfs/cnhlab001/cnh/projects/jrcasey/mse/test/mosek/9.1/toolbox/r2015a'	
addpath(genpath('~/mse_AMT/'));
addpath(genpath('~/mosek/'));
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
FileNames.StrMod_Path = 'data/GEM/StrMod.mat';
% List of strains to be analyzed
FileNames.strainList_Path = 'data/db/strainList.mat';
% Current OGTDat
FileNames.OGTDat_Path = 'data/db/OGTDat.csv';
% Cruise data
FileNames.CruiseDB_filename = 'data/envData/AMT13_Gridded2.csv';
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
Gridding.stationsVec =[{'A13_03'},{'A13_06'},{'A13_08'},{'A13_11'}, ...
    {'A13_13'},{'A13_15'},{'A13_16'},{'A13_19'},{'A13_22'},{'A13_25'}, ...
    {'A13_28'},{'A13_31'},{'A13_32'},{'A13_35'},{'A13_38'},{'A13_41'}, ...
    {'A13_42'},{'A13_45'},{'A13_48'},{'A13_51'},{'A13_54'},{'A13_57'}, ...
    {'A13_60'},{'A13_63'},{'A13_66'},{'A13_68'},{'A13_69'},{'A13_72'}, ...
    {'A13_75'},{'A13_76'},{'A13_77'},{'A13_78'}];
Gridding.stationsVec2 = [2 4 5 7 8 9 10 12 15 17 19 21 22 24 26 ...
    28 29 31 33 35 37 39 41 43 45 46 47 49 51 52 53 54];
Gridding.nStations = numel(Gridding.stationsVec);
% Depth (m)
Gridding.minZ = 10;
Gridding.maxZ = 200;
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
job_array_idx = str2num(getenv('SLURM_ARRAY_TASK_ID'));

% After running compile_AMT_subjobs.m, sometimes there are missing entries.
% Save those missing entry indices in a file called missingFileNo.mat and
% run those locally or back on the server from here: (comment out
% otherwise)
%  job_array_idx_temp = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%  load('data/output/missingFileNo.mat');
%  job_array_idx = missingFileNo(job_array_idx_temp);

% locate coordinates
[i,j,k] = ind2sub(size(idxMat),job_array_idx);

%% Run simulation
station = Gridding.stationsVec{i};
depth = Gridding.depthVec(j);
strName = Gridding.strNameVec{k};


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

% save local copy
%save(cell2str(FileNames.destination_fileName),'FullSolution');

% save on server
save(strcat('/nobackup1/jrcasey/','Solution_',num2str(job_array_idx),'.mat'),'Solution');

% change path to nobackup/jrcasey/ for saving the full output.