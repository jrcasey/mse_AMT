function [FullSolution] = get_AMT_Results_server(ResultsDirectory,PanGEM,FileNames,Gridding,CruiseData)
% Retrieves all individual solutions from the server and concatenates them
% into a single structure. Also saves Gridding and CruiseData into the same
% structure for analysis. 

% Inputs
% ResultsDirectory      -       String. Directory where the solutions mat
%                               files can be found. 
% FileNames             -       Structure. Should be saved from each run.
% Gridding              -       Structure. Should be saved from each run.
% CruiseData            -       Structure. Should be saved from each run.

% Outputs
% FullSolution          -       Structure. All solutions.

%% Preallocate structure matrices
FullSolution = struct;
%% Find index from SLURM
% preallocate a 3d matrix of dimensions nStations, nZ, nStr
idxMat = zeros(Gridding.nStations,Gridding.nZ,Gridding.nStr);


%% Load solutions and store in FullSolution.
% Import each file, identify its coordinates, and save output to the
% structure. 

files = dir(fullfile(ResultsDirectory, '*.mat'));
nFiles = numel(files);
idxMat = zeros(Gridding.nStations,Gridding.nZ,Gridding.nStr);

% Check that nFiles is consistent with Gridding
if nFiles ~= Gridding.nStations .* Gridding.nZ .* Gridding.nStr
    msg = 'number of files and gridding dimensions are inconsistent';
    error(msg)
end

% If there are files missing, find out which and add some blanks
expectedFiles = 1:(Gridding.nStations .* Gridding.nZ .* Gridding.nStr);
for a = 1:nFiles
    tempFile = files(a).name;
    startInd = regexp(tempFile,'_');
    endInd = regexp(tempFile,'.mat');
    fileNo(a) = str2num(tempFile(startInd+1:endInd-1));
end
fileNo2 = sort(fileNo);
missingFileNo = setdiff(expectedFiles,fileNo2);


%     
% 
% % Preallocate matrices to structure (need to load a solution to do this)
% load Solution_1
% for a = 1:Gridding.nStr
%     FullSolution.(Gridding.strNameVec{a}).Growth = zeros(Gridding.nStations,Gridding.nZ);
%     nFluxes = numel(Solution.Fluxes);
%     FullSolution.(Gridding.strNameVec{a}).Fluxes = zeros(Gridding.nStations,Gridding.nZ,nFluxes);
%     nBOF_coefs = numel(Solution.BOF_coefs);
%     FullSolution.(Gridding.strNameVec{a}).BOF_coefs = zeros(Gridding.nStations,Gridding.nZ,nBOF_coefs);
%     nTpOpt = numel(Solution.TpOpt);
%     FullSolution.(Gridding.strNameVec{a}).TpOpt = zeros(Gridding.nStations,Gridding.nZ,nTpOpt);
%     FullSolution.(Gridding.strNameVec{a}).r_opt = zeros(Gridding.nStations,Gridding.nZ);
%     npigAbs = numel(Solution.pigAbs);
%     FullSolution.(Gridding.strNameVec{a}).pigAbs = zeros(Gridding.nStations,Gridding.nZ,npigAbs);
%     nuptakeBounds = numel(Solution.uptakeBounds);
%     FullSolution.(Gridding.strNameVec{a}).uptakeBounds = zeros(Gridding.nStations,Gridding.nZ,nuptakeBounds);
%     FullSolution.(Gridding.strNameVec{a}).runtime = zeros(Gridding.nStations,Gridding.nZ);
% end
% 
% for a = 1:nFiles
%     load(strcat(ResultsDirectory,files(a).name)); % Will be called Solution
%     if isfield(Solution,'Fluxes')
%     % get subjob number
%     startChar = '_';
%     endChar = '\.';
%     startInd = regexp(files(a).name,startChar);
%     endInd = regexp(files(a).name,endChar);
%     job_array_idx = str2num(files(a).name(startInd+1:endInd-1));
%     
%     % get coordinates
%     [i,j,k] = ind2sub(size(idxMat),job_array_idx);
%     
%     % Assign data to FullSolution
%     FullSolution.(Solution.strName).Growth(i,j) = Solution.Growth;
%     FullSolution.(Solution.strName).Fluxes(i,j,:) = Solution.Fluxes;
%     FullSolution.(Solution.strName).BOF_coefs(i,j,:) = Solution.BOF_coefs;
%     FullSolution.(Solution.strName).TpOpt(i,j,:) = Solution.TpOpt;
%     FullSolution.(Solution.strName).r_opt(i,j) = Solution.r_opt;
%     FullSolution.(Solution.strName).pigAbs(i,j,:) = Solution.pigAbs;
%     FullSolution.(Solution.strName).uptakeBounds(i,j,:) = Solution.uptakeBounds;
%     FullSolution.(Solution.strName).runtime(i,j) = Solution.runtime;
%     
%     end
% end

% Preallocate matrices to structure (need to load a solution to do this)
load Solution_1
for a = 1:Gridding.nStr
    FullSolution.(Gridding.strNameVec{a}).Growth = zeros(Gridding.nZ,Gridding.nStations);
    nFluxes = numel(Solution.Fluxes);
    FullSolution.(Gridding.strNameVec{a}).Fluxes = zeros(Gridding.nZ,Gridding.nStations,nFluxes);
    nShadow = numel(Solution.Shadow);
    FullSolution.(Gridding.strNameVec{a}).Shadow = zeros(Gridding.nZ,Gridding.nStations,nShadow);
    nBOF_coefs = numel(Solution.BOF_coefs);
    FullSolution.(Gridding.strNameVec{a}).BOF_coefs = zeros(Gridding.nZ,Gridding.nStations,nBOF_coefs);
    nTpOpt = numel(Solution.TpOpt);
    FullSolution.(Gridding.strNameVec{a}).TpOpt = zeros(Gridding.nZ,Gridding.nStations,nTpOpt);
    FullSolution.(Gridding.strNameVec{a}).r_opt = zeros(Gridding.nZ,Gridding.nStations);
    npigAbs = numel(Solution.pigAbs);
    FullSolution.(Gridding.strNameVec{a}).pigAbs = zeros(Gridding.nZ,Gridding.nStations,npigAbs);
    nuptakeBounds = numel(Solution.uptakeBounds);
    FullSolution.(Gridding.strNameVec{a}).uptakeBounds = zeros(Gridding.nZ,Gridding.nStations,nuptakeBounds);
    FullSolution.(Gridding.strNameVec{a}).runtime = zeros(Gridding.nZ,Gridding.nStations);
    FullSolution.(Gridding.strNameVec{a}).StrMod1_growth = zeros(Gridding.nZ,Gridding.nStations);
    FullSolution.(Gridding.strNameVec{a}).StrMod2_growth = zeros(Gridding.nZ,Gridding.nStations);
    FullSolution.(Gridding.strNameVec{a}).StrMod3_growth = zeros(Gridding.nZ,Gridding.nStations);
    FullSolution.(Gridding.strNameVec{a}).StrMod4_growth = zeros(Gridding.nZ,Gridding.nStations);
    FullSolution.(Gridding.strNameVec{a}).StrMod5_growth = zeros(Gridding.nZ,Gridding.nStations);
end

for a = 1:nFiles
%for a = 17523:nFiles
    load(strcat(ResultsDirectory,files(a).name)); % Will be called Solution
    if isfield(Solution,'Fluxes')
    % get subjob number
    startChar = '_';
    endChar = '\.';
    startInd = regexp(files(a).name,startChar);
    endInd = regexp(files(a).name,endChar);
    job_array_idx = str2num(files(a).name(startInd+1:endInd-1));
    
    % get coordinates
    [i,j,k] = ind2sub(size(idxMat),job_array_idx);
    
    % Assign data to FullSolution
    FullSolution.(Solution.strName).Growth(j,i) = Solution.Growth;
    FullSolution.(Solution.strName).Fluxes(j,i,:) = Solution.Fluxes;
    FullSolution.(Solution.strName).Shadow(j,i,:) = Solution.Shadow;
    FullSolution.(Solution.strName).BOF_coefs(j,i,:) = Solution.BOF_coefs;
    FullSolution.(Solution.strName).TpOpt(j,i,:) = Solution.TpOpt;
    FullSolution.(Solution.strName).r_opt(j,i) = Solution.r_opt;
    FullSolution.(Solution.strName).pigAbs(j,i,:) = Solution.pigAbs;
    FullSolution.(Solution.strName).uptakeBounds(j,i,:) = Solution.uptakeBounds;
    FullSolution.(Solution.strName).runtime(j,i) = Solution.runtime;
    FullSolution.(Solution.strName).StrMod1_growth(j,i) = Solution.StrMod1_growth;
    FullSolution.(Solution.strName).StrMod2_growth(j,i) = Solution.StrMod2_growth;
    FullSolution.(Solution.strName).StrMod3_growth(j,i) = Solution.StrMod3_growth;
    if isempty(Solution.StrMod4_growth)
        FullSolution.(Solution.strName).StrMod4_growth(j,i) = NaN;
    else
        FullSolution.(Solution.strName).StrMod4_growth(j,i) = Solution.StrMod4_growth;
    end
    FullSolution.(Solution.strName).StrMod5_growth(j,i) = Solution.StrMod5_growth;
    end
end

%% Add gridding and cruise data
FullSolution.FileNames = FileNames;
FullSolution.Gridding = Gridding;
FullSolution.CruiseData = CruiseData;
FullSolution.dateAnalyzed = datestr(date,'yyyymmdd');
FullSolution.PanGEM = PanGEM;


end


