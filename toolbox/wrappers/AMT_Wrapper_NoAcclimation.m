%% MSE No acclimation
% This script re-runs each simulation using only the batch acclimated
% physiological state (no acclimation steps are done).

%% Load FullSolution_L2
load('data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;
Gridding = FullSolution.Gridding;
FileNames = FullSolution.FileNames;
CruiseData = FullSolution.CruiseData;

% load strain models
load(FileNames.StrMod_Path);

%% Load irradiance data
load(FileNames.IrrDat_fileName);
[IrrDat2] = standardizeIrr(IrrDat,Gridding.lambdaVec,Gridding.depthVec); %mmoles photons m-2 h-1 bandwidth*nm-1
IrrDat3 = reshape([IrrDat2{:}],numel(Gridding.depthVec),numel(Gridding.lambdaVec),numel(IrrDat2));

%% Pick a strain and a grid point to start with
strName = 'SB';
station_ind = 7;
depth_ind = 3;

Tp_init = squeeze(FullSolution.(strName).TpOpt(depth_ind,station_ind,:));
BOF_init = squeeze(FullSolution.(strName).BOF_coefs(depth_ind,station_ind,:))

%% Initialize StrMod
StrMod1 = StrMod.(strName);
StrMod1.id = strName;
StrMod1.description = strcat('Un-acclimated strain model with no specified constraints or compositions - Prochlorococcus strain',{' '},strName);

%% Prepare strain model for simulation
StrMod2 = prepEnvMod(StrMod1);
StrMod2.description = strcat('Un-acclimated strain model initialized with size and composition constraints - Prochlorococcus strain',{' '},strName);

for a = 1:Gridding.nZ
    for b = 1:Gridding.nStations
        depth_ind = a;
        station_ind = Gridding.stationsVec2(b);
        %% Run TpOpt without optimization
        % Load TpDat
        TpDat = readtable(FileNames.TpDat_fileName,'ReadVariableNames',true,'Delimiter',',');
        
        % Compile environmental data to retrieve and format TpOpt parameter table
        TpDat = getTpOptParams(CruiseData,StrMod2,TpDat,station_ind,depth_ind);
        
        % set transport bounds for StrMod
        [StrMod3, vOut] = updateModel_TpOpt(StrMod2,Tp_init,TpDat);
        
        % Save optimal cell radius
        r_opt = StrMod3.rInit;
        
        %% Run PhysOpt without optimization
        % Load physOpt and pigOpt constraints together
        Constraints = readtable(FileNames.PhysOptPigOptConstraints_fileName,'ReadVariableNames',true,'Delimiter',',');
        
        % parse entries and get crude fractions
        clb_ind = find(contains(Constraints.Constraint,'clb_'));
        cub_ind = find(contains(Constraints.Constraint,'cub_'));
        crudeFractions = strrep(Constraints.Constraint(clb_ind),'clb_','');
        
        % Prepare pigOpt
        % load PigDB
        PigDB = readtable(FileNames.PigDB_fileName,'ReadVariableNames',true);
        
        % Specify which pigments to include
        PigsIncluded = [{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'},{'alpha_Carotene'},{'Zeaxanthin'}];
        PigMW = [891.4731  905.4566  536.8726  568.8714];
        nPigs = numel(PigsIncluded);
        
        % Get reformated pigment specific absorption spectra
        PigDat = getSpecificAbsorption(PigDB,PigsIncluded,PigMW,Gridding.lambdaVec); %m2 mmol-1 bandwidth(nm)^-^1
        
        % Compute pigment absorption total
        PigAbs = nansum(PigDat .* repmat(squeeze(IrrDat3(depth_ind,:,station_ind))',1,nPigs) .* Gridding.bandwidth); % mmol photons [mmol pig]-1 h-1
        
        [StrMod4] = updateModel_pigPhysOpt(StrMod3,BOF_init,PigAbs,crudeFractions);
        StrMod4.description = strcat('Acclimated strain model Prochlorococcus strain',{' '},strName);
        
        tempSol3 = solveLP(StrMod3,1);
        tempSol4 = solveLP(StrMod4,1);
        if tempSol3.stat ==1
            mu3(a,b) = tempSol3.f;
        else
            mu3(a,b) = 0;
        end
        if tempSol4.stat ==1
            mu4(a,b) = tempSol4.f;
        else
            mu4(a,b) = 0;
        end
    end
end

    