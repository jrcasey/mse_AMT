function [Solution] = AMT_Simulation(strName, station, depth, Gridding, FileNames, Options)

% Run MSE pipeline for one strain on a sample from the MESO-SCOPE cruise. 

% Inputs
% strName       -       String. StrainID (e.g., MED4)
% station       -       Double. Station number (4-15)
% depth         -       Double. Depth from depthVec in Gridding (m)
% Gridding      -       Struct. Contains the ranges and resolutions for
%                       stations, depths, and wavelengths. 
% FileNames     -       Struct. Contains a list of filenames to be used for
%                       loading strain models, cruise data, and databases.
% Options       -       Struct. Contains various options like 'save',
%                       'maxIter', whether to save the acclimated StrMod's. 

% Outputs
% Solution      -       Struct. Solution contains growth rates, fluxes, BOF
%                       coefficients, uptake rate bounds, pigment specific
%                       absorption, transporter numbers, and cell size.


%% Load strain model
% StrainModelWrapper.m should be run for each new PanGEM version. Just
% update the PanGEM path and the StrMods path.
load(FileNames.StrMod_Path);
StrMod1 = StrMod.(strName);
StrMod1.id = strName;
StrMod1.description = strcat('Un-acclimated strain model with no specified constraints or compositions - Prochlorococcus strain',{' '},strName);

%% Prepare strain model for simulation
StrMod2 = prepEnvMod(StrMod1);
StrMod2.description = strcat('Un-acclimated strain model initialized with size and composition constraints - Prochlorococcus strain',{' '},strName);

%% Import Cruise Data
% Load and format nutrient concentrations
CruiseDB = readtable(FileNames.CruiseDB_filename,'Delimiter',',','ReadVariableNames',true);
CruiseData = getCruiseData(CruiseDB,Gridding.depthVec);

%% Import spectral irradiance
% load and format a HyperPro profile
load(FileNames.IrrDat_fileName);
[IrrDat2] = standardizeIrr(IrrDat,Gridding.lambdaVec,Gridding.depthVec); %mmoles photons m-2 h-1 bandwidth*nm-1
IrrDat3 = reshape([IrrDat2{:}],numel(Gridding.depthVec),numel(Gridding.lambdaVec),numel(IrrDat2));
CruiseData.PAR = squeeze(nansum(IrrDat3,2))';
% % Temporary!! Normalize to local noon. For now just use the maximum Irr
% NoonScalar = sum(CruiseData.PAR,2)./max(sum(CruiseData.PAR,2));
% CruiseData.PAR = CruiseData.PAR./repmat(NoonScalar,1,Gridding.nZ);
% for i = 1:Gridding.nZ
%     IrrDat3(i,:,:) = squeeze(IrrDat3a(i,:,:))./repmat(NoonScalar',numel(Gridding.lambdaVec),1);
% end

clear IrrDat IrrDat2

%% Choose a station and depth
%station_ind = find(station==Gridding.stationsVec); % Station index
station_ind = find(strcmp(station,CruiseData.Stations))
%station_name = strcat('Station_',num2str(station));
[junk, depth_ind] = min(abs(depth-Gridding.depthVec));

%% Optimal Transport
% Load TpDat
TpDat = readtable(FileNames.TpDat_fileName,'ReadVariableNames',true,'Delimiter',',');

% Compile environmental data to retrieve and format TpOpt parameter table
TpDat = getTpOptParams(CruiseData,StrMod2,TpDat,station_ind,depth_ind);

% Set up LP
fmax = 0.085;
[TpOptLP] = getTpOptLP(StrMod2, fmax, TpDat, Options.maxIter_TpOpt);

%%%%%%%%%%%%% TESTING LP PARAMS %%%%%%%%%%%%%
TpOptLP.options.StepTolerance = 1e-8;
% TpOptLP.options.OptimalityTolerance = 1e-6;
% TpOptLP.options.ConstraintTolerance = 1e-5;
% TpOptLP.options.ObjectiveLimit = -1000;

% Solve LP
tic
[n_opt,fval,exitflag,output,lambdaTemp,grad,hessian] = fmincon(TpOptLP);
toc

% Save S_star
T1 = CruiseData.T(station_ind,depth_ind);
S_star = calc_S_star(StrMod2,n_opt(1),TpDat,T1);

% Update model with TpRxns constraints
[StrMod3, vOut] = updateModel_TpOpt(StrMod2,n_opt,TpDat);

% Save optimal cell radius
r_opt = StrMod3.rInit;

%% MM Acclimation and Photoacclimation

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
a = nansum(PigDat .* repmat(squeeze(IrrDat3(depth_ind,:,station_ind))',1,nPigs) .* Gridding.bandwidth); % mmol photons [mmol pig]-1 h-1

% Set up LP
[physOptLP] = getPhysOptLP(StrMod3,Constraints, a, PigsIncluded, Options.maxIter_physOpt);

%%%%%%%%%%%%% TESTING LP PARAMS %%%%%%%%%%%%%
%physOptLP.options.StepTolerance = 1e-5;
%physOptLP.options.OptimalityTolerance = 1e-6;
%physOptLP.options.ConstraintTolerance = 1e-5;

% Solve LP
tic
[xOut,fval,exitflag,output] = fmincon(physOptLP);
toc
% Save absorption fluxes
pigAbs = a.*xOut(12:15); % mmol photons gDW-1 h-1

% Update model
[StrMod4] = updateModel_pigPhysOpt(StrMod3,xOut,a,crudeFractions);
StrMod4.description = strcat('Acclimated strain model Prochlorococcus strain',{' '},strName);

%% Temperature correction
OGTDat = readtable(FileNames.OGTDat_Path,'Delimiter',',','ReadVariableNames',true);
T = [273.15 + CruiseData.T(station_ind,depth_ind)];
TCorr = TemperatureCorr(OGTDat,strName,T);

%% Solve FBA and temporarily store output
[sol hsSolOut] = solveLP(StrMod4,1);
if sol.stat==1
    fluxes = sol.x.*TCorr(1);
    growth = -sol.f.*TCorr(1);
    shadow = hsSolOut.y;
else
    fluxes = zeros(numel(StrMod4.rxns),1);
    growth = 0;
    shadow = zeros(numel(StrMod4.mets),1);
end
BOF = xOut;

%% Store and save output
% replace zero growth entries with nans
if growth < 1e-3
growth = NaN;
fluxes = nan(numel(fluxes),1);
BOF = nan(numel(BOF),1);
n_opt = nan(numel(n_opt),1);
r_opt = NaN;
pigAbs = nan(numel(pigAbs),1);
vOut = nan(numel(vOut),1);
S_star = nan(numel(S_star),1);
end

% Store in results structure
Solution = struct;
Solution.Growth = growth;
Solution.Fluxes = fluxes;
Solution.Shadow = shadow;
Solution.BOF_coefs = BOF;
Solution.TpOpt = n_opt;
Solution.r_opt = r_opt;
Solution.pigAbs = pigAbs;
Solution.uptakeBounds = vOut;
Solution.S_star = S_star;
if Options.saveStrMod
    Solution.StrMod = StrMod4;
end

% Store growth rates for each model version
sol1 = solveLP(StrMod1,1); % Unacclimated, unconditioned
if sol1.stat==1
    Solution.StrMod1_growth = -sol1.f;
else
    Solution.StrMod1_growth = NaN;
end
sol2 = solveLP(StrMod2,1); % Unacclimated, conditioned
if sol2.stat==1
    Solution.StrMod2_growth = -sol2.f;
else
    Solution.StrMod2_growth = NaN;
end
sol3 = solveLP(StrMod3,1); % Tp acclimated, conditioned
if sol1.stat==1
    Solution.StrMod3_growth = -sol3.f;
else
    Solution.StrMod3_growth = NaN;
end
sol4 = solveLP(StrMod4,1); % Tp and MM acclimated, conditioned
if sol1.stat==1
    Solution.StrMod4_growth = -sol4.f;
else
    Solution.StrMod4_growth = NaN;
end
sol5 = solveLP(StrMod4,1); % Tp and MM acclimated, Temperature corrected, conditioned
if sol5.stat==1
    Solution.StrMod5_growth = -sol5.f.*TCorr;
else
    Solution.StrMod5_growth = NaN;
end


% Save
if Options.saveSolution
save('CBIOMES/Data/Environmental_Data/Cruises/MESO-SCOPE/Solution_20191122.mat','-struct','Solution');
end

end


