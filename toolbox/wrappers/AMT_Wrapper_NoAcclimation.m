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
% needs to be one with all Tp's... SB works. Also need to set Tp n at nG
strName = 'SB';
station_ind = 7;
depth_ind = 3;

Tp_init = squeeze(FullSolution.(strName).TpOpt(depth_ind,station_ind,:));
Tp_init(5) = Tp_init(4);
Tp_init(2) = Tp_init(4);
BOF_init = squeeze(FullSolution.(strName).BOF_coefs(depth_ind,station_ind,:))


for c = 1:Gridding.nStr
    strName = Gridding.strNameVec{c};
    
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
                mu3(a,b,c) = tempSol3.f;
            else
                mu3(a,b,c) = 0;
            end
            if tempSol4.stat ==1
                mu4(a,b,c) = tempSol4.f;
            else
                mu4(a,b,c) = 0;
            end
        end
    end
    
end

%% Histograms of all absolute and relative differences
dim1 = ceil(sqrt(Gridding.nStr));

% Fitness gains (absolute)
for a = 1:Gridding.nStr
    strName = Gridding.strNameVec{a};
    mu3_delta = FullSolution.(strName).StrMod3_growth - abs(mu3(:,:,a));
    mu4_delta = FullSolution.(strName).StrMod4_growth - abs(mu4(:,:,a));
    mu3_delta(find(mu3_delta==Inf)) = NaN; mu3_delta(find(mu3_delta==-Inf)) = NaN; mu3_delta(find(isnan(mu3_delta))) = NaN;
    mu4_delta(find(mu4_delta==Inf)) = NaN; mu4_delta(find(mu4_delta==-Inf)) = NaN; mu4_delta(find(isnan(mu4_delta))) = NaN;

    mu34_delta_abs(:,:,a) = [mu3_delta(:),mu4_delta(:)];
end

mu34_delta_abs2 = reshape(permute(mu34_delta_abs,[1 3 2]),size(mu34_delta_abs,1).*size(mu34_delta_abs,3),size(mu34_delta_abs,2))
mu34_delta_abs2(find(mu34_delta_abs2<=0.001)) = NaN;

% Fitness gains (relative increase)
for a = 1:Gridding.nStr
    strName = Gridding.strNameVec{a};
    mu3_delta = (FullSolution.(strName).StrMod3_growth - abs(mu3(:,:,a))) ./ abs(mu3(:,:,a));
    mu4_delta = (FullSolution.(strName).StrMod4_growth - abs(mu4(:,:,a))) ./ abs(mu4(:,:,a));
    mu3_delta(find(mu3_delta==Inf)) = NaN; mu3_delta(find(mu3_delta==-Inf)) = NaN; mu3_delta(find(isnan(mu3_delta))) = NaN;
    mu4_delta(find(mu4_delta==Inf)) = NaN; mu4_delta(find(mu4_delta==-Inf)) = NaN; mu4_delta(find(isnan(mu4_delta))) = NaN;

    mu34_delta_rel(:,:,a) = [mu3_delta(:),mu4_delta(:)];
end

mu34_delta_rel2 = reshape(permute(mu34_delta_rel,[1 3 2]),size(mu34_delta_rel,1).*size(mu34_delta_rel,3),size(mu34_delta_rel,2))
mu34_delta_rel2(find(mu34_delta_rel2<=0.001)) = NaN;
mu34_delta_rel2(find(mu34_delta_rel2>1)) = NaN;


figure
subplot(2,1,1)
histogram(mu34_delta_abs2(:,1),100);
hold on
histogram(mu34_delta_abs2(:,2),100);
legend('OptTrans','PhysOpt');
xlabel('fitness gain [h^-^1]')
set(gca,'FontSize',20)
subplot(2,1,2)
histogram(100*mu34_delta_rel2(:,1),100);
hold on
histogram(100*mu34_delta_rel2(:,2),100);
legend('OptTrans','PhysOpt');
xlabel('relative fitness gain [%]')
set(gca,'FontSize',20)



