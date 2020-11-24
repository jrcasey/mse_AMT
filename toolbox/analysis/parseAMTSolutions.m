function [StrainSolution, EcotypeSolution, PopulationSolution] = parseAMTSolutions(FullSolution)
%% parseAMTSolutions
% Generate structures from FullSolution for each strain, each ecotype (fittest strain), and the abundance weighted population

Gridding = FullSolution.Gridding;
PanGEM = FullSolution.PanGEM;
CruiseData = FullSolution.CruiseData;

%% Assign strains to ecotypes
orgDatabase = readtable('GitHub/mse_AMT/data/db/orgDatabase.csv','Delimiter',',','ReadVariableNames',true);
ecotypeList = [{'HLI'},{'HLII'},{'LLI'},{'LLII_LLIII'},{'LLIV'}];
ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];

for a = 1:Gridding.nStr
    strainInd = find(strcmp(Gridding.strNameVec{a},orgDatabase.StrainName));
    ecotype{a} = orgDatabase.Ecotype{strainInd};
    ecotypeInd(a) = find(strcmp(ecotype{a},ecotypeList));
end

for a = 1:numel(ecotypeList)
    strEco_idx{a} = find(ecotypeInd == a);
end
ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];


%% All strains
StrainSolution = struct;
for a = 1:Gridding.nStr
    StrainSolution.Growth(:,:,a) = FullSolution.(Gridding.strNameVec{a}).Growth;
    StrainSolution.Fluxes(:,:,:,a) = FullSolution.(Gridding.strNameVec{a}).Fluxes;
    StrainSolution.Shadow(:,:,:,a) = FullSolution.(Gridding.strNameVec{a}).Shadow;
    StrainSolution.BOF_coefs(:,:,:,a) = FullSolution.(Gridding.strNameVec{a}).BOF_coefs;
    StrainSolution.TpOpt(:,:,:,a) = FullSolution.(Gridding.strNameVec{a}).TpOpt;
    StrainSolution.r_opt(:,:,a) = FullSolution.(Gridding.strNameVec{a}).r_opt;
    StrainSolution.pigAbs(:,:,:,a) = FullSolution.(Gridding.strNameVec{a}).pigAbs;
    StrainSolution.S_star(:,:,:,a) = FullSolution.(Gridding.strNameVec{a}).S_star;
    StrainSolution.uptakeBounds(:,:,:,a) = FullSolution.(Gridding.strNameVec{a}).uptakeBounds;
    StrainSolution.StrMod1_growth(:,:,a) = FullSolution.(Gridding.strNameVec{a}).StrMod1_growth;
    StrainSolution.StrMod2_growth(:,:,a) = FullSolution.(Gridding.strNameVec{a}).StrMod2_growth;
    StrainSolution.StrMod3_growth(:,:,a) = FullSolution.(Gridding.strNameVec{a}).StrMod3_growth;
    StrainSolution.StrMod4_growth(:,:,a) = FullSolution.(Gridding.strNameVec{a}).StrMod4_growth;
    StrainSolution.StrMod5_growth(:,:,a) = FullSolution.(Gridding.strNameVec{a}).StrMod5_growth;
end

%% Winning strain from each ecotype
% Find index of winning strain
for a = 1:Gridding.nZ
    for b = 1:Gridding.nStations
        for k = 1:numel(ecotypeList2)
            [junk, growthMaxStrain_idx{k}(a,b)] = nanmax(StrainSolution.Growth(a,b,find(ecotypeInd==k)),[],3);
        end
    end
end


for a = 1:numel(ecotypeList2)
    growthMaxStrainEco_idx{a} = strEco_idx{a}(growthMaxStrain_idx{a});
end


EcotypeSolution = struct;
for a = 1:Gridding.nZ
    for b = 1:Gridding.nStations
        for k = 1:numel(ecotypeList2)
            EcotypeSolution.(ecotypeList{k}).Growth(a,b) = StrainSolution.Growth(a,b,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).Fluxes(a,b,:) = StrainSolution.Fluxes(a,b,:,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).Shadow(a,b,:) = StrainSolution.Shadow(a,b,:,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).BOF_coefs(a,b,:) = StrainSolution.BOF_coefs(a,b,:,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).TpOpt(a,b,:) = StrainSolution.TpOpt(a,b,:,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).r_opt(a,b) = StrainSolution.r_opt(a,b,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).pigAbs(a,b,:) = StrainSolution.pigAbs(a,b,:,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).S_star(a,b,:) = StrainSolution.S_star(a,b,:,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).uptakeBounds(a,b,:) = StrainSolution.uptakeBounds(a,b,:,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).StrMod1_growth(a,b) = StrainSolution.StrMod1_growth(a,b,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).StrMod2_growth(a,b) = StrainSolution.StrMod2_growth(a,b,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).StrMod3_growth(a,b) = StrainSolution.StrMod3_growth(a,b,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).StrMod4_growth(a,b) = StrainSolution.StrMod4_growth(a,b,growthMaxStrainEco_idx{k}(a,b));
            EcotypeSolution.(ecotypeList{k}).StrMod5_growth(a,b) = StrainSolution.StrMod5_growth(a,b,growthMaxStrainEco_idx{k}(a,b));
        end
    end
end

%% Weighted population 
PopulationSolution = struct;
PopulationSolution.Growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.Fluxes = zeros(Gridding.nZ,Gridding.nStations,numel(PanGEM.rxns));
PopulationSolution.Shadow = zeros(Gridding.nZ,Gridding.nStations,numel(PanGEM.rxns));
PopulationSolution.BOF_coefs = zeros(Gridding.nZ,Gridding.nStations,size(StrainSolution.BOF_coefs,3));
PopulationSolution.TpOpt = zeros(Gridding.nZ,Gridding.nStations,size(StrainSolution.TpOpt,3));
PopulationSolution.r_opt = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.pigAbs = zeros(Gridding.nZ,Gridding.nStations,size(StrainSolution.pigAbs,3));
PopulationSolution.S_star = zeros(Gridding.nZ,Gridding.nStations,size(StrainSolution.S_star,3));
PopulationSolution.uptakeBounds = zeros(Gridding.nZ,Gridding.nStations,size(StrainSolution.uptakeBounds,3));
PopulationSolution.StrMod1_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod2_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod3_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod4_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod5_growth = zeros(Gridding.nZ,Gridding.nStations);




massConversion = (4/3) .* pi() .* 280e-15 .* 2; % g cell-1
% calculate relative abudance of each ecotype. When doing the weighting,
% there's an issue with nan's, so let's zero out the relative abundance of
% nan growth instances

ecoAb = zeros(Gridding.nZ,Gridding.nStations,numel(ecotypeList));
for a = 1:numel(ecotypeList)
    ecoAb(:,:,a) = CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)'; % cells ml-1
    %[nanRowIdx, nanColIdx] = find(isnan(EcotypeSolution.(ecotypeList{i}).Growth)); % find row and column indices of nan's in the ecotype solutions
    %ecoAb(nanRowIdx,nanColIdx,i) = 0; % assign nan growth indices to be zero abundance
end
totalAb = nansum(ecoAb,3); % cells ml-1
for a = 1:numel(ecotypeList)
    relAb(:,:,a) = ecoAb(:,:,a) ./ totalAb; % relative abundance
end

% Caluclate weighted averages or absolute totals for all fields
for a = 1:numel(ecotypeList) 
    tempPopGrowth(:,:,a) = EcotypeSolution.(ecotypeList{a}).Growth  .* relAb(:,:,a);
    tempPopropt(:,:,a) = EcotypeSolution.(ecotypeList{a}).r_opt .* relAb(:,:,a);
    tempPopFluxes(:,:,:,a) = repmat(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)',1,1,numel(PanGEM.rxns)) .* EcotypeSolution.(ecotypeList{a}).Fluxes .* repmat(massConversion .* EcotypeSolution.(ecotypeList{a}).r_opt.^3,1,1,numel(PanGEM.rxns));
    tempPopShadow(:,:,:,a) = EcotypeSolution.(ecotypeList{a}).Shadow;
    tempPopBOF(:,:,:,a) = repmat(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.BOF_coefs,3)) .* EcotypeSolution.(ecotypeList{a}).BOF_coefs .* repmat(massConversion .* EcotypeSolution.(ecotypeList{a}).r_opt.^3,1,1,size(StrainSolution.BOF_coefs,3));
    %tempPopTpOpt(:,:,:,i) = repmat(CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)',1,1,size(tempTpOpt,3)) .* EcotypeSolution.(ecotypeList{i}).TpOpt .* repmat((4.*pi()) ./ (1e-6.*EcotypeSolution.(ecotypeList{i}).r_opt).^2,1,1,size(tempTpOpt,3));
    tempPopTpOpt(:,:,:,a) = repmat(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.TpOpt,3)) .* EcotypeSolution.(ecotypeList{a}).TpOpt;
    tempPoppigAbs(:,:,:,a) = repmat(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.pigAbs,3)) .* EcotypeSolution.(ecotypeList{a}).pigAbs .* repmat(massConversion .* EcotypeSolution.(ecotypeList{a}).r_opt.^3,1,1,size(StrainSolution.pigAbs,3));
    tempPopS_star(:,:,:,a) = EcotypeSolution.(ecotypeList{a}).S_star .* relAb(:,:,a);
    tempPopuptakeBounds(:,:,:,a) = repmat(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.uptakeBounds,3)) .* EcotypeSolution.(ecotypeList{a}).uptakeBounds .* repmat(massConversion .* EcotypeSolution.(ecotypeList{a}).r_opt.^3,1,1,size(StrainSolution.uptakeBounds,3));
    tempPopStrMod1_growth(:,:,a) = EcotypeSolution.(ecotypeList{a}).StrMod1_growth  .* relAb(:,:,a);
    tempPopStrMod2_growth(:,:,a) = EcotypeSolution.(ecotypeList{a}).StrMod2_growth  .* relAb(:,:,a);
    tempPopStrMod3_growth(:,:,a) = EcotypeSolution.(ecotypeList{a}).StrMod3_growth  .* relAb(:,:,a);
    tempPopStrMod4_growth(:,:,a) = EcotypeSolution.(ecotypeList{a}).StrMod4_growth  .* relAb(:,:,a);
    tempPopStrMod5_growth(:,:,a) = EcotypeSolution.(ecotypeList{a}).StrMod5_growth  .* relAb(:,:,a);
end


PopulationSolution.Growth = nansum(tempPopGrowth,3); % h-1 (weighted average)
PopulationSolution.Fluxes = nansum(tempPopFluxes,4); % mmol ml-1 h-1
PopulationSolution.Shadow = nansum(tempPopShadow,4); % ND
PopulationSolution.r_opt = nansum(tempPopropt,3);
PopulationSolution.BOF_coefs = nansum(tempPopBOF,4); % g ml-1
PopulationSolution.TpOpt = nansum(tempPopTpOpt,4); % Tp ml-1
PopulationSolution.pigAbs = nansum(tempPoppigAbs,4); % mmol ml-1 h-1
PopulationSolution.S_star = nansum(tempPopS_star,4); % mol m-3
PopulationSolution.uptakeBounds = nansum(tempPopuptakeBounds,4); % mmol ml-1 h-1
PopulationSolution.StrMod1_growth = nansum(tempPopStrMod1_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod2_growth = nansum(tempPopStrMod2_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod3_growth = nansum(tempPopStrMod3_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod4_growth = nansum(tempPopStrMod4_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod5_growth = nansum(tempPopStrMod5_growth,3); % h-1 (weighted average)

end
