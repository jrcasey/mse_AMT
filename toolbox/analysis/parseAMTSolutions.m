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

for i = 1:Gridding.nStr
    strainInd = find(strcmp(Gridding.strNameVec{i},orgDatabase.StrainName));
    ecotype{i} = orgDatabase.Ecotype{strainInd};
    ecotypeInd(i) = find(strcmp(ecotype{i},ecotypeList));
end

for i = 1:numel(ecotypeList)
    strEco_idx{i} = find(ecotypeInd == i);
end
ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];


%% All strains in each ecotype
StrainSolution = struct;
for i = 1:Gridding.nStr
    StrainSolution.Growth(:,:,i) = FullSolution.(Gridding.strNameVec{i}).Growth;
    StrainSolution.Fluxes(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).Fluxes;
    StrainSolution.Shadow(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).Shadow;
    StrainSolution.BOF_coefs(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).BOF_coefs;
    StrainSolution.TpOpt(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).TpOpt;
    StrainSolution.r_opt(:,:,i) = FullSolution.(Gridding.strNameVec{i}).r_opt;
    StrainSolution.pigAbs(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).pigAbs;
    StrainSolution.uptakeBounds(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).uptakeBounds;
    StrainSolution.StrMod1_growth(:,:,i) = FullSolution.(Gridding.strNameVec{i}).StrMod1_growth;
    StrainSolution.StrMod2_growth(:,:,i) = FullSolution.(Gridding.strNameVec{i}).StrMod2_growth;
    StrainSolution.StrMod3_growth(:,:,i) = FullSolution.(Gridding.strNameVec{i}).StrMod3_growth;
    StrainSolution.StrMod4_growth(:,:,i) = FullSolution.(Gridding.strNameVec{i}).StrMod4_growth;
    StrainSolution.StrMod5_growth(:,:,i) = FullSolution.(Gridding.strNameVec{i}).StrMod5_growth;
end

%% Winning strain from each ecotype
% Find index of winning strain
for i = 1:Gridding.nZ
    for j = 1:Gridding.nStations
        for k = 1:numel(ecotypeList2)
            [junk, growthMaxStrain_idx{k}(i,j)] = nanmax(StrainSolution.Growth(i,j,find(ecotypeInd==k)),[],3);
        end
    end
end


for i = 1:numel(ecotypeList2)
    growthMaxStrainEco_idx{i} = strEco_idx{i}(growthMaxStrain_idx{i});
end


EcotypeSolution = struct;
for i = 1:Gridding.nZ
    for j = 1:Gridding.nStations
        for k = 1:numel(ecotypeList2)
            EcotypeSolution.(ecotypeList{k}).Growth(i,j) = StrainSolution.Growth(i,j,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).Fluxes(i,j,:) = StrainSolution.Fluxes(i,j,:,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).Shadow(i,j,:) = StrainSolution.Shadow(i,j,:,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).BOF_coefs(i,j,:) = StrainSolution.BOF_coefs(i,j,:,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).TpOpt(i,j,:) = StrainSolution.TpOpt(i,j,:,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).r_opt(i,j) = StrainSolution.r_opt(i,j,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).pigAbs(i,j,:) = StrainSolution.pigAbs(i,j,:,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).uptakeBounds(i,j,:) = StrainSolution.uptakeBounds(i,j,:,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).StrMod1_growth(i,j) = StrainSolution.StrMod1_growth(i,j,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).StrMod2_growth(i,j) = StrainSolution.StrMod2_growth(i,j,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).StrMod3_growth(i,j) = StrainSolution.StrMod3_growth(i,j,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).StrMod4_growth(i,j) = StrainSolution.StrMod4_growth(i,j,growthMaxStrainEco_idx{k}(i,j));
            EcotypeSolution.(ecotypeList{k}).StrMod5_growth(i,j) = StrainSolution.StrMod5_growth(i,j,growthMaxStrainEco_idx{k}(i,j));
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
PopulationSolution.uptakeBounds = zeros(Gridding.nZ,Gridding.nStations,size(StrainSolution.uptakeBounds,3));
PopulationSolution.StrMod1_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod2_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod3_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod4_growth = zeros(Gridding.nZ,Gridding.nStations);
PopulationSolution.StrMod5_growth = zeros(Gridding.nZ,Gridding.nStations);




massConversion = (4/3) .* pi() .* 280e-15; % g cell-1
% calculate relative abudance of each ecotype for weighted averages
CruiseData = CruiseData;
relAb = zeros(Gridding.nZ,Gridding.nStations,numel(ecotypeList));
for i = 1:numel(ecotypeList)
    relAb(:,:,i) = CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)';
end
totalAb = nansum(relAb,3);
for i = 1:numel(ecotypeList)
    relAb2(:,:,i) = relAb(:,:,i) ./ totalAb;
end



for i = 1:numel(ecotypeList) 
    tempPopGrowth(:,:,i) = EcotypeSolution.(ecotypeList{i}).Growth  .* relAb2(:,:,i);
    tempPopropt(:,:,i) = EcotypeSolution.(ecotypeList{i}).r_opt .* relAb2(:,:,i);
    tempPopFluxes(:,:,:,i) = repmat(CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)',1,1,numel(PanGEM.rxns)) .* EcotypeSolution.(ecotypeList{i}).Fluxes .* repmat(massConversion .* EcotypeSolution.(ecotypeList{i}).r_opt.^3,1,1,numel(PanGEM.rxns));
    tempPopShadow(:,:,:,i) = EcotypeSolution.(ecotypeList{i}).Shadow;
    tempPopBOF(:,:,:,i) = repmat(CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.BOF_coefs,3)) .* EcotypeSolution.(ecotypeList{i}).BOF_coefs .* repmat(massConversion .* EcotypeSolution.(ecotypeList{i}).r_opt.^3,1,1,size(StrainSolution.BOF_coefs,3));
    %tempPopTpOpt(:,:,:,i) = repmat(CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)',1,1,size(tempTpOpt,3)) .* EcotypeSolution.(ecotypeList{i}).TpOpt .* repmat((4.*pi()) ./ (1e-6.*EcotypeSolution.(ecotypeList{i}).r_opt).^2,1,1,size(tempTpOpt,3));
    tempPopTpOpt(:,:,:,i) = repmat(CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.TpOpt,3)) .* EcotypeSolution.(ecotypeList{i}).TpOpt;
    tempPoppigAbs(:,:,:,i) = repmat(CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.pigAbs,3)) .* EcotypeSolution.(ecotypeList{i}).pigAbs .* repmat(massConversion .* EcotypeSolution.(ecotypeList{i}).r_opt.^3,1,1,size(StrainSolution.pigAbs,3));
    tempPopuptakeBounds(:,:,:,i) = repmat(CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)',1,1,size(StrainSolution.uptakeBounds,3)) .* EcotypeSolution.(ecotypeList{i}).uptakeBounds .* repmat(massConversion .* EcotypeSolution.(ecotypeList{i}).r_opt.^3,1,1,size(StrainSolution.uptakeBounds,3));
    tempPopStrMod1_growth(:,:,i) = EcotypeSolution.(ecotypeList{i}).StrMod1_growth  .* relAb2(:,:,i);
    tempPopStrMod2_growth(:,:,i) = EcotypeSolution.(ecotypeList{i}).StrMod2_growth  .* relAb2(:,:,i);
    tempPopStrMod3_growth(:,:,i) = EcotypeSolution.(ecotypeList{i}).StrMod3_growth  .* relAb2(:,:,i);
    tempPopStrMod4_growth(:,:,i) = EcotypeSolution.(ecotypeList{i}).StrMod4_growth  .* relAb2(:,:,i);
    tempPopStrMod5_growth(:,:,i) = EcotypeSolution.(ecotypeList{i}).StrMod5_growth  .* relAb2(:,:,i);
    
end


PopulationSolution.Growth = nansum(tempPopGrowth,3); % h-1 (weighted average)
PopulationSolution.Fluxes = nansum(tempPopFluxes,4); % mmol ml-1 h-1
PopulationSolution.Shadow = nansum(tempPopShadow,4); % ND
PopulationSolution.r_opt = nansum(tempPopropt,3);
PopulationSolution.BOF_coefs = nansum(tempPopBOF,4); % g ml-1
PopulationSolution.TpOpt = nansum(tempPopTpOpt,4); % Tp ml-1
PopulationSolution.pigAbs = nansum(tempPoppigAbs,4); % mmol ml-1 h-1
PopulationSolution.uptakeBounds = nansum(tempPopuptakeBounds,4); % mmol ml-1 h-1
PopulationSolution.StrMod1_growth = nansum(tempPopStrMod1_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod2_growth = nansum(tempPopStrMod2_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod3_growth = nansum(tempPopStrMod3_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod4_growth = nansum(tempPopStrMod4_growth,3); % h-1 (weighted average)
PopulationSolution.StrMod5_growth = nansum(tempPopStrMod5_growth,3); % h-1 (weighted average)

end
