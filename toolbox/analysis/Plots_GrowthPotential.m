%% Plots GrowthPotential
% Growth capabilities of each strain on sole C,N,P, and S substrates

%% Load data
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;

%% Parse out Gridding, CruiseData, FileNames, and PanGEM from FullSolution
Gridding = FullSolution.Gridding;
CruiseData = FullSolution.CruiseData;
FileNames = FullSolution.FileNames;

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

%% Load StrMod

load(FileNames.StrMod_Path)
StrMod_StrList = fieldnames(StrMod);
tempSol = solveLP(FullSolution.PanGEM,1);
transIdx = find(strcmp('Transport',FullSolution.PanGEM.subSystems))
auxo_idx = find(tempSol.x(transIdx) ~=0);
auxotrophyRxns = FullSolution.PanGEM.rxns(transIdx(auxo_idx));
auxotrophyRxns = {};

for a = 1:numel(StrMod_StrList)
    tempStrMod = StrMod.(StrMod_StrList{a});
    [soleSources_N{a}] = getSoleElementSources(tempStrMod, 'N', auxotrophyRxns);
    [soleSources_P{a}] = getSoleElementSources(tempStrMod, 'P', auxotrophyRxns);
end
[soleSources_N{end+1}] = getSoleElementSources(FullSolution.PanGEM, 'N', auxotrophyRxns);
[soleSources_P{end+1}] = getSoleElementSources(FullSolution.PanGEM, 'P', auxotrophyRxns);

% select only the subsets for which the PanGEM can grow
select_TpRxns_N = soleSources_N{end}.TpRxns(find(abs(soleSources_N{end}.Growth) > 1e-4));
select_TpRxns_P = soleSources_P{end}.TpRxns(find(abs(soleSources_P{end}.Growth) > 1e-4));

for a = 1:numel(soleSources_N)
    for b = 1:numel(select_TpRxns_N)
        tpIdx = find(strcmp(select_TpRxns_N{b},soleSources_N{a}.TpRxns));
        if ~isempty(tpIdx) & abs(soleSources_N{a}.Growth(tpIdx)) > 1e-3;
            Growth_N(b,a) = 1;
        else
            Growth_N(b,a) = 0;
        end
    end
end

for a = 1:numel(soleSources_P)
    for b = 1:numel(select_TpRxns_P)
        tpIdx = find(strcmp(select_TpRxns_P{b},soleSources_P{a}.TpRxns));
        if ~isempty(tpIdx) & abs(soleSources_P{a}.Growth(tpIdx)) > 1e-3;
            Growth_P(b,a) = 1;
        else
            Growth_P(b,a) = 0;
        end
    end
end

% concatenate
soleSources_NP = [select_TpRxns_N; select_TpRxns_P]
Growth_NP = [Growth_N;Growth_P];

% re-order strains
ecoOrder = [70 horzcat(strEco_idx{:})];

% condense and re-order substrates
sourcesList = [{'L-alanine'},{'L-arginine'},{'L-aspartate'},{'L-cysteine'},{'L-glutamine'},{'L-glycine'},{'L-serine'},...
    {'Betaine'},{'Choline'},{'Cyanate'},{'L-ornithine'},{'Urea'},...
    {'Ammonia'},{'Nitrite'},{'Nitrate'},...
    {'Orthophosphate'},{'Phosphoglycerate'},{'Phosphite'},{'Methylphosphonate'},{'Phosphonoacetate'},{'Aminoethylphosphonate'}]';
sourcesOrder = [6,8,10,16,21,24,33,11,14,15,31,34,7,30,29,38,40,42,37,41,36];



% xticks
xtickLabels = [StrMod_StrList;{'PanGEM'}];

ecotypeBackground = [ones(numel(sourcesList),1), repmat(2,numel(sourcesList),numel(strEco_idx{1})), repmat(3,numel(sourcesList),numel(strEco_idx{2})) repmat(4,numel(sourcesList),numel(strEco_idx{3})) repmat(5,numel(sourcesList),numel(strEco_idx{4})) repmat(6,numel(sourcesList),numel(strEco_idx{5}))]
X = Growth_NP(sourcesOrder,ecoOrder);
X(find(~X)) = 0.1
map = [0.5 0 0.5; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0.5 0]


figure
h=imagesc(ecotypeBackground)
colormap(map)
set(h,'AlphaData',X)
set(gca,'XTick',1:size(Growth_NP,2),'XTickLabels',xtickLabels(ecoOrder),'YTick',1:numel(sourcesList), 'YTickLabels',sourcesList);
xtickangle(90);
set(gca,'FontSize',16);


