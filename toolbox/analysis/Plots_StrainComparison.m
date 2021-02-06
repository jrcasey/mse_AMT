%% Plots StrainComparison
% Plots for the comparison of SB and MIT9314

%% Load data
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;

%% Parse out Gridding, CruiseData, FileNames, and PanGEM from FullSolution
Gridding = FullSolution.Gridding;
CruiseData = FullSolution.CruiseData;
PanGEM = FullSolution.PanGEM;

%% Retrieve solutions for strains, ecotypes, and population
[StrainSolution, EcotypeSolution, PopulationSolution] = parseAMTSolutions(FullSolution);

%% Gridded domain
x = CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;

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

%% Load strMod
load(FullSolution.FileNames.StrMod_Path);


%% Jaccard distance networks
minPanGEM = PanGEM;
minBOFind = find(strcmp('BIOMASSMinimal',minPanGEM.rxns));
minPanGEM.c = zeros(numel(minPanGEM.rxns),1);
minPanGEM.c(minBOFind) = 1;

[essentialRxns, results] = solveCS(minPanGEM);

for a = 1:Gridding.nStr
    for b = 1:Gridding.nStr
        S1 = full(StrMod.(Gridding.strNameVec{a}).S);
        S2 = full(StrMod.(Gridding.strNameVec{b}).S);
        S1(:,essentialRxns) = [];
        S2(:,essentialRxns) = [];
        S1(find(S1)) = 1;
        S2(find(S2)) = 1;
        similarity_nonEssential(a,b) = jaccard(S1,S2);
    end
end

% sort by ecotype
newOrder = horzcat(strEco_idx{:});

% sort by similarity within ecotypes
newOrder2 = ones(Gridding.nStr,1);
b=1;

for a = 1:numel(ecotypeList)
    %[val, newOrder2a{a}] = sort(sum(similarity_nonEssential(strEco_idx{a},:),2),'descend');

    [val, newOrder2a{a}] = sort(sum(similarity_nonEssential(strEco_idx{a},strEco_idx{a})),'descend');
    newOrder2(b:b+numel(strEco_idx{a})-1) = strEco_idx{a}(newOrder2a{a});
    b=b+numel(strEco_idx{a});
end

similarityEco_nonEssential = similarity_nonEssential(newOrder2,newOrder2);

    
figure
imagesc(similarityEco_nonEssential)
set(gca,'XTick',1:Gridding.nStr,'XTickLabels',Gridding.strNameVec(newOrder2),'YTick',1:Gridding.nStr,'YTickLabels',Gridding.strNameVec(newOrder2))
xtickangle(90)
hb = colorbar
ylabel(hb,'Jaccard Similarity')
set(hb,'FontSize',20)
colormap('jet')

%% growth comparison for SB and MIT-9314
str1 = 'SB';
str2 = 'MIT9314'; % MIT0604?

strIdx(1) = find(strcmp(str1,Gridding.strNameVec));
strIdx(2) = find(strcmp(str2,Gridding.strNameVec));


figure
z = 100*(StrainSolution.Growth(:,:,strIdx(1)) - 1.2*StrainSolution.Growth(:,:,strIdx(2))) ./ nanmean(StrainSolution.Growth(:,:,strIdx),3);
%z = 100*(StrainSolution.StrMod4_growth(:,:,strIdx(1)) - StrainSolution.StrMod4_growth(:,:,strIdx(2))) ./ nanmean(StrainSolution.StrMod4_growth(:,:,strIdx),3);
z(:,28) = NaN;
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
cmap = [1 1 1; jet(20)];
colormap(cmap)
%caxis([0 0.1])
hb = colorbar
ylabel(hb,'Relative growth rate %')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)
caxis([-100 150])


% Set differences

str1_Mod = StrMod.(str1);
str2_Mod = StrMod.(str2);

str1_present = find(sum(abs(str1_Mod.S),1));
str2_present = find(sum(abs(str2_Mod.S),1));

str1_not_str2 = setdiff(str1_present,str2_present);
str2_not_str1 = setdiff(str2_present,str1_present);

str1_not_str2_fluxes = StrainSolution.Fluxes(:,:,str1_not_str2,strIdx(1));
str2_not_str1_fluxes = StrainSolution.Fluxes(:,:,str2_not_str1,strIdx(2));

dim = ceil(sqrt(numel(str1_not_str2)));
figure
for a = 1:numel(str1_not_str2)
    subplot(dim,dim,a)
    imagesc(x,y,str1_not_str2_fluxes(:,:,a))
    colorbar
    title(PanGEM.rxns(str1_not_str2(a)))
end

dim = ceil(sqrt(numel(str2_not_str1)));
figure
for a = 1:numel(str2_not_str1)
    subplot(dim,dim,a)
    imagesc(x,y,str2_not_str1_fluxes(:,:,a))
    colorbar
    title(PanGEM.rxns(str2_not_str1(a)))
end


%% Difference between all fluxes

for a = 1:numel(PanGEM.rxns)
    z1 = StrainSolution.Fluxes(:,:,a,strIdx(1));
    z2 = StrainSolution.Fluxes(:,:,a,strIdx(2));
    z1_sum = nansum(nansum(abs(z1)));
    z2_sum = nansum(nansum(abs(z2)));
    % check if they're different
    if z1_sum > 1 || z2_sum > 1
        if z1_sum./z2_sum > 1.5 || z1_sum./z2_sum < 0.75
            isDifferent(a) = 1;
        end
    else isDifferent(a) = 0;
    end
end

different_idx = find(isDifferent);

dim1 = ceil(sqrt(numel(find(isDifferent))));
figure
for a = 1:numel(find(isDifferent))
    subplot(dim1,dim1,a)
    z = abs(StrainSolution.Fluxes(:,:,different_idx(a),strIdx(1))) - abs(StrainSolution.Fluxes(:,:,different_idx(a),strIdx(2)));
    imagesc(x,y,z)
    colorbar
    title(PanGEM.rxns(different_idx(a)))
end

%% Flux through all non-essential reactions
nonEssentialRxns = 1:numel(PanGEM.rxns);
nonEssentialRxns(essentialRxns) = [];


for a = 1:Gridding.nStr
    fluxesSum(:,a) = squeeze(nansum(nansum(abs(StrainSolution.Fluxes(:,:,nonEssentialRxns,a)),1),2));
end

nonEssential_nzFluxes = find(sum(fluxesSum,2)>1e-6);
nonEssential_nzRxns = PanGEM.rxns(nonEssentialRxns(nonEssential_nzFluxes));


dim1 = ceil(sqrt(numel(nonEssential_nzFluxes)));
figure
for a = 1:numel(nonEssential_nzFluxes)
    subplot(dim1,dim1,a)
    z = StrainSolution.Fluxes(:,:,nonEssentialRxns(nonEssential_nzFluxes(a)),strIdx(1)) - StrainSolution.Fluxes(:,:,nonEssentialRxns(nonEssential_nzFluxes(a)),strIdx(2));
    imagesc(x,y,z)
    colorbar
    title(nonEssential_nzRxns(a))
end

%% Sum of fluxes through PEPCKase

PEPCKase_idx = find(strcmp('R00345',PanGEM.rxns))
for a = 1:Gridding.nStr
    PEPCKase(a) = nansum(nansum(abs(StrainSolution.Fluxes(:,:,PEPCKase_idx,a))./StrainSolution.Growth(:,:,a)));
end
    

%% Compare average limitation for two strains
% Shadow prices as indicators of limitation

limitingMets = [{'Exciton'},{'L_Glutamate'}];
for i = 1:numel(limitingMets)
    limitingMets_idx(i) = find(strcmp(limitingMets{i},PanGEM.mets));
end


% N and E Shadow prices

figure
subplot(1,2,1)
z_raw_E = StrainSolution.Shadow(:,:,limitingMets_idx(1),strIdx(1));
z_norm_E = ( z_raw_E - nanmin(nanmin(z_raw_E)) ) ./ ( nanmax(nanmax(z_raw_E)) - nanmin(nanmin(z_raw_E)) );
z_raw_N = StrainSolution.Shadow(:,:,limitingMets_idx(2),strIdx(1));
z_norm_N = ( z_raw_N - nanmin(nanmin(z_raw_N)) ) ./ ( nanmax(nanmax(z_raw_N)) - nanmin(nanmin(z_raw_N)) );
z = z_norm_N - z_norm_E;
z(find(abs(z) < 1e-9)) = NaN;
imagesc(x,y,z)
caxis([-1 1])
cmap = [1 1 1; jet(1024)];
colormap(cmap);
cbh = colorbar
cbh.Ticks = linspace(-1,1,2)
cbh.TickLabels = [{}]
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)

subplot(1,2,2)
z_raw_E = StrainSolution.Shadow(:,:,limitingMets_idx(1),strIdx(2));
z_norm_E = ( z_raw_E - nanmin(nanmin(z_raw_E)) ) ./ ( nanmax(nanmax(z_raw_E)) - nanmin(nanmin(z_raw_E)) );
z_raw_N = StrainSolution.Shadow(:,:,limitingMets_idx(2),strIdx(2));
z_norm_N = ( z_raw_N - nanmin(nanmin(z_raw_N)) ) ./ ( nanmax(nanmax(z_raw_N)) - nanmin(nanmin(z_raw_N)) );
z = z_norm_N - z_norm_E;
z(find(abs(z) < 1e-9)) = NaN;
imagesc(x,y,z)
caxis([-1 1])
cmap = [1 1 1; jet(1024)];
colormap(cmap);
cbh = colorbar
cbh.Ticks = linspace(-1,1,2)
cbh.TickLabels = [{}]
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)



% a1 = axes;
% h1 = imagesc(x,y,z_norm_N);
% colormap(a1,'pink');
% c1 = colorbar;
% axis image
% 
% hold on
% a2 = axes;
% h2 = imagesc(x,y,z_norm_E);
% colormap(a2,'gray');
% c2 = colorbar;
% axis image




