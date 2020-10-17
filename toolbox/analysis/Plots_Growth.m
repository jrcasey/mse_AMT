%% Plot Growth rates
% Contour plots of growth rates for each strain, ecotype,
% and the full population. Also scatter plots of growth versus abundance
% for ecotypes and population

%% Load data
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;

%% Parse out Gridding, CruiseData, FileNames, and PanGEM from FullSolution
Gridding = FullSolution.Gridding;
CruiseData = FullSolution.CruiseData;

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

%% Strain contours
for a = 1:Gridding.nStr
    var = Gridding.strNameVec{a};
    varName = var;
    varUnits = 'h^-^1';
    saveMe = 1;
    fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Strains/','AMT13_Contour_Growth_',var,'.eps');
    z = StrainSolution.Growth(:,:,a);
    contour_AMT(x, y, z, var, varName, varUnits, saveMe, fileName);
end

%% Ecotype contours
for a = 1:numel(ecotypeList)
    var = ecotypeList{a};
    varName = var;
    varUnits = 'h^-^1';
    saveMe = 1;
    fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Ecotypes/','AMT13_Contour_Growth_',var,'.eps');
    z = EcotypeSolution.(var).Growth;
    contour_AMT(x, y, z, var, varName, varUnits, saveMe, fileName);
end

%% Population contour
var = 'Population';
varName = var;
varUnits = 'h^-^1';
saveMe = 1;
fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Population/','AMT13_Contour_Growth_',var,'.eps');
z = PopulationSolution.Growth;
contour_AMT(x, y, z, var, varName, varUnits, saveMe, fileName);

%% Ecotypes versus Zinser abundance

for a = 1:numel(ecotypeList)
    var = ecotypeList{a};
    varName = var;
    varUnits = 'h^-^1';
    fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Ecotypes/','AMT13_Growth_Abundance_',var,'.eps');
    x = EcotypeSolution.(var).Growth;
    y = CruiseData.(var)(Gridding.stationsVec2,:)'
    fig = figure('visible','off')
    plot(x,y,'.k','MarkerSize',20)
    xlabel('Growth rate [h^-^1]')
    ylabel('Abundance [cells ml^-^1]')
    set(gca,'FontSize',20)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    title(ecotypeList2{a})
    saveas(fig,fileName,'epsc')
end

%% Population versus abundance
var = 'Population';
varName = var;
varUnits = 'h^-^1';
fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Population/','AMT13_Growth_Abundance_',var,'.eps');
x = PopulationSolution.Growth;
y = CruiseData.Pro(Gridding.stationsVec2,:)'
fig = figure('visible','off')
plot(x,y,'.k','MarkerSize',20)
xlabel('Growth rate [h^-^1]')
ylabel('Abundance [cells ml^-^1]')
set(gca,'FontSize',20)
set(gca,'YScale','log')
title('Population')
saveas(fig,fileName,'epsc')


%% Compare with in situ observations
% These are the most directly comparable data we have found so far...
% Vaulot: In situ depth profiles of growth rates and abundances of Prochlorococcus
% were quantified in the eastern Equatorial Pacific in April (casts 1 and
% 2) and October (casts 3, 4, and 5) of 1992.
% Liu et al., 1995 (Landry) - Kanamycin method at ALOHA
% Liu et al., 1997 - cell cycle analysis in W and E equatorial Pac, and two
% cruises at ALOHA
% Ribalet - Gradients 1 cruise SeaFlow using MPM

fileName = '/Users/jrcasey/Documents/MATLAB/CBIOMES/Data/Environmental_Data/Growth_vs_abundance/insitu_growth_abundance.csv';
insitu = readtable(fileName,'Delimiter',',','ReadVariableNames',true);

% Vaulot subset
vaulot_idx = 1:59;

% find Liu/landry subset
landry_idx = 60:68;

% find Liu1997 subset
liu_idx = 68:104;

% find ribalet data
ribalet_idx = 105:209;


% assign variables
x1 = PopulationSolution.Growth;
y1 = CruiseData.Pro(Gridding.stationsVec2,:)'
x2 = insitu.growthRate(vaulot_idx) ./ 10;
y2 = insitu.Abundance(vaulot_idx);
x3 = insitu.growthRate(landry_idx)./10;
y3 = insitu.Abundance(landry_idx);
x4 = insitu.growthRate(liu_idx) ./ 10;
y4 = insitu.Abundance(liu_idx);
x5 = insitu.growthRate(ribalet_idx) ./ 10;
y5 = insitu.Abundance(ribalet_idx);

% remove zero growth
zero_mu_idx = find(x1==0);
y1(zero_mu_idx) = NaN;


fig = figure
h1 = plot(x1,y1,'.k','MarkerSize',20)
hold on
h2 = plot(x2,y2,'.r','MarkerSize',20)
h3 = plot(x3,y3,'.g','MarkerSize',20)
h4 = plot(x4,y4,'.b','MarkerSize',20)
%h5 = plot(x5,y5,'.c','MarkerSize',20)
xlabel('Growth rate [h^-^1]')
ylabel('Abundance [cells ml^-^1]')
set(gca,'FontSize',20)
set(gca,'YScale','log')
title('Population')
%legend([h1(1) h2(1) h3(1) h4(1) h5(1)],'MSE vs. Johnson et al., 2006','Vaulot et al., 1995','Liu et al., 1995','Liu et al., 1997','Ribalet','Location','SouthEast')
legend([h1(1) h2(1) h3(1) h4(1)],'MSE vs. Johnson et al., 2006','Vaulot et al., 1995','Liu et al., 1995','Liu et al., 1997','Location','SouthEast')

fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Population/','AMT13_Growth_Abundance_vs_literature',var,'.eps');

saveas(fig,fileName,'epsc')





