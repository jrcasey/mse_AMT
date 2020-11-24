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
    xlim([0.001 0.1])
    set(gca,'FontSize',20)
    %set(gca,'XScale','log')
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

% regression of the Holling 2 curve:
x1_vec = reshape(x1,size(x1,1).*size(x1,2),1)
y1_vec = reshape(y1,size(y1,1).*size(y1,2),1)

% Holling 2 fit
Holling2 = @(x,xdata) (xdata.*x(1))./(x(2) - xdata);
Holling2_fit = nlinfit(x1_vec, y1_vec, Holling2, [1e5 0.1])

% log fit
y1_vec(find(isnan(y1_vec))) = 0;
y1_vec2 = log10(y1_vec);
y1_vec2(find(y1_vec2==-Inf)) = 0.01;
logFit = fitlm(x1_vec,log10(y1_vec2))


fig = figure
h1 = plot(x1,y1,'.k','MarkerSize',20)
hold on
h2 = plot(x2,y2,'.r','MarkerSize',20)
%h3 = plot(x3,y3,'.g','MarkerSize',20)
h4 = plot(x4,y4,'.b','MarkerSize',20)
%h5 = plot(x5,y5,'.c','MarkerSize',20)
%x1_vec_lin = linspace(0,max(max(x1)),1000);
%h6 = plot(x1_vec_lin,(x1_vec_lin.*Holling2_fit(1))./(Holling2_fit(2)-x1_vec_lin),'-c','LineWidth',3)
xlabel('Growth rate [h^-^1]')
ylabel('Abundance [cells ml^-^1]')
set(gca,'FontSize',20)
set(gca,'YScale','log')
%set(gca,'XScale','log')
%title('Population')
%legend([h1(1) h2(1) h3(1) h4(1) h5(1)],'MSE vs. Johnson et al., 2006','Vaulot et al., 1995','Liu et al., 1995','Liu et al., 1997','Ribalet','Location','SouthEast')
legend([h1(1) h2(1) h4(1)],'MSE vs. Johnson et al., 2006','Vaulot et al., 1995','Liu et al., 1997','Location','SouthEast')

fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Population/','AMT13_Growth_Abundance_vs_literature.eps');

saveas(fig,fileName,'epsc')


% another coloring the depths
fig = figure
zColors = varycolor(Gridding.nZ);
for a = 1:Gridding.nZ
    plot(x1(a,:),y1(a,:),'.','MarkerSize',20,'MarkerFaceColor',zColors(a,:),'MarkerEdgeColor',zColors(a,:))
    hold on
    depthStr{a} = strcat(mat2str(Gridding.depthVec(a)),' m');
end
xlabel('Growth rate [h^-^1]')
ylabel('Abundance [cells ml^-^1]')
set(gca,'FontSize',20)
set(gca,'YScale','log')
title('Population')
legend(depthStr)

fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Population/AMT13_Growth_Abundance_Population_Depth.eps');
saveas(fig,fileName,'epsc');

%% Biomass versus production, depth integrated

for a = 1:Gridding.nStations
    zVec = Gridding.depthVec;
    prodVec = 1e3 .* PopulationSolution.Fluxes(:,a,find(strcmp('R00024',FullSolution.PanGEM.rxns))); % mol m-3 h-1
    biomassVec = 1e6 * CruiseData.Pro(Gridding.stationsVec2(a),:)' .* ( (4/3) .* pi() .* ((PopulationSolution.r_opt(:,a)).^3) .* 280 .* 1e-15 .* (1/12.011)); % mol m-3
    prodInt(a) = trapz(zVec,prodVec);
    biomassInt(a) = trapz(zVec,biomassVec);
end

fig = figure
plot(1000 * prodInt,1000 * biomassInt,'.k','MarkerSize',20)
xlabel('GPP [mmol C m^-^2 h^-^1]')
ylabel('Biomass [mmol C m^-^2]')
set(gca,'FontSize',20);

fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Growth/Population/AMT13_Prod_Biomass_Integrated.eps');
saveas(fig,fileName,'epsc');

% implied growth rate versus latitude
fig = figure
plot(CruiseData.Lat(Gridding.stationsVec2),prodInt./biomassInt,'.k','MarkerSize',20)
xlabel('Latitude')
ylabel('Biomass turnover rate [h^-^1]')
ylim([0 0.07])
set(gca,'FontSize',20);


%% Multi-color contours, abundance and growth
ecotypeColors = varycolor(numel(ecotypeList));

x = CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;

    
ecotypeColors2 = [220 97 26; 223 194 125; 208 28 139; 128 128 255; 1 133 113] ./ 255;

% Growth
% get mask for each ecotype
for a = 1:numel(ecotypeList)
    notnanIdx = find(~isnan(EcotypeSolution.(ecotypeList{a}).Growth));
    EcotypeMask{a} = zeros(numel(y),numel(x));
    EcotypeMask{a}(notnanIdx) = a;
end
fig = figure
colormap(ecotypeColors2)
for a = 1:numel(ecotypeList)
    alphaDat = EcotypeSolution.(ecotypeList{a}).Growth;
    %alphaDat = EcotypeSolution.(ecotypeList{a}).Growth .* CruiseData.Pro(Gridding.stationsVec2,:)';
    alphaDat(find(alphaDat==Inf)) = 0;
    h(a) = image(x,y,EcotypeMask{a});
    set(h(a),'AlphaData',(alphaDat - nanmin(nanmin(alphaDat)) ) ./ ( nanmax(nanmax(alphaDat)) - nanmin(nanmin(alphaDat)) ) )
    hold on
end
hc = colorbar('TickLabels',[ecotypeList2,{''}])
%ylabel(hc, 'Relative Growth Rate')
xlabel('Latitude')
ylabel('Depth [m]')
set(gca,'FontSize',20)


% Abundance
% get mask for each ecotype
for a = 1:numel(ecotypeList)
    notnanIdx = find(~isnan(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)'));
    EcotypeMask{a} = zeros(numel(y),numel(x));
    EcotypeMask{a}(notnanIdx) = a;
end
fig = figure
colormap(ecotypeColors2)
for a = 1:numel(ecotypeList)
    h(a) = image(x,y,EcotypeMask{a});
    set(h(a),'AlphaData',CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)' ./ nanmax(nanmax(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)')))
    hold on
end
hc = colorbar('TickLabels',[ecotypeList2,{''}])
%ylabel(hc, 'Relative Abundance')
xlabel('Latitude')
ylabel('Depth [m]')
set(gca,'FontSize',20)


%% Johnson-like figure (contours of growth and abundance)

fig = figure
for a = 1:numel(ecotypeList)
    subplot(2,5,a)
    z = log10(CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)');
    imagesc(x,y,z)
    caxis([0 5])
    if a==4
        caxis([0 3])
    end
    xlabel('Latitude')
    if a == 1
    ylabel('Depth')
    end
    colormap('jet')
    hc = colorbar
    %ylabel(hc, strcat(ecotypeList2{a},' [cells ml^-^1]'))
    set(gca,'FontSize',20)
end



for a = 1:numel(ecotypeList)
    subplot(2,5,a+5)
    z = EcotypeSolution.(ecotypeList{a}).Growth;
    imagesc(x,y,z)
    caxis([0 0.1])
    xlabel('Latitude')
    if a == 1
    ylabel('Depth')
    end
    colormap('jet')
    hc = colorbar
    %ylabel(hc, strcat(ecotypeList2{a},' [h^-^1]'))
    set(gca,'FontSize',20)
end


% for presentation

fig = figure
aVec = [2 1 3 5];
ecoVec = [{'eMIT9312'},{'eMED4'},{'eNATL2A'},{'eMIT9313'}]
for a = 1:4
    subplot(4,1,a)
    z = EcotypeSolution.(ecotypeList{aVec(a)}).Growth;
    imagesc(x,y,z)
    ht = text(min(x),40,ecoVec{a});
    caxis([0 0.1])
    set(ht,'FontSize',20,'color','w')
    if a<4
        set(gca,'XTick',[])
    end
    
    if a==4
    xlabel('Latitude')
    end
    ylabel('Depth')
    colormap('jet')
    hc = colorbar
    %ylabel(hc, strcat(ecotypeList2{a},' [h^-^1]'))
    set(gca,'FontSize',20)
end


%% Growth with abundance contours

fig = figure

for a = 1:numel(ecotypeList)
    
z1 = EcotypeSolution.(ecotypeList{a}).Growth;
z2 = CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)';
ecotype_levels = [1e-3 1e-4];
subplot(5,1,a)
h1 = imagesc(x,y,z1)
hold on
[c,h2] = contour(x,y,z2./1e6,ecotype_levels);
set(h2,'LineWidth',3,'ShowText','on','LineColor','w')
clabel(c,h2,'FontSize',20,'Color','w')
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Growth rate [h^-^1]')
set(gca,'FontSize',20)
end





