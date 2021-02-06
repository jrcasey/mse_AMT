%% Plot BOF compositions
% Contour plots of each BOF component for ecotypes and the population. 

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

%% BOF components, compositions, and elemental quotas
% list components
BOF_components = [{'DNA'},{'Lipid'},{'Carbohydrate'},{'Protein'},{'VitaCofactors'},{'RNA'},{'NB'},{'Ions'},{'BioPool'},{'Pigments'},{'CellWall'},{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'},{'alpha_Carotene'},{'Zeaxanthin'}];

%% Figure for paper: area plot for one strain
% Shows the average BOF composition as a function of density, with
% individual data points overlaid as points. 

y = CruiseData.Density(Gridding.stationsVec2,:)';
x = squeeze(StrainSolution.BOF_coefs(:,:,:,find(strcmp('MIT9313',Gridding.strNameVec))));
dMets = nansum(x(:,:,[5 7 8 9 10]),3);
pigs = nansum(x(:,:,[12:15]),3);
xOther = x(:,:,[1 2 3 4 6 11]);
x2 = cat(3,xOther, dMets, pigs)
xMean = squeeze(nanmean(x2,2));
xMean_norm = xMean ./ repmat(sum(xMean,2),1,size(xMean,2))
BOF_components2 = [{'DNA'},{'Lipid'},{'Carbohydrate'},{'Protein'},{'RNA'},{'CellWall'},{'dMets'},{'Pigments'}]

figure
plot(xMean_norm,repmat(Gridding.depthVec',1,size(xMean,2)),'-.','MarkerSize',20,'LineWidth',3);
set(gca,'YDir','reverse')
xlim([0 0.35]);
legend(BOF_components2)
set(gca,'FontSize',20)

%% Plot a BOF component against a cruise variable
ecotypeColors = varycolor(numel(ecotypeList));
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        %y = StrainSolution.BOF_coefs(:,:,15,strEco_idx{i}(j)) ./ ( StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j)) + StrainSolution.BOF_coefs(:,:,13,strEco_idx{i}(j)) + StrainSolution.BOF_coefs(:,:,14,strEco_idx{i}(j)) );
        y = StrainSolution.BOF_coefs(:,:,3,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x =  CruiseData.PAR(Gridding.stationsVec2,:)';
    %y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,15) ./ ( EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12) + EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) + EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,14) );
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,3);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x =  CruiseData.PAR(Gridding.stationsVec2,:)';
%y = PopulationSolution.BOF_coefs(:,:,15) ./ ( PopulationSolution.BOF_coefs(:,:,12) + PopulationSolution.BOF_coefs(:,:,13) + PopulationSolution.BOF_coefs(:,:,14) );
y = PopulationSolution.BOF_coefs(:,:,3);

%plot(x,y,'.','MarkerSize',15,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Protein')
%xlabel('E_d(474) / E_d(450)')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)




%% Determine elemental quotas and enthalpy for each strain, ecotype, and population for
% each element. 

%%%%%%%%%%%%%%%% Takes about 45 minutes. %%%%%%%%%%%%%%%%%%

% load EnthalpyDB
fileName = '/Users/jrcasey/Documents/MATLAB/CBIOMES/Data/Cheminformatics/Enthalpy_BOF/Enthalpy_BOF.csv'; 
EnthalpyDB = readtable(fileName,'Delimiter',',','ReadVariableNames',true);

% all strains
for i = 1:Gridding.nZ
    for j = 1:Gridding.nStations
        for k = 1:Gridding.nStr
            BOF_coefs = squeeze(StrainSolution.BOF_coefs(i,j,1:11,k));
            [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
            Quota(i,j,:,k) = sum(MMComposition.DW,1); % mmol element gDW-1
            BOF_coefs2 = squeeze(StrainSolution.BOF_coefs(i,j,:,k));
            Enthalpy(i,j,k) = getBiomassEnthalpy(PanGEM,BOF_coefs2,EnthalpyDB); % KJ gDW-1
        end
    end
end
% ecotype winners
for i = 1:Gridding.nZ
    for j = 1:Gridding.nStations
        for k = 1:numel(ecotypeList)
            BOF_coefs = squeeze(EcotypeSolution.(ecotypeList{k}).BOF_coefs(i,j,1:11));
            [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
            QuotaEco(i,j,:,k) = sum(MMComposition.DW,1); % mmol element gDW-1
            BOF_coefs2 = squeeze(EcotypeSolution.(ecotypeList{k}).BOF_coefs(i,j,:));
            EnthalpyEco(i,j,k) = getBiomassEnthalpy(PanGEM,BOF_coefs2,EnthalpyDB); % KJ gDW-1
        end
    end
end
% population
for i = 1:Gridding.nZ
    for j = 1:Gridding.nStations
            BOF_coefs = squeeze(PopulationSolution.BOF_coefs(i,j,1:11));
            [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
            QuotaPop(i,j,:) = sum(MMComposition.DW,1); % mmol element ml-1
            BOF_coefs2 = squeeze(PopulationSolution.BOF_coefs(i,j,:));
            EnthalpyPop(i,j) = getBiomassEnthalpy(PanGEM,BOF_coefs2,EnthalpyDB); % KJ ml-1
    end
end
    

%% Plot population C:N against DIN

% C:N vs DIN
x = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
y = QuotaPop(:,:,1) ./ QuotaPop(:,:,3);
xVec = reshape(x,size(x,1)*size(x,2),1);
yVec = reshape(y,size(y,1)*size(y,2),1);
[xVecAsc, xVecIdx] = sort(xVec,'ascend');
yVecAsc = yVec(xVecIdx);
xVecBinned = logspace(log10(min(xVec)),log10(max(xVec)),10);
for a = 1:numel(xVecBinned)-1
    inBinIdx= find(xVecAsc >= xVecBinned(a) & xVecAsc < xVecBinned(a+1));
    yMean(a) = nanmean(yVecAsc(inBinIdx));
    yStd(a) = nanstd(yVecAsc(inBinIdx));
end

fig = figure
plot(x,y,'.k','MarkerSize',15)
hold on
plot(xVecBinned(2:end),yMean,'-k','LineWidth',3);
%hErr = errorbar(xVecBinned(2:end),yMean,yStd,'-r');
hErr2 = fill([xVecBinned(2:end) fliplr(xVecBinned(2:end))],[(yMean + yStd) fliplr(yMean - yStd)],'g');
%set(hErr,'LineWidth',3)
set(hErr2,'LineStyle','none')
hErr2.FaceAlpha = 0.2;
set(gca,'XScale','log')
ylabel('C:N')
%xlabel('E_d(474) / E_d(450)')
xlabel('DIN [nM]')
set(gca,'FontSize',20)


saveas(fig,'/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/BOF/CN_vs_DIN.eps','epsc');

%% C:N vs DIN again, but with all strains, ecotypes, and population
ecotypeColors = varycolor(numel(ecotypeList));
fig = figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        %x = CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
        x = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
        y = Quota(:,:,1,strEco_idx{i}(j)) ./ Quota(:,:,3,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
    y = QuotaEco(:,:,1,i) ./ QuotaEco(:,:,3,i)
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
y = QuotaPop(:,:,1) ./ QuotaPop(:,:,3);
plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('C:N')
%xlabel('E_d(474) / E_d(450)')
xlabel('DIN [nM]')
set(gca,'FontSize',20)


saveas(fig,'/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/BOF/CN_vs_DIN_ALL.eps','epsc');



%% Integrated C:N:P

for a = 1:Gridding.nStations
    C = QuotaPop(:,a,1) .* 1e6 .* 1e-3; % mol m-3
    C_int(a) = trapz(Gridding.depthVec,C); % mol m-2
    N = QuotaPop(:,a,3) .* 1e6 .* 1e-3; % mol m-3
    N_int(a) = trapz(Gridding.depthVec,N); % mol m-2
    P = QuotaPop(:,a,5) .* 1e6 .* 1e-3; % mol m-3
    P_int(a) = trapz(Gridding.depthVec,P); % mol m-2
end

cmap = [1 1 1; jet(256)];
figure

% C:N
subplot(3,1,1)
imagesc(x,y,QuotaPop(:,:,1) ./ QuotaPop(:,:,3))
hc = colorbar
colormap(cmap)
xlabel('Latitude')
ylabel('Depth')
ylabel(hc,'C:N [mol C (mol N)^-^1]')
set(gca,'FontSize',20)

% C:P
subplot(3,1,2)
imagesc(x,y,QuotaPop(:,:,1) ./ QuotaPop(:,:,5))
hc = colorbar
colormap(cmap)
xlabel('Latitude')
ylabel('Depth')
ylabel(hc,'C:P [mol C (mol P)^-^1]')
set(gca,'FontSize',20)
subplot(3,1,3)
imagesc(x,y,QuotaPop(:,:,3) ./ QuotaPop(:,:,5))
hc = colorbar
colormap(cmap)
xlabel('Latitude')
ylabel('Depth')
ylabel(hc,'N:P [mol N (mol P)^-^1]')
set(gca,'FontSize',20)
% C:Chl
figure
z = (2.*1e-3.*12.011.* QuotaPop(:,:,1)) ./ (PopulationSolution.BOF_coefs(:,:,12) + PopulationSolution.BOF_coefs(:,:,13)); 
imagesc(x,y,z)
hc = colorbar
colormap(cmap)
xlabel('Latitude')
ylabel('Depth')
ylabel(hc,'C:Chl [g C (g dv-Chl)^-^1]')
set(gca,'FontSize',20)


%% Enthalpy
% ecotypes
fig = figure
for a = 1:numel(ecotypeList)
    subplot(5,1,a)
    imagesc(x,y,EnthalpyEco(:,:,a))
    set(gca,'clim',[75 85])
    hc = colorbar
    ylabel(hc,strcat(ecotypeList2{a},' \DeltaH_c (KJ gDW^-^1)'))
    colormap('jet')
    xlabel('Latitude')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

Enthalpy2 = Enthalpy
Enthalpy2(find(Enthalpy2==0)) = NaN
nanmin(nanmin(nanmin(Enthalpy2)))
nanmax(nanmax(nanmax(Enthalpy2)))


%% Carb vs lipid
z = (EcotypeSolution.HLI.BOF_coefs(:,:,2) ) ./ EcotypeSolution.HLI.BOF_coefs(:,:,3);

figure
for a = 1:numel(ecotypeList)
    subplot(5,1,a)
    z = (EcotypeSolution.(ecotypeList{a}).BOF_coefs(:,:,2) ) ./ EcotypeSolution.(ecotypeList{a}).BOF_coefs(:,:,3);
    imagesc(x,y,z)
    colorbar
    colormap('jet')
    xlabel('Latitude')
    ylabel('Depth')
    set(gca,'FontSize',20)
    
end


%% C:N and C:P contours

figure
subplot(2,1,1)
z = QuotaPop(:,:,1) ./ QuotaPop(:,:,3);
imagesc(x,y,z)
hc=colorbar
ylabel(hc,'C:N [mol C (mol N)^-^1]')
colormap('jet')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)
subplot(2,1,2)
z = QuotaPop(:,:,1) ./ QuotaPop(:,:,5);
imagesc(x,y,z)
hc=colorbar
ylabel(hc,'C:P [mol C (mol P)^-^1]')
colormap('jet')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)



%% Cell Size

figure
z = (4/3) .* pi() .* PopulationSolution.r_opt.^3 .* 280; % fg C cell-1
imagesc(x,y,z)
hc=colorbar
ylabel(hc,'Cell size (fg C cell^-^1')
cmap = [1 1 1; jet(256)];
colormap(cmap)
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)


%% Chl b/a
ecotypeColors = varycolor(numel(ecotypeList));
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/envData/IrrDat.mat');
aPeak_lambda = 450;
bPeak_lambda = 474;
aPeak_ind = find(IrrDat.Lambda == aPeak_lambda);
bPeak_ind = find(IrrDat.Lambda == bPeak_lambda);
for i = 1:Gridding.nStations
    Ed_ratio(:,i) = IrrDat.Data{Gridding.stationsVec2(i)}(:,bPeak_ind) ./ IrrDat.Data{Gridding.stationsVec2(i)}(:,aPeak_ind);
end

figure
subplot(2,1,1)
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x1 =  Ed_ratio;
        y1 = StrainSolution.BOF_coefs(:,:,13,strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
        plot(x1,y1,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x1 = Ed_ratio;
    y1 = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
    h{i} = plot(x1,y1,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
x1 =  Ed_ratio;
y1 = PopulationSolution.BOF_coefs(:,:,13) ./ PopulationSolution.BOF_coefs(:,:,12);
h{6} = plot(x1,y1,'.k','MarkerSize',15)
ylabel('Divinylchlorophyll b/a ratio')
xlabel('E_d(474) / E_d(450)')
xlim([1 3])
set(gca,'FontSize',20)
legend([h{1}(1) h{2}(1) h{3}(1) h{4}(1) h{5}(1) h{6}(a)],[ecotypeList2,{'Population'}])

subplot(2,1,2)
z = PopulationSolution.BOF_coefs(:,:,13) ./ PopulationSolution.BOF_coefs(:,:,12)
imagesc(x,y,z)
hc=colorbar
ylabel(hc,'Divinylchlorophyll b/a ratio')
colormap('jet')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)






