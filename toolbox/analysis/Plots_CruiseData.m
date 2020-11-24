%% Plot Cruise Data
% Contour plots and line plots of various cruise observations

%% Load data
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;

%% Parse out Gridding, CruiseData, FileNames, and PanGEM from FullSolution
Gridding = FullSolution.Gridding;
CruiseData = FullSolution.CruiseData;

%% Load variable names and units
varMeta = readtable('data/EnvData/Variables_AMT13.csv','ReadVariableNames',true,'Delimiter',',');

%% Gridded domain
x = CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;


%% Contours of each variable
for a = 1:numel(varMeta.varID)
    var = varMeta.varID{a};
    varName = varMeta.varName{a};
    varUnits = varMeta.varUnits{a};
    saveMe = 1;
    fileName = strcat('/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/CruiseVariables/','AMT13_Contour_',var,'.eps');
    z = CruiseData.(var)(Gridding.stationsVec2,:)';
    contour_AMT(x, y, z, var, varName, varUnits, saveMe, fileName);
end

%% Some specific contours

fig = figure
subplot(2,2,1)
var = 'T';
varName = 'Temperature';
varUnits = 'C';
z = CruiseData.(var)(Gridding.stationsVec2,:)';
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
title(varName)
set(gca,'FontSize',20)

subplot(2,2,2)
var = 'PAR';
varName = 'log PAR';
varUnits = 'mmol photons m^-^2 s^-^1';
z = log10(CruiseData.(var)(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
title(var)
set(gca,'FontSize',20)

subplot(2,2,3)
var = 'DIN';
varName = 'DIN';
varUnits = 'nM';
z = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
title(varName)
set(gca,'FontSize',20)

subplot(2,2,4)
var = 'DIP';
varName = 'DIP';
varUnits = 'nM';
z = CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
title(varName)
set(gca,'FontSize',20)



fig = figure
subplot(2,2,1)
var = 'Ammonia';
varName = 'log_1_0 Ammonia';
varUnits = 'nM';
z = log10(CruiseData.(var)(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
%title(varName)
set(gca,'FontSize',20)

subplot(2,2,2)
var = 'Nitrite';
varName = 'log_1_0 Nitrite';
varUnits = 'nM';
z = log10(CruiseData.(var)(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
%title(varName)
set(gca,'FontSize',20)

subplot(2,2,3)
var = 'Nitrate';
varName = 'log_1_0 Nitrate';
varUnits = 'nM';
z = log10(CruiseData.(var)(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
%title(varName)
set(gca,'FontSize',20)

subplot(2,2,4)
var = 'Orthophosphate';
varName = 'log_1_0 Phosphate';
varUnits = 'nM';
z = log10(CruiseData.(var)(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
%title(varName)
set(gca,'FontSize',20)


%% A selected depth profile

station_idx = 14;
y = Gridding.depthVec;
x1 = CruiseData.Ammonia(Gridding.stationsVec2(station_idx),:)';
x2 = CruiseData.Nitrite(Gridding.stationsVec2(station_idx),:)';
x3 = CruiseData.Nitrate(Gridding.stationsVec2(station_idx),:)';
x4 = CruiseData.Orthophosphate(Gridding.stationsVec2(station_idx),:)';
x5 = CruiseData.PAR(Gridding.stationsVec2(station_idx),:)';
x6 = CruiseData.T(Gridding.stationsVec2(station_idx),:)';

figure
subplot(2,2,1)
h1 = plot(x1,y,'-k','LineWidth',3)
ylabel('Depth [m]')
set(gca,'YDir','reverse')
set(gca,'XTick',[])

subplot(2,2,2)
h2 = plot(x2,y,'-g','LineWidth',3)
hold on
h3 = plot(x3,y,'-b','LineWidth',3)
h4 = plot(x4,y,'-c','LineWidth',3)
set(gca,'YDir','reverse')
set(gca,'XTick',[],'YTick',[])
subplot(2,2,3)
h2 = plot(x5,y,'-y','LineWidth',3)
set(gca,'YDir','reverse')
set(gca,'XTick',[],'YTick',[])
subplot(2,2,4)
h2 = plot(x6,y,'-r','LineWidth',3)
set(gca,'YDir','reverse')
set(gca,'XTick',[],'YTick',[])


%% Spectral irradiance profile

% load irradiance data
load(FullSolution.FileNames.IrrDat_fileName);
[IrrDat2] = standardizeIrr(IrrDat,Gridding.lambdaVec,Gridding.depthVec); %mmoles photons m-2 h-1 bandwidth*nm-1
IrrDat3 = reshape([IrrDat2{:}],numel(Gridding.depthVec),numel(Gridding.lambdaVec),numel(IrrDat2));

% load pigment absorption data
PigDB = readtable(FullSolution.FileNames.PigDB_fileName,'ReadVariableNames',true);


fig = figure
subplot(3,1,1)
h1 = plot(PigDB.lambda,PigDB.Divinylchlorophyll_a,'-g','LineWidth',3);
hold on
h2 = plot(PigDB.lambda,PigDB.Divinylchlorophyll_b,'-b','LineWidth',3);
h3 = plot(PigDB.lambda,PigDB.alpha_Carotene,'-r','LineWidth',3);
h4 = plot(PigDB.lambda,PigDB.Zeaxanthin,'-c','LineWidth',3);
%xlabel('Wavelength [nm]')
ylabel('Absorption [m^2 mg^-^1 nm^-^1]')
legend([h1(1) h2(1) h3(1) h4(1)],'dv-chlorophyll a','dv-chlorophyll b','\alpha-carotene','zeaxanthin')
set(gca,'FontSize',20)
set(gca,'XTick',[])
subplot(3,1,2:3)
imagesc(Gridding.lambdaVec,y,2*IrrDat3(:,:,35))
xlim([400 700])
hc = colorbar
xlabel('Wavelength [nm]')
ylabel('Depth [m]')
colormap('jet')
hc = colorbar
ylabel(hc, 'Irradiance [mmol photons m^-^2 s^-^1 nm^-^1]')
set(gca,'FontSize',20)




%% DIN with T overlay

z1 = CruiseData.Ammonia(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Nitrate(Gridding.stationsVec2,:)'
z2 = CruiseData.T(Gridding.stationsVec2,:)'
T_levels = [10 20 25];

fig = figure
h1 = imagesc(x,y,z1./1000)
hold on
[c,h2] = contour(x,y,z2,T_levels);
set(h2,'LineWidth',3,'ShowText','on','LineColor','w')
clabel(c,h2,'FontSize',20,'Color','w')
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'DIN [\muM]')
set(gca,'FontSize',20)


%% Prochlorococcus with PAR overlay
z1 = CruiseData.Pro(Gridding.stationsVec2,:)'
z2 = log10(CruiseData.PAR(Gridding.stationsVec2,:)')
PAR_levels = [1 2 3];
z2 = CruiseData.PAR(Gridding.stationsVec2,:)'
PAR_levels = [1 10 100 1000];

fig = figure
h1 = imagesc(x,y,z1)
hold on
[c,h2] = contour(x,y,z2,PAR_levels);
set(h2,'LineWidth',3,'ShowText','on','LineColor','w')
clabel(c,h2,'FontSize',20,'Color','w')
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Prochlorococcus [cells ml^-^1]')
set(gca,'FontSize',20)


%% Supplementary section figures

% Physics (4 panel contours)
figure
subplot(2,2,1)
z = CruiseData.T(Gridding.stationsVec2,:)';
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Temperature [degrees C]')
set(gca,'FontSize',20)
subplot(2,2,2)
z = CruiseData.Salinity(Gridding.stationsVec2,:)';
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Salinity [psu]')
set(gca,'FontSize',20)
subplot(2,2,3)
z = 1000 + CruiseData.Density(Gridding.stationsVec2,:)';
imagesc(x,y,z)
set(gca,'clim',[1022 1027])
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Density [kg m^-^3]')
set(gca,'FontSize',20)
subplot(2,2,4)
z = log10(CruiseData.PAR(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'log_1_0 PAR [\mumol photons m^-^2 s^-^1]')
set(gca,'FontSize',20)

% Nutrients
figure
subplot(2,2,1)
z = log10(CruiseData.Ammonia(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'log_1_0 Ammonia [nM]')
set(gca,'FontSize',20)
subplot(2,2,2)
z = log10(CruiseData.Nitrite(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'log_1_0 Nitrite [nM]')
set(gca,'FontSize',20)
subplot(2,2,3)
z = log10(CruiseData.Nitrate(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'log_1_0 Nitrate [nM]')
set(gca,'FontSize',20)
subplot(2,2,4)
z = log10(CruiseData.Orthophosphate(Gridding.stationsVec2,:)');
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'log_1_0 Phosphate [nM]')
set(gca,'FontSize',20)

% Biology
figure
subplot(4,2,1)
z = CruiseData.Fluor_Chl_a(Gridding.stationsVec2,:)';
imagesc(x,y,z)
set(gca,'clim',[0 0.7])
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Chlorophyll a [mg m^-^3]')
set(gca,'FontSize',20)
subplot(4,2,2)
z = CruiseData.Pro(Gridding.stationsVec2,:)';imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Prochlorococcus [cells ml^-^1]')
set(gca,'FontSize',20)
subplot(4,2,3)
z = CruiseData.Syn(Gridding.stationsVec2,:)';imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Synechococcus [cells ml^-^1]')
set(gca,'FontSize',20)
subplot(4,2,4)
z = CruiseData.Peuk(Gridding.stationsVec2,:)';imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Picoeukaryotes [cells ml^-^1]')
set(gca,'FontSize',20)
subplot(4,2,5)
z = CruiseData.Neuk(Gridding.stationsVec2,:)';imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Nanoeukaryotes [cells ml^-^1]')
set(gca,'FontSize',20)
subplot(4,2,6)
z = CruiseData.Cocolithophores(Gridding.stationsVec2,:)';imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Cocolithophores [cells ml^-^1]')
set(gca,'FontSize',20)
subplot(4,2,7)
z = CruiseData.Cyrptophytes(Gridding.stationsVec2,:)';imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Cryptophytes [cells ml^-^1]')
set(gca,'FontSize',20)
subplot(4,2,8)
z = CruiseData.Hbac1(Gridding.stationsVec2,:)';imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, 'Heterotrophic bacteria [cells ml^-^1]')
set(gca,'FontSize',20)









