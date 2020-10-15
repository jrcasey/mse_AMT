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




