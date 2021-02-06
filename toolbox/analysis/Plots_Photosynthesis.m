%% Plots photosynthesis
% Selected plots of fluxes for strains, ecotypes,
% or full populations. 

% QY, PvsI, PBmax, alpha, Ek, C:chl

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

for a = 1:Gridding.nStr
    strainInd = find(strcmp(Gridding.strNameVec{a},orgDatabase.StrainName));
    ecotype{a} = orgDatabase.Ecotype{strainInd};
    ecotypeInd(a) = find(strcmp(ecotype{a},ecotypeList));
end

for a = 1:numel(ecotypeList)
    strEco_idx{a} = find(ecotypeInd == a);
end
ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];

% Color palate for ecotypes
ecotypeColors = varycolor(numel(ecotypeList));


% get indices for Rubisco, pig absorption, CEF, Mehler, PSI (R09542pc),
% PSII (R09503), and PR (R01334)
rnxNames = [{'R00024'},{'DVChla_abs'},{'DVChlb_abs'},{'aCar_abs'},{'Zeax_DVChla_abs'}, {'Zeax_DVChlb_abs'},{'CEF'},{'MehlerPSI'},{'R09503'},{'R09542pc'},{'R01334'}];
for a = 1:numel(rnxNames)
    rxn_idx(a) = find(strcmp(rnxNames{a},FullSolution.PanGEM.rxns));
end

%% QY
fig=figure('position',[0 0 1000 500])
for a = 1:numel(ecotypeList)
    for b = 1:numel(strEco_idx{a})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{a}(b)) ./ ( StrainSolution.Fluxes(:,:,rxn_idx(2),strEco_idx{a}(b)) + StrainSolution.Fluxes(:,:,rxn_idx(3),strEco_idx{a}(b)) + StrainSolution.Fluxes(:,:,rxn_idx(4),strEco_idx{a}(b)) ); 
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',5);
        hold on
    end
end
for a = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(1)) ./ ( EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(2)) + EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(3)) + EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(4)) );
    h{a} = plot(x,y,'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15)
    hold on
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ ( PopulationSolution.Fluxes(:,:,rxn_idx(2))  + PopulationSolution.Fluxes(:,:,rxn_idx(3))  + PopulationSolution.Fluxes(:,:,rxn_idx(4))  ); 
   
h{6} = plot(x,y,'.k','MarkerSize',15)
set(gca,'XScale','log')
ylabel('Quantum Yield [mol C (mol photons)^-^1]')
xlabel('PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)
legend([h{1}(1), h{2}(1),h{3}(1),h{4}(1),h{5}(1), h{6}(1)],[ecotypeList2,{'Population'}])


%% PvsI

fig=figure('position',[0 0 1000 500])
for a = 1:numel(ecotypeList)
    for b = 1:numel(strEco_idx{a})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{a}(b)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{a}(b));
        plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',5);
        hold on
    end
end
for a = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(1)) ./ EcotypeSolution.(ecotypeList{a}).BOF_coefs(:,:,12);
    h{a} = plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15);
    hold on        
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ PopulationSolution.BOF_coefs(:,:,12);
h{6} = plot(x,y.*12.011.*(1/1000),'.k','MarkerSize',15)

set(gca,'XScale','log')
%set(gca,'YScale','log')
ylabel('P^B [mg C (mg Chl a)^-^1 h^-^1]')
xlabel('PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)
%legend([h{1}(1), h{2}(1),h{3}(1),h{4}(1),h{5}(1), h{6}(1)],[ecotypeList2,{'Population'}])


%% PvsI params
% Calculate PBm, alphaB, and Ek
% PB = PBm * tanh(alphaB*E/PBm)
% Ek = PBm/alphaB
% x2(1) is PBm, x2(2) is alphaB

for a = 1:numel(ecotypeList)
    for b = 1:numel(strEco_idx{a})
        for c = 1:Gridding.nStations
            x = CruiseData.PAR(Gridding.stationsVec2(c),:)';
            y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{a}(b)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{a}(b));
            x0 = [2 0.05];
            Ed = x;
            PB = y(:,c).*12.011.*(1/1000);
            PB(find(isnan(PB))) = 0;
            PvsIfun = @(x2)(x2(1).*tanh((x2(2).*Ed)./x2(1)))-PB;
            yPred = lsqnonlin(PvsIfun,x0);
            PBm{a}(b,c) = abs(yPred(1));
            alphaB{a}(b,c) = abs(yPred(2));
            Ek{a}(b,c) = abs(yPred(1)) ./ abs(yPred(2));
        end
    end
end

% Ek violins
for a = 1:numel(ecotypeList)
    Ek2(:,a) = nanmean(Ek{a},1);
    PBm2(:,a) = nanmean(PBm{a},1);
    alphaB2(:,a) = nanmean(alphaB{a},1);
end

fig=figure('position',[0 0 300 1000])
subplot(2,1,1) % PBm
h = boxplot(PBm2,'Notch','on','Labels',ecotypeList2,'Whisker',1,'Colors',ecotypeColors)
ylabel('P^B_m [mg C (mg Chl a)^-^1 h^-^1]')
set(h,{'linew'},{2})
set(gca,'XTick',[])
set(gca,'FontSize',20)
%ylim([0 350])
subplot(2,1,2) % Ek
h = boxplot(Ek2,'Notch','on','Labels',ecotypeList2,'Whisker',1,'Colors',ecotypeColors)
ylabel('E_k [\mu mol quanta m^-^2 s^-^1]')
set(h,{'linew'},{2})
set(gca,'FontSize',20)
set(gca,'XTick',[])
ylim([0 350])


%% Spectral irradiance and absorption
% load irradiance data
load(FullSolution.FileNames.IrrDat_fileName);
[IrrDat2] = standardizeIrr(IrrDat,Gridding.lambdaVec,Gridding.depthVec); %mmoles photons m-2 h-1 bandwidth*nm-1
IrrDat3 = reshape([IrrDat2{:}],numel(Gridding.depthVec),numel(Gridding.lambdaVec),numel(IrrDat2));

% load pigment absorption data
PigDB = readtable(FullSolution.FileNames.PigDB_fileName,'ReadVariableNames',true);

fig=figure
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
imagesc(Gridding.lambdaVec,Gridding.depthVec,2*IrrDat3(:,:,35))
xlim([400 700])
hc = colorbar
xlabel('Wavelength [nm]')
ylabel('Depth [m]')
colormap('jet')
hc = colorbar
ylabel(hc, 'Irradiance [\mumol photons m^-^2 s^-^1 nm^-^1]')
set(gca,'FontSize',20)


%% All together
fig = figure

subplot(2,4,1)
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

subplot(2,4,5)
imagesc(Gridding.lambdaVec,Gridding.depthVec,2*IrrDat3(:,:,35))
xlim([400 700])
hc = colorbar
xlabel('Wavelength [nm]')
ylabel('Depth [m]')
colormap('jet')
hc = colorbar
ylabel(hc, 'Irradiance [\mumol photons m^-^2 s^-^1 nm^-^1]')
set(gca,'FontSize',20)

subplot(2,4,2:3)
for a = 1:numel(ecotypeList)
    for b = 1:numel(strEco_idx{a})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{a}(b)) ./ ( StrainSolution.Fluxes(:,:,rxn_idx(2),strEco_idx{a}(b)) + StrainSolution.Fluxes(:,:,rxn_idx(3),strEco_idx{a}(b)) + StrainSolution.Fluxes(:,:,rxn_idx(4),strEco_idx{a}(b)) ); 
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',5);
        hold on
    end
end
for a = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(1)) ./ ( EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(2)) + EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(3)) + EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(4)) );
    h{a} = plot(x,y,'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15)
    hold on
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ ( PopulationSolution.Fluxes(:,:,rxn_idx(2))  + PopulationSolution.Fluxes(:,:,rxn_idx(3))  + PopulationSolution.Fluxes(:,:,rxn_idx(4))  ); 
   
h{6} = plot(x,y,'.k','MarkerSize',15)
set(gca,'XScale','log')
ylabel('Quantum Yield [mol C (mol photons)^-^1]')
%xlabel('PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'XTick',[])
set(gca,'FontSize',20)
legend([h{1}(1), h{2}(1),h{3}(1),h{4}(1),h{5}(1), h{6}(1)],[ecotypeList2,{'Population'}])

subplot(2,4,6:7)
for a = 1:numel(ecotypeList)
    for b = 1:numel(strEco_idx{a})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{a}(b)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{a}(b));
        plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',5);
        hold on
    end
end
for a = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(1)) ./ EcotypeSolution.(ecotypeList{a}).BOF_coefs(:,:,12);
    h{a} = plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15);
    hold on        
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ PopulationSolution.BOF_coefs(:,:,12);
h{6} = plot(x,y.*12.011.*(1/1000),'.k','MarkerSize',15)

set(gca,'XScale','log')
%set(gca,'YScale','log')
ylabel('P^B [mg C (mg Chl a)^-^1 h^-^1]')
xlabel('PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)


subplot(2,4,4) % PBm
hb = boxplot(PBm2,'Notch','on','Labels',ecotypeList2,'Whisker',1,'Colors',ecotypeColors)
ylabel('P^B_m [mg C (mg Chl a)^-^1 h^-^1]')
set(hb,{'linew'},{2})
set(gca,'XTick',[])
set(gca,'FontSize',20)
%ylim([0 350])
subplot(2,4,8) % Ek
hb = boxplot(Ek2,'Notch','on','Labels',ecotypeList2,'Whisker',1,'Colors',ecotypeColors)
ylabel('E_k [\mu mol quanta m^-^2 s^-^1]')
set(hb,{'linew'},{2})
set(gca,'FontSize',20)
%set(gca,'XTick',[])
ylim([0 350])


%% QY, PvsI, PmB, Ek

fig = figure

subplot(2,3,1:2)
for a = 1:numel(ecotypeList)
    for b = 1:numel(strEco_idx{a})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{a}(b)) ./ ( StrainSolution.Fluxes(:,:,rxn_idx(2),strEco_idx{a}(b)) + StrainSolution.Fluxes(:,:,rxn_idx(3),strEco_idx{a}(b)) + StrainSolution.Fluxes(:,:,rxn_idx(4),strEco_idx{a}(b)) ); 
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',5);
        hold on
    end
end
for a = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(1)) ./ ( EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(2)) + EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(3)) + EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(4)) );
    h{a} = plot(x,y,'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15)
    hold on
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ ( PopulationSolution.Fluxes(:,:,rxn_idx(2))  + PopulationSolution.Fluxes(:,:,rxn_idx(3))  + PopulationSolution.Fluxes(:,:,rxn_idx(4))  ); 
   
h{6} = plot(x,y,'.k','MarkerSize',15)
set(gca,'XScale','log')
ylabel('Quantum Yield [mol C (mol photons)^-^1]')
%xlabel('PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'XTick',[])
set(gca,'FontSize',20)
legend([h{1}(1), h{2}(1),h{3}(1),h{4}(1),h{5}(1), h{6}(1)],[ecotypeList2,{'Population'}])

subplot(2,3,4:5)
for a = 1:numel(ecotypeList)
    for b = 1:numel(strEco_idx{a})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{a}(b)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{a}(b));
        plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',5);
        hold on
    end
end
for a = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxn_idx(1)) ./ EcotypeSolution.(ecotypeList{a}).BOF_coefs(:,:,12);
    h{a} = plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15);
    hold on        
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ PopulationSolution.BOF_coefs(:,:,12);
h{6} = plot(x,y.*12.011.*(1/1000),'.k','MarkerSize',15)

set(gca,'XScale','log')
%set(gca,'YScale','log')
ylabel('P^B [mg C (mg Chl a)^-^1 h^-^1]')
xlabel('PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)


subplot(2,3,3) % PBm
hb = boxplot(PBm2,'Notch','on','Labels',ecotypeList2,'Whisker',1,'Colors',ecotypeColors)
ylabel('P^B_m [mg C (mg Chl a)^-^1 h^-^1]')
set(hb,{'linew'},{2})
set(gca,'XTick',[])
set(gca,'FontSize',20)
%ylim([0 350])
subplot(2,3,6) % Ek
hb = boxplot(Ek2,'Notch','on','Labels',ecotypeList2,'Whisker',1,'Colors',ecotypeColors)
ylabel('E_k [\mu mol quanta m^-^2 s^-^1]')
set(hb,{'linew'},{2})
set(gca,'FontSize',20)
%set(gca,'XTick',[])
ylim([0 350])





%% Figures for the cartoon

% irr and abs
figure
y1 = 1.2*( squeeze(IrrDat3(15,:,31)) - min(squeeze(IrrDat3(15,:,31))) ) ./ ( max(squeeze(IrrDat3(15,:,31))) - min(squeeze(IrrDat3(15,:,31))) );
y2 = ( PigDB.Divinylchlorophyll_a - min(PigDB.Divinylchlorophyll_a)) ./ ( max(PigDB.Divinylchlorophyll_a) - min(PigDB.Divinylchlorophyll_a));
y3 = ( PigDB.Divinylchlorophyll_b - min(PigDB.Divinylchlorophyll_b)) ./ ( max(PigDB.Divinylchlorophyll_b) - min(PigDB.Divinylchlorophyll_b));
y4 = ( PigDB.alpha_Carotene - min(PigDB.alpha_Carotene)) ./ ( max(PigDB.alpha_Carotene) - min(PigDB.alpha_Carotene));
y5 = ( PigDB.Zeaxanthin - min(PigDB.Zeaxanthin)) ./ ( max(PigDB.Zeaxanthin) - min(PigDB.Zeaxanthin));
h1 = plot(Gridding.lambdaVec,y1,'-k','LineWidth',3);
hold on
h2 = plot(PigDB.lambda,y2,'-g','LineWidth',3);
h3 = plot(PigDB.lambda,y3,'-b','LineWidth',3);
h4 = plot(PigDB.lambda,y4,'-r','LineWidth',3);
h5 = plot(PigDB.lambda,y5,'-c','LineWidth',3);

set(gca,'YTick',[]);
xlabel('Wavelength [nm]')
set(gca,'FontSize',20);
legend('Irradiance','Abs. DV-Chl a','Abs. DV-Chl b','Abs. \alpha-Car','Abs. Zeax')



