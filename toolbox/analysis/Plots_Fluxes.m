%% Plots Fluxes
% Selected plots of fluxes for strains, ecotypes,
% or full populations. 
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

% Color palate for ecotypes
ecotypeColors = varycolor(numel(ecotypeList));

%% Primary production

% 14C
fileName = 'CBIOMES/Data/Environmental_Data/Cruises/AMT13/PP_14C.csv';
PP_dat = readtable(fileName,'Delimiter',',','ReadVariableNames',true);
% PP_14C is in units of mmol C m-3 d-1, with roughly 12 hour incubations.
PP_stations = unique(PP_dat.Station);
for a = 1:numel(PP_stations)
    stIdx(a) = find(strcmp(CruiseData.Stations,PP_stations{a}));
    PP_lat(a) = CruiseData.Lat(stIdx(a));
    PP_lon(a) = CruiseData.Lon(stIdx(a));
    PP_stIdx{a} = find(strcmp(PP_dat.Station,PP_stations{a}));
end

% Extrapolate (same value at 55% PAR and surface, zero at 200 m)
for a = 1:numel(PP_stations)
    clear PPvec PPvec2 zVec zVec2
    zVec = PP_dat.Depth(PP_stIdx{a});
    zVec2 = [0 zVec' max(zVec)+10 200];
    PPvec = PP_dat.Total(PP_stIdx{a});
    PPvec2 = [PPvec(1) PPvec' 0 0];
    PP_14C_int(a) = trapz(zVec2,PPvec2);
end

% Winkler GPP
Winkler = CruiseData.GP;
Winkler(find(isnan(Winkler))) = 0;
Winkler(find(Winkler<0)) = 0;
for a = 1:size(Winkler,1)
    Winkler_GPP_int(a) = trapz(Gridding.depthVec,Winkler(a,:));
end
Winkler_GPP_int(find(Winkler_GPP_int == 0)) = NaN;

% Model NPP
CO2EX_idx = find(strcmp('CO2EX',FullSolution.PanGEM.rxns));
HCO3EX_idx = find(strcmp('HCO3EX',FullSolution.PanGEM.rxns));
NPP = PopulationSolution.Fluxes(:,:,CO2EX_idx) + PopulationSolution.Fluxes(:,:,HCO3EX_idx); % mmol ml-1 h-1
NPP2 = NPP .* 1e6 .* 12
for a = 1:Gridding.nStations
    Pro_NPPint(a) = trapz(Gridding.depthVec,-NPP2(:,a));
end

% Model GPP
Rubisco_idx = find(strcmp('R00024',FullSolution.PanGEM.rxns));
GPP = PopulationSolution.Fluxes(:,:,Rubisco_idx); % mmol ml-1 h-1
GPP2 = GPP .* 1e6 .* 12; % umol l-1 d-1
% need PQ to calculate O2 units
O2EX_idx = find(strcmp('O2EX',FullSolution.PanGEM.rxns));
PQ = PopulationSolution.Fluxes(:,:,O2EX_idx) ./ (PopulationSolution.Fluxes(:,:,HCO3EX_idx) + PopulationSolution.Fluxes(:,:,CO2EX_idx) );
GPP3 = GPP2 .* -PQ;
GPP3(find(isnan(GPP3))) = 0;
for a = 1:Gridding.nStations
    Pro_GPPint(a) = trapz(Gridding.depthVec,GPP3(:,a));
end


figure
subplot(2,1,1)
h1 = plot(PP_lat,PP_14C_int,'.k','MarkerSize',20);
hold on
h2 = plot(CruiseData.Lat,Winkler_GPP_int,'.r','MarkerSize',20);
h3 = plot(x,Pro_NPPint,'ok')
h4 = plot(x,Pro_GPPint,'or')
xlabel('Latitude')
ylabel('mmol m^-^2 d^-^1')
set(gca,'FontSize',20)
legend('^1^4C-PP (C)','Winkler-GPP (O_2)','Pro-NPP (C)','Pro-GPP (O_2)');

subplot(2,1,2)
imagesc(x,y,1000*GPP3)
set(gca,'clim',[0 500])
colormap('jet')
hc = colorbar
ylabel(hc,'Pro-GPP [nmol O_2 L^-^1 d^-^1]')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)



% delineate chlorophyll regions
Chl2 = CruiseData.Fluor_Chl_a(Gridding.stationsVec2,:)'; % mg m-3
Chl2(find(isnan(Chl2))) = 0;
for a = 1:Gridding.nStations
    Chl_int(a) = trapz(Gridding.depthVec,Chl2(:,a));
end

% let's say 15 mg m-2 is the threshold for low versus high
low_idx = find(Chl_int < 15);
high_idx = find(Chl_int >=15);

lowChl_Pro_GPP = mean(Pro_GPPint(low_idx));
highChl_Pro_GPP = mean(Pro_GPPint(high_idx));
lowChl_Winkler = nanmean(Winkler_GPP_int(Gridding.stationsVec2(low_idx)));
highChl_Winkler = nanmean(Winkler_GPP_int(Gridding.stationsVec2(high_idx)));

lowChl_ProContribution = lowChl_Pro_GPP ./ lowChl_Winkler;
highChl_ProContribution = highChl_Pro_GPP ./ highChl_Winkler;

lowChl_ProContribution2 = nanmean(Pro_GPPint(low_idx) ./ Winkler_GPP_int(Gridding.stationsVec2(low_idx)) )
lowChl_ProContribution2_sd = nanstd(Pro_GPPint(low_idx) ./ Winkler_GPP_int(Gridding.stationsVec2(low_idx)) )

highChl_ProContribution2 = nanmean(Pro_GPPint(high_idx) ./ Winkler_GPP_int(Gridding.stationsVec2(high_idx)) )
highChl_ProContribution2_sd = nanstd(Pro_GPPint(high_idx) ./ Winkler_GPP_int(Gridding.stationsVec2(high_idx)) )

%% f ratio
nh3_idx = find(strcmp('AmmoniaEX',PanGEM.rxns));
no2_idx = find(strcmp('NitriteEX',PanGEM.rxns));
no3_idx = find(strcmp('NitrateEX',PanGEM.rxns));
z = PopulationSolution.Fluxes(:,:,no3_idx) ./ nansum(PopulationSolution.Fluxes(:,:,[nh3_idx no2_idx no3_idx]),3)
figure
imagesc(x,y,z)
hc = colorbar
ylabel(hc,'f-ratio')
colormap('jet')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)




%% Photosynthesis reactions

rxnIdx1 = find(strcmp('R09503',FullSolution_L2.PanGEM.rxns));
rxnIdx2 = find(strcmp('R09542pc',FullSolution_L2.PanGEM.rxns));
fig=figure
for a = 1:numel(ecotypeList)
    subplot(2,3,a)
    z = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxnIdx1);
    %z = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxnIdx1) ./ EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxnIdx2)
    imagesc(x,y,z)
    xlabel('Latitude')
    ylabel('Depth');
    colorbar
    colormap('jet')
    set(gca,'FontSize',20)
end

rxnIdx1 = find(strcmp('LightTRANS',FullSolution_L2.PanGEM.rxns));
fig = figure
for a = 1:numel(ecotypeList)
    z = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxnIdx1);
    z = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxnIdx1) ./ EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,rxnIdx2)
    plot(CruiseData.PAR(Gridding.stationsVec2,:)',z,'.','MarkerSize',15,'MarkerFaceColor',ecotypeColors(a,:),'MarkerEdgeColor',ecotypeColors(a,:));
    hold on
    
end
set(gca,'XScale','log')

%% alpha-car absorption fraction
for a = 1:Gridding.nStr
    aCar_fraction(:,:,a) = FullSolution_L2.(Gridding.strNameVec{a}).pigAbs(:,:,3) ./ (FullSolution_L2.(Gridding.strNameVec{a}).pigAbs(:,:,1) + FullSolution_L2.(Gridding.strNameVec{a}).pigAbs(:,:,2) + FullSolution_L2.(Gridding.strNameVec{a}).pigAbs(:,:,3));
end

for a = 1:numel(strEco_idx)
    aCar_fraction_ecoMean(:,:,a) = nanmean(aCar_fraction(:,:,strEco_idx{a}),3);
end



