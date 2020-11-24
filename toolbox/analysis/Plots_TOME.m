%% Plots of temperature and TOME

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

%% Moore et al., fits to Arrhenius

% data
Tex = [11 13 14 15 17 18 19 21 22 23 24 25 26 27 28 29 30 31 33];
MED4_mu = [0.09 0.2 0.24 0.23 0.3 0.3 0.32 0.35 0.37 0.4 0.39 0.4 0.38 0.34 0.02 0 0 0 0];
MIT9515_mu = [0.1 0.19 0.21 0.24 0.25 0.27 0.3 0.34 0.39 0.42 0.46 0.47 0.43 0.28 0.01 0 0 0 0];
MIT9215_mu = [0 0 0 0.01 0.16 0.17 0.22 0.29 0.36 0.4 0.44 0.5 0.43 0.47 0.43 0.43 0.45 0.43 0.01];
MIT9312_mu = [0 0 0 0.01 0.28 0.29 0.32 0.36 0.49 0.54 0.58 0.62 0.58 0.61 0.47 0.33 0.01 0 0];
NATL2A_mu = [0 0.01 0.18 0.21 0.27 0.29 0.31 0.36 0.43 0.47 0.48 0.47 0.42 0.36 0.01 0 0 0 0];
MIT9313_mu = [0 0 0 0 0 0 0 0 0.37 0.39 0.43 0.44 0.45 0.47 0.38 0.01 0 0 0];
MED4_mu = [0.09 0.2 0.24 0.23 0.3 0.3 0.32 NaN 0.37 0.4 0.39 0.4 0.38 0.34 0.02 NaN NaN NaN NaN];
MIT9515_mu = [0.1 0.19 0.21 0.24 NaN 0.27 0.3 NaN NaN NaN 0.46 0.47 0.43 0.28 0.01 NaN NaN NaN NaN];
MIT9215_mu = [NaN NaN NaN 0.01 0.16 0.17 0.22 NaN 0.36 0.4 0.44 0.5 0.43 0.47 0.43 0.43 0.45 0.43 0.01];
MIT9312_mu = [NaN NaN NaN 0.01 0.28 0.29 0.32 0.36 0.49 0.54 0.58 0.62 0.58 0.61 0.47 0.33 0.01 NaN NaN];
NATL2A_mu = [NaN 0.01 0.18 NaN 0.27 NaN 0.31 NaN 0.43 0.47 0.48 0.47 0.42 0.36 0.01 NaN NaN NaN NaN];
MIT9313_mu = [NaN NaN NaN NaN NaN NaN NaN NaN 0.37 0.39 0.43 0.44 0.45 0.47 0.38 0.01 NaN NaN NaN];

% consolidate
muDat = [Tex; MED4_mu; MIT9515_mu; MIT9215_mu; MIT9312_mu; NATL2A_mu; MIT9313_mu];
rngDat = [10 29; 13 27; 17 31; 17 29; 17 29; 14 27; 14 28]

% average by ecotype and consolidate again
HLI_mu = nanmean(muDat([2 3],:));
HLII_mu = nanmean(muDat([4 5],:));
LLI_mu = muDat(6,:);
LLII_III_mu = nanmean(muDat([6 7],:));
LLIV_mu = muDat(7,:);
muDatEcotype = [Tex; HLI_mu; HLII_mu; LLI_mu; LLII_III_mu; LLIV_mu];
rngDatEcotype = [10 29; 17 31; 14 27; 14 29; 14 28]

% Let's lump all the data together, and get rid of the data that falls outside the linear domain.
% We need to eyeball this...
% First, get rid of the low values... from the plot this looks like all
% values log(mu) < -2.5
muDatEcotype(find(log(muDatEcotype) < -2)) = NaN;
% Now let's set the cutoff for the upper limit to be 1/T > 0.04;
muDatEcotype2 = muDatEcotype;
muDatEcotype2(:,find(1./Tex < 0.038)) = NaN;

Tex2 = repmat(Tex,1,5) + 273.15;
muDatEcotype3 = reshape(muDatEcotype2(2:end,:)',5.*size(muDatEcotype2,2),1);

% Cool now let's get the regression parameters
mdl1 = LinearModel.fit(1./Tex2,log(muDatEcotype3));
A = exp(mdl1.Coefficients.Estimate(1)); 
R = 8.314;
Ea = -mdl1.Coefficients.Estimate(2).*R; % J mol-1

fig = figure
plot(Tex2(1:numel(Tex)), A.*exp(-Ea./(R.*(Tex2(1:numel(Tex))))),'-k','LineWidth',3)
hold on
plot(Tex2, muDatEcotype3,'.k','MarkerSize',20);
xlabel('Temperature (K)')
ylabel('\mu [d^-^1]');
xlim([285 300])
set(gca,'FontSize',20)
legend('Arrhenius (E_a = 53 KJ mol^-^1)','Data','Location','NorthWest');

%fileName = '/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/Temperature_TOME/Arrenius_fit.eps';
%saveas(fig,fileName,'epsc')

%% TOME results
OGTDat = readtable('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/db/OGTDat.csv','Delimiter',',','ReadVariableNames',true);
T_vec = linspace(10,35,100);
StrainName = 'MIT9312';
strInd = find(strcmp(StrainName,OGTDat.Strain));
OGT = 273.15 + OGTDat.OGT_old(strInd); % K

Ea = 5.2733e4; %J mol-1
R = 8.314;
A = 1;
OGT_rate = exp(-Ea./(R.*(OGT)));
T = 273.15 + T_vec
InSitu_rate = exp(-Ea./(R.*(T))).*(1-exp(T-(OGT+2)));

% Lower bound linear correction (optional, comment out for just the
% Arrhenius function). Currently using 10 degrees less than OGT as the
% upper bound of this linear function, and 15 degrees less  than OGT as the
% lower bound (zero rate).
for a = 1:numel(T)
    if T(a) < OGT - 10
        OGT_minus10_rate = exp(-Ea./(R.*(OGT-10)));
        LB_slope = OGT_minus10_rate / 5;
        dT = OGT - 10 - T(a);
        InSitu_rate(a) = InSitu_rate(a) - dT.*LB_slope;
    end
    % Check that we have a positive rate
    if InSitu_rate(a) < 0
        InSitu_rate(a) = 0;
    end
end



% temperature correction
TCorr = InSitu_rate ./ OGT_rate;

% store each strain separately
TCorr1 = TCorr;
TCorr2 = TCorr;


%% Violins TOME by ecotype

for a = 1:numel(ecotypeList)
    strList = Gridding.strNameVec(strEco_idx{a})
    for b = 1:numel(strList)
        strIdx = find(strcmp(strList{b},OGTDat.Strain));
        strOGT{a}(b) = OGTDat.OGT_old(strIdx(1));
    end
end

% violin plot
figure
violin(strOGT)
set(gca,'XTickLabels',ecotypeList2)
ylabel('OGT (degrees C)')
set(gca,'FontSize',20)


% Plot against observed data
fig = figure
subplot(1,3,1)
plot(Tex2(1:numel(Tex)), A.*exp(-Ea./(R.*(Tex2(1:numel(Tex))))),'-k','LineWidth',3)
hold on
plot(Tex2, muDatEcotype3,'.k','MarkerSize',20);
xlabel('Temperature (C)')
ylabel('Growth rate [d^-^1]');
xlim([285 300])
set(gca,'FontSize',20)
legend('Arrhenius (E_a = 53 KJ mol^-^1)','Data','Location','NorthWest');
subplot(1,3,2)
violin(strOGT)
set(gca,'XTickLabels',ecotypeList2)
ylabel('OGT (degrees C)')
set(gca,'FontSize',20)
subplot(1,3,3)
h1 = plot(T,TCorr1./max(TCorr1),'-r','LineWidth',3)
hold on
h2 = plot(T,TCorr2./max(TCorr2),'-k','LineWidth',3)
h3 = plot(Tex+273.15,MED4_mu ./ nanmax(MED4_mu),'.r','MarkerSize',20)
h4 = plot(Tex+273.15,MIT9312_mu ./ nanmax(MIT9312_mu),'.k','MarkerSize',20)
xlabel('Temperature [C]')
ylabel('Relative growth rate')
set(gca,'FontSize',20)
legend('MED4 Predicted','MIT-9312 Predicted','MED4 Observed','MIT9312 Observed','Location','NorthWest')


