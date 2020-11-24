%% Plot nutrient uptake and transporters
% Contour plots of S*. 
% Analyze n in terms of n_max and n_opt to see where SA limitation could have occurred. 
% Contour plot of f-ratio


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

%% Contours of S_star for population


figure
subplot(2,1,1)
imagesc(x,y,CruiseData.Ammonia(Gridding.stationsVec2,:)')
set(gca,'clim',[0 60])
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)
hc = colorbar
ylabel(hc,'Ambient [nM]')
colormap('jet')
subplot(2,1,2)
z = PopulationSolution.S_star(:,:,1);
z(find(~z)) = NaN;
imagesc(x,y,1e6*z)
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)
hc = colorbar
ylabel(hc,'S^* [nM]')
colormap('jet')




% histograms of S-star and ambient concentrations for ammonia
s = PopulationSolution.S_star(:,:,1);
s = s(:);
s(find(~s)) = NaN;
s = 1e6.*s; %nM

c = CruiseData.Ammonia(Gridding.stationsVec2,:)'; %nM

binRange = [0:5:100];

hcs = histcounts(s,[binRange Inf]);
hcc = histcounts(c,[binRange Inf]);

figure
hb = bar(binRange,[hcc;hcs]',0.8)
xlabel('Ammonia [nM]')
set(gca,'FontSize',20)
legend('Ambient','S^*')


%% Diffusion limitation

figure
subplot(2,2,1) % ammonia
z = CruiseData.Ammonia(Gridding.stationsVec2,:)' - (1e6*PopulationSolution.S_star(:,:,1));
diffLim_idx = find(z<0);
z2 = zeros(size(z,1),size(z,2));
z2(diffLim_idx) = 1;
imagesc(x,y,z2)
xlabel('Latitude')
ylabel('Depth')
title('Ammonia')
set(gca,'FontSize',20)
colormap('jet')

subplot(2,2,2) % Nitrite
z = CruiseData.Nitrite(Gridding.stationsVec2,:)' - (1e6*PopulationSolution.S_star(:,:,4));
diffLim_idx = find(z<0);
z2 = zeros(size(z,1),size(z,2));
z2(diffLim_idx) = 1;
imagesc(x,y,z2)
xlabel('Latitude')
ylabel('Depth')
title('Nitrite')
set(gca,'FontSize',20)
colormap('jet')

subplot(2,2,3) % Nitrate
z = CruiseData.Nitrate(Gridding.stationsVec2,:)' - (1e6*PopulationSolution.S_star(:,:,3));
diffLim_idx = find(z<0);
z2 = zeros(size(z,1),size(z,2));
z2(diffLim_idx) = 1;
imagesc(x,y,z2)
xlabel('Latitude')
ylabel('Depth')
title('Nitrate')
set(gca,'FontSize',20)
colormap('jet')

subplot(2,2,4) % Orthophosphate
z = CruiseData.Orthophosphate(Gridding.stationsVec2,:)' - (1e6*PopulationSolution.S_star(:,:,2));
diffLim_idx = find(z<0);
z2 = zeros(size(z,1),size(z,2));
z2(diffLim_idx) = 1;
imagesc(x,y,z2)
xlabel('Latitude')
ylabel('Depth')
title('Phosphate')
set(gca,'FontSize',20)
colormap('jet')

%% Plot n opt and n for one strain
TpDat = readtable(FullSolution.FileNames.TpDat_fileName,'ReadVariableNames',true,'Delimiter',',');
S_vec = logspace(-6,-1,100);
f_max = 0.085; % Maximum fraction of the cell surface that can be covered by transporters... wonkiest term in here. See text for a discussion of n_max.

% Nitrate
D = getDiffusivity(TpDat.TpMets_MW(3),25);
kcat = TpDat.kcat(3);
r = 0.35e-6;
A = TpDat.A(3);
mtc = D./r;
% capture probability
alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless
ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3
nG = 49;
n_star_D = 4*pi()*r*D.*S_vec ./ kcat;
n_star_G = ( ks_p + S_vec ) ./ ( (S_vec ./ nG) - (kcat./(D.*r)) )

% get min set
n_max = f_max.*(4.*pi().*(r.^2))./A; % cell-1
n_max2 = repmat(n_max,1,numel(S_vec));
n_star_G(find(n_star_G<0)) = NaN;
n_star_Nitrate = nanmin([n_star_G;n_star_D;n_max2]);


% Nitrite
D = getDiffusivity(TpDat.TpMets_MW(4),25);
kcat = TpDat.kcat(4);
r = 0.35e-6;
A = TpDat.A(4);
mtc = D./r;
% capture probability
alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless
ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3
nG = 49;
n_star_D = 4*pi()*r*D.*S_vec ./ kcat;
n_star_G = ( ks_p + S_vec ) ./ ( (S_vec ./ nG) - (kcat./(D.*r)) )

% get min set
n_max = f_max.*(4.*pi().*(r.^2))./A; % cell-1
n_max2 = repmat(n_max,1,numel(S_vec));
n_star_G(find(n_star_G<0)) = NaN;
n_star_Nitrite = nanmin([n_star_G;n_star_D;n_max2]);

% Ammonia
D = getDiffusivity(TpDat.TpMets_MW(1),25);
kcat = TpDat.kcat(1);
r = 0.35e-6;
A = TpDat.A(1);
mtc = D./r;
% capture probability
alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless
ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3
nG = 80;
n_star_D = 4*pi()*r*D.*S_vec ./ kcat;
n_star_G = ( ks_p + S_vec ) ./ ( (S_vec ./ nG) - (kcat./(D.*r)) )

% get min set
n_max = f_max.*(4.*pi().*(r.^2))./A; % cell-1
n_max2 = repmat(n_max,1,numel(S_vec));
n_star_G(find(n_star_G<0)) = NaN;
n_star_Ammonia = nanmin([n_star_G;n_star_D;n_max2]);

% Orthophosphate
D = getDiffusivity(TpDat.TpMets_MW(2),25);
kcat = TpDat.kcat(2);
r = 0.35e-6;
A = TpDat.A(2);
mtc = D./r;
% capture probability
alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless
ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3
nG = 55;
n_star_D = 4*pi()*r*D.*S_vec ./ kcat;
n_star_G = ( ks_p + S_vec ) ./ ( (S_vec ./ nG) - (kcat./(D.*r)) )

% get min set
n_max = f_max.*(4.*pi().*(r.^2))./A; % cell-1
n_max2 = repmat(n_max,1,numel(S_vec));
n_star_G(find(n_star_G<0)) = NaN;
n_star_Orthophosphate = nanmin([n_star_G;n_star_D;n_max2]);




strName = 'MIT9312'
strIdx = find(strcmp(strName,Gridding.strNameVec));
S_vec2 = logspace(0,5,100); %nM

fig = figure
h1 = plot(S_vec2,n_star_Ammonia,'-g','LineWidth',3);
hold on
h2 = plot(S_vec2,n_star_Nitrite,'-k','LineWidth',3);
h3 = plot(S_vec2,n_star_Nitrate,'-r','LineWidth',3);
h4 = plot(S_vec2,n_star_Orthophosphate,'-b','LineWidth',3);
plot(CruiseData.Ammonia(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,2,strIdx),'.g','MarkerSize',20);
plot(CruiseData.Nitrite(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,5,strIdx),'.k','MarkerSize',20);
plot(CruiseData.Nitrate(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,4,strIdx),'.r','MarkerSize',20);
plot(CruiseData.Orthophosphate(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,3,strIdx),'.b','MarkerSize',20);
set(gca,'XScale','log')
xlim([0 100000])
xlabel('Nutrient concentration [nM]')
ylabel('Transporters per cell')
legend([h1(1) h2(1) h3(1) h4(1)],'Ammonia','Nitrite','Nitrate','Phosphate','Location','NorthEast')
set(gca,'FontSize',20)

% contour of vmax
figure
z = 1e15 .* 3600 .* (PopulationSolution.TpOpt(:,:,4) ./ CruiseData.Pro(Gridding.stationsVec2,:)') .* TpDat.kcat(4);
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)
colormap('jet')
hc = colorbar
ylabel(hc,'v_m_a_x [fmol cell^-^1 h^-^1]')
set(gca,'clim',[0 0.25])


% just ammonia n

D = getDiffusivity(TpDat.TpMets_MW(1),25);
kcat = TpDat.kcat(1);
r = 0.35e-6;
A = TpDat.A(1);
mtc = D./r;
% capture probability
alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless
ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3
nG = 80;
n_star_D = 4*pi()*r*D.*S_vec ./ kcat;
n_star_G = ( ks_p + S_vec ) ./ ( (S_vec ./ nG) - (kcat./(D.*r)) )

% get min set
n_max = f_max.*(4.*pi().*(r.^2))./A; % cell-1
n_max2 = repmat(n_max,1,numel(S_vec));
n_star_G(find(n_star_G<0)) = NaN;
n_star_Ammonia = nanmin([n_star_G;n_star_D;n_max2]);

% histograms of S-star and ambient concentrations for ammonia

fig = figure % MIT9312
strName = 'MIT9312'
strIdx = find(strcmp(strName,Gridding.strNameVec));

subplot(2,1,1)
s = StrainSolution.S_star(:,:,1,strIdx);
s = s(:);
%s(find(~s)) = NaN;
s = 1e6.*s; %nM
c = CruiseData.Ammonia(Gridding.stationsVec2,:)'; %nM
binRange = [0:5:100];
hcs = histcounts(s,[binRange Inf]);
hcc = histcounts(c,[binRange Inf]);

hb = bar(binRange,[hcc;hcs]',0.8)
xlabel('Ammonium [nM]')
set(gca,'FontSize',20)
legend('Ambient','S^*')
subplot(2,1,2)
h1 = plot(S_vec2,n_star_Ammonia,'-k','LineWidth',3);
hold on
plot(CruiseData.Ammonia(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,2,strIdx),'.k','MarkerSize',20);
xlim([0 100])
xlabel('Ammonium [nM]')
ylabel('Transporters per cell')
set(gca,'FontSize',20)
%legend([h1(1) h2(1) h3(1)],'Ammonium','Nitrite','Nitrate')


fig = figure % NATL1A
strName = 'NATL1A'
strIdx = find(strcmp(strName,Gridding.strNameVec));

subplot(2,1,1)
s = StrainSolution.S_star(:,:,1,strIdx);
s = s(:);
%s(find(~s)) = NaN;
s = 1e6.*s; %nM
c = CruiseData.Ammonia(Gridding.stationsVec2,:)'; %nM
binRange = [0:5:100];
hcs = histcounts(s,[binRange Inf]);
hcc = histcounts(c,[binRange Inf]);

hb = bar(binRange,[hcc;hcs]',0.8)
xlabel('Ammonium [nM]')
set(gca,'FontSize',20)
legend('Ambient','S^*')
subplot(2,1,2)
h1 = plot(S_vec2,n_star_Ammonia,'-k','LineWidth',3);
hold on
h2 = plot(S_vec2,n_star_Nitrite,'-r','LineWidth',3);
plot(CruiseData.Ammonia(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,2,strIdx),'.k','MarkerSize',20);
plot(CruiseData.Nitrite(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,5,strIdx),'.r','MarkerSize',20);
%strName = 'SB'
%strIdx = find(strcmp(strName,Gridding.strNameVec));
%plot(CruiseData.Ammonia(Gridding.stationsVec2,:)',StrainSolution.TpOpt(:,:,2,strIdx),'.g','MarkerSize',20);
%set(gca,'XScale','log')
xlim([0 100])
xlabel('Nutrient concentration [nM]')
ylabel('Transporters per cell')
set(gca,'FontSize',20)
legend([h1(1) h2(1) ],'Ammonium','Nitrite')

%% Loop through each strain, each station, each depth, and each exomet
% load strMod
load(FullSolution.FileNames.StrMod_Path);

typeStrVec = [{'MED4'},{'MIT9312'},{'SS120'},{'MIT9211'},{'MIT9313'}]

for a = 1:numel(typeStrVec)
    strName = typeStrVec{a};
    StrMod1 = StrMod.(strName);
    StrMod1.id = strName;
    StrMod1.description = strcat('Un-acclimated strain model with no specified constraints or compositions - Prochlorococcus strain',{' '},strName);
    StrMod2 = prepEnvMod(StrMod1);
    StrMod2.description = strcat('Un-acclimated strain model initialized with size and composition constraints - Prochlorococcus strain',{' '},strName);
    
    for b = 1:Gridding.nZ
        depth_ind = b;
        for c = 1:Gridding.nStations
            station_ind = Gridding.stationsVec2(c);
            S_env = 1e-6*CruiseData.Ammonia(Gridding.stationsVec2(c),depth_ind); % mol m-3
            
            %% Load Tp data
            % Load TpDat
            TpDat = readtable(FullSolution.FileNames.TpDat_fileName,'ReadVariableNames',true,'Delimiter',',');
            
            % Compile environmental data to retrieve and format TpOpt parameter table
            TpDat = getTpOptParams(CruiseData,StrMod2,TpDat,station_ind,depth_ind);
            
            %% Compile parameters for Pro growth on ammonia
            
            % find VmaxG
            tempSol = solveLP(StrMod2,1);
            TpMet = 'Ammonia';
            TpRxn = 'AmmoniaTRANS';
            TpRxn_idx = find(strcmp(TpRxn,FullSolution.PanGEM.rxns));
            TpRxn_flux = abs(tempSol.x(TpRxn_idx)); % mmol gDW-1 h-1
            ro = StrainSolution.r_opt(depth_ind,c,a); % diameter (um);
            if isnan(ro)
                ro = 0.3135;
            end
            DW = 4/3.*pi().*(ro.^3) .* 280; % fg cell-1
            cell_n = 1e15/DW;
            VmaxG = TpRxn_flux * (1/3600) .* (1/1000) .* (1./cell_n); % mol cell-1 s-1
            nG = VmaxG ./ 1.37e-22;
            r = 1e-6*ro
            %VmaxG = 1.77e-18; % Maximum uptake rate in replete batch (mol cell-1 s-1), derived from FBA
            %nG = 2775; % Transporter abundance in replete batch (cell-1), measured by Schmidt et al., 2016
            % Cell radius (m), measured by Schmidt et al., 2016
            A = 3.91e-17; % PtsG catchment area (m2), measured using modeled protein structure.
            T = CruiseData.T(station_ind,depth_ind); % Temperature (C), measured by Schmidt et al., 2016
            if isnan(T)
                T = CruiseData.T(station_ind,depth_ind-1);
            end
            MW = 17; % Molecular weight of ammonia (Da). Be sure that your molecule is hydrated before calculating MW
            u = 0; % Advective velocity (m s-1)
            
            %% define a vector of substrate concentrations (easier to see in log space)
            
            S_res = 500; % resolution in the S domain
            S = logspace(-6,-2,S_res); % Nutrient concentration range (mole m-3). You may need to adjust this for another transport process
            [junk, S_env_idx] = min(abs(S-S_env));
            %% Compute molecular diffusivity
            
            D = getDiffusivity(MW,T); % Hydrated molecular diffusivity (m2 s-1)
            
            %% Compute transporter in vivo catalytic rate
            
            kcat = VmaxG./nG; % catalytic constant (molecules tp-1 s-1)
            
            %% Compute maximum number of transporters
            
            f_max = 0.085; % Maximum fraction of the cell surface that can be covered by transporters... wonkiest term in here. See text for a discussion of n_max.
            n_max = f_max.*(4.*pi().*(r.^2))./A; % cell-1
            
            %% Get glucose-replete batch-acclimated parameters
            
            [ks_d, ks_p] = getK(VmaxG, kcat, nG, r, A, D, u); % porter limit and diffusive limit of half saturation concentration (mol m-3)
            
            %% Calculate various critical concentrations
            
            % Lower bound to growth limited domain of n-star
            S_G_lb = (nG.*kcat) ./ (D.*r); % mol m-3
            S_G_lb_idx = min(find(S-S_G_lb >= 0)); % nearest concentration index in S
            
            % Calculate S-star. This is the positive root of the polynomial
            % intersection of n_optP and n_optD.
            kprime = (kcat./(r .* D));
            a1 = 1./(kprime.*nG);
            b1 = -2;
            c1 = -ks_p;
            p = [a1 b1 c1];
            y0 = roots(p);
            S_star = y0(1); % mol m-3
            S_star_idx = min(find(S-S_star >= 0)); % nearest concentration index in S
            
            % Lower and upper bounds of surface area limitation domain (true if optimal
            % n excedes n_max
            max_n_opt = (ks_p + S(S_star_idx)) ./ ((S(S_star_idx)./nG) - (kcat./(D.*r)));
            if max_n_opt > n_max
                S_SA_lb = (n_max .* kcat) ./ (r .* D); % mol m-3
                S_SA_lb_idx = min(find(S-S_SA_lb >= 0)); % nearest concentration index in S
                S_SA_ub = ( ((n_max .* kcat)./(D.*r)) + ks_p ) ./ ( (n_max./nG) -1 ); % mol m-3
                S_SA_ub_idx = min(find(S-S_SA_ub >= 0)); % nearest concentration index in S
                SA_transition = true;
            else
                SA_transition = false;
            end
            
            %% Compute optimal number of transporters
            
            % In the growth limit range, such that v = vmax
            n_optG = (ks_p + S) ./ ((S./nG) - (kcat./(D.*r))); % cell-1
            % Constrain to the feasible domain (nans outside)
            n_optG(1:S_G_lb_idx-1) = NaN;
            
            % In the diffusive domain, such that v = vD
            n_optD = (S.*r.*D)./kcat; % cell-1
            
            % In the surface area limited domain, such that v = n_max*kcat*S/(ks_p + ks_d + S),
            % where ks_d is evaluated at n_max
            n_optSA = nan(1,numel(S));
            if SA_transition
                n_optSA(S_SA_lb_idx:S_SA_ub_idx) = n_max; % cell-1
            end
            
            % Minimum of all three sets
            n_opt_minDP = nanmin([n_optG;n_optSA;n_optD]); % cell-1
            
            %% Save results for each exomet
            S_star2(b,c,a) = S_star;
            n_opt2(b,c,a) = n_opt_minDP(S_env_idx);
            n2(b,c,a) = StrainSolution.TpOpt(b,c,2,a);
            r2(b,c,a) = ro;
        end
    end
end




%% Compute uptake for optimal n

for a = 1:numel(S)
    Vmax(a) = n_opt_minDP(a).*kcat; % mole cell-1 s-1
    [ks_d(a), junk] = getK(Vmax(a), kcat, n_opt_minDP(a), r, A, D, u); % mol m-3
    [v(a)] = getUptake(S(a),Vmax(a),ks_d(a), ks_p); % mole cell-1 s-1
end
ks = ks_d + ks_p; % mol m-3
v_P = v.*1e15.*3600; % fmol cell-1 h-1

%% Michaelis-Menten kinetics for glucose-replete acclimated cells, just for comparison

[ks_d3, junk] = getK(VmaxG, kcat, nG, r, A, D, u); % mol m-3
for a = 1:numel(S)
    [v3(a)] = getUptake(S(a),VmaxG,ks_d3, ks_p); % mole cell-1 s-1
end
ks3 = ks_d3 + ks_p; % mol m-3
v_P3 = v3.*1e15.*3600; % mole cell-1 s-1

%% Compute maximum diffusive flux toward the cell
% Based on Berg and Purcell, 1977
vD = 4*pi()*r .* D .* S .* ( (n_optD.*A) ./ (4.*pi().*(r.^2)) ); % mole cell-1 s-1

%% Uptake in the n vs S plane
n_res = 100; % resolution in the n domain
nLim = 1.2*max_n_opt; % choose some maximum value for the domain... you'll need to adjust this to zoom in on the region you're interested in.
n_Plane = linspace(1,nLim,n_res);
% loop through the plane and compute uptake rate at each pixel
for a = 1:numel(n_Plane)
    Vmax_Plane = n_Plane(a).*kcat; % mole cell-1 s-1
    [ks_d_Plane, ks_p_Plane] = getK(Vmax_Plane, kcat, n_Plane(a), r, A, D, u); % mol m-3
    for b = 1:numel(S)
        [v_Plane(a,b)] = getUptake(S(b), Vmax_Plane, ks_d_Plane, ks_p_Plane); % mole cell-1 s-1
    end
end
v_Plane2 = v_Plane .* 1e15 .* 3600; % fmol cell-1 h-1

%% Add experimental data from Schmidt et al., 2016
% 
% nE = [2381 3919 6753 9402]; % transporters per cell corresponding to 0.12, 0.20, 0.35, and 0.50 h-1
% rE = 1e-7.*[7.683 7.919 8.306 8.628]; % radius (m)
% muE = [0.12 0.20 0.35 0.50]; % growth rate (h-1)
% mu_max = 1.80; % maximum growth rate (h-1) using Schmidt's LB batch growth rate
% 
% % calculate steady-state glucose concentrations in the chemostats
% for a = 1:numel(nE)    
%     VmaxE(a) = kcat.*nE(a); % mole cell-1 s-1
%     [ks_d_E(a), ks_p_E(a)] = getK(VmaxE(a), kcat, nE(a), rE(a), A, D, u); % mol m-3 
% end
% ks_E = ks_p_E + ks_d_E; % mol m-3
% S_E = (muE .* ks_E) ./ (mu_max - muE); % mol m-3
% for a = 1:numel(nE)
%     [v_E(a)] = getUptake(S_E(a), VmaxE(a), ks_d_E(a), ks_p_E(a)); % mole cell-1 s-1
% end
% v_E2 = v_E .* 1e15 .* 3600; % fmol cell-1 h-1

%% Plot results

figure
subplot(2,1,1) % plot uptake in Tp vs S space. Overlay data
h1 = contourf(S,n_Plane,v_Plane2,20,'LineStyle','none');
hold on
h2 = plot(S,n_opt_minDP,'-r','LineWidth',4);
if SA_transition
    h5 = plot([S_SA_lb S_SA_lb],[0 nLim],'--g','LineWidth',3)
    h6 = plot([S_SA_ub S_SA_ub],[0 nLim],'--g','LineWidth',3)
else
    h7 = plot([S_star S_star],[0 nLim],'--g','LineWidth',3)
end
h3 = plot([min(S) max(S)],[n_max n_max],'--k','LineWidth',3)
%h4 = plot([S_E 33],[nE nG],'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',15,'LineWidth',3);
xlabel('S_\infty [mol m^-^3]');
ylabel('PtsG [cell^-^1]');
xlim([1e-6 1e-2])
ylim([0 nLim])
h = colorbar;
colormap('jet')
ylabel(h, 'Uptake rate [fmol cell^-^1 h^-^1]','FontSize',20);
set(gca,'FontSize',20);
set(gca,'XScale','log');
if max(n_opt_minDP) > nLim % Ugly AF but works
    if SA_transition
        hl = legend([ h2, h3, h5, h6],[{'n^*'},{'n_m_a_x'},{'S_S_A^l^b'},{'S_S_A^u^b'}])
    else
        hl = legend([ h2, h3, h7],[{'n^*'},{'n_m_a_x'},{'S^*'}])
    end
else
    if SA_transition
        hl = legend([h2, h5, h6],[{'n^*'},{'S_S_A^l^b'},{'S_S_A^u^b'}])
    else
        hl = legend([h2, h7],[{'n^*'},{'S^*'}])
    end
end
set(hl,'FontSize',20,'EdgeColor','w')

subplot(2,1,2) % Uptake versus concentration
% plot([S_E 33],[v_E2 VmaxG.*1e15.*3600],'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',15,'LineWidth',3);
plot(S,v_P,'-r','LineWidth',4);
hold on
plot(S,v_P3,'-c','LineWidth',4);
plot(S,vD.*1e15.*3600,'--b','LineWidth',3);
plot([S_star S_star],[0 2.*VmaxG.*1e15.*3600],'--g','LineWidth',3)
plot([min(S) max(S)],[VmaxG.*1e15.*3600 VmaxG.*1e15.*3600],'--','Color',[253, 187, 132]./256,'LineWidth',3);
set(gca,'XScale','log');
xlim([1e-6 1e-2])
ylim([0 1.2.*VmaxG.*1e15.*3600])
xlabel('S_\infty [mol m^-^3]')
ylabel('Uptake Rate [fmol cell^-^1 h^-^1]')
hl = legend('Optimal','Batch-acclimated','Maximum Diffusive Flux','S^*','V_m_a_x^G','Location','SouthEast')
set(hl,'EdgeColor','w')
set(gca,'FontSize',20)


%% SA limitation
% all strains
for a = 1:Gridding.nStr
    SA = 4 .* pi() .* (1e-6.*StrainSolution.r_opt(:,:,a)).^2;
    for b = 1:4
        nTp = StrainSolution.TpOpt(:,:,b+1,a);
        SA_Tp(:,:,b) = nTp .* TpDat.A(b);
    end
    SA_Tp_total = nansum(SA_Tp,3);
    pct_available_SA(:,:,a) = SA_Tp_total ./ (0.085.*SA);
end

% ecotypes
for a = 1:numel(ecotypeList)
    SA = 4 .* pi() .* (1e-6.*EcotypeSolution.(ecotypeList{a}).r_opt).^2;
    for b = 1:4
        nTp = EcotypeSolution.(ecotypeList{a}).TpOpt(:,:,b+1);
        SA_Tp(:,:,b) = nTp .* TpDat.A(b);
    end
    SA_Tp_total = nansum(SA_Tp,3);
    pct_available_SA_eco(:,:,a) = SA_Tp_total ./ (0.085.*SA);
end

figure
str = 'MIT9312';
str_idx = find(strcmp(str,Gridding.strNameVec));

imagesc(x,y,pct_available_SA(:,:,str_idx) - 0.25)
hc = colorbar
ylabel(hc,'Fraction of available surface occupied')
colormap('jet')
set(gca,'clim',[0 1])
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)


figure
imagesc(x,y,pct_available_SA_eco(:,:,1) - 0.25)
hc = colorbar
ylabel(hc,'Fraction of available surface occupied')
colormap('jet')
set(gca,'clim',[0 1])
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)

% population
SA = 4 .* pi() .* (1e-6.*PopulationSolution.r_opt).^2;
for b = 1:4
    nTp = PopulationSolution.TpOpt(:,:,b+1) ./ CruiseData.Pro(Gridding.stationsVec2,:)';
    SA_Tp(:,:,b) = nTp .* TpDat.A(b);
end
SA_Tp_total = nansum(SA_Tp,3);
pct_available_SA_pop = SA_Tp_total ./ (0.085.*SA);

figure
imagesc(x,y,pct_available_SA_pop)
hc = colorbar
ylabel(hc,'Fraction of available surface occupied')
colormap('jet')
set(gca,'clim',[0 1])
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)

