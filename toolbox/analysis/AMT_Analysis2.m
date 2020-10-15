%% AMT Analysis 2
% A bit more pointed analysis and plots for the paper...
% Sections:
% Growth versus abundance
% Photosynthesis and photoacclimation
% Macromolecular acclimation
% Nutrient uptake and transporters
% Metabolic fluxes

% Data processing: one idea is to include in line plots a few different
% categories: all strains (thin lines), top strain from each ecotype
% (thicker lines), and the abundance weighted population totals (thickest
% line). 

%% Load Level 2 FullSolution
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;

%% Parse out Gridding, CruiseData, FileNames, and PanGEM from FullSolution
Gridding = FullSolution.Gridding;
FileNames = FullSolution.FileNames;
CruiseData = FullSolution.CruiseData;
PanGEM = FullSolution.PanGEM;

%% Retrieve Strain, Ecotype, and Population solution structures

[StrainSolution, EcotypeSolution, PopulationSolution] = parseAMTSolutions(FullSolution)

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


%% Plot color scheme
ecotypeColors = varycolor(numel(ecotypeList));
strainEcotypeColors = [ecotypeColors repmat(0.05,size(ecotypeColors,1),1)];


%% Cruise Variables
figure
subplot(2,2,1)
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,CruiseData.T(Gridding.stationsVec2,:)')
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
colorbar
title('Temperature')
set(gca,'FontSize',20)
subplot(2,2,2)
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(CruiseData.PAR(Gridding.stationsVec2,:)'))
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
colorbar
title('log PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)
subplot(2,2,3)
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)'))
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
colorbar
title('log DIN [nM]')
set(gca,'FontSize',20)
subplot(2,2,4)
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(CruiseData.Orthophosphate(Gridding.stationsVec2,:)'))
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
colorbar
title('log Phosphate [nM]')
set(gca,'FontSize',20)

%% Cell size
figure
PopDW = massConversion .* (PopulationSolution.r_opt.^3);
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,1e15*PopDW)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
colorbar
set(gca,'FontSize',20)
caxis([30 80])

% 
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = StrainSolution.r_opt(:,:,strEco_idx{i}(j));
        y = CruiseData.PAR(Gridding.stationsVec2,:)';
        %y = StrainSolution.Growth(:,:,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    %subplot(3,2,i)
    x = EcotypeSolution.(ecotypeList{i}).r_opt;
    y = CruiseData.PAR(Gridding.stationsVec2,:)';
    %y = EcotypeSolution.(ecotypeList{i}).Growth;
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',10);
    set(gca,'YScale','log');
    hold on
end
x = PopulationSolution.r_opt;
y = CruiseData.PAR(Gridding.stationsVec2,:)';
plot(x,y,'.k','MarkerSize',20);
%set(gca,'YScale','log')
%set(gca,'XScale','log')
xlabel('Cell Size')
ylabel('PAR')
set(gca,'FontSize',20)
%% Growth versus abundance
% Thinking of two plots here: the first would go in the SI, or maybe as an
% inset to the second figure, just a direct comparison between ecotype
% abundance and growth rate. The second could go in the paper, a figure
% showing population weighted growth rates against total cell numbers.
% Would be nice to include some bounds with a population growth model, and
% indicate the carrying capacity on the same plot.

% First figure (ecotype growth v abundance
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        %plot(CruiseData.T(Gridding.stationsVec2,:)',squeeze(StrainSolution.Growth(:,:,strEco_idx{i}(j))),'-','Color',strainEcotypeColors(i,:),'LineWidth',0.5)
        %p = plot(CruiseData.Pro(Gridding.stationsVec2,:)',squeeze(StrainSolution.Growth(:,:,strEco_idx{i}(j))),'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        x = CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)';
        y = StrainSolution.Growth(:,:,strEco_idx{i}(j))
        p = plot(y,x,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        
        %p.Color(4) = 0.2;
        hold on
    end
end
for i = 1:numel(ecotypeList)
    %subplot(3,2,i)
    x = CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Growth;
    plot(y,x,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',10);
    %set(gca,'YScale','log');
    hold on
end
% plot(PopulationSolution.Growth,CruiseData.Pro(Gridding.stationsVec2,:)','.k','MarkerSize',15);
set(gca,'YScale','log')
%set(gca,'XScale','log')
ylabel('Abundance [cells ml^-^1]')
xlabel('Growth rate [h^-^1]')
set(gca,'FontSize',20)

% Same but in separate panels
figure
for a = 1:numel(ecotype)
    subplot(1,5,a)
    for b = 1:numel(strEco_idx{a})
        x = CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)';
        y = StrainSolution.Growth(:,:,strEco_idx{a}(b))
        p = plot(y,x,'.','MarkerFaceColor',ecotypeColors(a,:),'MarkerEdgeColor',ecotypeColors(a,:),'MarkerSize',5);
        hold on
    end
    x = CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{a}).Growth;
    plot(y,x,'.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',10);
end
set(gca,'YScale','log')
ylabel('Abundance [cells ml^-^1]')
xlabel('Growth rate [h^-^1]')
set(gca,'FontSize',20)


% Second figure (
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        %plot(CruiseData.T(Gridding.stationsVec2,:)',squeeze(StrainSolution.Growth(:,:,strEco_idx{i}(j))),'-','Color',strainEcotypeColors(i,:),'LineWidth',0.5)
        %p = plot(CruiseData.Pro(Gridding.stationsVec2,:)',squeeze(StrainSolution.Growth(:,:,strEco_idx{i}(j))),'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        p = plot(squeeze(StrainSolution.Growth(:,:,strEco_idx{i}(j))),CruiseData.Pro(Gridding.stationsVec2,:)','.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        
        %p.Color(4) = 0.2;
        hold on
    end
end
for i = 1:numel(ecotypeList)
    %plot(CruiseData.T(Gridding.stationsVec2,:)',EcotypeSolution.(ecotypeList{i}).Growth,'-','Color',ecotypeColors(i,:),'LineWidth',1)
    plot(EcotypeSolution.(ecotypeList{i}).Growth,CruiseData.Pro(Gridding.stationsVec2,:)','.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
%plot(CruiseData.T(Gridding.stationsVec2,:)',PopulationSolution.Growth,'-k','LineWidth',3)
plot(PopulationSolution.Growth,CruiseData.Pro(Gridding.stationsVec2,:)','.k','MarkerSize',15);
set(gca,'YScale','log')
%set(gca,'XScale','log')
ylabel('Total abundance [cells ml^-^1]')
xlabel('Growth rate [h^-^1]')
set(gca,'FontSize',20)

% 
% muVec = linspace(0,0.12,100);
% muMax = 0.12
% M = 3e5;
% k1 = 0.1
% k2 = 4
% P = 2;
% 
% n1 = (M.*(muVec./muMax).^P)./(k1.^P+(muVec./muMax).^P);
% n2 = (M.*(muVec./muMax).^P)./(k2.^P+(muVec./muMax).^P);
% %figure
% plot(muVec,n1,'-k','LineWidth',3);
% hold on
% plot(muVec,n2,'-r','LineWidth',3);

%set(gca,'YScale','log')
%ylim([100 1e6])


% Just the total population (remove zero-growth rate data)
zeroMu = find(PopulationSolution.Growth ==0)
abundance_nonZeroMu = CruiseData.Pro(Gridding.stationsVec2,:)'
abundance_nonZeroMu(zeroMu) = 0;
figure
plot(PopulationSolution.Growth,abundance_nonZeroMu,'.k','MarkerSize',15);
set(gca,'YScale','log')
%set(gca,'XScale','log')
ylabel('Abundance [cells ml^-^1]')
xlabel('Growth rate [h^-^1]')
set(gca,'FontSize',20)



% regress a couple functions to the abundance vs mu relationship
mu = reshape(PopulationSolution.Growth,size(PopulationSolution.Growth,1).*size(PopulationSolution.Growth,2),1)
B = reshape(abundance_nonZeroMu,size(abundance_nonZeroMu,1).*size(abundance_nonZeroMu,2),1)

% Holling 2 fit
MM_fct = @(x,xdata) (xdata.*x(1))./(x(2) - xdata);
B_MM = nlinfit(mu, B, MM_fct, [1e5 0.1])
% log fit
logB = log(B);
logB(find(logB==-inf)) = NaN;
mdl = LinearModel.fit(mu,logB);

figure
plot(PopulationSolution.Growth,abundance_nonZeroMu,'.k','MarkerSize',15);
hold on
plot(mu,(mu.*B_MM(1))./(B_MM(2)-mu),'.g','MarkerSize',10);
%plot(mu,log(mdl.Coefficients.Estimate(2).*mu)+10000,'-r')
set(gca,'YScale','log')
%set(gca,'XScale','log')
ylabel('Abundance [cells ml^-^1]')
xlabel('Growth rate [h^-^1]')
set(gca,'FontSize',20)


%% Johnson et al., 2006, Science vs MSE
% Growth anomaly
figure
subplot(2,6,1)
z = PopulationSolution.Growth;
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
colormap(jet)
caxis([0 0.1])
colorbar
title('Total Population [h^-^1]')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',10)
for i = 1:numel(ecotypeList)
    subplot(2,6,i+1)
    z = EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    caxis([0 0.1])
    colorbar
    title(ecotypeList{i})
    xlabel('Latitude')
    ylabel('Depth')
    set(gca,'FontSize',10)
end
subplot(2,6,7)
z = CruiseData.Pro(Gridding.stationsVec2,:)';
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(z));
colormap(jet)
caxis([0 5.2])
colorbar
title('Total Population [log_1_0 cells ml^-^1]')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',10)
for i = 1:numel(ecotypeList)
    subplot(2,6,i+7)
    z = CruiseData.(ecotypeList{i})(Gridding.stationsVec2,:)';
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(z));
    colormap(jet)
    caxis([0 5.2])
    colorbar
    title(ecotypeList{i})
    xlabel('Latitude')
    ylabel('Depth')
    set(gca,'FontSize',10)
end

%% MSE acclimation steps

% Absolute change in growth rates at each step
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = StrainSolution.StrMod3_growth(:,:,strEco_idx{i}(j));
        y = StrainSolution.StrMod4_growth(:,:,strEco_idx{i}(j));
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = EcotypeSolution.(ecotypeList{i}).StrMod3_growth;
    y = EcotypeSolution.(ecotypeList{i}).StrMod4_growth;
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
%x = CruiseData.PAR(Gridding.stationsVec2,:)';
%y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ ( PopulationSolution.Fluxes(:,:,rxn_idx(2))  + PopulationSolution.Fluxes(:,:,rxn_idx(3))  + PopulationSolution.Fluxes(:,:,rxn_idx(4))  ); 
   
%plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
%set(gca,'XScale','log')
ylabel('After PhysOpt')
xlabel('After TpOpt')
set(gca,'FontSize',20)

% Violin plots of growth rates at each step of the MSE, by depth
for a = 1:numel(ecotypeList)
    x_s3(a,:) = nanmean(EcotypeSolution.(ecotypeList{a}).StrMod3_growth,2);
    x_s4(a,:) = nanmean(EcotypeSolution.(ecotypeList{a}).StrMod4_growth,2);
end

PhysOpt_delta = x_s4./x_s3;
PhysOpt_delta(find(PhysOpt_delta <1)) = 1;
figure
plot(100*(PhysOpt_delta-1),Gridding.depthVec,'.-','MarkerSize',20,'LineWidth',3)
xlabel('% growth rate increase')
ylabel('Depth')
set(gca,'FontSize',20)
set(gca,'YDir','reverse')
legend


%% Growth limitation regions for the total population
% Shadow prices as indicators of limitation
limitingMets = [{'Exciton'},{'Ammonia'},{'Nitrite'},{'Nitrate'},{'Orthophosphate'},{'L_Glutamate'}];
for i = 1:numel(limitingMets)
    limitingMets_idx(i) = find(strcmp(limitingMets{i},PanGEM.mets));
end

% find the minimum depth index of each shadow price depth profile
for a = 1:numel(ecotypeList)
    for b = 1:Gridding.nStations
        for c = 1:numel(limitingMets)
            shadowVec = EcotypeSolution.(ecotypeList{a}).Shadow(:,b,limitingMets_idx(c));
            N = numel(shadowVec);
            ySEM = std(shadowVec)/sqrt(N);      
            CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
            yCI95 = bsxfun(@times, ySEM, CI95(:));  
            maxShadow_idx2Vec = find(shadowVec > yCI95(2));
            if isempty(maxShadow_idx2Vec)
                maxShadow_idx2(a,b,c) = Gridding.nZ;
            else
                maxShadow_idx2(a,b,c) = min(maxShadow_idx2Vec);
            end
            
            [maxShadow1(a,b,c) maxShadow_idx1(a,b,c)] = nanmax(shadowVec);
            
        end
    end
end

% let's just start with light and nitrogen (exciton and glut)
for a = 1:numel(ecotypeList)
    limMat = zeros(Gridding.nZ,Gridding.nStations);
    for b = 1:Gridding.nStations       
        limMat(maxShadow_idx1(a,b,6):end,b) = 1;
        limMat(maxShadow_idx1(a,b,1):end,b) = limMat(maxShadow_idx2(a,b,1),b) + 2;
    end
    limMat2(:,:,a) = limMat;
end
figure
for a = 1:numel(ecotypeList)
    subplot(1,5,a)
    imagesc(Gridding.stationsVec2,Gridding.depthVec,limMat2(:,:,a));
    colormap('jet')
    colorbar
    set(gca,'FontSize',20)
    xlabel('Lat')
    ylabel('Depth')
end


% another approach, get a binary vector for all values above the 95%CI and
% add them together
binaryShadowMat = zeros(Gridding.nZ,Gridding.nStations,numel(ecotypeList),numel(limitingMets));

for a = 1:numel(ecotypeList)
    for b = 1:Gridding.nStations
        for c = 1:numel(limitingMets)
            shadowVec = EcotypeSolution.(ecotypeList{a}).Shadow(:,b,limitingMets_idx(c));
            N = numel(shadowVec);
            ySEM = std(shadowVec)/sqrt(N);      
            CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
            yCI95 = bsxfun(@times, ySEM, CI95(:));  
            maxShadow_idx2Vec = find(shadowVec > yCI95(2));
            binaryShadowMat(maxShadow_idx2Vec,b,a,c) = 1;
        end
    end
end




    
% Transporter space limitation
TpDat = readtable('GitHub/mse_AMT/data/db/TpDat.csv','Delimiter',',','ReadVariableNames',true);


for a = 1:numel(ecotypeList)
    SA_cell(:,:,a) = 4 * pi() * (1e-6 * EcotypeSolution.(ecotypeList{a}).r_opt).^2;
    for b = 1:numel(TpDat.A)
        SA_Tp(:,:,b) = EcotypeSolution.(ecotypeList{a}).TpOpt(:,:,b+1) .* TpDat.A(b);
    end
    SA_Tp_sum(:,:,a) = nansum(SA_Tp,3);
end

% find indices for each profile of where SA is maxed out
binarySALim = zeros(Gridding.nZ,Gridding.nStations,numel(ecotypeList));
for a = 1:numel(ecotypeList)
    SAMat = SA_Tp_sum(:,:,a) ./ (SA_cell(:,:,a) .* 0.1049);
    maxIdx = find(SAMat > 0.99);
    [maxIdx_row maxIdx_col] = ind2sub(size(SAMat),maxIdx);
    for b = 1:numel(maxIdx);
        binarySALim (maxIdx_row(b),maxIdx_col(b),a) = 1;
    end
end

        
figure
for a = 1:numel(ecotypeList)
    subplot(1,5,a)
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,binaryShadowMat(:,:,a,6) + 2*binaryShadowMat(:,:,a,1) + 4*binarySALim(:,:,a));
    colormap('jet')
    colorbar
    title(strcat(ecotypeList2{a},' Limiting resource'))
    xlabel('Latitude')
    ylabel('Depth')
    set(gca,'FontSize',20)

end



Gridding.strNameVec

% Same but just pick a single strain this time
% MIT_9314 is 32, SB is 36
strInd = 36;
binaryShadowMat = zeros(Gridding.nZ,Gridding.nStations,numel(limitingMets));

for b = 1:Gridding.nStations
    for c = 1:numel(limitingMets)
        shadowVec = StrainSolution.Shadow(:,b,limitingMets_idx(c),strInd);
        N = numel(shadowVec);
        ySEM = std(shadowVec)/sqrt(N);
        CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:));
        maxShadow_idx2Vec = find(shadowVec > yCI95(2));
        binaryShadowMat(maxShadow_idx2Vec,b,c) = 1;
    end
end

% Transporter space limitation


SA_cell = 4 * pi() * (1e-6 * StrainSolution.r_opt(:,:,strInd)).^2;
for b = 1:numel(TpDat.A)
    SA_Tp(:,:,b) = StrainSolution.TpOpt(:,:,b+1,strInd) .* TpDat.A(b);
end
SA_Tp_sum = nansum(SA_Tp,3);


% find indices for each profile of where SA is maxed out
binarySALim = zeros(Gridding.nZ,Gridding.nStations);
SAMat = SA_Tp_sum ./ (SA_cell .* 0.1049);
maxIdx = find(SAMat > 0.99);
[maxIdx_row maxIdx_col] = ind2sub(size(SAMat),maxIdx);
for b = 1:numel(maxIdx)
    binarySALim (maxIdx_row(b),maxIdx_col(b)) = 1;
end



figure

imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,binaryShadowMat(:,:,6) + 2*binaryShadowMat(:,:,1));
colormap('jet')
title(strcat(Gridding.strNameVec{strInd},' Limiting resource'))
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)



% Compare average limitation for two strains
% MIT_9314 is 32, SB is 36
strInd = 32;
binaryShadowMat = zeros(Gridding.nZ,Gridding.nStations,numel(limitingMets));

for b = 1:Gridding.nStations
    for c = 1:numel(limitingMets)
        shadowVec = StrainSolution.Shadow(:,b,limitingMets_idx(c),strInd);
        N = numel(shadowVec);
        ySEM = std(shadowVec)/sqrt(N);
        CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:));
        maxShadow_idx2Vec = find(shadowVec > yCI95(2));
        binaryShadowMat(maxShadow_idx2Vec,b,c) = 1;
    end
end

% average limitation vs depth
% nitrogen
x = mean(binaryShadowMat(:,:,6),2);
x_norm = (x - min(x))./(max(x) - min(x));
% light
y = mean(binaryShadowMat(:,:,1),2);
y_norm = (y - min(y))./(max(y) - min(y));
Nlim(:,1) = x;
Llim(:,1) = y;

% transition depth
r = binaryShadowMat(:,:,6) + 2*binaryShadowMat(:,:,1);
for a = 1:Gridding.nStations
    xSt = r(:,a);
    xDiff = diff(xSt);
    trans_idx = find(xDiff == 2);
    if ~isempty(trans_idx)
    trans_z(a) = Gridding.depthVec(trans_idx(end));
    end
end
trans_z2(:,1) = trans_z;

figure

imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,binaryShadowMat(:,:,6) + 2*binaryShadowMat(:,:,1));
colormap('jet')
title(strcat(Gridding.strNameVec{strInd},' Limiting resource'))
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)


strInd = 36;
binaryShadowMat = zeros(Gridding.nZ,Gridding.nStations,numel(limitingMets));

for b = 1:Gridding.nStations
    for c = 1:numel(limitingMets)
        shadowVec = StrainSolution.Shadow(:,b,limitingMets_idx(c),strInd);
        N = numel(shadowVec);
        ySEM = std(shadowVec)/sqrt(N);
        CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:));
        maxShadow_idx2Vec = find(shadowVec > yCI95(2));
        binaryShadowMat(maxShadow_idx2Vec,b,c) = 1;
    end
end

% average limitation vs depth
% nitrogen
x = mean(binaryShadowMat(:,:,6),2);
x_norm = (x - min(x))./(max(x) - min(x));
% light
y = mean(binaryShadowMat(:,:,1),2);
y_norm = (y - min(y))./(max(y) - min(y));

Nlim(:,2) = x;
Llim(:,2) = y;

% transition depth
r = binaryShadowMat(:,:,6) + 2*binaryShadowMat(:,:,1);
for a = 1:Gridding.nStations
    xSt = r(:,a);
    xDiff = diff(xSt);
    trans_idx = find(xDiff == 2);
    if ~isempty(trans_idx)
    trans_z(a) = Gridding.depthVec(trans_idx(1));
    end
end
trans_z2(:,2) = trans_z;

figure
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,binaryShadowMat(:,:,6) + 2*binaryShadowMat(:,:,1));
colormap('jet')
title(strcat(Gridding.strNameVec{strInd},' Limiting resource'))
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)

figure
subplot(1,2,1)
plot(Nlim(:,1),Gridding.depthVec,'-k','LineWidth',3);
hold on;
plot(Nlim(:,2),Gridding.depthVec,'-r','LineWidth',3);
xlabel('Limiting factor')
ylabel('Depth')
title(Gridding.strNameVec{strInd});
set(gca,'FontSize',20);
title(Gridding.strNameVec{strInd});
subplot(1,2,2)
plot(Llim(:,1),Gridding.depthVec,'-k','LineWidth',3);
hold on;
plot(Llim(:,2),Gridding.depthVec,'-r','LineWidth',3);
xlabel('Limiting factor')
ylabel('Depth')
title(Gridding.strNameVec{strInd});
set(gca,'FontSize',20);
title(Gridding.strNameVec{strInd});


figure
%s = pcolor(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)'))
s = pcolor(CruiseData.Lat(Gridding.stationsVec2(1:end-1)),Gridding.depthVec,log10(CruiseData.Nitrate(Gridding.stationsVec2(1:end-1),:)' + CruiseData.Nitrite(Gridding.stationsVec2(1:end-1),:)'))

hold on
h1 = plot(CruiseData.Lat(Gridding.stationsVec2(1:end-1)),trans_z2((1:end-1),1),'.k','MarkerSize',40);
hold on;
h2 = plot(CruiseData.Lat(Gridding.stationsVec2(1:end-1)),trans_z2((1:end-1),2),'.r','MarkerSize',30);
s.LineStyle = 'none';
%s.FaceColor = 'interp';
set(gca,'YDir','reverse')
xlabel('Latitude')
ylabel('Depth')
colormap('parula')
h = colorbar;
ylabel(h, 'log_1_0 Nitrate + Nitrite [nM]')
title('Transition from N-lim to Light-lim')
legend([h1, h2],'MIT-9314','SB')
set(gca,'FontSize',20)

% Growth rate comparison
figure
z = 100*(StrainSolution.Growth(:,:,36) - StrainSolution.Growth(:,:,32)) ./ nanmean(StrainSolution.Growth(:,:,[32 36]),3);
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
colormap('jet')
%caxis([0 0.1])
colorbar
title('SB versus MIT-9314 growth rate [%]')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)

% QY

rnxNames = [{'R00024'},{'DVChla_abs'},{'DVChlb_abs'},{'aCar_abs'},{'CEF'},{'MehlerPSI'},{'R09503'},{'R09542pc'},{'R01334'}];
for i = 1:numel(rnxNames)
    rxn_idx(i) = find(strcmp(rnxNames{i},PanGEM.rxns));
end

figure 
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y1 = StrainSolution.Fluxes(:,:,rxn_idx(1),32) ./ ( StrainSolution.Fluxes(:,:,rxn_idx(2),32) + StrainSolution.Fluxes(:,:,rxn_idx(3),32) + StrainSolution.Fluxes(:,:,rxn_idx(4),32) ); 
y2 = StrainSolution.Fluxes(:,:,rxn_idx(1),36) ./ ( StrainSolution.Fluxes(:,:,rxn_idx(2),36) + StrainSolution.Fluxes(:,:,rxn_idx(3),36) + StrainSolution.Fluxes(:,:,rxn_idx(4),36) ); 
h1 = plot(x,y1,'.k','MarkerSize',40);
hold on
h2 = plot(x,y2,'.r','MarkerSize',30);
set(gca,'XScale','log')
ylabel('Quantum Yield')
xlabel('PAR []')
legend([h1(1), h2(1)],[{'MIT-9314'},{'SB'}])
set(gca,'FontSize',20)

% fluxes at some station
x1 = (abs(squeeze(StrainSolution.Fluxes(:,:,:,32) ) ./ squeeze(StrainSolution.Growth(:,:,32))));
x2 = (abs(squeeze((StrainSolution.Fluxes(:,:,:,36))) ./ squeeze(StrainSolution.Growth(:,:,36))));
PanGEM.rxnNames(find(abs((x1-x2))./nanmean([x1, x2],2) > 1))

% fluxes across the entire domain
x1 = (abs(squeeze(StrainSolution.Fluxes(:,:,:,32) ) ./ squeeze(StrainSolution.Growth(:,:,32))));
x2 = (abs(squeeze((StrainSolution.Fluxes(:,:,:,36))) ./ squeeze(StrainSolution.Growth(:,:,36))));
cumDelta_idx = [];
for a = 1:Gridding.nZ
    for b = 1:Gridding.nStations
        delta_idx = find(abs( x1(a,b,:) - x2(a,b,:) ) ./ nanmean([x1(a,b,:),x2(a,b,:)],2) > 1.1);
        cumDelta_idx = [cumDelta_idx; delta_idx];
    end
end
uniqueDelta_idx = unique(cumDelta_idx);
for a = 1:numel(uniqueDelta_idx)
    freqDelta_idx(a) = numel(find(cumDelta_idx == uniqueDelta_idx(a))) ./ numel(uniqueDelta_idx);
end
[freqDelta_idx_sorted, freqDelta_idx_sorted_idx] = sort(freqDelta_idx,'descend');
PanGEM.rxnNames(uniqueDelta_idx(freqDelta_idx_sorted_idx(1:100)))
PanGEM.rxnNames(diffDelta_idx)

figure
plot(1:numel(x1),x1-x2,'ok');
hold on;
text(1:numel(x1),x1-x2,PanGEM.rxns);
set(gca,'XScale','log','YScale','log')

% get difference between networks
load('GitHub/mse_AMT/data/GEM/targStrMod3.mat');
SB_notMIT9314_idx = find(sum(abs(targStrMod3.SB.S),1) - sum(abs(targStrMod3.MIT_9314.S),1) >0)
MIT9314_notSB_idx = find(sum(abs(targStrMod3.SB.S),1) - sum(abs(targStrMod3.MIT_9314.S),1) <0)

PanGEM.rxnNames(SB_notMIT9314_idx)
PanGEM.rxnNames(MIT9314_notSB_idx)

SB_Rxns_idx = find(sum(abs(targStrMod3.SB.S),1));
MIT9314_Rxns_idx = find(sum(abs(targStrMod3.MIT_9314.S),1));
shared = intersect(SB_Rxns_idx,MIT9314_Rxns_idx);
shared_nRxns = numel(shared);
%% Photosynthesis and photoacclimation
% A few plots here, several multipanel
% QY, carbohydrate content, Chl/C, b/a, zeax/other abs,
% glycolysis/gluconeogenesis, TCA/oPPP, Mehler, CEF, others

rnxNames = [{'R00024'},{'DVChla_abs'},{'DVChlb_abs'},{'aCar_abs'},{'CEF'},{'MehlerPSI'},{'R09503'},{'R09542pc'},{'R01334'}];
for i = 1:numel(rnxNames)
    rxn_idx(i) = find(strcmp(rnxNames{i},PanGEM.rxns));
end

% QY
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{i}(j)) ./ ( StrainSolution.Fluxes(:,:,rxn_idx(2),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,rxn_idx(3),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,rxn_idx(4),strEco_idx{i}(j)) ); 
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(1)) ./ ( EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(2)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(3)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(4)) );
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ ( PopulationSolution.Fluxes(:,:,rxn_idx(2))  + PopulationSolution.Fluxes(:,:,rxn_idx(3))  + PopulationSolution.Fluxes(:,:,rxn_idx(4))  ); 
   
plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Quantum Yield')
xlabel('PAR [\mumol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)

figure
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ ( PopulationSolution.Fluxes(:,:,rxn_idx(2))  + PopulationSolution.Fluxes(:,:,rxn_idx(3))  + PopulationSolution.Fluxes(:,:,rxn_idx(4))  );
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,y)
colorbar
xlabel('Latitude')
ylabel('Depth [m]')
set(gca,'FontSize',20)

% Carbohydrate content (BOF index 3)
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.BOF_coefs(:,:,3,strEco_idx{i}(j));
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,3);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end

%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Carbohydrate content [g gDW^-^1]')
xlabel('PSII Flux [mmol quanta gDW^-^1 h^-^1]')
set(gca,'FontSize',20)

figure % whole population had different units, let's plot that separately
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.BOF_coefs(:,:,3);
   
plot(x,y,'.k','MarkerSize',15)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Carbohydrate content [g ml^-^1]')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)

% b/a ratio
% for this we'll need IrrDat to calculate the peak absorption ratio (450
% for dvchla and 474 for dvchlb
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/envData/IrrDat.mat');
aPeak_lambda = 450;
bPeak_lambda = 474;
aPeak_ind = find(IrrDat.Lambda == aPeak_lambda);
bPeak_ind = find(IrrDat.Lambda == bPeak_lambda);
for i = 1:Gridding.nStations
    Ed_ratio(:,i) = IrrDat.Data{Gridding.stationsVec2(i)}(:,bPeak_ind) ./ IrrDat.Data{Gridding.stationsVec2(i)}(:,aPeak_ind);
end

figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x =  Ed_ratio;
        y = StrainSolution.BOF_coefs(:,:,13,strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = Ed_ratio;
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
%x =  Ed_ratio;
%y = PopulationSolution.BOF_coefs(:,:,13) ./ PopulationSolution.BOF_coefs(:,:,12);
%plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
%set(gca,'XScale','log')
ylabel('DV chl b : DV chl a ratio')
xlabel('E_d(474) / E_d(450)')
%xlabel('PAR []')
set(gca,'FontSize',20)
legend([h{1}(1) h{2}(1) h{3}(1) h{4}(1) h{5}(1)],ecotypeList2)

% b/a just against PAR
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x =  CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.BOF_coefs(:,:,13,strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x =  CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x =  CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.BOF_coefs(:,:,13) ./ PopulationSolution.BOF_coefs(:,:,12);
plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('DV chl b : DV chl a ratio')
%xlabel('E_d(474) / E_d(450)')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)


% Chl : c
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x =  CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
%x =  CruiseData.PAR(Gridding.stationsVec2,:)';
%y = PopulationSolution.BOF_coefs(:,:,12);
%plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Chlorophyll : C ratio')
%xlabel('E_d(474) / E_d(450)')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)



% Zeax : other pigs
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.BOF_coefs(:,:,15,strEco_idx{i}(j)) ./ ( StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j)) + StrainSolution.BOF_coefs(:,:,13,strEco_idx{i}(j)) + StrainSolution.BOF_coefs(:,:,14,strEco_idx{i}(j)) );
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x =  CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,15) ./ ( EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12) + EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) + EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,14) );
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x =  CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.BOF_coefs(:,:,15) ./ ( PopulationSolution.BOF_coefs(:,:,12) + PopulationSolution.BOF_coefs(:,:,13) + PopulationSolution.BOF_coefs(:,:,14) );
%plot(x,y,'.','MarkerSize',15,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Zeax : LHP')
%xlabel('E_d(474) / E_d(450)')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)

% CEF
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(5),strEco_idx{i}(j));
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(5));
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(5));
   
%plot(x,y,'.k','MarkerSize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Cyclic electron flow (PSI) [mmol gDW^-^1 h^-^1]')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)

figure % population separately again
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(5));
plot(x,y,'.k','MarkerSize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Cyclic electron flow (PSI) [mmol ml^-^1 h^-^1]')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)

% Mehler
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(6),strEco_idx{i}(j));
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(6));
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(5));
   
%plot(x,y,'.k','MarkerSize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Mehler reaction [mmol gDW^-^1 h^-^1]')
xlabel('PAR []')
set(gca,'FontSize',20)
ylim([0.01 100])

figure % population separately again
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(6));
plot(x,y,'.k','MarkerSize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Cyclic electron flow (PSI) [mmol ml^-^1 h^-^1]')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)

% P vs I curve
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{i}(j));
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(1));
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
%x = CruiseData.PAR(Gridding.stationsVec2,:)';
%y = PopulationSolution.Fluxes(:,:,rxn_idx(5));
   
%plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Carbon fixation rate [mmol gDW^-^1 h^-^1]')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)
%ylim([0.01 100])


% PB vs I curve (photosynthesis normalized to chlorophyll (and fit a curve)

figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
        plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(1)) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
    plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on        
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ PopulationSolution.BOF_coefs(:,:,12);
p = plot(x,y.*12.011.*(1/1000),'.k','MarkerSize',15)



set(gca,'XScale','log')
%set(gca,'YScale','log')
ylabel('P^B [mg C (mg Chl a)^-^1 h^-^1]')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)
%ylim([0.01 100])

% Same but linear
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.PAR(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
        plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.PAR(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,rxn_idx(1)) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
    plot(x,y.*12.011.*(1/1000),'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on        
end
x = CruiseData.PAR(Gridding.stationsVec2,:)';
y = PopulationSolution.Fluxes(:,:,rxn_idx(1)) ./ PopulationSolution.BOF_coefs(:,:,12);
plot(x,y.*12.011.*(1/1000),'.k','MarkerSize',15)



%set(gca,'XScale','log')
%set(gca,'YScale','log')
ylabel('P^B [mg C (mg Chl a)^-^1 h^-^1]')
xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20) 

% Calculate PBm, alphaB, and Ek
% PB = PBm * tanh(alphaB*E/PBm)
% Ek = PBm/alphaB
% x2(1) is PBm, x2(2) is alphaB

for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        for k = 1:Gridding.nStations
            x = CruiseData.PAR(Gridding.stationsVec2(k),:)';
            y = StrainSolution.Fluxes(:,:,rxn_idx(1),strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
            x0 = [2 0.05];
            Ed = x;
            PB = y(:,k).*12.011.*(1/1000);
            PB(find(isnan(PB))) = 0;
            PvsIfun = @(x2)(x2(1).*tanh((x2(2).*Ed)./x2(1)))-PB;
            yPred = lsqnonlin(PvsIfun,x0);
            PBm{i}(j,k) = abs(yPred(1));
            alphaB{i}(j,k) = abs(yPred(2));
            Ek{i}(j,k) = abs(yPred(1)) ./ abs(yPred(2));
        end
    end
end

EdVec = logspace(-1,3.5,100);
% percentile method
CIFcn = @(x3,p)prctile(x3,abs([0,100]-(100-p)/2));
p = 75;
% calculate mean and CI's for each ecotype across all stations, plot.
for i = 1:numel(ecotypeList)
    PBm_mean(i) = nanmean(nanmean(PBm{i}));
    PBm_CI(i,:) = CIFcn(reshape(PBm{i},size(PBm{i},1)*size(PBm{i},2),1),p);
    alphaB_mean(i) = nanmean(nanmean(alphaB{i}));
    alphaB_CI(i,:) = CIFcn(reshape(alphaB{i},size(alphaB{i},1)*size(alphaB{i},2),1),p);
    Ek_mean(i) = nanmean(nanmean(Ek{i}));
    Ek_mean_station(i,:) = nanmean(Ek{i},1)
    Ek_CI(i,:) = CIFcn(reshape(Ek{i},size(Ek{i},1)*size(Ek{i},2),1),p);
end

figure
boxplot(Ek_mean_station(:,1:end-1)','Notch','on','Labels',ecotypeList2,'Whisker',1)
ylabel('E_k [\mu mol quanta m^-^2 s^-^1]')
set(gca,'FontSize',20)
ylim([0 350])


for i = 1:numel(ecotypeList)
    plot(EdVec,PBm_mean(i).*tanh((alphaB_mean(i).*EdVec)./PBm_mean(i)),'-','Color',ecotypeColors(i,:),'LineWidth',3)
    hold on
    plot(EdVec,PBm_CI(i,1).*tanh((alphaB_CI(i,2).*EdVec)./PBm_CI(i,1)),'-','Color',strainEcotypeColors(i,:),'LineWidth',3)
    plot(EdVec,PBm_CI(i,2).*tanh((alphaB_CI(i,1).*EdVec)./PBm_CI(i,2)),'-','Color',strainEcotypeColors(i,:),'LineWidth',3)

end
set(gca,'XScale','log')

% PBm vs alphaB
figure
for i = 1:numel(ecotypeList)
    plot(alphaB{i},PBm{i},'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end
xlim([0 0.1])
ylim([0 5])
ylabel('P^B_{max}')
xlabel('\alpha^B')
set(gca,'FontSize',20)
%set(gca,'XScale','log','YScale','log')






%% Macromolecular acclimation
% Let's look at relationships between macromolecular composition, elemental
% stoichiometry, and nutrient availability
crudeFractions2 = {'DNA';'Lipid';'Carbohydrate';'Protein';'VitaCofactors';'RNA';'NB';'Ions';'BioPool';'Pigments';'CellWall';'Divinylchlorophyll_a';'Divinylchlorophyll_b';'alpha_Carotene';'Zeaxanthin'};
% compute cellular elemental stoichiometry

% all strains
for i = 1:Gridding.nZ
    for j = 1:Gridding.nStations
        for k = 1:Gridding.nStr
            BOF_coefs = squeeze(StrainSolution.BOF_coefs(i,j,1:11,k));
            [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
            Quota(i,j,:,k) = sum(MMComposition.DW,1); % mmol element gDW-1
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
        end
    end
end
% population
for i = 1:Gridding.nZ
    for j = 1:Gridding.nStations
            BOF_coefs = squeeze(PopulationSolution.BOF_coefs(i,j,1:11));
            [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
            QuotaPop(i,j,:) = sum(MMComposition.DW,1); % mmol element ml-1
    end
end
% Phosphorus vs C:P
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
        y = Quota(:,:,1,strEco_idx{i}(j)) ./ Quota(:,:,5,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x =  CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
    y = QuotaEco(:,:,1,i) ./ QuotaEco(:,:,5,i);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x =  CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
y = QuotaPop(:,:,1) ./ QuotaPop(:,:,5);
plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('C:P')
%xlabel('E_d(474) / E_d(450)')
xlabel('Orthophosphate [nM]')
set(gca,'FontSize',20)

% DIN vs C:N
figure
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

% Same but just the population ratio, add the average line

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

figure
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

figure
plot(x,y,'.k','MarkerSize',15)
set(gca,'XScale','log')
ylabel('C:N')
%xlabel('E_d(474) / E_d(450)')
xlabel('DIN [nM]')
set(gca,'FontSize',20)


% DIN : DIP vs N:P
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = ( CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)' ) ./ CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
        y = Quota(:,:,3,strEco_idx{i}(j)) ./ Quota(:,:,5,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = ( CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)' ) ./ CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
    y = QuotaEco(:,:,3,i) ./ QuotaEco(:,:,5,i)
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x = ( CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)' ) ./ CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
y = QuotaPop(:,:,3) ./ QuotaPop(:,:,5);
plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('N:P')
%xlabel('E_d(474) / E_d(450)')
xlabel('DIN : DIP')
set(gca,'FontSize',20)

% Droop style (nitrogen). Growth rate vs N quota
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = 1e15.*((1/1000) .* Quota(:,:,3,strEco_idx{i}(j)) ).*( 1e-15 .* (4/3) .* pi() .* (StrainSolution.r_opt(:,:,strEco_idx{i}(j))).^3 .* 150 );
        y = StrainSolution.Growth(:,:,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = 1e15.*((1/1000) .* QuotaEco(:,:,3,i) ).*( 1e-15 .* (4/3) .* pi() .* (EcotypeSolution.(ecotypeList{i}).r_opt).^3 .* 150 );
    y = EcotypeSolution.(ecotypeList{i}).Growth;
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end

%set(gca,'YScale','log')
%set(gca,'XScale','log')
ylabel('Growth Rate')
%xlabel('E_d(474) / E_d(450)')
xlabel('N composition')
set(gca,'FontSize',20)

% population Droop
figure
x = ( 1e15.*((1/1000) .* QuotaPop(:,:,3) ).*( 1e-15 .* (4/3) .* pi() .* (PopulationSolution.r_opt).^3 .* 150 ) ) ./ CruiseData.Pro(Gridding.stationsVec2,:)';
y = PopulationSolution.Growth;
plot(x,y,'.k','MarkerSize',15)
ylabel('Growth Rate')
xlabel('N composition')
set(gca,'FontSize',20)
%set(gca,'XScale','log')


% Protein : Carb vs N
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = ( CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)' );
        y = StrainSolution.BOF_coefs(:,:,4,strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,3,strEco_idx{i}(j));
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = ( CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)' );
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,4) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,3);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
x = ( CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)' );
y = PopulationSolution.BOF_coefs(:,:,4) ./ PopulationSolution.BOF_coefs(:,:,3);
plot(x,y,'.k','MarkerSize',15)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Protein : Carb')
%xlabel('E_d(474) / E_d(450)')
xlabel('DIN  [nM]')
set(gca,'FontSize',20)

% Contours of CP and CN
figure
for a = 1:numel(ecotypeList)
    subplot(1,5,a)
    z = QuotaEco(:,:,1,a) ./ QuotaEco(:,:,3,a);
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap('jet')
    colorbar
    xlabel('Latitude')
    ylabel('Depth')
    set(gca,'FontSize',20)
    caxis([5 7.5])
end


figure
z = QuotaPop(:,:,1) ./ QuotaPop(:,:,3);
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
colormap('jet')
colorbar
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)

figure
z = PopulationSolution.BOF_coefs(:,:,4) ./ PopulationSolution.BOF_coefs(:,:,3);
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
colormap('jet')
colorbar
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)

figure
z = PopulationSolution.BOF_coefs(:,:,1) ./ QuotaPop(:,:,1);
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
colormap('jet')
colorbar
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)


%% Figures for C:N ratio for Mick's Tallk (2020)

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

figure
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

figure
plot(x,y,'.k','MarkerSize',15)
set(gca,'XScale','log')
ylabel('C:N')
%xlabel('E_d(474) / E_d(450)')
xlabel('DIN [nM]')
set(gca,'FontSize',20)

% pie's at the high and low end
z = PopulationSolution.BOF_coefs(:,:,4) ./ PopulationSolution.BOF_coefs(:,:,3); 
[highDIN_rowIdx, highDIN_colIdx]  = find(z < 1.2);
y_s1 = PopulationSolution.BOF_coefs(highDIN_rowIdx,highDIN_colIdx,:);
yMean_s1 = squeeze(mean(nanmean(y_s1,1)));
yMean_s1_norm = yMean_s1 ./ sum(yMean_s1)

[lowDIN_rowIdx, lowDIN_colIdx]  = find(z > 3);
y_s2 = PopulationSolution.BOF_coefs(lowDIN_rowIdx,lowDIN_colIdx,:);
yMean_s2 = squeeze(mean(nanmean(y_s2,1)));
yMean_s2_norm = yMean_s2 ./ sum(yMean_s2)

figure
labels = repmat({''},size(yMean_s1_norm));%option 1
h = pie(yMean_s1_norm,crudeFractions2)

figure
h = pie(yMean_s2_norm,crudeFractions2)

% plot nitrite transporters verusus nitrite
x = CruiseData.Nitrite(Gridding.stationsVec2,:)';
y = StrainSolution.TpOpt(:,:,5,find(strcmp('MIT9313',Gridding.strNameVec)));

figure
plot(x,y,'.k','MarkerSize',15)
set(gca,'XScale','log')
set(gca,'FontSize',20)
xlabel('Nitrite [nM]')
ylabel('NitA [transporters cell^-^1]')

% Carbohydrate content (BOF index 3) vs absorption
figure
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = StrainSolution.Fluxes(:,:,find(strcmp('R00086th',PanGEM.rxns)),strEco_idx{i}(j));
        y = StrainSolution.BOF_coefs(:,:,3,strEco_idx{i}(j));
        plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,find(strcmp('R00086th',PanGEM.rxns)));
    y = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,3);
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
    hold on
end

%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Carbohydrate content [g gDW^-^1]')
xlabel('ATP Synthase Flux [mmol gDW^-^1 h^-^1]')
set(gca,'FontSize',20)



% 
% % Pick subsets and plot Prt:Carb ratio of those against light
% [highDIN_rowIdx, highDIN_colIdx]  = find(x>200);
% x_s1 = CruiseData.PAR(Gridding.stationsVec2(highDIN_colIdx),highDIN_rowIdx)';
% y_s1 = PopulationSolution.BOF_coefs(highDIN_rowIdx,highDIN_colIdx,4) ./ PopulationSolution.BOF_coefs(highDIN_rowIdx,highDIN_colIdx,3);
% 
% [lowDIN_rowIdx, lowDIN_colIdx]  = find(x<200);
% x_s2 = CruiseData.PAR(Gridding.stationsVec2(lowDIN_colIdx),lowDIN_rowIdx)';
% y_s2 = PopulationSolution.BOF_coefs(lowDIN_rowIdx,lowDIN_colIdx,4) ./ PopulationSolution.BOF_coefs(lowDIN_rowIdx,lowDIN_colIdx,3);
% 
% % x = x_s1;
% % y = y_s1;
% % xVec = reshape(x,size(x,1)*size(x,2),1);
% % yVec = reshape(y,size(y,1)*size(y,2),1);
% % [xVecAsc, xVecIdx] = sort(xVec,'ascend');
% % yVecAsc = yVec(xVecIdx);
% % xVecBinned = logspace(log10(min(xVec)),log10(max(xVec)),20);
% % for a = 1:numel(xVecBinned)-1
% %     inBinIdx= find(xVecAsc >= xVecBinned(a) & xVecAsc < xVecBinned(a+1));
% %     yMean(a) = nanmean(yVecAsc(inBinIdx));
% %     yStd(a) = nanstd(yVecAsc(inBinIdx));
% % end
% 
% figure
% plot(x_s1,y_s1,'.r','MarkerSize',10);
% hold on;
% plot(x_s2,y_s2,'.b','MarkerSize',10);
% 
% %hErr2 = fill([xVecBinned(2:end) fliplr(xVecBinned(2:end))],[(yMean + yStd) fliplr(yMean - yStd)],'g');
% %set(hErr,'LineWidth',3)
% %set(hErr2,'LineStyle','none')
% %hErr2.FaceAlpha = 0.2;
% set(gca,'FontSize',20)
% set(gca,'XScale','log')



%% Nutrient uptake and transporters
figure
% Orthophosphate vs Tp
subplot(4,1,1)

for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
        y = 1e-12.*StrainSolution.TpOpt(:,:,3,strEco_idx{i}(j)) ./ ( 4*pi() * ( (1e-6) .* StrainSolution.r_opt(:,:,strEco_idx{i}(j)) ).^2 );
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
    y = 1e-12.*EcotypeSolution.(ecotypeList{i}).TpOpt(:,:,3) ./ ( 4*pi() * ( (1e-6) .* EcotypeSolution.(ecotypeList{i}).r_opt ).^2 );
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
% x =  CruiseData.PAR(Gridding.stationsVec2,:)';
% y = PopulationSolution.BOF_coefs(:,:,6)
% plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Transporter density [\mum^-^2]')
%xlabel('E_d(474) / E_d(450)')
xlabel('Orthophosphate [nM]')
set(gca,'FontSize',10)
xlim([1 30000])

% Ammonia vs Tp
subplot(4,1,2)

for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.Ammonia(Gridding.stationsVec2,:)';
        y = 1e-12.*StrainSolution.TpOpt(:,:,2,strEco_idx{i}(j)) ./ ( 4*pi() * ( (1e-6) .* StrainSolution.r_opt(:,:,strEco_idx{i}(j)) ).^2 );
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.Ammonia(Gridding.stationsVec2,:)';
    y = 1e-12.*EcotypeSolution.(ecotypeList{i}).TpOpt(:,:,2) ./ ( 4*pi() * ( (1e-6) .* EcotypeSolution.(ecotypeList{i}).r_opt ).^2 );
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
% x =  CruiseData.PAR(Gridding.stationsVec2,:)';
% y = PopulationSolution.BOF_coefs(:,:,6)
% plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Transporter density [\mum^-^2]')
%xlabel('E_d(474) / E_d(450)')
xlabel('Ammonia [nM]')
set(gca,'FontSize',10)
xlim([1 30000])



% Nitrate vs Tp
subplot(4,1,3)
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.Nitrate(Gridding.stationsVec2,:)';
        y = 1e-12.*StrainSolution.TpOpt(:,:,4,strEco_idx{i}(j)) ./ ( 4*pi() * ( (1e-6) .* StrainSolution.r_opt(:,:,strEco_idx{i}(j)) ).^2 );
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.Nitrate(Gridding.stationsVec2,:)';
    y = 1e-12.*EcotypeSolution.(ecotypeList{i}).TpOpt(:,:,4) ./ ( 4*pi() * ( (1e-6) .* EcotypeSolution.(ecotypeList{i}).r_opt ).^2 );
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
% x =  CruiseData.PAR(Gridding.stationsVec2,:)';
% y = PopulationSolution.BOF_coefs(:,:,6)
% plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Transporter density [\mum^-^2]')
%xlabel('E_d(474) / E_d(450)')
xlabel('Nitrate [nM]')
set(gca,'FontSize',10)
xlim([1 30000])

% Nitrite vs Tp
subplot(4,1,4)
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        x = CruiseData.Nitrite(Gridding.stationsVec2,:)';
        y = 1e-12.*StrainSolution.TpOpt(:,:,5,strEco_idx{i}(j)) ./ ( 4*pi() * ( (1e-6) .* StrainSolution.r_opt(:,:,strEco_idx{i}(j)) ).^2 );
        plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    x = CruiseData.Nitrite(Gridding.stationsVec2,:)';
    y = 1e-12.*EcotypeSolution.(ecotypeList{i}).TpOpt(:,:,5) ./ ( 4*pi() * ( (1e-6) .* EcotypeSolution.(ecotypeList{i}).r_opt ).^2 );
    plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
% x =  CruiseData.PAR(Gridding.stationsVec2,:)';
% y = PopulationSolution.BOF_coefs(:,:,6)
% plot(x,y,'.k','MarkerSize',20)
%set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('NitA Transporter density [\mum^-^2]')
%xlabel('E_d(474) / E_d(450)')
xlabel('Nitrite [nM]')
set(gca,'FontSize',10)
xlim([1 30000])
tightfig


% Population New production
EX_rxns = [{'NitrateEX'},{'NitriteEX'},{'AmmoniaEX'}];
for i = 1:numel(EX_rxns)
    EX_rxns_idx(i) = find(strcmp(EX_rxns{i},PanGEM.rxns));
end
x =  CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;
z = PopulationSolution.Fluxes(:,:,EX_rxns_idx(1)) ./ (PopulationSolution.Fluxes(:,:,EX_rxns_idx(1)) + PopulationSolution.Fluxes(:,:,EX_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,EX_rxns_idx(3)) );
figure
imagesc(x,y,z);
colormap('jet');
colorbar
ylabel('Depth')
xlabel('Latitude')

set(gca,'FontSize',20)

figure
zInt = sum(PopulationSolution.Fluxes(:,:,EX_rxns_idx(1)),1) ./ (sum(PopulationSolution.Fluxes(:,:,EX_rxns_idx(1)),1) + sum(PopulationSolution.Fluxes(:,:,EX_rxns_idx(2)),1) + sum(PopulationSolution.Fluxes(:,:,EX_rxns_idx(3)),1) );
plot(x,zInt,'ok');

% uptake vs conc ratios
Uptake_rxns = [{'NitrateEX'},{'NitriteEX'},{'AmmoniaEX'},{'OrthophosphateEX'}];
for i = 1:numel(Uptake_rxns)
    Uptake_rxns_idx(i) = find(strcmp(Uptake_rxns{i},PanGEM.rxns));
end
figure

% Nitrate : Nitrite
subplot(3,1,1)
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        %x = StrainSolution.Fluxes(:,:,Uptake_rxns_idx(1),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,Uptake_rxns_idx(2),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,Uptake_rxns_idx(3),strEco_idx{i}(j));
        x = CruiseData.Nitrate(Gridding.stationsVec2,:)' ./ CruiseData.Nitrite(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,Uptake_rxns_idx(1),strEco_idx{i}(j)) ./ StrainSolution.Fluxes(:,:,Uptake_rxns_idx(2),strEco_idx{i}(j));
        plot(abs(x),abs(y),'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    %x = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(1)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(2)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(3));
    x = CruiseData.Nitrate(Gridding.stationsVec2,:)' ./ CruiseData.Nitrite(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(1)) ./ EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(2));
    plot(abs(x),abs(y),'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
plot(x,x,'-k','LineWidth',3);

%x = PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(1)) + PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(3));
%y = PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(4));
%plot(x,y,'.k','MarkerSize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Nitrate : Nitrite Uptake')
%xlabel('E_d(474) / E_d(450)')
xlabel('Nitrate : Nitrite Concentration')
set(gca,'FontSize',20)
%xlim([1e-2 1])

% Nitrite : ammonia
subplot(3,1,2)
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        %x = StrainSolution.Fluxes(:,:,Uptake_rxns_idx(1),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,Uptake_rxns_idx(2),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,Uptake_rxns_idx(3),strEco_idx{i}(j));
        x = CruiseData.Nitrite(Gridding.stationsVec2,:)' ./ CruiseData.Ammonia(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,Uptake_rxns_idx(2),strEco_idx{i}(j)) ./ StrainSolution.Fluxes(:,:,Uptake_rxns_idx(3),strEco_idx{i}(j));
        plot(abs(x),abs(y),'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    %x = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(1)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(2)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(3));
    x = CruiseData.Nitrite(Gridding.stationsVec2,:)' ./ CruiseData.Ammonia(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(2)) ./ EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(3));
    plot(abs(x),abs(y),'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
plot(x,x,'-k','LineWidth',3);

%x = PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(1)) + PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(3));
%y = PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(4));
%plot(x,y,'.k','MarkerSize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Nitrite : Ammonia Uptake')
%xlabel('E_d(474) / E_d(450)')
xlabel('Nitrite : Ammonia Concentration')
set(gca,'FontSize',20)
%xlim([0 1])

% Nitrate : ammonia
subplot(3,1,3)
for i = 1:numel(ecotypeList)
    for j = 1:numel(strEco_idx{i})
        %x = StrainSolution.Fluxes(:,:,Uptake_rxns_idx(1),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,Uptake_rxns_idx(2),strEco_idx{i}(j)) + StrainSolution.Fluxes(:,:,Uptake_rxns_idx(3),strEco_idx{i}(j));
        x = CruiseData.Nitrate(Gridding.stationsVec2,:)' ./ CruiseData.Ammonia(Gridding.stationsVec2,:)';
        y = StrainSolution.Fluxes(:,:,Uptake_rxns_idx(1),strEco_idx{i}(j)) ./ StrainSolution.Fluxes(:,:,Uptake_rxns_idx(3),strEco_idx{i}(j));
        plot(abs(x),abs(y),'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
        hold on
    end
end
for i = 1:numel(ecotypeList)
    %x = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(1)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(2)) + EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(3));
    x = CruiseData.Nitrate(Gridding.stationsVec2,:)' ./ CruiseData.Ammonia(Gridding.stationsVec2,:)';
    y = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(1)) ./ EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,Uptake_rxns_idx(3));
    plot(abs(x),abs(y),'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
    hold on
end
plot(x,x,'-k','LineWidth',3);

%x = PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(1)) + PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(3));
%y = PopulationSolution.Fluxes(:,:,Uptake_rxns_idx(4));
%plot(x,y,'.k','MarkerSize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('Nitrate : Ammonia Uptake')
%xlabel('E_d(474) / E_d(450)')
xlabel('Nitrate : Ammonia Concentration')
set(gca,'FontSize',20)
%xlim([0 1])

%% Metabolic fluxes

% PQ
PQ_rxns = [{'O2EX'},{'HCO3EX'},{'CO2EX'}];
for i = 1:numel(PQ_rxns)
    PQ_rxns_idx(i) = find(strcmp(PQ_rxns{i},PanGEM.rxns));
end
x =  CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;
z = PopulationSolution.Fluxes(:,:,PQ_rxns_idx(1)) ./ (PopulationSolution.Fluxes(:,:,PQ_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,PQ_rxns_idx(3)) );
figure
imagesc(x,y,z);
colormap('jet');
colorbar
ylabel('Depth')
xlabel('Latitude')
set(gca,'FontSize',20)


% NP vs PQ
NP = PopulationSolution.Fluxes(:,:,EX_rxns_idx(1)) ./ ( PopulationSolution.Fluxes(:,:,EX_rxns_idx(1)) + PopulationSolution.Fluxes(:,:,EX_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,EX_rxns_idx(3)) );
PQ = PopulationSolution.Fluxes(:,:,PQ_rxns_idx(1)) ./ (PopulationSolution.Fluxes(:,:,PQ_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,PQ_rxns_idx(3)) );
plot(abs(NP),abs(PQ),'ok')

% NPP/GPP
PP_rxns = [{'R00024'},{'HCO3EX'},{'CO2EX'}];
for i = 1:numel(PP_rxns)
    PP_rxns_idx(i) = find(strcmp(PP_rxns{i},PanGEM.rxns));
end


figure
x =  CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;
z = -(PopulationSolution.Fluxes(:,:,PP_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,PP_rxns_idx(3)) ) ./ PopulationSolution.Fluxes(:,:,PP_rxns_idx(1));
imagesc(x,y,z);
colormap('jet');
colorbar
ylabel('Depth')
xlabel('Latitude')
set(gca,'FontSize',20)

% Just NPP
figure
% NPP population contour
subplot(2,2,1)
x =  CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;
z = PopulationSolution.Fluxes(:,:,PP_rxns_idx(1)) - ( PopulationSolution.Fluxes(:,:,PP_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,PP_rxns_idx(3)) );
imagesc(x,y,log10(1e9*z));
colormap('jet');
colorbar
ylabel('Depth')
xlabel('Latitude')
title('NPP [nmol L^-^1 h^-^1]')
set(gca,'FontSize',20)

% NPP integrated, population and ecotypes
subplot(2,2,3)
zPopInt = trapz(Gridding.depthVec,z*1e6);
plot(x,zPopInt,'.k','MarkerSize',20);
hold on
for a = 1:numel(ecotypeList)
    ecoNPP = EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,PP_rxns_idx(1)) - ( EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,PP_rxns_idx(2)) + EcotypeSolution.(ecotypeList{a}).Fluxes(:,:,PP_rxns_idx(3)) );
    ecoNPP_Total = ecoNPP .* CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)' .* massConversion .* EcotypeSolution.(ecotypeList{a}).r_opt.^3;
    ecoNPP_Total(find(isnan(ecoNPP_Total))) = 0;
    zEcoInt = trapz(Gridding.depthVec,ecoNPP_Total*1e6);
    zEcoInt2(:,a) = zEcoInt; 
    plot(x,zEcoInt,'-.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15,'Color',ecotypeColors(a,:),'LineWidth',3)
    hold on
end
ylim([1e-5 10])
set(gca,'FontSize',20)
set(gca,'YScale','log')
xlabel('Latitude')
ylabel('Int. NPP (10-200m) [mmol m-2 h-1]')

% Biomass integrated, population and ecotypes versus NPP
subplot(2,2,2)
biomassPop = CruiseData.Pro(Gridding.stationsVec2,:)' .* massConversion .* PopulationSolution.r_opt.^3; % g ml-1
biomassPopInt = trapz(Gridding.depthVec,biomassPop*1e9*(1/12.011)); %mmol m-2
plot(biomassPopInt, zPopInt,'.k','MarkerSize',20);
hold on
for a = 1:numel(ecotypeList)
    biomassEco = CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)' .* massConversion .* EcotypeSolution.(ecotypeList{a}).r_opt.^3;
    biomassEco(find(isnan(biomassEco))) = 0;
    biomassEcoInt = trapz(Gridding.depthVec,biomassEco*1e9*(1/12.011)); %mmol m-2
    h{a} = plot(biomassEcoInt,zEcoInt2(:,a), '.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15);
    hold on
end
xlim([1e-4 100])
ylim([1e-5 10])
set(gca,'FontSize',20)
set(gca,'YScale','log','XScale','log')
xlabel('Int. Biomass (10-200m) [mmol m-2]')
ylabel('Int. NPP (10-200m) [mmol m-2 h-1]')
legend([h{1}(1) h{2}(1) h{3}(1) h{4}(1) h{5}(1)],ecotypeList2)

subplot(2,2,4)
cellsPop = CruiseData.Pro(Gridding.stationsVec2,:)'; % cells ml-1
cellsPopInt = trapz(Gridding.depthVec,cellsPop*1e9); %cells m-2
plot(cellsPopInt, zPopInt,'.k','MarkerSize',20);
hold on
for a = 1:numel(ecotypeList)
    cellsEco = CruiseData.(ecotypeList{a})(Gridding.stationsVec2,:)';
    cellsEco(find(isnan(cellsEco))) = 0;
    cellsEcoInt = trapz(Gridding.depthVec,cellsEco*1e9); %cells m-2
    h{a} = plot(cellsEcoInt,zEcoInt2(:,a), '.','MarkerEdgeColor',ecotypeColors(a,:),'MarkerFaceColor',ecotypeColors(a,:),'MarkerSize',15);
    hold on
end
%xlim([0 100])
ylim([1e-5 10])
set(gca,'FontSize',20)
set(gca,'YScale','log','XScale','log')
xlabel('Int. Abundance (10-200m) [cells m-2]')
ylabel('Int. NPP (10-200m) [mmol m-2 h-1]')


% NPP vs 14C PP
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
    PPint{a} = trapz(zVec2,PPvec2);
end

figure
plot(PP_lat,cell2mat(PPint),'-k','LineWidth',3)
hold on
plot(CruiseData.Lat(Gridding.stationsVec2),8*biomassEcoInt,'-r','LineWidth',3)




% ATP/NAD(P)H
tempMets = [{'ATP'},{'NADH'},{'NADPH'}];
for i = 1:numel(tempMets)
    tempMets_idx(i) = find(strcmp(tempMets{i},PanGEM.mets));
    tempMets_rxns{i} = find(PanGEM.S(tempMets_idx(i),:));
end
% compute flux sums for each
ATP_NAD_ratio = nansum(abs(PopulationSolution.Fluxes(:,:,tempMets_rxns{1})),3) ./ ( nansum(abs(PopulationSolution.Fluxes(:,:,tempMets_rxns{2})),3) + nansum(abs(PopulationSolution.Fluxes(:,:,tempMets_rxns{3})),3) );
figure
imagesc(x,y,ATP_NAD_ratio);
colormap('jet');
colorbar
ylabel('Depth')
xlabel('Latitude')
set(gca,'FontSize',20)
caxis([0.25 1.25])

% Oxygen CO2 and bicarb
GasEx_rxns = [{'O2EX'},{'CO2EX'},{'HCO3EX'}];
for i = 1:numel(GasEx_rxns)
    GasEx_rxns_idx(i) = find(strcmp(GasEx_rxns{i},PanGEM.rxns));
end
x =  CruiseData.Lat(Gridding.stationsVec2);
y = Gridding.depthVec;
z = -(PopulationSolution.Fluxes(:,:,GasEx_rxns_idx(2)) + PopulationSolution.Fluxes(:,:,GasEx_rxns_idx(3)) ) ./ PopulationSolution.Fluxes(:,:,GasEx_rxns_idx(1));

figure
imagesc(x,y,z);
colormap('jet');
colorbar
ylabel('Depth')
xlabel('Latitude')
set(gca,'FontSize',20)




