%% AMT Analysis
% Retrieves results from the MIT cluster, compiles them into a single structure, and then does a bunch of analyses

%% Load dependencies

load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/Gridding.mat');
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FileNames.mat');
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/CruiseData.mat');
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/GEM/PanGEM.mat');

%% Compile results (this is done on the server now with compile_MESO_SCOPE_subjobs and get_MESO_SCOPE_Results_server)
% ResultsDirectory = 'CBIOMES/Data/Environmental_Data/Cruises/MESO-SCOPE/mse_Results/';
% [FullSolution] = get_MESO_SCOPE_Results(ResultsDirectory,PanGEM,FileNames,Gridding,CruiseData);
% 
% save('CBIOMES/Data/Environmental_Data/Cruises/MESO-SCOPE/FullSolution.mat','FullSolution');

%% Load compiled results
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution.mat');
%% Assign strains to ecotypes
orgDatabase = readtable('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/db/orgDatabase.csv','Delimiter',',','ReadVariableNames',true);
ecotypeList = [{'HLI'},{'HLII'},{'LLI'},{'LLII_LLIII'},{'LLIV'}];
for i = 1:Gridding.nStr
    strainInd = find(strcmp(Gridding.strNameVec{i},orgDatabase.StrainName));
    ecotype{i} = orgDatabase.Ecotype{strainInd};
    ecotypeInd(i) = find(strcmp(ecotype{i},ecotypeList));
end


%% Clean up bad data
% This is sketchy af but it seems to be TpOpt that is failing. Bad data
% seem to show up in the phosphate transporters as being pegged to the
% upper bound. So let's filter those out and interpolate
for i = 1:Gridding.nStr
    maxTp = nanmax(nanmax(FullSolution.(Gridding.strNameVec{i}).TpOpt(:,:,3)));
    delta = 0.02; % allow some room here
    badInd = find(FullSolution.(Gridding.strNameVec{i}).TpOpt(:,:,3) > maxTp*(1-delta));
    [badInd_z, badInd_station] = ind2sub(size(FullSolution.(Gridding.strNameVec{i}).TpOpt(:,:,3)),badInd);
    % now replace with nans
    for j = 1:numel(badInd)
        FullSolution.(Gridding.strNameVec{i}).Growth(badInd_z(j),badInd_station(j)) = NaN;
        FullSolution.(Gridding.strNameVec{i}).Fluxes(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{i}).Shadow(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{i}).BOF_coefs(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{i}).TpOpt(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{i}).r_opt(badInd_z(j),badInd_station(j)) = NaN;
        FullSolution.(Gridding.strNameVec{i}).pigAbs(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{i}).uptakeBounds(badInd_z(j),badInd_station(j),:) = NaN;
    end
    clear badInd badInd_station badInd_z
    % now find the nans and interpolate
    for j = 1:Gridding.nStations
        % Growth
        x=FullSolution.(Gridding.strNameVec{i}).Growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{i}).Growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % Fluxes
        for k = 1:size(FullSolution.(Gridding.strNameVec{i}).Fluxes(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{i}).Fluxes(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{i}).Fluxes(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % Shadow prices
        for k = 1:size(FullSolution.(Gridding.strNameVec{i}).Shadow(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{i}).Shadow(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{i}).Shadow(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % BOF
        for k = 1:size(FullSolution.(Gridding.strNameVec{i}).BOF_coefs(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{i}).BOF_coefs(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{i}).BOF_coefs(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % TpOpt
        for k = 1:size(FullSolution.(Gridding.strNameVec{i}).TpOpt(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{i}).TpOpt(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{i}).TpOpt(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % r_opt
        x=FullSolution.(Gridding.strNameVec{i}).r_opt(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{i}).r_opt(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % pigAbs
        for k = 1:size(FullSolution.(Gridding.strNameVec{i}).pigAbs(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{i}).pigAbs(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{i}).pigAbs(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % uptakeBounds
        for k = 1:size(FullSolution.(Gridding.strNameVec{i}).uptakeBounds(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{i}).uptakeBounds(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{i}).uptakeBounds(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % StrMod1 growth
        x=FullSolution.(Gridding.strNameVec{i}).StrMod1_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{i}).StrMod1_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod2 growth
        x=FullSolution.(Gridding.strNameVec{i}).StrMod2_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{i}).StrMod2_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod3 growth
        x=FullSolution.(Gridding.strNameVec{i}).StrMod3_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{i}).StrMod3_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod4 growth
        x=FullSolution.(Gridding.strNameVec{i}).StrMod4_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{i}).StrMod4_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod5 growth
        x=FullSolution.(Gridding.strNameVec{i}).StrMod5_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{i}).StrMod5_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
    end
    
end
FullSolution_L2 = FullSolution;
save('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat','FullSolution_L2');
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;
%% Compute elemental stoichiometry
% for k = 1:nStr
%     for i = 1:2
%         for j = 1:numel(Gridding.depthVec)
%             tempMod = FullSolution.(strNameVec{k}).StrMod{i}{j};
%             [MMComposition] = getMMElementalStoichiometry(tempMod);
%             MMCompositionPartialTotal = sum(MMComposition.DW,1); % mmol element gDW-1
%             
%             % add in the pigments
%             PigsIncluded = [{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'},{'alpha_Carotene'},{'Zeaxanthin'}];
%             elementAbbrevs = [{'C'},{'H'},{'N'},{'O'},{'P'},{'S'}];
%             for q=1:numel(PigsIncluded)
%                 PigsInd = find(strcmp(PigsIncluded{q},tempMod.mets));
%                 PigsFormulas = tempMod.metFormulas(PigsInd);
%                 [elements, useMat, exitFlag, MW] = parseFormulas(PigsFormulas);
%                 [junk, junk2, elementInd] = intersect(elementAbbrevs,elements.abbrevs);
%                 PigsComposition = useMat(elementInd);% mol element [mol pig]-1
%                 PigsCoef = squeeze(FullSolution.(strNameVec{k}).BOF(i,j,11+q)); % g gDW-1
%                 PigsContribution(q,:) = PigsCoef .* (1./MW) .* PigsComposition .* 1000; % mmol element gDW-1
%             end
%             PigsContributionTotal = nansum(PigsContribution,1);
%             
%             % Compute total
%             MMCompositionTotal(k,i,j,:) = MMCompositionPartialTotal + PigsContributionTotal;
%             % check for sum
%             TotalDW(k,i,j) = 1e-3.* sum(squeeze(MMCompositionTotal(k,i,j,:))'.*[12.0107 1.00794 14.0067 15.9994 30.973762 32.065]);
%         end
%     end
% end
% %% Trap data
% fileName_traps = 'CBIOMES/Data/Environmental_Data/Cruises/MESO-SCOPE/Other_cruise_data/Traps.csv';
% TrapData = readtable(fileName_traps,'ReadVariableNames',true,'Delimiter',',');
% 
% % Cruise Tracks
% figure
% plot(CruiseData.Lon,CruiseData.Lat,'.k','MarkerSize',20);
% hold on
% ht = text(CruiseData.Lon,CruiseData.Lat,num2cell(CruiseData.Stations));
% plot(TrapData.Lon,TrapData.Lat,'.r','MarkerSize',20);
% grid on
% xlabel('Longitude')
% ylabel('Latitude')
% legend('Cruise Stations','Trap Deployments')
% set(gca,'FontSize',20)
% set(ht,'FontSize',20)

%% Ecotype averages
for i = 1:Gridding.nStr
    tempGrowth(:,:,i) = FullSolution.(Gridding.strNameVec{i}).Growth;
    tempFluxes(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).Fluxes;
    tempShadow(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).Shadow;
    tempBOF(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).BOF_coefs;
    tempTpOpt(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).TpOpt;
    tempr_opt(:,:,i) = FullSolution.(Gridding.strNameVec{i}).r_opt;
    temppigAbs(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).pigAbs;
    tempuptakeBounds(:,:,:,i) = FullSolution.(Gridding.strNameVec{i}).uptakeBounds;
end
% Mean
EcotypeSolution = struct;
for i = 1:numel(ecotypeList)
    EcotypeSolution.(ecotypeList{i}).Growth = nanmean(tempGrowth(:,:,find(ecotypeInd==i)),3);
    EcotypeSolution.(ecotypeList{i}).Fluxes = nanmean(tempFluxes(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).Shadow = nanmean(tempShadow(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).BOF_coefs = nanmean(tempBOF(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).TpOpt = nanmean(tempTpOpt(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).r_opt = nanmean(tempr_opt(:,:,find(ecotypeInd==i)),3);
    EcotypeSolution.(ecotypeList{i}).pigAbs = nanmean(temppigAbs(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).uptakeBounds = nanmean(tempuptakeBounds(:,:,:,find(ecotypeInd==i)),4);
end

% Max
EcotypeSolution = struct;
for i = 1:numel(ecotypeList)
    EcotypeSolution.(ecotypeList{i}).Growth = nanmax(tempGrowth(:,:,find(ecotypeInd==i)),[],3);
    EcotypeSolution.(ecotypeList{i}).Fluxes = nanmax(tempFluxes(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).Shadow = nanmax(tempShadow(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).BOF_coefs = nanmax(tempBOF(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).TpOpt = nanmax(tempTpOpt(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).r_opt = nanmax(tempr_opt(:,:,find(ecotypeInd==i)),[],3);
    EcotypeSolution.(ecotypeList{i}).pigAbs = nanmax(temppigAbs(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).uptakeBounds = nanmax(tempuptakeBounds(:,:,:,find(ecotypeInd==i)),[],4);
end
%% Contour plots (ecotypes)
ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];


% Growth
figure
for i = 1:numel(ecotypeList)
    subplot(2,3,i)
    z = EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    caxis([0 0.1])
    colorbar
    title(ecotypeList2{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% Highest growth ecotype
for i = 1:numel(ecotypeList)
    allGrowthEcotype(:,:,i) = EcotypeSolution.(ecotypeList{i}).Growth;
end
[maxGrowthVal, maxGrowthInd] = nanmax(allGrowthEcotype,[],3);
meanGrowthEcotype = nanmean(allGrowthEcotype,3);
figure
colors = varycolor(numel(ecotypeList));
z = maxGrowthInd;
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
colormap(colors)
colorbar('Ticks',[1 2 3 4 5],'TickLabels',ecotypeList)
title('Dominant ecotype')
xlabel('Station')
ylabel('Depth')
set(gca,'FontSize',20)

% Highest growth rate (all strains)
maxGrowthVal2 = nanmax(tempGrowth,[],3);
figure
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,maxGrowthVal2);
colormap('jet')
colorbar
title('Max growth')
xlabel('Station')
ylabel('Depth')
set(gca,'FontSize',20)

% Growth anomaly
figure
for i = 1:numel(ecotypeList)
    subplot(2,3,i)
    z = EcotypeSolution.(ecotypeList{i}).Growth - meanGrowthEcotype;
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    caxis([-0.03 0.03])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% A BOF coeficients
crudeFractions2 = {'DNA';'Lipid';'Carbohydrate';'Protein';'VitaCofactors';'RNA';'NB';'Ions';'BioPool';'Pigments';'CellWall';'Divinylchlorophyll_a';'Divinylchlorophyll_b';'alpha_Carotene';'Zeaxanthin'};
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,3);
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    %caxis([-0.05 0.05])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end


% Cell Size
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = EcotypeSolution.(ecotypeList{i}).r_opt;
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    %caxis([-0.05 0.05])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end
figure
z = EcotypeSolution.LLI.r_opt;
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
colormap(jet)
%caxis([-0.05 0.05])
colorbar
title('LLI cell size')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)


% b/a ratio
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    %caxis([-0.05 0.05])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% chlorophyll?
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = (EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) + EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12)) .* EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(CruiseData.Lat(Gridding.stationsVec),Gridding.depthVec,z');
    colormap(jet)
    %caxis([-0.05 0.05])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% chlorophyll total?

for i = 1:numel(ecotypeList)
    z_temp = (EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) + EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12)) .* EcotypeSolution.(ecotypeList{i}).Growth;
end
z = sum(z_temp,3);
figure
subplot(1,2,1)
imagesc(Gridding.stationsVec,Gridding.depthVec,z');
colormap(jet)
%caxis([-0.05 0.05])
colorbar
title('chlorophyll synthesis rate?')
xlabel('Station')
ylabel('Depth')
set(gca,'FontSize',20)
subplot(1,2,2)
imagesc(Gridding.stationsVec,Gridding.depthVec,CruiseData.Chl_GFF');
colormap(jet)
%caxis([-0.05 0.05])
colorbar
title('Observed chlorophyll')
xlabel('Station')
ylabel('Depth')
set(gca,'FontSize',20)    


% a flux
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    fluxInd = find(strcmp('R00369',PanGEM.rxns))
    %z = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd);
    z = EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd) ./ EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    %caxis([0 500])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% a transporter 
% 2: ammonia, 3: orthophosphate, 4: nitrate, 5. nitrite
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = EcotypeSolution.(ecotypeList{i}).TpOpt(:,:,3);
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    %caxis([-0.05 0.05])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end


% a pig absorption
% 1 - dvchla 2 - dvchlb 3 - acar 4 - zeax
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = EcotypeSolution.(ecotypeList{i}).pigAbs(:,:,4);
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(z));
    colormap(jet)
    %caxis([-0.05 0.05])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% ratio of two fluxes
% PSII/PSI R09503/R09542pc (multiply PSII by 2)
% 

figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    fluxInd1 = find(strcmp('R03140',PanGEM.rxns));
    fluxInd2 = find(strcmp('R00024',PanGEM.rxns));
    z =  2.*EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd1) ./ EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd2);
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    %caxis([-0.05 0.05])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% ratio of two fluxes
% Quantum yield as moles C mole photons absorbed

figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    fluxInd1 = find(strcmp('R00024',PanGEM.rxns));
    fluxInd2 = find(strcmp('DVChla_abs',PanGEM.rxns));
    fluxInd3 = find(strcmp('DVChlb_abs',PanGEM.rxns));
    fluxInd4 = find(strcmp('aCar_abs',PanGEM.rxns));
    z =  2.*EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd1) ./ ...
        ( EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd2) + ...
        EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd3) + ...
        EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,fluxInd4) )
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    caxis([0 0.125])
    colorbar
    title(ecotypeList2{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% GS-GOGAT
% R00248 
% R00093 
% R00114 
% R00355 
% R00253 (GS) 
% R00372 (Photoresp?) 
% R00021 (ferredoxin mediated)

% A sum of fluxes through a metabolite
FS_met = 'Oxygen_th';
FS_metInd = find(strcmp(FS_met,PanGEM.mets));
FS_metRxns_ind = find(full(PanGEM.S(FS_metInd,:)));

figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = squeeze(nansum(abs(EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,FS_metRxns_ind)),3));    
    %z = squeeze(nansum(abs(EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,FS_metRxns_ind)),3)) ./ EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(jet)
    %caxis([0 150])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end


% Ratio of two metabolite turnover rates
FS_met1 = 'Oxygen';
FS_metInd1 = find(strcmp(FS_met1,PanGEM.mets));
FS_metRxns_ind1 = find(full(PanGEM.S(FS_metInd1,:)));
FS_met2 = 'Photon_e';
FS_metInd2 = find(strcmp(FS_met2,PanGEM.mets));
FS_metRxns_ind2 = find(full(PanGEM.S(FS_metInd2,:)));
% FS_met3 = 'AMP';
% FS_metInd3 = find(strcmp(FS_met3,PanGEM.mets));
% FS_metRxns_ind3 = find(full(PanGEM.S(FS_metInd3,:)));
figure
for i = 1:numel(ecotypeList)
    subplot(1,5,i)
    z = squeeze(nansum(abs(EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,FS_metRxns_ind1)),3)) ./ ...
        squeeze(nansum(abs(EcotypeSolution.(ecotypeList{i}).Fluxes(:,:,FS_metRxns_ind2)),3));
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,(z));
    colormap(jet)
    %caxis([0 1])
    colorbar
    title(ecotypeList{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end



% A cruise variable vs growth
colors = varycolor(numel(ecotypeList));
figure
for i = 1:numel(ecotypeList)
    z =  EcotypeSolution.(ecotypeList{i}).Growth;
    h{i}=plot(CruiseData.Pro(Gridding.stationsVec2,:),z','.','MarkerSize',30,'MarkerFaceColor',colors(i,:));
    hold on
    xlabel('Variable')
    ylabel('Growth rate')
    hold on
end
%set(gca,'XScale','log','YScale','log')
set(gca,'XScale','log');
% legend([h{1} h{2} h{3} h{4} h{5}],ecotypeList)
set(gca,'FontSize',20)

%% Density surface
sigmaT = CruiseData.Density(Gridding.stationsVec2,:)';
targSigmaT = 26.5;
[junk, deltaT_idx] = nanmin(abs(sigmaT - targSigmaT));
figure
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,sigmaT);
hold on
plot(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec(deltaT_idx),'-k','LineWidth',4)
colorbar;
caxis([20 27]);
colormap('jet')

% growth along isopycnal
for i = 1:Gridding.nStations
    maxMu(i) = maxGrowthVal2(deltaT_idx(i),i);
end

figure
plot(CruiseData.Lat(Gridding.stationsVec2),maxMu,'ok')

%% A nice one for the website
% Growth
figure
for i = 1:numel(ecotypeList)
    subplot(2,3,i)
    z = EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z);
    colormap(gca,jet)
    caxis([0 0.1])
    colorbar
    title(ecotypeList2{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

% Highest growth ecotype
for i = 1:numel(ecotypeList)
    allGrowthEcotype(:,:,i) = EcotypeSolution.(ecotypeList{i}).Growth;
end
[maxGrowthVal, maxGrowthInd] = nanmax(allGrowthEcotype,[],3);
meanGrowthEcotype = nanmean(allGrowthEcotype,3);
subplot(2,3,6)
colors = varycolor(numel(ecotypeList));
z = maxGrowthInd;
imagesc(CruiseData.Lat(Gridding.stationsVec).depthVec,z');
colormap(gca,colors)
%colormap(jet)
colorbar('Ticks',[1 2 3 4 5],'TickLabels',ecotypeList2)
title('Dominant ecotype')
xlabel('Station')
ylabel('Depth')
set(gca,'FontSize',20)

%% Another nice one for Mick
figure
subplot(2,1,1)
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,maxGrowthVal2);
colormap('jet')
colorbar
title('Max growth rate, all strains [h^-^1]')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)
subplot(2,1,2)
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,CruiseData.Pro(Gridding.stationsVec2,:)');
colormap('jet')
colorbar
title('Cell abundance [ml^-^1]')
xlabel('Latitude')
ylabel('Depth')
set(gca,'FontSize',20)







%% Plot Cruise data

vars = [{'T'},{'Nitrate'},{'Orthophosphate'},{'PAR'},{'Nitrite'},{'Pro'}];
figure
for i = 1:numel(vars)
    subplot(2,3,i)
    z = CruiseData.(vars{i});
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z');
    colormap(gca,jet)
    %caxis([0 0.1])
    colorbar
    title(vars{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

load('GitHub/mse_AMT/data/envData/IrrDat.mat')
figure
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(IrrDat.PAR(:,Gridding.stationsVec2)))
colormap('jet')
caxis([-1 2])
colorbar

%% Plot Ecotype abundances
vars = [{'HLI'},{'HLII'},{'LLI'},{'LLII_LLIII'},{'LLIV'},{'Pro'}];

% linear
figure
for i = 1:numel(vars)
    subplot(2,3,i)
    z = CruiseData.(vars{i});
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,z');
    colormap(gca,jet)
    %caxis([0 0.1])
    colorbar
    title(vars{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end




% log
figure
for i = 1:numel(vars)
    subplot(2,3,i)
    z = CruiseData.(vars{i});
    imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(z(Gridding.stationsVec2,:))');
    colormap(gca,jet)
    caxis([0 6])
    colorbar
    title(vars{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
    zTotal(:,:,i) = z;
end
subplot(2,3,6)
z2 = squeeze(nansum(zTotal,3));
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,log10(z2(Gridding.stationsVec2,:))')
colormap(gca,jet)
caxis([0 6])
colorbar
title('Total')
xlabel('Station')
ylabel('Depth')
set(gca,'FontSize',20)


% compare ecotype abundance to growth rate
vars = [{'HLI'},{'HLII'},{'LLI'},{'LLII_LLIII'},{'LLIV'}];
figure
for i=1:numel(vars)
    subplot(2,3,i)
    x = CruiseData.(vars{i});
    y = EcotypeSolution.(vars{i}).Growth;
    plot(x(Gridding.stationsVec2,:),y','ok');
    xlabel('Abundance')
    ylabel('Growth rate')
    set(gca,'FontSize',20)
    set(gca,'XScale','log')
end

figure
plot(CruiseData.Pro(Gridding.stationsVec2,:),maxGrowthVal2','ok')
%set(gca,'XScale','log','YScale','log');
set(gca,'XScale','log')
xlabel('Abundance')
ylabel('Growth rate')
set(gca,'FontSize',20)

% Ricker Model
% N_t+1


%% Range of growth rates for all strains by ecotype
EcotypeSolution = struct;
for i = 1:numel(ecotypeList)
    EcotypeSolution.(ecotypeList{i}).Growth = nanmax(tempGrowth(:,:,find(ecotypeInd==i)),[],3) - nanmin(tempGrowth(:,:,find(ecotypeInd==i)),[],3);
    EcotypeSolution.(ecotypeList{i}).Fluxes = nanmax(tempFluxes(:,:,:,find(ecotypeInd==i)),[],4) - nanmin(tempFluxes(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).BOF_coefs = nanmax(tempBOF(:,:,:,find(ecotypeInd==i)),[],4) - nanmin(tempBOF(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).TpOpt = nanmax(tempTpOpt(:,:,:,find(ecotypeInd==i)),[],4) - nanmin(tempTpOpt(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).r_opt = nanmax(tempr_opt(:,:,find(ecotypeInd==i)),[],3) - nanmin(tempr_opt(:,:,find(ecotypeInd==i)),[],3);
    EcotypeSolution.(ecotypeList{i}).pigAbs = nanmax(temppigAbs(:,:,:,find(ecotypeInd==i)),[],4) - nanmin(temppigAbs(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).uptakeBounds = nanmax(tempuptakeBounds(:,:,:,find(ecotypeInd==i)),[],4) - nanmin(tempuptakeBounds(:,:,:,find(ecotypeInd==i)),[],4);
end

ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];

% Growth
figure
for i = 1:numel(ecotypeList)
    subplot(2,3,i)
    z = EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(Gridding.stationsVec,Gridding.depthVec,z');
    colormap(jet)
    caxis([0 0.06])
    colorbar
    title(ecotypeList2{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

%% CV of growth rates for all strains by ecotype
EcotypeSolution = struct;
for i = 1:numel(ecotypeList)
    EcotypeSolution.(ecotypeList{i}).Growth = nanstd(tempGrowth(:,:,find(ecotypeInd==i)),[],3) ./ nanmean(tempGrowth(:,:,find(ecotypeInd==i)),3);
    EcotypeSolution.(ecotypeList{i}).Fluxes = nanstd(tempFluxes(:,:,:,find(ecotypeInd==i)),[],4) - nanmean(tempFluxes(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).BOF_coefs = nanstd(tempBOF(:,:,:,find(ecotypeInd==i)),[],4) - nanmean(tempBOF(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).TpOpt = nanstd(tempTpOpt(:,:,:,find(ecotypeInd==i)),[],4) - nanmean(tempTpOpt(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).r_opt = nanstd(tempr_opt(:,:,find(ecotypeInd==i)),[],3) - nanmean(tempr_opt(:,:,find(ecotypeInd==i)),3);
    EcotypeSolution.(ecotypeList{i}).pigAbs = nanstd(temppigAbs(:,:,:,find(ecotypeInd==i)),[],4) - nanmean(temppigAbs(:,:,:,find(ecotypeInd==i)),4);
    EcotypeSolution.(ecotypeList{i}).uptakeBounds = nanstd(tempuptakeBounds(:,:,:,find(ecotypeInd==i)),[],4) - nanmean(tempuptakeBounds(:,:,:,find(ecotypeInd==i)),4);
end

ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];

% Growth
figure
for i = 1:numel(ecotypeList)
    subplot(2,3,i)
    z = 100 .* EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(Gridding.stationsVec,Gridding.depthVec,z');
    colormap(jet)
    caxis([0 50])
    colorbar
    title(ecotypeList2{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

%% Std of growth rates for all strains by ecotype
EcotypeSolution = struct;
for i = 1:numel(ecotypeList)
    EcotypeSolution.(ecotypeList{i}).Growth = nanstd(tempGrowth(:,:,find(ecotypeInd==i)),[],3);
    EcotypeSolution.(ecotypeList{i}).Fluxes = nanstd(tempFluxes(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).BOF_coefs = nanstd(tempBOF(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).TpOpt = nanstd(tempTpOpt(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).r_opt = nanstd(tempr_opt(:,:,find(ecotypeInd==i)),[],3);
    EcotypeSolution.(ecotypeList{i}).pigAbs = nanstd(temppigAbs(:,:,:,find(ecotypeInd==i)),[],4);
    EcotypeSolution.(ecotypeList{i}).uptakeBounds = nanstd(tempuptakeBounds(:,:,:,find(ecotypeInd==i)),[],4);
end


ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];

% Growth
figure
for i = 1:numel(ecotypeList)
    subplot(2,3,i)
    z = 100 .* EcotypeSolution.(ecotypeList{i}).Growth;
    imagesc(Gridding.stationsVec,Gridding.depthVec,z');
    colormap(jet)
    %caxis([0 0.05])
    colorbar
    title(ecotypeList2{i})
    xlabel('Station')
    ylabel('Depth')
    set(gca,'FontSize',20)
end

%% Messing around with HLI

T = 273.15 + CruiseData.T(Gridding.stationsVec2,:);
OGT_MED4 = 273.15 + 28.97;
%OGT_MED4 = 273.15 + 28.97;
Ea = 5.2733e4; %J mol-1
R = 8.314;
A = 1;

for i = 1:Gridding.nStations
    for j = 1:Gridding.nZ
        OGT_rate = exp(-Ea./(R.*(OGT_MED4)));
        InSitu_rate = exp(-Ea./(R.*(T(i,j)))).*(1-exp(T(i,j)-(OGT_MED4+2)));
        if InSitu_rate <0
            InSitu_rate = 0;
        end
        
        % temperature correction
        TCorr(i,j) = InSitu_rate ./ OGT_rate;
        
    end
end

figure
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,TCorr');
colormap('jet')
colorbar

TCorr1 = TCorr;
TCorr2 = TCorr;
TCorr3 = TCorr1 - TCorr2;
TCorr3(find(TCorr3<0))=0;

figure
imagesc(CruiseData.Lat(Gridding.stationsVec2),Gridding.depthVec,TCorr3');
colormap('jet')
colorbar

%% Contour Plots (all strains)
clear temp
for i = 1:Gridding.nStr
    temp(:,:,i) = FullSolution.(Gridding.strNameVec{i}).Growth;
end
meanGrowth = nanmean(temp,3);


dim1 = ceil(sqrt(Gridding.nStr));
figure
for i = 1:Gridding.nStr
    subplot(dim1,dim1,i)
    z = FullSolution.(Gridding.strNameVec{i}).Growth';
    % deal with nans
    [r, c] = find(isnan(z));
    imagesc(Gridding.stationsVec2,Gridding.depthVec,z);
    mycolormap = [1 1 1; jet(48)];
    colormap(mycolormap)
    %caxis([-0.05 0.05])
    caxis([0 0.12])
    hold on
    text(1,1,Gridding.strNameVec{i})
    colorbar
end

% A flux
dim1 = ceil(sqrt(Gridding.nStr));
figure
for i = 1:Gridding.nStr
    subplot(dim1,dim1,i)
    fluxInd = find(strcmp('PyruvateEX',PanGEM.rxns));
    z = FullSolution.(Gridding.strNameVec{i}).Fluxes(:,:,fluxInd)';
    imagesc(Gridding.stationsVec2,Gridding.depthVec,z);
    text(1,1,Gridding.strNameVec{i})
    colorbar
end   
    
% Ratio of 2 fluxes
dim1 = ceil(sqrt(Gridding.nStr));
figure
for i = 1:Gridding.nStr
    subplot(dim1,dim1,i)
    fluxInd1 = find(strcmp('R01070',PanGEM.rxns));
    fluxInd2 = find(strcmp('R01068',PanGEM.rxns));
    z =  2.*FullSolution.(Gridding.strNameVec{i}).Fluxes(:,:,fluxInd2)' ./ FullSolution.(Gridding.strNameVec{i}).Fluxes(:,:,fluxInd1)';
    imagesc(Gridding.stationsVec,Gridding.depthVec,z);
    text(1,1,Gridding.strNameVec{i})
    colorbar
    %caxis([1.1 1.4])

end   

% A BOF component
dim1 = ceil(sqrt(Gridding.nStr));
figure
for i = 1:Gridding.nStr
    subplot(dim1,dim1,i)
    z = FullSolution.(Gridding.strNameVec{i}).BOF_coefs(:,:,12)';
    imagesc(Gridding.stationsVec,Gridding.depthVec,z);
    text(1,1,Gridding.strNameVec{i})
    colorbar
end   

% chl b/a
dim1 = ceil(sqrt(Gridding.nStr));
figure
for i = 1:Gridding.nStr
    subplot(dim1,dim1,i)
    z = FullSolution.(Gridding.strNameVec{i}).BOF_coefs(:,:,13)' ./ FullSolution.(Gridding.strNameVec{i}).BOF_coefs(:,:,12)';
    imagesc(Gridding.stationsVec,Gridding.depthVec,z);
    text(1,1,Gridding.strNameVec{i})
    colorbar
end   

% a transporter
dim1 = ceil(sqrt(Gridding.nStr));
figure
for i = 1:Gridding.nStr
    subplot(dim1,dim1,i)
    z = FullSolution.(Gridding.strNameVec{i}).TpOpt(:,:,3)';
    imagesc(Gridding.stationsVec,Gridding.depthVec,z);
    text(1,1,Gridding.strNameVec{i})
    colorbar
end   


% cell size
dim1 = ceil(sqrt(Gridding.nStr));
figure
for i = 1:Gridding.nStr
    subplot(dim1,dim1,i)
    z = FullSolution.(Gridding.strNameVec{i}).r_opt';
    imagesc(Gridding.stationsVec,Gridding.depthVec,z);
    text(1,1,Gridding.strNameVec{i})
    colorbar
end   
%% N supply integrated all together
dim1 = ceil(sqrt(Gridding.nStr));
for i = 1:Gridding.nStr
    for j = 1:Gridding.nStations
        fluxInd1 = find(strcmp('NitrateEX',PanGEM.rxns));
        fluxInd1b = find(strcmp('NitriteEX',PanGEM.rxns));
        fluxInd2 = find(strcmp('AmmoniaEX',PanGEM.rxns));
        flux1NotNaNInd = find(~isnan(FullSolution.(Gridding.strNameVec{i}).Fluxes(j,:,fluxInd1)));
        flux1bNotNaNInd = find(~isnan(FullSolution.(Gridding.strNameVec{i}).Fluxes(j,:,fluxInd1b)));
        flux2NotNaNInd = find(~isnan(FullSolution.(Gridding.strNameVec{i}).Fluxes(j,:,fluxInd2)));
        flux1Int = trapz(Gridding.depthVec(flux1NotNaNInd),FullSolution.(Gridding.strNameVec{i}).Fluxes(j,flux1NotNaNInd,fluxInd1));
        flux1bInt = trapz(Gridding.depthVec(flux1bNotNaNInd),FullSolution.(Gridding.strNameVec{i}).Fluxes(j,flux1bNotNaNInd,fluxInd1b));
        flux2Int = trapz(Gridding.depthVec(flux2NotNaNInd),FullSolution.(Gridding.strNameVec{i}).Fluxes(j,flux2NotNaNInd,fluxInd2));
        N_supply(i,j) = (flux1Int+flux1bInt) ./ (flux2Int+flux1Int+flux1bInt);
    end
end


for j = 1:Gridding.nStations
    intNO23(j) = trapz(Gridding.depthVec,1e-6*(CruiseData.Nitrate(j,:)+CruiseData.Nitrite(j,:)));
end

figure
plot(CruiseData.SLA(1:end-1),N_supply(:,1:end-1),'.k','MarkerSize',30);
xlabel('SLA (cm)')
ylabel('Fraction of N supplied by NO2+NO3')
set(gca,'FontSize',20)

figure
plot(intNO23(1:end-1),N_supply(:,1:end-1),'.k','MarkerSize',30);
xlabel('Integrated Nitrate 0-180m')
ylabel('Fraction of N supplied by NO2+NO3')
set(gca,'FontSize',20)

figure
plot(intNO23(1:end-1),CruiseData.SLA(1:end-1),'.k','MarkerSize',30);
xlabel('Integrated Nitrate 0-180m')
ylabel('SLA (cm)')
set(gca,'FontSize',20)
%% Plots
% Exclude solutions where the LP didn't converge. For now just looking at
% the runtimes.
for k = 1:Gridding.nStr
    for i = 1:Gridding.nStations
        qualInd{k}{i} = find(FullSolution.(Gridding.strNameVec{k}).runtime(i,:)>10 & FullSolution.(Gridding.strNameVec{k}).runtime(i,:)<2000);
    end
end
lineVec = [{'-'},{'--'}];
colorVec = [{'r'},{'g'},{'b'},{'m'},{'k'}];
markerVec = [{'.'},{'*'}]
% Growth rate
figure
for k = 1:Gridding.nStr
    for i = 1:Gridding.nStations
        legendEntry{k}{i} = strcat(Gridding.strNameVec{k},'-Station-',num2str(Gridding.stationsVec(i)));
        %plot(FullSolution.(strNameVec{k}).Growth(i,qualInd{k}{i}),-Gridding.depthVec(qualInd{k}{i}),'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        y = -Gridding.depthVec(qualInd{k}{i});
        x = FullSolution.(Gridding.strNameVec{k}).Growth(i,qualInd{k}{i});
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        %plot(FullSolution.(strNameVec{k}).Growth(i,qualInd{k}{i}),-Gridding.depthVec(qualInd{k}{i}),'Marker',markerVec{i},'Color',colorVec{k},'MarkerSize',20,'LineStyle','none');
        %plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        plot(x,y,'-')
        hold on; 
    end
end
xlabel('Growth rate [h^-^1]')
ylabel('Depth [m]')
legendEntries = [legendEntry{:}];
%legend('Station 6 MED4','Station 12 MED4','Station 6 MIT9313','Station 12 MIT9313','Station 6 SB','Station 12 SB','Location','NorthWest')
%legend(legendEntries)
set(gca,'FontSize',20)
clear nanInd

% Protein
figure
for k = 1:Gridding.nStr
    for i = 1:Gridding.nStations
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(Gridding.strNameVec{k}).BOF_coefs(i,qualInd{k}{i},4));
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        %plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        plot(x,y,'-')
        hold on; 
    end
end
xlabel('Protein content [g gDW^-^1]')
ylabel('Depth [m]')
%legend(legendEntries)
set(gca,'FontSize',20)
clear nanInd

% Carbohydrate
figure
for k = 1:nStr
    for i = 1:2
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(strNameVec{k}).BOF(i,qualInd{k}{i},3));
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        hold on; 
    end
end
xlabel('Carbohydrate content [g gDW^-^1]')
ylabel('Depth [m]')
legend(legendEntries)
set(gca,'FontSize',20)
clear nanInd

% Chlorophyll (total)
figure
for k = 1:nStr
    for i = 1:2
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(strNameVec{k}).BOF(i,qualInd{k}{i},13))+squeeze(FullSolution.(strNameVec{k}).BOF(i,qualInd{k}{i},12));
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        hold on; 
    end
end
xlabel('Total Divinlychlorophyll content [g gDW^-^1]')
ylabel('Depth [m]')
legend(legendEntries)
set(gca,'FontSize',20)
clear nanInd

% Chlorophyll b/a ratio
figure
for k = 1:Gridding.nStr
    for i = 1:Gridding.nStations
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(Gridding.strNameVec{k}).BOF_coefs(i,qualInd{k}{i},13))./squeeze(FullSolution.(Gridding.strNameVec{k}).BOF_coefs(i,qualInd{k}{i},12));
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        %plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        plot(x,y,'-')
        hold on; 
    end
end
xlabel('Divinylhlorophyll b/a ratio [ND]')
ylabel('Depth [m]')
%legend(legendEntries)
set(gca,'FontSize',20)
clear nanInd

% A transporter
figure
for k = 1:nStr
    for i = 1:2
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(strNameVec{k}).TpOpt(i,qualInd{k}{i},4));
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        hold on; 
    end
end
xlabel('Transporters per cell')
ylabel('Depth [m]')
legend(legendEntries)
set(gca,'FontSize',20)
clear nanInd

% A flux
fluxInd = find(strcmp('NADPHDHth',PanGEM.rxns));
figure
for k = 1:Gridding.nStr
    for i = 1:Gridding.nStations
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(Gridding.strNameVec{k}).Fluxes(i,qualInd{k}{i},fluxInd));
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        %plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        plot(x,y,'-')
        hold on; 
    end
end
xlabel('Flux [mmol gDW^-^1 h^-^1]')
ylabel('Depth [m]')
%legend(legendEntries)
set(gca,'FontSize',20)
title(PanGEM.rxns(fluxInd),'FontSize',20)
clear nanInd

% GS-GOGAT
% Glutamine synthetase: R00253
% GOGAT: R00114
% GOGAT temporary: R00248
fluxInd1 = find(strcmp('R00253',FullSolution.(strNameVec{1}).StrMod{1}{1}.rxns)); % 0.054
fluxInd2 = find(strcmp('R00248',FullSolution.(strNameVec{1}).StrMod{1}{1}.rxns)); % 0.226
fluxInd3 = find(strcmp('R04173',FullSolution.(strNameVec{1}).StrMod{1}{1}.rxns));
figure
for k = 1:nStr
    for i = 1:2
        subplot(1,3,1)
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(strNameVec{k}).Fluxes(i,qualInd{k}{i},fluxInd1))
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        hold on; 
        xlabel('Flux [mmol gDW^-^1 h^-^1]')
        ylabel('Depth [m]')
        set(gca,'FontSize',20)
        title('Glutamine synthetase','FontSize',20)
        clear nanInd
        subplot(1,3,2)
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(strNameVec{k}).Fluxes(i,qualInd{k}{i},fluxInd2))
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        hold on; 
        xlabel('Flux [mmol gDW^-^1 h^-^1]')
        ylabel('Depth [m]')
        set(gca,'FontSize',20)
        title('Glutamate dehydrogenase','FontSize',20)
        clear nanInd
        subplot(1,3,3)
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(strNameVec{k}).Fluxes(i,qualInd{k}{i},fluxInd3))
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        hold on; 
        xlabel('Flux [mmol gDW^-^1 h^-^1]')
        ylabel('Depth [m]')
        set(gca,'FontSize',20)
        title('Phosphoserine transaminase','FontSize',20)
        clear nanInd
    end
end

legend(legendEntries)

% PQ
fluxInd_x1 = find(strcmp('CO2EX',FullSolution.(strNameVec{1}).StrMod{1}{1}.rxns));
fluxInd_x2 = find(strcmp('HCO3EX',FullSolution.(strNameVec{1}).StrMod{1}{1}.rxns));
fluxInd_x3 = find(strcmp('O2EX',FullSolution.(strNameVec{1}).StrMod{1}{1}.rxns));

figure
for k = 1:nStr
    for i = 1:2
        y = -Gridding.depthVec(qualInd{k}{i});
        x = squeeze(FullSolution.(strNameVec{k}).Fluxes(i,qualInd{k}{i},fluxInd_x3)) ./ ( squeeze(FullSolution.(strNameVec{k}).Fluxes(i,qualInd{k}{i},fluxInd_x1)) + squeeze(FullSolution.(strNameVec{k}).Fluxes(i,qualInd{k}{i},fluxInd_x2)) );
        nanInd = find(isnan(x));
        x(nanInd) = [];
        y(nanInd) = [];
        plot(x,y,'Marker',markerVec{i},'Color',colorVec{k},'LineStyle',lineVec{i},'MarkerSize',20,'LineWidth',3);
        hold on; 
    end
end
xlabel('PQ [mol O2 (mol CO2)^-^1]')
ylabel('Depth [m]')
legend(legendEntries)
set(gca,'FontSize',20)
title(FullSolution.(strNameVec{1}).StrMod{1}{1}.rxns(fluxInd),'FontSize',20)
clear nanInd








% MM area fill


figure
n =1;
for k = 1:nStr
    for i = 1:2
        subplot(nStr,2,n)
        area(Gridding.depthVec(qualInd{k}{i}),squeeze(FullSolution.(strNameVec{k}).BOF(i,qualInd{k}{i},:)));
        n=n+1;
        title(strcat(strNameVec{k},'station_',num2str(i)))
    end
end







% Stoichiometry
figure
n=1
for k = 1:nStr
    for i = 1:2
        subplot(nStr,2,n)
        plot(squeeze(MMCompositionTotal(k,i,qualInd{k}{i},1))./squeeze(MMCompositionTotal(k,i,qualInd{k}{i},3)),Gridding.depthVec(qualInd{k}{i}))
        n=n+1;
        title(strcat(strNameVec{k},'station_',num2str(i)))
    end
end


%% Figure out why there is a big gap in the results for all strains, near the equator
load('GitHub/mse_AMT/data/envData/IrrDat.mat');
for i = 1:numel(IrrDat.Station)
PAR2 = 
figure
subplot(2,2,1)
imagesc(CruiseData.Lat(Gridding.stationsVec),Gridding.depthVec,CruiseData.T')
colorbar
subplot(2,2,2)
imagesc(CruiseData.Lat(Gridding.stationsVec),Gridding.depthVec,EcotypeSolution.HLI.Growth')
colorbar
subplot(2,2,3)
imagesc(CruiseData.Lat(Gridding.stationsVec),Gridding.depthVec,CruiseData.Orthophosphate')
colorbar
subplot(2,2,4)
imagesc(CruiseData.Lat(Gridding.stationsVec),Gridding.depthVec,EcotypeSolution.HLI.Growth')
colorbar