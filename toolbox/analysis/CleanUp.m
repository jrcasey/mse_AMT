%% Post Processing
% OptTrans sometimes does not converge, so we'll clean up the L1
% FullSolution by removing bad data and interpolating wherever necessary.
% It seems that solutions that did not converge will pin the number of
% orthophosphate transporters to the upper bound (or thereabouts), so that
% quick-and-dirty approach is implemented here. 

% This script generates the Level 2 FullSolution. 

%% Load compiled results (from compile_AMT_subjobs.m)
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution.mat');
% parse out gridding
Gridding = FullSolution.Gridding;

% need to provide a bit of room below the upper bound
delta = 0.02; % allow some room here

% loop through
for i = 1:Gridding.nStr
    maxTp = nanmax(nanmax(FullSolution.(Gridding.strNameVec{i}).TpOpt(:,:,3)));
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

%% Save output
FullSolution_L2 = FullSolution;
save('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat','FullSolution_L2');

