%% Post Processing
% OptTrans sometimes does not converge, so we'll clean up the L1
% FullSolution by removing bad data and interpolating wherever necessary.
% It seems that solutions that did not converge will pin the number of
% orthophosphate transporters to the upper bound (or thereabouts), so that
% quick-and-dirty approach is implemented here. 

% This script generates the Level 2 FullSolution.

% Takes about 15 minutes

%% Load compiled results (from compile_AMT_subjobs.m)
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L1.mat');
FullSolution = FullSolution_L1;
clear FullSolution_L1
% parse out gridding
Gridding = FullSolution.Gridding;

% need to provide a bit of room below the upper bound
delta = 0.02; % allow some room here

% loop through
for a = 1:Gridding.nStr
    maxTp = nanmax(nanmax(FullSolution.(Gridding.strNameVec{a}).TpOpt(:,:,3)));
    badInd = find(FullSolution.(Gridding.strNameVec{a}).TpOpt(:,:,3) > maxTp*(1-delta));
    [badInd_z, badInd_station] = ind2sub(size(FullSolution.(Gridding.strNameVec{a}).TpOpt(:,:,3)),badInd);
    % now replace with nans
    for j = 1:numel(badInd)
        FullSolution.(Gridding.strNameVec{a}).Growth(badInd_z(j),badInd_station(j)) = NaN;
        FullSolution.(Gridding.strNameVec{a}).Fluxes(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{a}).Shadow(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{a}).BOF_coefs(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{a}).TpOpt(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{a}).r_opt(badInd_z(j),badInd_station(j)) = NaN;
        FullSolution.(Gridding.strNameVec{a}).pigAbs(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{a}).uptakeBounds(badInd_z(j),badInd_station(j),:) = NaN;
        FullSolution.(Gridding.strNameVec{a}).S_star(badInd_z(j),badInd_station(j),:) = NaN;        
    end
    clear badInd badInd_station badInd_z
    % now find the nans and interpolate
    for j = 1:Gridding.nStations
        % Growth
        x=FullSolution.(Gridding.strNameVec{a}).Growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{a}).Growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % Fluxes
        for k = 1:size(FullSolution.(Gridding.strNameVec{a}).Fluxes(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{a}).Fluxes(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{a}).Fluxes(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % Shadow prices
        for k = 1:size(FullSolution.(Gridding.strNameVec{a}).Shadow(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{a}).Shadow(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{a}).Shadow(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % BOF
        for k = 1:size(FullSolution.(Gridding.strNameVec{a}).BOF_coefs(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{a}).BOF_coefs(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{a}).BOF_coefs(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % TpOpt
        for k = 1:size(FullSolution.(Gridding.strNameVec{a}).TpOpt(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{a}).TpOpt(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{a}).TpOpt(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % r_opt
        x=FullSolution.(Gridding.strNameVec{a}).r_opt(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{a}).r_opt(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % pigAbs
        for k = 1:size(FullSolution.(Gridding.strNameVec{a}).pigAbs(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{a}).pigAbs(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{a}).pigAbs(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % uptakeBounds
        for k = 1:size(FullSolution.(Gridding.strNameVec{a}).uptakeBounds(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{a}).uptakeBounds(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{a}).uptakeBounds(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % S_star
        for k = 1:size(FullSolution.(Gridding.strNameVec{a}).S_star(:,j,:),3)
            x=FullSolution.(Gridding.strNameVec{a}).S_star(:,j,k);
            nanx = isnan(x);
            if sum(nanx) ~= numel(x) & sum(~nanx) > 1
                t = Gridding.depthVec;
                temp = interp1(t(~nanx),x(~nanx),t);
                FullSolution.(Gridding.strNameVec{a}).S_star(:,j,k) = interp1(t(~nanx),x(~nanx),t);
            end
        end
        % StrMod1 growth
        x=FullSolution.(Gridding.strNameVec{a}).StrMod1_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{a}).StrMod1_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod2 growth
        x=FullSolution.(Gridding.strNameVec{a}).StrMod2_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{a}).StrMod2_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod3 growth
        x=FullSolution.(Gridding.strNameVec{a}).StrMod3_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{a}).StrMod3_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod4 growth
        x=FullSolution.(Gridding.strNameVec{a}).StrMod4_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{a}).StrMod4_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
        % StrMod5 growth
        x=FullSolution.(Gridding.strNameVec{a}).StrMod5_growth(:,j);
        nanx = isnan(x);
        if sum(nanx) ~= numel(x) & sum(~nanx) > 1
            t = Gridding.depthVec;
            temp = interp1(t(~nanx),x(~nanx),t);
            FullSolution.(Gridding.strNameVec{a}).StrMod5_growth(:,j) = interp1(t(~nanx),x(~nanx),t);
        end
    end
    
end

%% Flux consistency
% Some issue has been found where non-zero fluxes are associated with
% reactions which are not actually present. Let's loop through and zero
% these out
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/GEM/StrMod.mat');
for a = 1:Gridding.nStr
    % get indices of zero entries in S
    z_idx = find(sum(abs(StrMod.(Gridding.strNameVec{a}).S),1)==0);
    % assign zero flux to corresponding fluxes
    FullSolution.(Gridding.strNameVec{a}).Fluxes(:,:,z_idx) = 0;
end

%% Save output
FullSolution_L2 = FullSolution;
save('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat','FullSolution_L2');

