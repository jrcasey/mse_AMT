function [TpDat2] = getTpOptParams(CruiseData,StrMod,TpDat,station_ind,depth_ind)

% Computes, formats and compiles S, A, u, D, kcat, SA, and T into a
% structure for input into TpOpt_MultiVar

%% Compile indices
TpMets = TpDat.TpMets;
nTpMets = numel(TpMets);
TpRxns = TpDat.TpRxns;
nTpRxns = numel(TpRxns);

% EnvFields = fieldnames(CruiseData);
% for i = 1:nTpMets
%     TpMets_ind(i) = find(strcmp(TpMets{i},EnvFields));

%% Compile environmental data
% Get temperature
TpDat.T = repmat(CruiseData.T(station_ind,depth_ind),nTpRxns,1); % degC
% Get substrate concentrations and hydrated molecular diffusivity
TpDat.S = zeros(nTpMets,1);
for i = 1:nTpMets
    TpDat.S(i) = 1e-6.*CruiseData.(TpMets{i})(station_ind,depth_ind); % mol m-3
    TpMets_MW = TpDat.TpMets_MW(i); % g mol-1
    TpDat.D(i) = getDiffusivity(TpMets_MW,TpDat.T(i)); % m2 s-1
end
TpDat2 = TpDat;

end

