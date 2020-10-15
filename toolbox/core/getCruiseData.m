function CruiseDat = getCruiseData(CruiseDB, depthVec)
%% Standardize AMT Data

Data = CruiseDB;
% target fields (nutrients and cell counts for now)
% varNames = [{'Nitrate'},{'Nitrite'},{'Ammonia'},{'Orthophosphate'},{'Fe2'},...
%     {'Zn2'},{'Cadmium'},{'Ni'},{'Copper'},{'Pb'},{'Mn'},{'Cobalt_ion'},...
%     {'Pro'},{'Syn'},{'Euk'},{'Hbac'}];

varNames = Data.Properties.VariableNames(6:end); % just the nutrients for now
% Find and index stations
uStations = unique(Data.Station);
nStations = numel(uStations);
for i = 1:nStations
    station_ind{i} = find(strcmp(uStations{i},Data.Station));
end

% Interpolate (linear)
CruiseDat = struct;
CruiseDat.Stations = uStations;
CruiseDat.Depth = depthVec;
for i = 1:nStations
    for j = 1:numel(varNames) 
        clear x v vq
        x = Data.Depth(station_ind{i});
        v = Data.(varNames{j})(station_ind{i});
        % deal with the nans
        nanInd = find(isnan(v));
        x(nanInd) = [];
        v(nanInd) = [];
        if numel(find(~isnan(v)))<3
            %nuts(i,j,:) = zeros(numel(depthVec),1);
            CruiseDat.(varNames{j})(i,:) = zeros(numel(depthVec),1);
        else
        vq = interp1(x,v,depthVec,'linear');
        %nuts(i,j,:) = vq;
        CruiseDat.(varNames{j})(i,:) = vq;
        end
    end
    CruiseDat.Lat(i) = Data.Lat(station_ind{i}(1));
    CruiseDat.Lon(i) = Data.Lon(station_ind{i}(1));
    CruiseDat.Date(i) = Data.Date(station_ind{i}(1));
    %CruiseDat.SLA(i) = Data.SLA_cm(station_ind{i}(1));
end

% Change some units to nM
CruiseDat.NitratePlusNitrite = CruiseDat.NitratePlusNitrite .* 1000;
CruiseDat.Nitrite = CruiseDat.Nitrite .* 1000;
CruiseDat.Orthophosphate = CruiseDat.Orthophosphate .* 1000;
CruiseDat.Ammonium = 7 + (CruiseDat.Ammonium .* 1000);
CruiseDat.Silicate = CruiseDat.Silicate .* 1000;
% Compute nitrate
CruiseDat.Nitrate = CruiseDat.NitratePlusNitrite - CruiseDat.Nitrite;
CruiseDat.Nitrate(find(CruiseDat.Nitrate <0)) = 0;

% Some name changes too
CruiseDat.T = CruiseDat.Temp;
CruiseDat.Ammonia = CruiseDat.Ammonium;
rmfield(CruiseDat,{'Temp','Ammonium'});




end