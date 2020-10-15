function [IrrDat2] = standardizeIrr(IrrDat,lambdaVec,depthVec)
%% Process HyperPro data from cast dimensions to standardized dimensions
% Transforms a hyperpro profile from native dimensions to standard
% dimensions

% Inputs
% IrrDat                -       Structure. Fields include Station, Depth,
%                               and Data as cell arrays. Downwelling 
%                               irradiance data in units of microW cm-2 nm-1 
% lambda                -       Double. Wavelength vector to transform to
% depth                 -       Double. Depth vector to transform to

% Outputs
% IrrDat2               -       Array. Reformated irradiance data.
%                               Downwelling irradiance data in units of 
%                               micromoles photons m-2 s-1

% 20191028
% John R. Casey

%% Import HyperPro data from ascii files from MESO-SCOPE Cruise (Comment Out)
% 
% fileList_path = 'CBIOMES/Data/Environmental_Data/HyperPro/MESO-SCOPE/ASCII/FileList.txt';
% fileList = readtable(fileList_path,'ReadVariableNames',false,'Delimiter','tab');
% Lambda = 354:2:796; %nm
% for i = 1:numel(fileList)
%     fileName = strcat('CBIOMES/Data/Environmental_Data/HyperPro/MESO-SCOPE/ASCII/',fileList.Var1{i});
%     Dat = readtable(fileName,'ReadVariableNames',true,'Delimiter','tab','HeaderLines',47);
%     Station(i) = i;
%     Depth{i} = Dat.x0_3;
%     Data{i} = Dat(:,9:230);
% end
% IrrDat.Lambda = Lambda;
% IrrDat.Station = Station;
% IrrDat.Depth = Depth;
% IrrDat.Data = Data;
% 
% save('CBIOMES/Data/Environmental_Data/HyperPro/MESO-SCOPE/ASCII/IrrDat.mat','IrrDat');
% 
% clear
% load('CBIOMES/Data/Environmental_Data/HyperPro/MESO-SCOPE/ASCII/IrrDat.mat')
% 



%% Match data to grid
% Interpolate to a grid with input dimensions
minLambda = min(lambdaVec);
maxLambda = max(lambdaVec);
bandWidth = (maxLambda - minLambda)/(numel(lambdaVec)-1);

minDepth = min(depthVec);
maxDepth = max(depthVec);
depthInterval = (maxDepth - minDepth)/(numel(depthVec)-1);

Data = IrrDat.Data;


for i = 1:numel(Data);
    clear X Xg Y Yg Xg2 Yg2 V Xq Yq Xq2 Yq2 F zi
    X = IrrDat.Depth;
    Xg = linspace(min(X),max(X),numel(X));
    Y = IrrDat.Lambda;
    Yg = linspace(min(Y),max(Y),numel(Y));
    [Xg2,Yg2] = ndgrid(Xg,Yg);
    V = Data{i};
    V(isnan(V)) = 0;
    Xq = linspace(minDepth,maxDepth,numel(depthVec));
    %Yq = linspace(348,804,numel(348:bandWidth:804));
    Yq = linspace(minLambda,maxLambda,numel(lambdaVec));

    [Xq2,Yq2] = ndgrid(Xq,Yq);
    F = griddedInterpolant(Xg2,Yg2,V,'linear','none');
    zi = F(Xq2,Yq2);
    Data2{i} = zi;
end

%% Convert W m-2 nm-1 to micromoles photons m-2 s-1
% solve photon energies (J mol-1)
c = 3.0e8; % speed of light (m s-1)
planck = 6.63e-34; % planck's constant (J s)
avogadro = 6.02e23; % Avogadro's number (mol-1)
Elambda = avogadro*c*planck./(Yq./1e9); % J mol-1
for i = 1:numel(Data2)
    % convert to mmol of photons m-2 h-1
    Data3{i} = (bandWidth .* 1e6 .* Data2{i} .* (1./1000) .* (3600)) ./ repmat(Elambda,size(Data2{i},1),1); % depth x wavelength
end

IrrDat2 = Data3; % mmol photons m-2 h-1;


