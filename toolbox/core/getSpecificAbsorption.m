function PigDB2 = getSpecificAbsorption(PigDB,PigsIncluded,PigMW,lambda)

%% Assemble pigment absorption from database
% Retrieves specific absorption data from the Pigment database (PigDB)
% compiled by Bricaud et al., 2004. Resizes to a specified wavelength range
% and bandwidth which is to be matched by the downwelling irradiance data
% for the environmental simulation. 
nLambda = numel(lambda);
nPigs = numel(PigsIncluded);

% fix the BS problem with Matlab2017b
PigDBPigs = PigDB.Properties.VariableNames;
if ~isempty(find(strcmp('x___lambda',PigDBPigs)))
    PigDB.lambda = PigDB.x___lambda;
end

maxLambdaInd = find(PigDB.lambda == max(lambda));
PigDB(maxLambdaInd+1:end,:) = [];
% get pigs indices
for i = 1:nPigs
    PigDBPigs_ind(i) = find(strcmp(PigsIncluded{i},PigDBPigs));
end

% reshape pig database to the wavelength vector size
PigDat = zeros(nLambda,nPigs);
StartInd = find(lambda==min(PigDB.lambda));
EndInd = find(lambda==max(PigDB.lambda));
PigDat(StartInd:EndInd,:) = table2array(PigDB(:,PigDBPigs_ind)); % m2 mg-1 nm-1

% convert units
PigDB2 = PigDat .* repmat(PigMW,nLambda,1); % m2 mmol-1 nm-1

end
