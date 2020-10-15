function [newMod, checkSum] = BOFadjust(model,crudeFraction,newValue)
% input
% model
% crudeFraction - String. name of crudeFraction met
% newValue - Double. g/gDW
%
% Output
% newMod - Struct. GEM with BOF adjusted to reflect changes in the
% composition 
% checkSum

% index of crudeBOF reaction
crudeBOFind = find(strcmp('BIOMASSCRUDE',model.rxns));

% index of crude fraction to change
fractionInd = find(strcmp(crudeFraction,model.mets));

% get all the crude fractions from crudeBOF
crudeBOFcoefs = full(model.S(:,crudeBOFind));
crudeBOFmetsInd = find(crudeBOFcoefs<0); % LHS metInds

% exclude ATP and H2O
ATPind = find(strcmp('ATP',model.mets));
H2Oind = find(strcmp('H2O',model.mets));
crudeBOFmetsInd(crudeBOFmetsInd == ATPind) = [];
crudeBOFmetsInd(crudeBOFmetsInd == H2Oind) = [];

% find the index of the crude fraction being changed within the others
changeFractionInd = find(crudeBOFmetsInd==fractionInd);

crudeBOFmets = model.mets(crudeBOFmetsInd); % LHS mets

% get the original value of the fraction being adjusted
oldValue = crudeBOFcoefs(fractionInd);

% get the change
change = abs(newValue) - abs(oldValue); % positive means other fractions need to be reduced

% get the adjustment to the other fractions
nonTargetCoefs = crudeBOFcoefs(crudeBOFmetsInd);
nonTargetCoefs(changeFractionInd) = [];
sumNonTarget = sum(nonTargetCoefs);
nonTargetProportion = nonTargetCoefs./sumNonTarget;
newCoefs = nonTargetCoefs + change.*nonTargetProportion;
newCoefs2 = [newCoefs(1:changeFractionInd-1); -abs(newValue); newCoefs(changeFractionInd:numel(newCoefs))];

checkSum = sum(newCoefs2);

newMod = model;
newMod.S(crudeBOFmetsInd,crudeBOFind) = newCoefs2;
