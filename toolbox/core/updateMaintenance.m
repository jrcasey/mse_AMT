function [GAM, NGAM, newModel] = updateMaintenance(model,cell_r,cellDW)
% Use the empirical relationship described by Lynch and Marinov, 2015. PNAS
% to compute the growth and non-growth associated maintenance costs.
% Updates the model.

% Inputs
% model             Stucture. GEM
% cell_r            Double. Cell radius (m cell-1)
% cellDW            Double. Cell dry weight (g cell-1)

% Outputs
% GAM               Double. Growth associated maintance (mmol ATP gDW-1
%                   h-1)
% NGAM              Double. Non-growth associated maintance (mmol ATP h-1)
% newModel          Structure. GEM with updated GAM and NGAM

%% Calculate volume and dry weight
V = (4/3).* pi() .* ((1e6.*cell_r).^3); % Verity conversion (microns^3)
DW = 1e15 .* cellDW .* 2; % Luria, 1960 (fg cell-1)

%% NGAM (and convert to SBML units)
Cm = 0.39 * V.^0.88; % 10^9 ATP molecules cell-1 h-1
NGAM = Cm .*(1e9) .* (1/6.022e23) .* (1000) .* (1./DW) .* (1e15); % mmol ATP gDW-1 h-1

%% GAM (and convert to SBML units)
Cg = 26.92.*V.^0.97;
preGAM = Cg .* (1e9) .* (1/6.022e23) .* (1000) .* (1./DW) .* (1e15);

% Remove protein synthesis costs and metabolic costs, which are
% already accounted for in the models. Calculated as predicted GAM minus 1/2 sum
% of ATP flux, normalized to the growth rate, minus protein synthesis ATP
% cost, minus NGAM
sol = solveLP(model,1);
ATPind = find(strcmp('ATP',model.mets));
ATPrxnInd = find(full(model.S(ATPind,:)));
ATPFluxSum = 1/2*sum(abs(sol.x(ATPrxnInd)));
ATPFluxSumNorm = ATPFluxSum ./ -sol.f;
prtSynthInd = find(strcmp('ProteinSynthesis',model.rxns));
prtATP = abs(full(model.S(ATPind,prtSynthInd)));

%GAM = preGAM - (ATPFluxSumNorm - prtATP - NGAM);
GAM = preGAM - (prtATP + NGAM);

%% Update model
newModel = model;
crudeBOFind = find(strcmp('BIOMASSCRUDE',model.rxns));
NGAMind = find(strcmp('Maintenance',model.rxns));
metsToChangeInd = getIndexes(model,[{'ATP'},{'ADP'},{'Orthophosphate'},{'H2O'}],'mets');
% ATP + H20 => ADP + Orthophosphate... assign the LHS as negative and RHS
% as positive
newModel.S(metsToChangeInd([1 4]),crudeBOFind) = -GAM;
newModel.S(metsToChangeInd([1 4]),NGAMind) = -NGAM;

newModel.S(metsToChangeInd([2 3]),crudeBOFind) = GAM;
newModel.S(metsToChangeInd([2 3]),NGAMind) = NGAM;

%fprintf(' GAM is now %s mmol gDW-1 h-1.',mat2str(round(GAM,1)));
%fprintf('\n NGAM is now %s mmol gDW-1 h-1.',mat2str(round(NGAM,2)));
%fprintf('\n Model updated \n');





end