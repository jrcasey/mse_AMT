function [mu] = pigPhysOpt(x,a,model)

% assign x to BOF mets
crudeBOFind = find(strcmp('BIOMASSCRUDE',model.rxns));
crudeFractions = [{'DNA'},{'Lipid'},{'Carbohydrate'},{'Protein'},...
    {'VitaCofactors'},{'RNA'},{'NB'},{'Ions'},{'BioPool'},...
    {'Pigments'},{'CellWall'},{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'},{'alpha_Carotene'},{'Zeaxanthin'}];
for i = 1:numel(crudeFractions)
    crudeFractionsInd(i) = find(strcmp(crudeFractions{i},model.mets));
end
xInit = full(model.S(crudeFractionsInd,crudeBOFind));
model.S(crudeFractionsInd,crudeBOFind) = -x;
weight = 1e-2;
penalty = weight .* sum((abs(x)-abs(xInit')).^2);

% compute pigment absorption
absRxns = [{'DVChla_abs'},{'DVChlb_abs'},{'aCar_abs'},{'Zeax_DVChla_abs'},{'Zeax_DVChlb_abs'}];
for i = 1:numel(absRxns)
    absRxns_ind(i) = find(strcmp(absRxns{i},model.rxns));
end

% DVChla
a_DVChla = x(12).*a(1);
model.lb(absRxns_ind(1)) = min([a_DVChla; 300]); % mmol photons gDW-1 h-1
model.ub(absRxns_ind(1)) = min([a_DVChla; 300]); % mmol photons gDW-1 h-1
% DVChlb
a_DVChlb = x(13).*a(2);
model.lb(absRxns_ind(2)) = min([a_DVChlb; 300]); % mmol photons gDW-1 h-1
model.ub(absRxns_ind(2)) = min([a_DVChlb; 300]); % mmol photons gDW-1 h-1
% aCar
a_aCar = x(14).*a(3);
model.lb(absRxns_ind(3)) = min([a_aCar; 300]); % mmol photons gDW-1 h-1
model.ub(absRxns_ind(3)) = min([a_aCar; 300]); % mmol photons gDW-1 h-1
% Zeax_DVChla
a_Zeax = x(15).*a(4);
a_Zeax_DVChla = a_Zeax .* (a_DVChla ./ (a_DVChla+a_DVChlb));
model.lb(absRxns_ind(4)) = min([a_Zeax_DVChla; 300]); % mmol photons gDW-1 h-1
model.ub(absRxns_ind(4)) = min([a_Zeax_DVChla; 300]); % mmol photons gDW-1 h-1
% Zeax_DVChlb
a_Zeax_DVChlb = a_Zeax .* (a_DVChlb ./ (a_DVChla+a_DVChlb));
model.lb(absRxns_ind(5)) = min([a_Zeax_DVChlb; 300]); % mmol photons gDW-1 h-1
model.ub(absRxns_ind(5)) = min([a_Zeax_DVChlb; 300]); % mmol photons gDW-1 h-1

%[sol hsSolOut] = solveLP(model,1);
sol = solveLP(model);
if sol.stat==1
    mu = sol.f.*(1 + penalty);
else
    mu = 0;
end

end
