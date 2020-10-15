function [StrMod_pigPhysOpt] = updateModel_pigPhysOpt(StrMod,xOut,a,crudeFractions)

% Update BOF coefficients
crudeBOFind = find(strcmp('BIOMASSCRUDE',StrMod.rxns)); % Crude BOF index
for i = 1:numel(crudeFractions)
    crudeFractionsInd(i) = find(strcmp(crudeFractions{i},StrMod.mets));
end
StrMod.S(crudeFractionsInd,crudeBOFind) = -abs(xOut);

% Update absorption reactions bounds
absRxns = [{'DVChla_abs'},{'DVChlb_abs'},{'aCar_abs'},{'Zeax_DVChla_abs'},{'Zeax_DVChlb_abs'}];
for i = 1:numel(absRxns)
    absRxns_ind(i) = find(strcmp(absRxns{i},StrMod.rxns));
end

% DVChla
a_DVChla = xOut(12).*a(1);
StrMod.lb(absRxns_ind(1)) = min([a_DVChla; 300]); % mmol photons gDW-1 h-1
StrMod.ub(absRxns_ind(1)) = min([a_DVChla; 300]); % mmol photons gDW-1 h-1
% DVChlb
a_DVChlb = xOut(13).*a(2);
StrMod.lb(absRxns_ind(2)) = min([a_DVChlb; 300]); % mmol photons gDW-1 h-1
StrMod.ub(absRxns_ind(2)) = min([a_DVChlb; 300]); % mmol photons gDW-1 h-1
% aCar
a_aCar = xOut(14).*a(3);
StrMod.lb(absRxns_ind(3)) = min([a_aCar; 300]); % mmol photons gDW-1 h-1
StrMod.ub(absRxns_ind(3)) = min([a_aCar; 300]); % mmol photons gDW-1 h-1
% Zeax_DVChla
a_Zeax = xOut(15).*a(4);
a_Zeax_DVChla = a_Zeax .* (a_DVChla ./ (a_DVChla+a_DVChlb));
StrMod.lb(absRxns_ind(4)) = min([a_Zeax_DVChla; 300]); % mmol photons gDW-1 h-1
StrMod.ub(absRxns_ind(4)) = min([a_Zeax_DVChla; 300]); % mmol photons gDW-1 h-1
% Zeax_DVChlb
a_Zeax_DVChlb = a_Zeax .* (a_DVChlb ./ (a_DVChla+a_DVChlb));
StrMod.lb(absRxns_ind(5)) = min([a_Zeax_DVChlb; 300]); % mmol photons gDW-1 h-1
StrMod.ub(absRxns_ind(5)) = min([a_Zeax_DVChlb; 300]); % mmol photons gDW-1 h-1

StrMod_pigPhysOpt = StrMod;

end
