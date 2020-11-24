function soleSources = getSoleElementSources(model, element, auxotrophyRxns)

%% get transporter reactions for each element
transRxnInd = find(strcmp('Transport',model.subSystems));

% get transport mets
[transMetInd, junk] = find(full(model.S(:,transRxnInd)));
% screen out the intracellular mets
extracellularInd = find(model.metComps == 3);
transMetInd2 = intersect(transMetInd,extracellularInd);
transMets = model.mets(transMetInd2);

%% Get auxotrophy reactions indices
excludedRxnsInd = find(contains(model.rxns,auxotrophyRxns));
% Also need to exclude the gasses
excludedRxns2 = [{'CO2TRANS'},{'HCO3TRANS'}];
excludedRxnsInd2 = find(contains(model.rxns,excludedRxns2));
excludedRxnsInd3 = union(excludedRxnsInd,excludedRxnsInd2);
%% parse out formulas
formulas = model.metFormulas(transMetInd2);
[elements, useMat, exitFlag, MW] = parseFormulas(formulas, false,false,true);

% get indx of element
bioElementsInd = find(strcmp(element,elements.abbrevs));

% Get elemental composition of each exogenous metabolite
CInd = find(useMat(:,bioElementsInd(1)));
CSources = transMets(CInd);
Catoms = useMat(CInd,bioElementsInd(1));

% now link back to their respective transporter reactions
for i = 1:numel(CSources)
    temp = find(strcmp(CSources{i},model.mets));
    CSourceRxnInd = find(full(model.S(temp,:)));
    nTemp = numel(CSourceRxnInd);
    % make sure it's a transporter
    CSourceRxnInd2{i} = intersect(CSourceRxnInd,transRxnInd);
    CSourceRxns{i} = model.rxns(CSourceRxnInd2{i});
    if nTemp ==1
        CatomRxns{i} = Catoms(i);
    else
    CatomRxns{i} = repmat(Catoms(i),nTemp-1,1);
    end
end
CSourceRxns2a = vertcat(CSourceRxns{:});
CSourceRxnInd3a = vertcat(CSourceRxnInd2{:});
CatomRxn2a = vertcat(CatomRxns{:});

%% Exclude required transport reactions reactions
[CSourceRxnInd3, ia] = setdiff(CSourceRxnInd3a,excludedRxnsInd3);
CSourceRxns2 = CSourceRxns2a(ia);
CatomRxn2 = CatomRxn2a(ia);

%% Make sure that the model cannot grow without any source of each element
sol2 = solveLP(model,1);

% Growth on no carbon source
ClimMod = model;
ClimMod.lb(CSourceRxnInd3) = 0;
noCsol = solveLP(ClimMod,1);

% all ok?
threshold = 1e-8;
if abs(noCsol.f) < threshold 
    magicGrowth = 'No'
else magicGrowth = 'Yes';
end

%% Growth on each substrate as sole element source
threshold = -1e-8;
for i = 1:numel(CSourceRxnInd3)
    tempMod = ClimMod;
    tempMod.lb(CSourceRxnInd3(i)) = -1;
    tempSol = solveLP(tempMod,1);
    if tempSol.f < threshold
        Cgrowth(i) = tempSol.f;
        Cyield(i) = tempSol.f./(tempSol.x(CSourceRxnInd3(i))./CatomRxn2(i));
    else Cgrowth(i) = NaN;
        Cyield(i) = NaN;
    end
end

%% Store output
soleSources = struct;
soleSources.TpRxns = CSourceRxns2;
soleSources.Element = element;
soleSources.Atoms = CatomRxn2;
soleSources.Growth = Cgrowth;
soleSources.Yield = Cyield;
soleSources.MagicGrowth = magicGrowth;

end