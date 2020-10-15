function [StrMod_TpOpt, vOut] = updateModel_TpOpt(StrMod,xOut,TpDat)

% Update cell size
StrMod.rInit = xOut(1);

% Get indices
TpMets = TpDat.TpMets;
nTpMets = numel(TpMets);
TpRxns = TpDat.TpRxns;
nTpRxns = numel(TpRxns);

for i = 1:nTpRxns
    TpRxns_ind(i) = find(strcmp(TpRxns{i},StrMod.rxns));
end

% Cell size and mass
r = 1e-6 .* xOut(1); % m cell-1
SA = 4.* pi() .* (r.^2);
q = 1e-15 .* 4/3.*pi().*(r.*1e6).^3 .* 280 .* 2; % g DW cell-1

% Tp abundance, size and mass
n = xOut(2:5);

% Params for paramsArmstrong
kcat = TpDat.kcat;
Vmax = n.*kcat;
A = TpDat.A;
D = TpDat.D;
u = 0;
S = TpDat.S;

% Solve for uptake rates
for i = 1:nTpRxns
    [ks_dl(i), ks_pl(i)] = paramsArmstrong(Vmax(i), kcat(i), n(i), r, A(i), D(i), u);
    v(i) = uptakeArmstrong(S(i),Vmax(i),ks_dl(i), ks_pl(i)); % mole cell-1 s-1
end

% Convert to FBA units
v_FBA = 1000 .* v .* (1 ./ q) .* 3600; % mmol gDW-1 h-1

% Set constraints
StrMod.lb(TpRxns_ind) = 0;
StrMod.ub(TpRxns_ind) = 0;
for i = 1:nTpRxns
    sense = TpDat.TpRxns_Sense(i);
    if sense == -1
        StrMod.lb(TpRxns_ind(i)) = -v_FBA(i);
        StrMod.ub(TpRxns_ind(i)) = 0;
    elseif sense == 1
        StrMod.lb(TpRxns_ind(i)) = 0;
        StrMod.ub(TpRxns_ind(i)) = v_FBA(i);
    end
end

StrMod_TpOpt = StrMod;
vOut = abs(v_FBA);
end