function [mu] = TpOpt(x,StrMod,TpDat)

% LP objective function for TpOpt

% Get indices
TpMets = TpDat.TpMets;
nTpMets = numel(TpMets);
TpRxns = TpDat.TpRxns;
nTpRxns = numel(TpRxns);

for i = 1:nTpRxns
    TpRxns_ind(i) = find(strcmp(TpRxns{i},StrMod.rxns));
end

% Cell size and mass
r = 1e-6 .* x(1); % m cell-1
SA = 4.* pi() .* (r.^2); % m2 cell-1
q_native = 1e3.* 4/3.*pi().*(r.^3) .* 280 .* 2; % g DW cell-1

% Tp abundance, size and mass
n = x(2:5); % Tp cell-1
TpMass = n.* TpDat.TpRxns_MW' .* (1/6.022e23); % g Tp cell-1
% add to quota 
q = q_native + 1000.*sum(TpMass); % g DW cell-1

% Params for paramsArmstrong
kcat = 0.1.*TpDat.kcat'; % mol Tp-1 s-1
Vmax = n.*kcat; % mol cell-1 s-1
A = TpDat.A;
D = 10.*TpDat.D;
u = 0;
S = TpDat.S;
%S2 = logspace(-7,-1,100)
% Solve for uptake rates
for i = 1:nTpRxns
    [ks_dl(i), ks_pl(i)] = paramsArmstrong(Vmax(i), kcat(i), n(i), r, A(i), D(i), u); % mol m-3
    v(i) = uptakeArmstrong(S(i),Vmax(i),ks_dl(i), ks_pl(i)); % mole cell-1 s-1
%     for j = 1:numel(S2)
%        v(i,j) = uptakeArmstrong(S2(j),Vmax(i),ks_dl(i), ks_pl(i)); % mole cell-1 s-1
%        v_FBA(i,j) = 1000 .* v(i,j) .* (1 ./ q) .* 3600;
%     end
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

% Solve FBA
%sol = solveLP(StrMod,1);
sol = solveLP(StrMod);
if sol.stat==1
    mu = sol.f;
else
    mu = 0;
end


end










