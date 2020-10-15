function S_star = calc_S_star(StrMod2,r_opt,TpDat,T)
%% Compute S_star for each exoMet after TpOpt


%% Parameters
% calculate nG for each TpRxn (uses the customized StrMod2 under optimal
% growth conditions)
nExoMetRxns = numel(TpDat.TpRxns);
StrMod2_sol = solveLP(StrMod2,1);

for a = 1:nExoMetRxns
    rxnIdx = find(strcmp(TpDat.TpRxns{a},StrMod2.rxns));
    flux = StrMod2_sol.x(rxnIdx); % mmol gDW-1 h-1
    DW = 1e-15 * 4/3.*pi().*(r_opt.^3) .* 280; % g DW cell-1
    flux2(a) = abs(flux) .* (1/1000) .* (DW) .* (1/3600);
    if a > 2
        nG(a) = flux2(1) ./ TpDat.kcat(a);
    else
        nG(a) = flux2(a) ./ TpDat.kcat(a);
    end
end



nExoMets = numel(TpDat.TpMets);
for a = 1:nExoMets
    MW = TpDat.TpMets_MW(a); % g mol-1
    D = getDiffusivity(MW,T); % m2 s-1
    kcat = TpDat.kcat(a); % mol Tp-1 s-1
    A = TpDat.A(a); % m2
    r = 1e-6*r_opt; % m
    u=0; % m s-1
    
    mtc = D./r + u/2; % m s-1
    alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless
    ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3
    
    kprime = (kcat./(r .* D));
    if nG(a) ~= 0
        a1 = 1./(kprime.*nG(a));
        b1 = -2;
        c1 = -ks_p;
        p1 = [a1 b1 c1];
        y0 = roots(p1);
        S_star(a) = y0(1); % mol m-3
    else
        S_star(a) = NaN;
    end
end

end
