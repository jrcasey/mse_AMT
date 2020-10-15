function [TpOptLP] = getTpOptLP(StrMod, fmax, TpDat, maxIter)

SA = 4.* pi() .* ((1e-6 .* StrMod.rInit).^2);

A = [0 TpDat.A']; % m2 transporter area 
b = fmax.*SA; 
x0 = [StrMod.rInit repmat(0.2,1,numel(TpDat.A))];
xlb = [StrMod.rInit.*0.9 zeros(1,numel(TpDat.A))];
xub = [StrMod.rInit.*1.1 (b ./ TpDat.A')];

options = optimoptions('fmincon','ConstraintTolerance',1e-3,'MaxIterations',maxIter);


% Compile TpOpt structure for the LP solver
TpOptLP = struct;
TpOptLP.objective = @(x)TpOpt(x,StrMod,TpDat);
TpOptLP.x0 = x0;
TpOptLP.Aineq = A;
TpOptLP.bineq = b;
TpOptLP.Aeq = [];
TpOptLP.beq = [];
TpOptLP.lb = xlb;
TpOptLP.ub = xub;
TpOptLP.nonlcon = [];
TpOptLP.solver = 'fmincon';
TpOptLP.options = options;

end