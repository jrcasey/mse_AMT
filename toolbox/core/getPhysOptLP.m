function [physOptLP] = getPhysOptLP(StrMod, Constraints, a, PigsIncluded, maxIter)

%% Set up physOpt and pigOpt LP


%% Indexing
% parse entries
clb_ind = find(contains(Constraints.Constraint,'clb_'));
cub_ind = find(contains(Constraints.Constraint,'cub_'));

% crude fraction names
crudeFractions = strrep(Constraints.Constraint(clb_ind),'clb_','');

% BOF index
crudeBOFind = find(strcmp('BIOMASSCRUDE',StrMod.rxns));

% Get crude fraction indices and inital values
for i = 1:numel(crudeFractions)
    crudeFractionsInd(i) = find(strcmp(crudeFractions{i},StrMod.mets));
    crudeFractionsOpt(i) = abs(full(StrMod.S(crudeFractionsInd(i),crudeBOFind)));
end

% Get ecotype index
strName = StrMod.id;
fileName = 'data/db/orgDatabase.csv';
orgDatabase = readtable(fileName,'ReadVariableNames',true,'Delimiter',',');
orgInd = find(contains(orgDatabase.StrainName,strName));
ecotype = orgDatabase.Ecotype{orgInd};
ecotype_ind = find(strcmp(ecotype,Constraints.Properties.VariableNames));

% Get pigment indices 
for i = 1:numel(PigsIncluded)
    A_PigInd(i) = find(strcmp(PigsIncluded{i},crudeFractions));
end
%% Set up LP
% Ax<b
A = zeros(4,numel(crudeFractions));
% Assign Plb and Pub vectors
A(1,A_PigInd) = 1;
A(2,A_PigInd) = -1;
b(1) = Constraints.(ecotype)(2); %Plb
b(2) = Constraints.(ecotype)(1); %Pub
b(3) = 0;
b(4) = 0;
% Assign Rlb and Rub vectors
A(3,A_PigInd(1)) = Constraints.(ecotype)(3); % Rlb
A(3,A_PigInd(2)) = -1;
A(4,A_PigInd(1)) = -Constraints.(ecotype)(4); % Rub
A(4,A_PigInd(2)) = 1;

% xlb and xub
xlb = Constraints.(ecotype)(clb_ind);
xub = Constraints.(ecotype)(cub_ind);

% x0 initial
x0 = crudeFractionsOpt;
x0(1) = xub(1);
% nl constraint (simplex)
nlinconst2 = @nlinconst;

options = optimoptions('fmincon','ConstraintTolerance',1e-6,'MaxIterations',maxIter);

%% Save output to structure
physOptLP = struct;
physOptLP.objective = @(x)pigPhysOpt(x,a,StrMod);
physOptLP.x0 = x0;
physOptLP.Aineq = A;
physOptLP.bineq = b;
physOptLP.Aeq = [];
physOptLP.beq = [];
physOptLP.lb = xlb;
physOptLP.ub = xub;
physOptLP.nonlcon = nlinconst2;
physOptLP.solver = 'fmincon';
physOptLP.options = options;


end




