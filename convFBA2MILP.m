function model = convFBA2MILP(model, ReactionDB)

bigM = 1e6;
% % save the original direction reversibilities
[num_mets_org, num_rxns] = size(model.S);

% formatting the metabolite and reaction names to remove brackets
for i=1:num_mets_org
    newmetname = model.mets{i};
    newmetname = strrep(newmetname,'[','_');
    newmetname = strrep(newmetname,']','');
    newmetname = strrep(newmetname,'(','_');
    newmetname = strrep(newmetname,')','');
    model.mets{i} = newmetname;
end

for i=1:num_rxns
    newrxnname = model.rxns{i};
    newrxnname = strrep(newrxnname,'(','_');
    newrxnname = strrep(newrxnname,')','');
    newrxnname = strrep(newrxnname,'[','_');
    newrxnname = strrep(newrxnname,']','');
    model.rxns{i} = newrxnname;
end

% if all reversible parameter is indicated then we assume all reactions are reversible first except biomass
% and create the Irrev version of the model

% create the A matrix using the S matrix first

[modelIrrev, ~, ~, ~] = convertToIrreversibleAgador(model);
model.A = modelIrrev.S;
[num_mets, num_vars] = size(modelIrrev.S);
objective = modelIrrev.rxns(find(modelIrrev.c));

% check that num of metabolites remain the same
if num_mets ~= num_mets_org
    error('number of metabolites do not match!')
end

% Initialize fields (in case of existing, they are erased)
model.constraintType = [];
model.constraintNames = [];
model.rhs = [];
model.varNames = [];
model.vartypes = [];

% create the constraint type vector first for mass balances
% and the rhs vector first for mass balances
for i = 1:num_mets
    model.constraintType{i,1} = '=';
    model.constraintNames{i,1} = strcat('M_',model.mets{i});
    model.rhs(i,1) = 0;
end

% create the variable type vector and setting their bounds
% and lower and upper bounds for the flux variables
% upperbounds should not be negative
model.var_lb = zeros(num_vars,1);
model.var_ub = zeros(num_vars,1);
for i = 1:num_vars
    if modelIrrev.ub(i) < 0
        modelIrrev.ub(i) = 0;
    end
    if modelIrrev.lb(i) < 0
        modelIrrev.lb(i) = 0;
    end
    model.var_ub(i) = modelIrrev.ub(i);
    model.var_lb(i) = modelIrrev.lb(i);
    model.varNames = modelIrrev.rxns;
    model.vartypes{i,1} = 'C';
end

for i = 1:num_rxns
    
    F_flux_index = find(ismember(model.varNames, strcat('F_', model.rxns{i})));
    R_flux_index = find(ismember(model.varNames, strcat('R_', model.rxns{i})));
    
    % We DON'T add thermodynamic constraints! We only add constraints
    % for the forward and reverse fluxes.
    fprintf('generating only use constraints for reaction %s\n', model.rxns{i});
    
    model = addNewVariableInTFA(model, strcat('FU_', model.rxns{i}), 'B', [0 1]);
    FU_index = size(model.varNames,1);
    model = addNewVariableInTFA(model, strcat('BU_', model.rxns{i}), 'B', [0 1]);
    BU_index = size(model.varNames,1);
    % create the prevent simultaneous use constraints
    % SU_rxn: FU_rxn + BU_rxn <= 1
    CLHS.varIDs    = [FU_index  BU_index];
    CLHS.varCoeffs = [+1        +1      ];
    model = addNewConstraintInTFA(model, strcat('SU_', model.rxns{i}), '<', CLHS, 1);
    % create constraints that control fluxes with their use variables
    % UF_rxn: F_rxn - 1000 FU_rxn < 0
    CLHS.varIDs    = [F_flux_index  FU_index];
    CLHS.varCoeffs = [+1            -bigM   ];
    model = addNewConstraintInTFA(model, strcat('UF_', model.rxns{i}), '<', CLHS, 0);
    % UR_rxn: R_rxn - 1000 RU_rxn < 0
    CLHS.varIDs    = [R_flux_index  BU_index];
    CLHS.varCoeffs = [+1            -bigM   ];
    model = addNewConstraintInTFA(model, strcat('UR_', model.rxns{i}), '<', CLHS, 0);
end



% creating the objective
model.f = zeros(size(model.A, 2), 1);


% if objective not found return error message
if ~isempty(find(ismember(model.varNames, objective)))
    model.f(find(ismember(model.varNames, objective))) = 1;
else
    disp('Objective not found');
end

model.objtype = -1; % 1: minimize, -1: maximize
model.types = ReactionDB.vartypes;

end