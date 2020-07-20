function model = convToTFBA_BNICE(model,ReactionDB)


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

% if the model has LMPD reactions, change the rev value of the LMPDs to 1,
% to have R_ reaction in the modelIrrev.
LMPD_ind = find(contains(model.rxns,'LMPD'));
if ~isempty(LMPD_ind)
    model.rev(LMPD_ind) = 1;
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

% creating the objective
model.f = zeros(num_vars,1);

% if objective not found return error message
if ~isempty(find(ismember(model.varNames, objective)))
    model.f(find(ismember(model.varNames, objective))) = 1;
else
    disp('Objective not found');
end

model.objtype = -1; % 1: minimize, -1: maximize
model.types = ReactionDB.vartypes;

end