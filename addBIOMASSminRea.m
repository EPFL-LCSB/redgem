function [model, ActRxnFormulas, sol, net_fluxes] = addBIOMASSminRea(model, Ind_OtherReacs, DB_AlbertyUpdate, minPercent, ImposeThermodynamics, CplexParameters)

%what minimum percentage of original maximal growth do we want to calculate
%extract the "subnetwork" with? default = 0.99
if ~exist('minPercent','var') || isempty(minPercent)
    minPercent = 0.9;
end

% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

for i=1:length(model.rxns)
    mets_rxns=model.mets(find(model.S(:,i)));
    if ismember('co2_e', mets_rxns)
        if length(mets_rxns)==1
            model.lb(i)=0;
            break;
        end
    end
end

% Converts, again, the model to TFBA, BUT without actually adding ANY
% thermodynamic information to any reaction
if strcmp(ImposeThermodynamics, 'yes')
    verboseFlag = false;
    model = convToTFA(model, DB_AlbertyUpdate, model.rxns, 'NO');
elseif strcmp(ImposeThermodynamics, 'no')
    fprintf('In the current case we do not impose thermodynamic constraints\n')
    model = convFBA2MILP(model, DB_AlbertyUpdate);
else
    error('Wrong option for ImposeThermodynamics')
end

sol = solveTFAmodelCplex(model,[],[],mipTolInt,emphPar,feasTol,scalPar,[]);
% Set lower biomass production
model.var_lb(find(model.f))=minPercent*sol.val;

% Reset all the objectives to zero
model.f=zeros(length(model.varNames),1);
% Get indices of user-variables (F/R)
forward  = getAllVar(model, {'F'});
forward  = forward(Ind_OtherReacs);
backward = getAllVar(model, {'R'});
backward = backward(Ind_OtherReacs);

[~, num_vars] = size(model.A);

for i = 1:length(Ind_OtherReacs)
    model.varNames(length(model.varNames)+1) = strcat('BFUSE_', model.rxns(Ind_OtherReacs(i)));
    model.var_ub(length(model.varNames))=1;
    model.var_lb(length(model.varNames))=0;
    model.vartypes(length(model.varNames)) = {'B'};
    model.f(length(model.varNames))=1;
end

if length(forward) == length(backward)
    for i = 1:length(forward)
        [num_constr,~] = size(model.A);
        model.rhs(num_constr+1,1) = 100;
        model.constraintNames{num_constr+1,1} = strcat('BFMER_',num2str(i));
        model.constraintType{num_constr+1,1} = '<';
        model.A(num_constr+1,forward(i)) = 1;
        model.A(num_constr+1,backward(i)) = 1;
        model.A(num_constr+1,num_vars+i) = 100;
    end
else
    error('ATTENTION: forward and backward reactions are not equal!!')
end

sol = solveTFAmodelCplex(model, [], [], mipTolInt, emphPar, feasTol, scalPar, []);

% Find all the BFUSE indices
Ind_BFUSE = getAllVar(model, {'BFUSE'});
bfuse = sol.x(Ind_BFUSE);
% Find which of the BFUSE variables are zero, i.e. which of the reactions
rxnsBFUSE = find(bfuse<1e-9);
% If you model crashes around here, it is probable that either it is
% infeasible or you have only the null solution
ActBVarNames = model.varNames(Ind_BFUSE(rxnsBFUSE));
ActRxnNames = strrep(ActBVarNames, 'BFUSE_', '');
[~, id_ActRxnNames] = ismember(ActRxnNames, model.rxns);
ActRxnFormulas = printRxnFormula(model, model.rxns(id_ActRxnNames));
ActRxnFormulas = [model.rxns(id_ActRxnNames) ActRxnFormulas];
F_rxns = strcat('F_', ActRxnFormulas(:,1));
R_rxns = strcat('R_', ActRxnFormulas(:,1));
[~, F_rxn_ids] = ismember(F_rxns, model.varNames);
[~, R_rxn_ids] = ismember(R_rxns, model.varNames);
net_fluxes = sol.x(F_rxn_ids) - sol.x(R_rxn_ids);

end

