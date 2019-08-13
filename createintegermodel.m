function [model] = createintegermodel(model)
 
 
[num_mets_org,num_rxns] = size(model.S);
 
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
 
modelIrrev = convertToIrreversibleAgador(model);
A = modelIrrev.S;
[num_mets,num_vars] = size(modelIrrev.S);
objective = modelIrrev.rxns(find(modelIrrev.c));
 
% check that num of metabolites remain the same
 
if num_mets ~= num_mets_org
  error('number of metabolites do not match!') 
end
 
% create the constraint type vector first for mass balances
% and the rhs vector first for mass balances
 
for i=1:num_mets
   constraintType{i,1} = '=';
   constraintNames{i,1} = strcat('M_',model.mets{i});
   rhs(i,1) = 0;
end
 
% create the variable type vector and setting their bounds
% and lower and upper bounds for the flux variables
% upperbounds should not be negative
vartypes = {};
var_lb = zeros(num_vars,1);
var_ub = zeros(num_vars,1);
 
for i=1:num_vars
 
   if modelIrrev.ub(i) < 0
       modelIrrev.ub(i) = 0;
   end
 
   if modelIrrev.lb(i) < 0
       modelIrrev.lb(i) = 0;
   end
 
   var_ub(i) = modelIrrev.ub(i);
   var_lb(i) = modelIrrev.lb(i);
   varNames = modelIrrev.rxns;
   vartypes{i,1} = 'C';
end
 
 
for i=1:num_rxns
 
   [num_constr,num_vars] = size(A);
 
   varNames{length(varNames)+1,1} = strcat('FU_',model.rxns{i});
   var_ub(length(varNames),1) = 1;
   var_lb(length(varNames),1) = 0;
   vartypes{length(varNames),1} = 'B';
 
   [num_constr,num_vars] = size(A);
 
   varNames{length(varNames)+1,1} = strcat('BU_',model.rxns{i});
   var_ub(length(varNames),1) = 1;
   var_lb(length(varNames),1) = 0;
   vartypes{length(varNames),1} = 'B';
 
   % create the prevent simultaneous use constraints
   % SU_rxn: FU_rxn + BU_rxn <= 1
   [num_constr,num_vars] = size(A);
   rhs(num_constr+1,1) =  1;
   constraintNames{num_constr+1,1} = strcat('SU_',model.rxns{i});   
   constraintType{num_constr+1,1} = '<';
 
   A(num_constr+1,length(varNames)-1) = 1;
   A(num_constr+1,length(varNames)) = 1; 
 
   % create constraints that control fluxes with their use variables
   % UF_rxn: F_rxn - 60 FU_rxn < 0
 
   [num_constr,num_vars] = size(A);
   rhs(num_constr+1,1) =  0;
   constraintNames{num_constr+1,1} = strcat('UF_',model.rxns{i});   
   constraintType{num_constr+1,1} = '<';
 
   F_flux_var = find(ismember(varNames,strcat('F_',model.rxns{i})));
   FU_var_index = find(ismember(varNames,strcat('FU_',model.rxns{i})));
   A(num_constr+1,F_flux_var) = 1;
   A(num_constr+1,FU_var_index) = -60;
 
   % UR_rxn: R_rxn - 60 RU_rxn < 0
   [num_constr,num_vars] = size(A);
   rhs(num_constr+1,1) =  0;
   constraintNames{num_constr+1,1} = strcat('UR_',model.rxns{i});   
   constraintType{num_constr+1,1} = '<';
 
   R_flux_var = find(ismember(varNames,strcat('R_',model.rxns{i})));
   BU_var_index = find(ismember(varNames,strcat('BU_',model.rxns{i})));
   A(num_constr+1,R_flux_var) = 1;
   A(num_constr+1,BU_var_index) = -60;
end
 
f=zeros(length(varNames),1);
[~, ind]=ismember(objective, varNames);
f(ind)=1;
% collecting the new thermodynamics model
model.A = A;
model.varNames = varNames;
model.vartypes = vartypes;
model.var_lb = var_lb;
model.var_ub = var_ub;
model.objtype = -1; % 1 - minimize, -1 - maximize
model.f = f; % objective vector for TFBA problem
model.constraintNames = constraintNames;
model.constraintType = constraintType;
model.rhs = rhs;
