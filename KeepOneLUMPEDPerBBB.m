function [model_out, MinLumpRxnNames, MinLumpRxnFormulas] = KeepOneLUMPEDPerBBB(model, indLumped, growthRate, MinFluxOfLumped)

% If there is no minimal flux imposed for every lumped, set it to a default
% value of 10^-7
if ~exist('MinFluxOfLumped','var') || isempty(MinFluxOfLumped)
    MinFluxOfLumped = 1e-07;
end

model_in = model;

% the lower bound for biomass. should be the maximum from the model with all Smin Lumps.
model.var_lb(find(model.f)) = growthRate;

% ADD ALL THE CONSTRAINTS TO THE MODEL:
%
flagLUMPED = 'AtMaxOneLUMPEDPerBBB';
% flagLUMPED = 'MinNumOfLUMPED';

if strcmp(flagLUMPED,'MinNumOfLUMPED')
    modelWithConst = AddConstraintsForLUMPEDRxns(model,indLumped,MinFluxOfLumped);
else
    % Extract two tokens from the Lumped reaction names:
    % - number of lumped
    % - metabolite name
    LMPDTokens = regexp(model.varNames(indLumped),'LMPD_(\d+)_([\w.\-]+)','tokens');
    LMPDTokens = [LMPDTokens{:}]';
    LMPDTokens = [LMPDTokens{:}];
    LMPDTokens = reshape(LMPDTokens,2,size(LMPDTokens,2)/2)';
    UniqueMetsOrdered = unique(LMPDTokens(:,2),'stable');
    modelWithConst = model;
    for i=1:size(UniqueMetsOrdered,1)
        indLumpedPerBBB = find_cell(UniqueMetsOrdered(i),LMPDTokens(:,2));
        modelWithConst = AddConstraintsForLUMPEDRxns(model,indLumped,MinFluxOfLumped);
    end
end

dmodel = modelWithConst;
sol = solveTFAmodelCplex(dmodel);
ind = getAllVar(dmodel, {'BFUSE'});
if isempty(sol.x)
    MinLumpRxnFormulas = {'NA'};
    MinLumpRxnNames = {'NA'};
else
    bfuse = sol.x(ind);
    rxnsBFUSE = find(bfuse==0);
    if isempty(rxnsBFUSE)
        MinLumpRxnFormulas = {'NA'};
        MinLumpRxnNames = {'NA'};
    else
        MinLumpRxnNames = dmodel.varNames(ind(rxnsBFUSE));
        MinLumpRxnNames = strrep(MinLumpRxnNames, 'BFUSE_', '');
        [~, ba] = ismember(MinLumpRxnNames, dmodel.rxns);
        MinLumpRxnFormulas = printRxnFormula(dmodel, dmodel.rxns(ba));
    end
end

% We would like to remove all the other lumped reactions apart from this
% minimal set. So we would like to remove :
LumpRxnNamesToRemove = setdiff(model_in.rxns(indLumped),MinLumpRxnNames);

model_out = removeRxns(model_in,LumpRxnNamesToRemove,false,false);
% Find those unused metabolites and keep the indices
UnusedMets = model_out.mets(any(sum(abs(model_out.S),2) == 0,2));
% UnusedMetIDs = find_cell(UnusedMets,model_out.mets);

if (~isempty(UnusedMets))
    model_out = removeMetabolites(model_out, UnusedMets, false);
end

end

function modelWithConst = AddConstraintsForLUMPEDRxns(model,indLumped,MinFluxOfLumped)
% INPUT
% - model: model without constraints
%
% OUTPUT
% - modelWithConst: model with constraints

% Get the forward and backward variable indices of the LUMPED reactions
forward = getAllVar(model, {'F'});
backward = getAllVar(model, {'R'});
forward  = forward(indLumped);
backward = backward(indLumped);

if length(forward)~=length(backward)
    error('Number of forward reactions is not the same as number of backward reactions!!')
else
    NLmpd = length(forward);
end

% Set the objective vector of the model to zeros:
model.f = zeros(length(model.f),1);
% Create for every reaction a binary variable and set it as an objective
for i=1:length(indLumped)
    model.varNames(length(model.varNames)+1) = strcat('BFUSE_', model.rxns(indLumped(i)));
    model.var_ub(length(model.varNames)) = 1;
    model.var_lb(length(model.varNames)) = 0;
    model.vartypes(length(model.varNames)) = {'B'};
    model.f(length(model.varNames)) = 1;
end

ind_bfuse = getAllVar(model, {'BFUSE'});

% Constraint: F + R + 60*BFUSE < 60
% BFUSE = 1 ---> Reaction is closed (F+R <= 0)
% BFUSE = 0 ---> Reaction is active (F+R < 60)
for i = 1:NLmpd
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1) =  60;
    model.constraintNames{num_constr+1,1} = strcat('BFMER_',num2str(i));
    model.constraintType{num_constr+1,1} = '<';
    model.A(num_constr+1,forward(i)) = 1;
    model.A(num_constr+1,backward(i)) = 1;
    model.A(num_constr+1,ind_bfuse(i)) = 60;
end

% Constraint: F + R + BFUSE > 1e-07
% BFUSE = 1 ---> Reaction is closed (F + R = 0)
% BFUSE = 0 ---> Reaction is active with residual flux (F + R = 0)
[num_constr,~] = size(model.A);
for i = 1:NLmpd
    model.rhs(num_constr+i,1) = MinFluxOfLumped;
    model.constraintNames{num_constr+i,1} = strcat('BFMER2_',num2str(i));
    model.constraintType{num_constr+i,1} = '>';
    model.A(num_constr+i,forward(i)) = 1;
    model.A(num_constr+i,backward(i)) = 1;
    model.A(num_constr+i,ind_bfuse(i)) = 1;
end

% Create two binary variables:
% one for the forward (FUSE):
for i = 1:NLmpd
    model.varNames(length(model.varNames)+1) = strcat('FUSE_', model.rxns(indLumped(i)));
    model.var_ub(length(model.varNames)) = 1;
    model.var_lb(length(model.varNames)) = 0;
    model.vartypes(length(model.varNames)) = {'B'};
    model.A(:, length(model.varNames)) = 0;
    model.f(length(model.varNames)) = 0;
end
% and one for the backward (BUSE) flux
for i = 1:NLmpd
    model.varNames(length(model.varNames)+1) = strcat('BUSE_', model.rxns(indLumped(i)));
    model.var_ub(length(model.varNames)) = 1;
    model.var_lb(length(model.varNames)) = 0;
    model.vartypes(length(model.varNames)) = {'B'};
    model.A(:, length(model.varNames)) = 0;
    model.f(length(model.varNames)) = 0;
end

% Extract the indices of these forward and backward binary variables:
ind_fu = getAllVar(model, {'FUSE'});
ind_bu = getAllVar(model, {'BUSE'});


% Constraint: F  < 60*FUSE
% FUSE = 0 ---> Forward reaction is closed
% FUSE = 1 ---> Forward reaction is active
[num_constr,~] = size(model.A);
for i = 1:NLmpd
    model.rhs(num_constr+i,1) =  0;
    model.constraintNames{num_constr+i,1} = strcat('BFMER3_',num2str(i));
    model.constraintType{num_constr+i,1} = '<';
    model.A(num_constr+i,forward(i)) = 1;
    model.A(num_constr+i, ind_fu(i)) = -60;
end

% Constraint: R  < 60*BUSE
% BUSE = 0 ---> Backward reaction is closed
% BUSE = 1 ---> Backward reaction is active
[num_constr,~] = size(model.A);
for i = 1:NLmpd
    model.rhs(num_constr+i,1) =  0;
    model.constraintNames{num_constr+i,1} = strcat('BFMER4_',num2str(i));
    model.constraintType{num_constr+i,1} = '<';
    model.A(num_constr+i,backward(i)) = 1;
    model.A(num_constr+i, ind_bu(i)) = -60;
end

% Constraint:  BUSE + FUSE + BFUSE  < 1
% At most one binary variable can be equal to one. Explanation: because of
% tFA, if FUSE=1, then BUSE=0 and vice versa. If any of the two is one,
% i.e. FUSE=1 or BUSE=1, then the reaction is active, so BFUSE=0. I.e.
% these two types of integer variables are exclusive
[num_constr,~] = size(model.A);
for i = 1:NLmpd
    model.rhs(num_constr+i,1) =  1;
    model.constraintNames{num_constr+i,1} = strcat('BFMER5_',num2str(i));
    model.constraintType{num_constr+i,1} = '<';
    model.A(num_constr+i,ind_bfuse(i)) = 1;
    model.A(num_constr+i, ind_fu(i)) = 1;
    model.A(num_constr+i, ind_bu(i)) = 1;
end

modelWithConst = model;

end


