function model = removefutiles(model, rxns, directions, checkgrowth, CplexParameters)

% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

mets = unique(rxns(:,2));

for i = 1:size(mets)
    [~, ba] = ismember(rxns(:, 2), mets(i));
    to_check = rxns(find(ba), 1);
    dir_to_check = directions(find(ba));
    [~, ba] = ismember(to_check, model.rxns);
    if length(find(ba>0))>1
        to_check = to_check(find(ba));
        dir_to_check = dir_to_check(find(ba));
        pairs = combnk(1:length(to_check),2);
        for j = 1:size(pairs, 1)
            model_orig = model;
            num_constr = size(model.A, 1);
            model.rhs(num_constr+1,1) = 1;
            
            model.constraintType{num_constr+1,1} = '<';
            if dir_to_check(pairs(j,1))==dir_to_check(pairs(j,2))
                rxn1 = strcat('FU_', to_check{pairs(j,1)});
                rxn2 = strcat('BU_', to_check{pairs(j,2)});
                rxn3 = strcat('BU_', to_check{pairs(j,1)});
                rxn4 = strcat('FU_', to_check{pairs(j,2)});
            else
                rxn1 = strcat('FU_', to_check{pairs(j,1)});
                rxn2 = strcat('FU_', to_check{pairs(j,2)});
                rxn3 = strcat('BU_', to_check{pairs(j,1)});
                rxn4 = strcat('BU_', to_check{pairs(j,2)});
            end
            
            model.constraintNames{num_constr+1,1} = strcat('Futile_',rxn1,'_',rxn2,'_',num2str(j));
            
            [~, ba] = ismember(rxn1, model.varNames);
            [~, bb] = ismember(rxn2, model.varNames);
            [~, bc] = ismember(rxn3, model.varNames);
            [~, bd] = ismember(rxn4, model.varNames);
            model.A(num_constr+1, ba) = 1;
            model.A(num_constr+1, bb) = 1;
            
            num_constr = size(model.A,1);
            model.rhs(num_constr+1,1) = 1;
            
            model.constraintNames{num_constr+1,1} = strcat('Futile_',rxn3,'_',rxn4,'_',num2str(j));
            
            model.constraintType{num_constr+1,1} = '<';
            
            model.A(num_constr+1, bc) = 1;
            model.A(num_constr+1, bd) = 1;
            
            if checkgrowth
                sol = solveTFAmodelCplex(model, 300, [], mipTolInt, emphPar, feasTol, scalPar, []);
                if isempty(sol.x)
                    model = model_orig;
                elseif sol.val<0.01
                    model = model_orig;
                end
            end
        end
    end
end

end