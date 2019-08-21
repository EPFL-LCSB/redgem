function [DPs,objectives,model] = findDP_YIELD_TimeLimit(model,maxNumEPs,sol,use_vectors_indices,TimeLimitSecs,SizeLUMPNetwork,CplexParameters)
% This function introduces integer cuts to find all DPs in a model

% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

num_sols = 0;

count=1;
while ((num_sols < maxNumEPs) && ~(isempty(sol)))
    
    [num_cons,num_vars] = size(model.A);
    
    if ~(isempty(sol.x))
        num_sols = num_sols + 1;
        objectives(num_sols,1) = sol.val;
        DPs(:,num_sols) = sol.x;
        
        % we find all the use vectors and formulate them into a new integer cut
        % constraint
        use_vectors_sol = ones(num_vars,1);
        use_vectors_sol(use_vectors_indices) = sol.x(use_vectors_indices);
        active_use_vectors = find(use_vectors_sol<1e-9);
        if count==1
            smallest_subnetwork=length(active_use_vectors);
        end
        new_constraint = zeros(num_vars,1);
        new_constraint(active_use_vectors) = 1;
        % 1*activeUseVar1 + 1*activeUseVar2 + 1*activeUseVar3 + ... > 0.5
        model.A(num_cons+1,:) = new_constraint;
        model.rhs(num_cons+1) = 0.5;
        model.constraintNames{num_cons+1} = ['CUT_' num2str(num_sols)];
        model.constraintType{num_cons+1} = '>';
        sol = solveTFAmodelCplex(model, TimeLimitSecs, [], mipTolInt, emphPar, feasTol, scalPar, []);

        if isempty(sol.x)
            break;
        end
        length_sol = sol.x(use_vectors_indices);
        length_sol = length(find(length_sol<1e-9));
        if (length_sol - smallest_subnetwork) > SizeLUMPNetwork
            break;
        end

        fprintf('Number of DPs:\t%d\n',num_sols);
        
    end
    count=count+1;
end


