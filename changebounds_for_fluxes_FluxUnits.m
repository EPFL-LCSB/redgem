function model = changebounds_for_fluxes_FluxUnits(model)
% This function finds all the forward and reverse ('F','R') variables, and
% checks wether any of them has values higher than 60 (6*10, where 6 is the
% number of carbons of glucose, and 10 is the glucose uptake rate. In
% principle, no flux should have higher value.
% 
% INPUTS
% - model: TFBA model
%
% OUTPUTS
% - model: TFBA model with constrained forward and reverse flux variables
% ('F','R')

ind_for=getAllVar(model, {'F'});
ind_back=getAllVar(model, {'R'});

bound_for=model.var_ub(ind_for);
bound_back=model.var_ub(ind_back);

for i=1:length(ind_for)
    
    if bound_for(i)>60
        model.var_ub(ind_for(i))=60;
    end
    
    if bound_back(i)>60
        model.var_ub(ind_back(i))=60;
    end
    
end
    
end