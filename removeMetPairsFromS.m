function Smatrix = removeMetPairsFromS(model,met_pairs,inorganic)
% function to remove cofactors pairs from reactions as well as small
% metabolites for graph search
%
% INPUTS:
% model refers to the FBA model structure
% met_pairs should be delimited by a colon met1:met2 and in terms of
% compound IDs
% inorganic are small metabolites that should be removed

Smatrix = model.S;

if isa(inorganic,'double') % if the input is directly the met indeces
    inorganic_indices = inorganic;
else 
    inorganic_indices = find(ismember(model.metSEEDID, inorganic));
end

Smatrix(inorganic_indices,:) = 0;

[~, num_rxns] = size(model.S);
all_mets = {};

for i = 1:length(met_pairs)
    mets = splitString(met_pairs{i},':');
    mets_to_add{1,1} = mets{1};
    mets_to_add{1,2} = mets{2};
    all_mets = [all_mets; mets_to_add];
end

for i = 1:num_rxns
    mets_match = {};
    mets = model.mets(find(model.S(:,i)));
    for j = 1:size(all_mets,1)
        [~, ba] = ismember(all_mets(j,:), mets);
        if length(find(ba)) == 2
            mets_match = [mets_match all_mets(j,:)];
        end
    end
    if exist('mets_match', 'var')
        [~, ba] = ismember(mets_match, model.mets);
        if length(mets_match) > 1
            Smatrix(ba,i) = 0;
        end
    end
    clear mets_match
end

end