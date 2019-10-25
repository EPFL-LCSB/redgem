function GEMmodel = preventBBBuptake(GEMmodel)

% Find biomass index. ATTENTION!!! MAKE SURE THAT GEMmodel.c IS SET TO
% OPTIMIZE THE BIOMASS!! WE NEED TO ADD A CHECK BEFOREHAND!!!
id_bmass = find(GEMmodel.c);
% Find biomass building blocks
bbbs = GEMmodel.mets(GEMmodel.S(:,id_bmass)<0);

bbbs_e = cellfun(@(x) [x(1:end-2),'_e'],bbbs,'UniformOutput' ,false);

bio_drains = [];
for i = 1:length(GEMmodel.rxns)
    no_of_mets = length(find(GEMmodel.S(:,i)));
    if no_of_mets == 1
        if ismember(GEMmodel.mets(find(GEMmodel.S(:,i))),bbbs_e)
            bio_drains=[bio_drains;GEMmodel.rxns(i)];
        end
    end
end

% Exclude water. We use the seedid to identify water to exchange rxns be
% sure that its is organism and GEM independent. 'hydroxide ion' is also
% assigned as cpd00001!

id_h2o = find_cell('h2o_e',GEMmodel.mets);
id_h2o_e_in_h2o = find_cell('e',GEMmodel.metCompSymbol(id_h2o));
id_h2o_e = id_h2o(id_h2o_e_in_h2o);
id_h2o_e_rxns = find(GEMmodel.S(id_h2o_e,:)~=0);
NumMets_in_h2o_e_rxns = sum(abs(GEMmodel.S(:,id_h2o_e_rxns)));

id_h2o_e_exchange_rxn = id_h2o_e_rxns(find(NumMets_in_h2o_e_rxns == 1));

if size(id_h2o_e_exchange_rxn,1)>1
    error('There exist more than 1 h2o exchange reactions in the model!!??')
end

h2o_e_exchange_rxn = GEMmodel.rxns(id_h2o_e_exchange_rxn);
bio_drains = setdiff(bio_drains, h2o_e_exchange_rxn);

[~, b] = ismember(bio_drains, GEMmodel.rxns);
% Set the lower bound of these to zero (they can only be secreted).
GEMmodel.lb(b) = 0;
end