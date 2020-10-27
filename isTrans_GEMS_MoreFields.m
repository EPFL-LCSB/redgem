function model = isTrans_GEMS_MoreFields(model, mets_core)

unique_comp = unique(model.metCompSymbol);
for im = 1:length(unique_comp)
    mets_core = strrep(mets_core, strcat('_', unique_comp{im}), '');
end
mets_core_all={};
for im = 1:length(unique_comp)
    mets_core_all = [mets_core_all, strcat(mets_core, '_', unique_comp{im})];
end
model.isCoreTrans = zeros(length(model.rxns),1);
model.extra_enzymatic = zeros(length(model.rxns),1);

for i=1:length(model.rxns)
    mets = model.mets(find(model.S(:, i)));
    comps = model.metCompSymbol(find(model.S(:, i)));
    mets = strrep(mets, '_c', '');
    mets = strrep(mets, '_p', '');
    mets = strrep(mets, '_e', '');
    mets = strrep(mets, '_g', '');
    mets = strrep(mets, '_l', '');
    mets = strrep(mets, '_m', '');
    mets = strrep(mets, '_n', '');
    mets = strrep(mets, '_r', '');
    mets = strrep(mets, '_s', '');
    mets = strrep(mets, '_x', '');
    mets = strrep(mets, '_j', '');
    mets = strrep(mets, '_o', '');
    mets = strrep(mets, '_q', '');
    mets = strrep(mets, '_t', '');
    mets = strrep(mets, '_v', '');
    mets = strrep(mets, '_w', '');
    unique_mets = unique(mets);
    unique_comp = unique(comps);
    if length(mets)~=length(unique_mets) || length(unique_comp)>1
        
        model.isTrans(i,1)=1;
        
        mets_to_check = model.mets(find(model.S(:, i)));
        [ab, bmets] = ismember(mets_to_check, mets_core_all);
        if length(mets_to_check)==length(find(bmets))
            model.isCoreTrans(i,1) = 1;
        end
        
        if ismember('e', unique_comp)
            model.isExchangeTrans(i,1) = 1;
        else
            model.isExchangeTrans(i,1) = 0;
        end
    else
        model.isTrans(i,1) = 0;
        model.isExchangeTrans(i,1) = 0;
    end
    if length(mets)>1 && length(unique_comp)==1
        if strcmp(unique_comp, 'e')
            model.extra_enzymatic(i,1) = 1;
        else
            model.extra_enzymatic(i,1) = 0;
        end
    end
end

drains = extract_drains(model);
[~,id] = ismember(drains, model.rxns);
model.isExchange = zeros(length(model.rxns),1);
model.isExchange(id) = 1;


end