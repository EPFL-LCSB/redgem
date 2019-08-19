function met_pairs_to_remove=preparemetpairs(model, metpairs)
m=1;
compartments=unique(model.metCompSymbol);
for i=1:length(metpairs(:,1))
    for j=1:length(compartments)
        met1=strcat(metpairs(i,1), '_', compartments(j));
        met2=strcat(metpairs(i,2), '_', compartments(j));
        [~, bmets]=ismember([met1;met2], model.mets);
        if ~ismember(0, bmets)
        met_pairs_to_remove(m,1)=strcat(metpairs(i,1), '_', compartments(j));
        met_pairs_to_remove(m,2)=strcat(metpairs(i,2), '_', compartments(j));
        m=m+1;
        end
    end
end
met_pairs_to_remove=strcat(met_pairs_to_remove(:,1), ':', met_pairs_to_remove(:,2));