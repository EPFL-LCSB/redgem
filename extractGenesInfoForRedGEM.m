function RedModelWGenes = extractGenesInfoForRedGEM(RedModel, GEM, idGEMRxnsInRedGEM, idLUMPEDRxnsInRedGEM)
% This function finds the gene information of the reduced model in the GEM,
% and creates the corresponding fields with the potentially reduced
% information. The lumped reactions have no gene reaction information
% (we leave it empty)

% initialize output model with gene information
RedModelWGenes = RedModel;

if size(idGEMRxnsInRedGEM,2) > size(idGEMRxnsInRedGEM,1)
    idGEMRxnsInRedGEM = idGEMRxnsInRedGEM';
end

if size(idLUMPEDRxnsInRedGEM,2) > size(idLUMPEDRxnsInRedGEM,1)
    idLUMPEDRxnsInRedGEM = idLUMPEDRxnsInRedGEM';
end

% Find the indices of the reactions of the reduced model in the GEM:
idRedGEMRxnsInGEM = find_cell(RedModelWGenes.rxns(idGEMRxnsInRedGEM), GEM.rxns);

% The following subset of gene-reaction rules are relevant to the reduced
% model:
RedModelWGenes.grRules = repmat({''}, size(RedModelWGenes.rxns,1), 1);
RedModelWGenes.grRules(1:size(idGEMRxnsInRedGEM,1)) = GEM.grRules(idRedGEMRxnsInGEM);

% Find the non-repeated gene names that are present in the reduced model
% Examples of "grRules" field:
% - E. coli: '((b0463 and b0462 and b3035) or (b0463 and b2470 and b3035))'
% Remove "and", and "or" operators, as well as right/left parentheses
PresentGenesInRedGEM = regexprep(RedModelWGenes.grRules,'(or)|(and)|\)|\(','');
% Some models have different notation for the logical opperators... this is
% a headache...
if any(~cellfun(@isempty, regexpi(PresentGenesInRedGEM, '( OR )|(  AND )')))
    PresentGenesInRedGEM = regexprep(PresentGenesInRedGEM,' OR ',' ');
    PresentGenesInRedGEM = regexprep(PresentGenesInRedGEM,' AND ',' ');
    PresentGenesInRedGEM = regexprep(PresentGenesInRedGEM,'\(|\)','');
end
% Remove spaces from both sides of the rule strings
PresentGenesInRedGEM = regexprep(PresentGenesInRedGEM,'^( +)|( +)$','');

% Split the gene names that are separated with spaces
[A, B] = regexp(PresentGenesInRedGEM,'( +)','split','tokens');
idNonEmpty = ~cellfun(@(x) isempty(x), B);
C = A(idNonEmpty);

% Get the non-repeated set of genes (excluding the emplty one)
PresentGenesInRedGEM = setdiff(unique([C{:}]),{''})';

% Set the non-repeated genes to the reduced model:
RedModelWGenes.genes = PresentGenesInRedGEM;

% Find the indices of the genes of the reduced model in the GEM:
idRedGEMGenesInGEM = find_cell(RedModelWGenes.genes, GEM.genes);

% Get the reaction gene matrix for the reduced model
RedModelWGenes.rxnGeneMat = zeros(size(RedModelWGenes.rxns,1), size(RedModelWGenes.genes,1));
RedModelWGenes.rxnGeneMat(1:size(idRedGEMRxnsInGEM,2),:) = GEM.rxnGeneMat(idRedGEMRxnsInGEM, idRedGEMGenesInGEM);


end