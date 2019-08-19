function redundant_table = ReacRedundancy(M)

% Does the model have any redundant reactions?
% Print-out to a variable all the reactions of the M
GEMRxns = printRxnFormula(M,M.rxns,0);
% Set to S the stoichiometric matrix of the M
S = M.S;
% Find all the unique columns-entries of S, (OR row entries of S'), and set
% them to x.
[x,~] = unique(S','rows');

% Find all the indices of the unique columns of the S'
[~,x2] = ismember(S',x,'rows');
% find the occurences of these unique columns
[a,~] = hist(x2,unique(x2));
% which of these unique entries have more than 1 occurencies? These are the
% reactions that are found at least two times in our model!
[c,~] = ismember(S',x(a>1,:),'rows');

redundant_table = [];
if ~isempty(find(c))
    warning('Attention: The following reactions appear to be exactly the same. Please make sure that this is intentional! E.g. regulated by different gene')
    if isfield(M,'rxnNames')
        redundant_table = [M.rxnNames(find(c)) M.rxns(find(c)) GEMRxns(find(c)) M.grRules(find(c)) num2cell([M.lb(find(c)) M.ub(find(c))])]
    else
        redundant_table = [M.rxns(find(c)) GEMRxns(find(c)) M.grRules(find(c)) num2cell([M.lb(find(c)) M.ub(find(c))])]
    end
end