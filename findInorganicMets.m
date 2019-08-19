function [InorgMetSEEDIDs,InorgMetDetails] = findInorganicMets(model,NoFormulaStr,ExtraFormulas,ToExclude)
% Simple function to find inorganic metabolites from COBRA compatible
% models that have a metFormula entry (model.metFormulas).
%
% INPUTS
% - model: A COBRA compatible model structire
% - NoFormulaStr: A cell array containing characteristic strings used in
%   model.metFormulas to indicate that no chemical formula is available for
%   a particular metabolite. This is model specific (default is: 'NA',
%   'None','NULL','').
% - ExtraFormulas: To the list of inorganic compounds, we might also want
%   to add some compounds that are not captured by the simplistic approach
%   of this function. If not set, by default we also add cyanate (CNO),
%   bicarbonate (CHO3), carbon dioxide (CO2), thiocyanate (CNS), and
%   hydrogen cyanide (CHN). 
% - FormulasToExclude: some metabolites whose formula is misleading, and might end
%   up in the inorganic list, such as:
%   S2Fe2: [2Fe-2S] iron-sulfur cluster, which is not inorganic
% 
% OUTPUTS
% - InorgMetSEEDIDs: SEEDIDs of all the inorganic metabolites with no
%   repetitions.
% - InorgMets: cell array of size n-by-5, where n is the number of
%   inorganic compounds found in the model. The 5 columns correspond to:
%   1) Index-number of metabolites in the model metabolite list
%   2).mets
%   3).metNames
%   4).metFormulas
%   5).metSEEDID
% 
% LCSB, EPFL, December 2014


if nargin<4
    FormulasToExclude = {'S2Fe2','S4Fe4'};
end

if nargin<3
    ExtraFormulas = {'CNO','CHO3','CO2','CNS','CHN'};
end

if nargin<2
    NoFormulaStr = {'NA','None','NULL',''};
end
    
% Make sure that the provided model contains a field with the chemical formulas
model_fields = fieldnames(model);

if ~ismember(model_fields,'metFormulas')
    error('Model File Does not contain Metabolite Formulas or the variable is not named model.metFormulas')
end

% Find the indices of the metabolites whose chemical formulas does not
% contain 'C' followed by another capital letter or any number. For example
% indices of elements such as 'Ca', 'Cd', etc will be found:
InorganicMetID = find(cellfun(@isempty,regexp(model.metFormulas,'(?-i)C[A-Z_0-9]')));
% Attention: are there any falsely identified organic metabolites?
FalseInorgMetIds = intersect(InorganicMetID,find_cell(FormulasToExclude,model.metFormulas));
if ~isempty(FalseInorgMetIds)
    InorganicMetID = setdiff(InorganicMetID,FalseInorgMetIds);
end
% From these indices (InorganicMetID) exclude the indices that correspond
% to "NoFormulaStr".
NoFormulaStrID = find(ismember(model.metFormulas,NoFormulaStr));
InorganicMetID = setdiff(InorganicMetID,NoFormulaStrID);
% To these indices (InorganicMetID) add the indices corresponding to the 
% "ExtraFormulas"
ExtraFormulasID = find(ismember(model.metFormulas,ExtraFormulas));
InorganicMetID = union(InorganicMetID,ExtraFormulasID);

if size(InorganicMetID,2) > size(InorganicMetID,1)
    InorganicMetID = InorganicMetID';
end

% Append all useful properties in one cell
InorgMetDetails = [num2cell(InorganicMetID) ...
                   model.mets(InorganicMetID) ...
                   model.metNames(InorganicMetID) ...
                   model.metFormulas(InorganicMetID) ...
                   model.metSEEDID(InorganicMetID)];

% Find the indices of the metabolites that have no SEEDID      
NotCpdID = find(cellfun(@isempty,regexp(model.metSEEDID(InorganicMetID),'cpd')));

% If there exist metabolites without SEEDID, take action accordingly:
if ~isempty(NotCpdID)
    warning('Some of the inorganic metabolites have no cpd-number:')
    InorgMetDetails(NotCpdID,:)
end

% Find all the SEEDIDs without repetition, excluding the strings that
% correspond to not existent cpd-number (e.g.: 'NA')
if ~isempty(InorgMetDetails)
    InorgMetSEEDIDs = setdiff(unique(InorgMetDetails(:,5)),unique(InorgMetDetails(NotCpdID,5)));
else
    % Sometimes, especially for toy models it could be the case that there
    % exist no inorganic metabolites. This is not per se an error!
    warning('There appear to be NO inorganic metabolites in your model!')
    InorgMetSEEDIDs = [];
    InorgMetDetails = [];
end

      

end