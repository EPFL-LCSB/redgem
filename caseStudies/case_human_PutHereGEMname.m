function [OriginalGEM, GEMmodel, core_ss, Biomass_rxns, met_pairs_to_remove, InorgMetSEEDIDs, BBBsToExclude, ExtraCellSubsystem, OxPhosSubsystem] = ...
    case_human(GEMname, ZeroZeroGEMbounds, ListForInorganicMets, ListForCofactorPairs, SelectedSubsystems, AddExtracellularSubsystem, DB_AlbertyUpdate)

% Model Specific Settings:
fprintf('Loading the GEM for human...\n')
GEM_filename = GEMname;
%GetModelFromGITresources(GEM_filename)
GEMmodel = load(['./GEMs/',GEM_filename]);
GEMmodel = GEMmodel.model;
GEMmodel.rxns(find(ismember(GEMmodel.rxns,'DM_atp_c'))) = {'ATPM'};
GEMmodel.varNames(find(ismember(GEMmodel.varNames,'F_DM_atp_c'))) = {'F_ATPM'};
GEMmodel.varNames(find(ismember(GEMmodel.varNames,'R_DM_atp_c'))) = {'R_ATPM'};
OriginalGEM = GEMmodel;

OxPhosSubsystem = 'Oxidative phosphorylation';

fprintf('For this model:\n')
fprintf('- The core subsystems  are:\n')
if iscell(SelectedSubsystems)
    core_ss = SelectedSubsystems
elseif ischar(SelectedSubsystems) && strcmp(SelectedSubsystems, 'default')
    core_ss = {'Citric acid cycle'
    'Pentose phosphate pathway'
    'Glycolysis/gluconeogenesis'
    'Oxidative phosphorylation'}
else
    error('Error defining the subsystems!')
end  

fprintf('- The biomass reaction is:\n')
Biomass_rxns = {'biomass'}
ObjRxn = {'biomass'};

%%    Model Preliminary Checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does the model have any redundant reactions?
% Print-out to a variable all the reactions of the GEMmodel
GEMRxns = printRxnFormula(GEMmodel,GEMmodel.rxns,0);
% Set to S the stoichiometric matrix of the GEMmodel
S = GEMmodel.S;
% Find all the unique columns-entries of S, (OR row entries of S'), and set
% them to x.
[x,~] = unique(S','rows');
% S(:,y) % this is x'
% Are there any redundant/same reactions?
% size(S,2)==size(x,1);
% Find all the indices of the unique columns of the S'
[~,x2] = ismember(S',x,'rows');
% find the occurences of these unique columns
[a,~] = hist(x2,unique(x2));
% which of these unique entries have more than 1 occurencies? These are the
% reactions that are found at least two times in our model!
[c,~] = ismember(S',x(a>1,:),'rows');
% S(:,c);
if ~isempty(find(c))
    warning('Attention: The following reactions appear to be exactly the same. Please make sure that this is intentional! E.g. regulated by different gene')
    [GEMmodel.rxnNames(find(c)) GEMmodel.rxns(find(c)) GEMRxns(find(c)) GEMmodel.grRules(find(c)) num2cell([GEMmodel.lb(find(c)) GEMmodel.ub(find(c))])]
end

% If the flux bounds are too high, set them to +/- 100
% ATTENTION!!! We should make this autmmatic based on the carbon
% molecules of the carbon source!!!
GEMmodel.lb(find(GEMmodel.lb <-100)) = -100;
GEMmodel.ub(find(GEMmodel.ub >+100)) = +100;

% Some reactions in the GEM model can have both lower and upper
% bounds set to zero.
idZeroZeroGEMbounds = find(GEMmodel.lb==0 & GEMmodel.ub==0);
% Would you like to proceed with these bounds,
% or would you like to open these bounds ([-100 +100])?
if ~isempty(idZeroZeroGEMbounds)
    fprintf('The following reactions appear to have zero-zero GEM bounds:\n')
    GEMmodel.rxns(idZeroZeroGEMbounds)
    if strcmp(ZeroZeroGEMbounds,'Original')
        % Do nothing, just proceed!
        fprintf('We keep the zero-zero GEM bounds!\n')
    elseif strcmp(ZeroZeroGEMbounds,'DefineCustom')
        % Create a field with the changes introduced here!!
        error('Not implemented yet')
    elseif strcmp(ZeroZeroGEMbounds,'OpenTo100')
        fprintf('We choose to open the zero-zero GEM bounds to [-100 +100]\n')
        GEMmodel.lb(idZeroZeroGEMbounds) = -100;
        GEMmodel.ub(idZeroZeroGEMbounds) = +100;
    else
        error('Wrong option!')
    end
end

BBBsToExclude = []; % there are not BBBs that we want to exclude for lumpGEM
ExtraCellSubsystem = []; % no extracellular subsystem for this human case

%%    Inorganic Metabolites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ListForInorganicMets,'curated')
    % load curated list from mat-file
    InorgMetSEEDIDs_curated = load('inorganic150109.mat');
    InorgMetSEEDIDs_curated = InorgMetSEEDIDs_curated.inorganic;
    InorgMetSEEDIDs = InorgMetSEEDIDs_curated;
elseif strcmp(ListForInorganicMets,'automatic')
    % Call function to calculate the inorganic metabolite list
    InorgMetSEEDIDs = findInorganicMets(GEMmodel);
else
    error('Wrong option!')
end

%%    Metabolite Pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ListForCofactorPairs,'curated')
    % load curated list from mat-file
    metpairsHUMAN = load('metpairsHUMAN.mat');
    metpairsHUMAN = metpairsHUMAN.metpairsHUMAN;
    met_pairs_to_remove = preparemetpairs(GEMmodel, metpairsHUMAN);
elseif strcmp(ListForCofactorPairs,'automatic')
    % Call function to calculate the lost of cofactor pairs
    % List of SEEDIDs of small metabolites:
    List_small_mets = {'h','h2o','co2','o2','pi','ppi','nh4','h2o2'};
    % Based on this list, find the corresponding SEEDIDs
    small_metsSEEDIDs = unique(GEMmodel.metSEEDID(~cellfun(@isempty,regexp(GEMmodel.mets, cell2mat(cellfun(@(x) ['^',x,'_|'],List_small_mets,'UniformOutput',false))))));
    all_met_pairs = countMetPairs_alex(GEMmodel,small_metsSEEDIDs);
    % Set the Minimum number of occurences of each metabolite pair,
    % acting as a threshold for considering a pair to be a valid cofactor pair.
    MinNumOccurence = 5;
    met_pairs_to_remove = findCofactorPairs_alex(all_met_pairs,MinNumOccurence);
else
    error('Wrong option!')
end

% GEMmodel=prepModelforTFA(GEMmodel, DB_AlbertyUpdate, GEMmodel.CompartmentData);

end