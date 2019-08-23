function [OriginalGEM, GEMmodel, core_ss, Biomass_rxns, met_pairs_to_remove, InorgMetSEEDIDs, BBBsToExclude, ExtraCellSubsystem, OxPhosSubsystem] = ...
    case_yeast_iMM904bmDrain(GEMname, ZeroZeroGEMbounds, ListForInorganicMets, ...
    ListForCofactorPairs, SelectedSubsystems, AddExtracellularSubsystem, DB_AlbertyUpdate)

% Model Specific Settings:
fprintf('Loading the GEM for Yeast iMM904 ...\n')
GEM_filename = 'iMM904bmDrain.mat';
% GetModelFromGITresources(GEM_filename)
GEMmodel = load(['./GEMs/',GEM_filename]);
GEMmodel = GEMmodel.iMM904bmDrain;
OriginalGEM = GEMmodel;

OxPhosSubsystem = 'Oxidative Phosphorylation';

fprintf('For this model:\n')
fprintf('- The core subsystems  are:\n')
if iscell(SelectedSubsystems)
    core_ss = SelectedSubsystems
elseif ischar(SelectedSubsystems) && strcmp(SelectedSubsystems, 'default')
    core_ss = {'Citric Acid Cycle'
    'Pentose Phosphate Pathway'
    'Glycolysis/Gluconeogenesis'
    'Pyruvate Metabolism'
    'Oxidative Phosphorylation'}
else
    error('Error defining the subsystems!')
end  

if iscell(AddExtracellularSubsystem)
    fprintf('- The extracellular subsystm includes the following reactions:\n')
    ExtraCellSubsystem = AddExtracellularSubsystem;
elseif ischar(AddExtracellularSubsystem) && strcmp(AddExtracellularSubsystem,'default')
    fprintf('- The extracellular subsystm includes the following reactions:\n')
    ExtraCellSubsystem = {'EX_succ(e)';
        'EX_ac(e)';
        'EX_etoh(e)';
        'EX_glyc(e)';
        'EX_lac-D(e)';
        'EX_akg(e)';
        'EX_for(e)';
        'EX_pyr(e)';
        'EX_co2(e)';
        'EX_mal-L(e)'}
elseif ischar(AddExtracellularSubsystem) && strcmp(AddExtracellularSubsystem,'no')
    fprintf('- NO extracellular subsystem defined!')
    ExtraCellSubsystem = [];
elseif ischar(AddExtracellularSubsystem) && strcmp(AddExtracellularSubsystem,'automatic')
    fprintf('Extracting the medium automatically from the GEM\n')
    error('Pierre has started working on this, but needs to be tested!')
    medium = get_medium(GEMmodel);
    ExtraCellSubsystem = table2cell(medium(:,1));
else
    error('Error defining the ExtracellularSubsystem!')
end

fprintf('- The biomass reactions are:\n')
load iMM904_pools.mat

Biomass_rxns=iMM904_pools;
% To be modified!!! Quick fix for robustYeastMeeting!!!
%Biomass_rxns = GEMmodel.rxns(find_cell(Biomass_rxns,GEMmodel.rxnNames));
% ObjRxn = {'yeast 5 biomass pseudoreaction'};
% % To be modified!!! Quick fix for robustYeastMeeting!!!
% ObjRxn = GEMmodel.rxns(find_cell(ObjRxn,GEMmodel.rxnNames));

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

if strcmp(FluxUnits,'mmol')
    % Do nothing
elseif strcmp(FluxUnits,'mumol') && strcmp(Organism,'ecoli')
    GEMmodel.ub = UnitFactor*GEMmodel.ub;
    GEMmodel.lb = UnitFactor*GEMmodel.lb;
    % ATTENTION, THIS IS ONLY FOR THIS ECOLI GEM MODEL!!!!!
    GEMmodel.S(:,[7 8]) = GEMmodel.S(:,[7 8])*UnitFactor;
else
    error('Wrong option!')
end



%% If metabolites appear in the biomass equation in more than one
%% compartments (e.g. cytoplasmic and periplasmic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the biomass there are some metabolites that appear at the same
% time in the cytoplasm and in the periplasm. This causes problems
% because provided that the transports do not belong to the
% non-core-reactions that are used for the lumping, we risk of
% having exactly the same reaction for two different BBBs. This
% creates problems of singularity later on in the non-linear
% models. Find these  metabolites:


bbbz = GEMmodel.mets(find(GEMmodel.S(:,find(GEMmodel.c))<0));
bbbz_no_comp = cellfun(@(x) x(1:end-2),bbbz,'UniformOutput',false);
[unique_bbbz_no_comp,order_unique ]= unique(bbbz_no_comp);
[~,b] = ismember(bbbz_no_comp,unique_bbbz_no_comp);
[c,d] = hist(b,unique(b));
fprintf('The following metabolites appear as bbbs in more than one compartments in the biomass equation')
CP_bbbz_no_comp = unique_bbbz_no_comp(d(find(c>1)))
% Initialize BBBs to exclude:
BBBsToExclude = [];
for i=1:size(CP_bbbz_no_comp)
    SameBBBsDifferentCs = bbbz(find_cell(CP_bbbz_no_comp(i),bbbz_no_comp));
    DifferentCs = cellfun(@(x) {x(end)},SameBBBsDifferentCs);
    if sum(ismember(DifferentCs,{'c','p'}),1)==2
        % Check if there is a transport reaction between these metabilites
        TransRxn_CP = GEMmodel.rxns(find(sum(abs(full(GEMmodel.S(find_cell(SameBBBsDifferentCs,GEMmodel.mets),:))),1)==2));
        if ~isempty(TransRxn_CP)
            % From this pair we keep track of the periplasmic BBB
            % so that we can exclude it from the bbb-list. It will
            % be produced by the transport.
            BBBsToExclude = [BBBsToExclude SameBBBsDifferentCs(find_cell({'p'},DifferentCs))];
            fprintf('There exists a transport reaction between the BBBs %s & %s, reaction: %s.\n Therefore we will not generate lumped reaction for %s\n',...
                SameBBBsDifferentCs{1},SameBBBsDifferentCs{2},TransRxn_CP{1},SameBBBsDifferentCs{find_cell({'p'},DifferentCs)})
        end
    else
        error('Are you sure your biomass equation is assigned correctly?')
    end
end

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
    met_pairs_to_remove = load('met_pairs_to_remove_Yeast.mat');
    met_pairs_to_remove = met_pairs_to_remove.met_pairs_to_remove;
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

% GEMmodel = prepModelforTFA(GEMmodel,DB_AlbertyUpdate,GEMmodel.CompartmentData, false, false);

end