function [OriginalGEM, GEMmodel, core_ss, Biomass_rxns, met_pairs_to_remove, InorgMetSEEDIDs, BBBsToExclude, ExtraCellSubsystem, OxPhosSubsystem] = ...
    case_putida_iJN746(GEMname, ZeroZeroGEMbounds, ListForInorganicMets, ListForCofactorPairs, SelectedSubsystems, AddExtracellularSubsystem, DB_AlbertyUpdate)

% Model Specific Settings:
fprintf('Loading the GEM for Putida...\n')
GEM_filename = 'putida_GEM_20150812.mat';
GetModelFromGITresources(GEM_filename)
GEMmodel = load(['./GEMs/',GEM_filename]);
GEMmodel = GEMmodel.ttmodel;
OriginalGEM = GEMmodel;

OxPhosSubsystem = 'Oxidative phosphorylation';

fprintf('For this model:\n')
fprintf('- The core subsystems  are:\n')
if iscell(SelectedSubsystems)
    core_ss = SelectedSubsystems
elseif ischar(SelectedSubsystems) && strcmp(SelectedSubsystems, 'default')
    core_ss = {'TCA Cycle'
    'Pentose Phosphate Pathway'
    'Gluconeogenesis'
    'Glycolysis'
    'Pyruvate Met'}
else
    error('Error defining the subsystems!')
end  

if iscell(AddExtracellularSubsystem)
    fprintf('- The extracellular subsystm includes the following reactions:\n')
    ExtraCellSubsystem = AddExtracellularSubsystem;
elseif ischar(AddExtracellularSubsystem) && strcmp(AddExtracellularSubsystem,'default')
    fprintf('- The extracellular subsystm includes the following reactions:\n')
    ExtraCellSubsystem = {'DM_succ_e'
        'DM_ac_e'
        'DM_glyc_e'
        'DM_lac_D_e'
        'DM_akg_e'
        'DM_mal_L_e'}
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

% 'Oxidative Phosphorylation'} We do not add oxidative
% phosphorylation, because its reactions are replaced later by
% ETC reactions
fprintf('- The biomass reactions are:\n')
Biomass_rxns = {'GROWTH'}
ObjRxn = {'GROWTH'};


% If the flux bounds are too high, set them to +/- 100
% ATTENTION!!! We should make this autmmatic based on the carbon
% molecules of the carbon source!!!
GEMmodel.lb(find(GEMmodel.lb <-100)) = -100;
GEMmodel.ub(find(GEMmodel.ub >+100)) = +100;

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
    met_pairs_to_remove = load('met_pairs_to_remove160127.mat');
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
% GEMmodel=prepModelforTFA(GEMmodel, DB_AlbertyUpdate, GEMmodel.CompartmentData);
end