function [RedModel, activeRxns, LumpedRxnFormulas, bbbNames, DPsAll, rxns] = redGEM(varargin)
% This function generates a reduced model of a GEM (GEnomescale Model).
%
% INPUTS
% - L: is an adjustable parameter indicating how far from each subsystem is
%   the algorithm allowed to search in order to find connections among
%   the subsystems. More precicely, it indicates how distant, in terms of
%   reactions, can the metabolites of any two core sybsystems be from any
%   potential interconnecting metabolite.
% - D: given that the algorithm has searched paths of up to length L, D
%   specifies how many lengths to include in the new model. if not
%   specified, this is equal to L.
% - startFromMin: default = 'yes'. if 'yes', indicates that the L's that are
%   included in the model start at the shortest distance between subsystem
%   pairs. if zero, the model is built with the first D=L levels regardless
%   of the presence of paths.
% - viewStats: options to plot statistics on the connecting pathways.
%   default is 'no'
% - OnlyConnectExclusiveMets: if this is 'yes', than only metabolites that
%   are exlusive to each subsystem will be connected. if 'no', even if a
%   metabolite is shared, it will still be connected. default = 'yes'
% - ConnectIntracellularSubsystems: if this is 'yes', than metabolites that are
%   within a subsystem will be connected to eachother. if 'no', they will
%   not be connected. default is 'no';
%   ----------------------------------------------------------------------
%   Example: L=10 D=4 startFromMin= 'yes'
%   the adjacency matrix will be calculated up to length 10, but only the
%   first 4 lengths that connect subsystems will be included. If subsystem
%   A and B have a D=1 connection, then the paths of length 1:4 will be
%   included in the model. However, the shortest distance between
%   subsystems A and C is 3 so all paths of length 3:7 will be included.
%   ----------------------------------------------------------------------
% - stats: options to plot statistics on the connecting pathways. default=0

% OUTPUTS
% - RedModel: Reduced model
% - rxns_total:
% - activeRxns:
% - overall_am:
% - overall_am_names_second:
% - DPs_all:

tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTING UP THE OPTIONS OF THE REDUCTION  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If there is no input to the function, assign a string name to varargin
if isempty(varargin)
    varargin = {'NoInput'};
end

% Set all non-specified flags to empty. In this way, whatever has not been
% specified by the "varargin", will be later on asked with individual
% questions. All the variables that have been specified in "varargin" will
% overwrite these empty settings.
[paramRedGEM, paramNames] = redGEMQuestionnaire(varargin{1});


% Assign the value names to all the parameters of the current redGEM run
eval(cell2mat(cellfun(@(x) [x,'= getfield(paramRedGEM,''',x,''');'],paramNames,'UniformOutput', false)'));

if ~exist([output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname],'dir')
    mkdir([output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname])
end

if ~exist([output_PATH,'/UserOutputs/RunTimes'],'dir')
    mkdir([output_PATH,'/UserOutputs/RunTimes'])
end

connectingpaths_folder = [output_PATH,'/UserOutputs/ConnectingPaths'];
if ~exist(connectingpaths_folder,'dir')
    mkdir(connectingpaths_folder)
end

if ~exist([output_PATH,'/UserOutputs/Models/',Organism],'dir')
    mkdir([output_PATH,'/UserOutputs/Models/',Organism])
end

% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET NECESSARY THE PATHS AND CHECK THE VERSIONS USED  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restore the default path
restoredefaultpath

% Add required paths:
% - Cobra and the tFBA path
addpath(genpath(TFA_PATH))

% - Cplex path
addpath(genpath(CPLEX_PATH))

% - redGEM path
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(path))

if isempty(which('cplex.p'))
    error('You need to add CPLEX to the Matlab-path!!')
end

% Set the solver
changeCobraSolver('cplex_direct','LP')

%% %%%%%%%%%%%%%%%%%%%%%%
%% LOAD NECESSARY FILES %
%% %%%%%%%%%%%%%%%%%%%%%%

% Load the thermodynamic database used
load(thermo_data_PATH)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the properties of the GEM that is reduced %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['[OriginalGEM, GEMmodel, core_ss, Biomass_rxns, met_pairs_to_remove, InorgMetSEEDIDs, BBBsToExclude, ExtraCellSubsystem, OxPhosSubsystem] = '...
    case_filename,'(GEMname, ZeroZeroGEMbounds, ListForInorganicMets, ListForCofactorPairs, SelectedSubsystems, AddExtracellularSubsystem, DB_AlbertyUpdate);'])

%%
% Model Reaction Redundancy Check
redundant_table = ReacRedundancy(GEMmodel);


%% Prevent (close) the uptakes of biomass building blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The reason why we have to do this is because otherwise the reduction
% might skip the generation of lumped reactions for some biomass building
% blocks.
GEMmodel_orig=GEMmodel;
if strcmp(PreventBBBuptake, 'yes')
    GEMmodel = preventBBBuptake(GEMmodel);
end

%% Adding ETC as a subsystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(AddETCAsSubsystem,'yes')
    fprintf('The ETC has been added to the core subsystems\n')
    GEMmodel = generate_ETC_SubSyst(GEMmodel,OxPhosSubsystem);
    core_ss = [core_ss; {'ETC_Rxns'}];
elseif strcmp(AddETCAsSubsystem,'no')
    fprintf('The ETC has not been added to the core subsystems\n')
else
    error('Wrong option!')
end

%% Adding extracellular metabolites as a subsystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(ExtraCellSubsystem)
    [~, ba] = ismember(ExtraCellSubsystem, GEMmodel.rxns);
    GEMmodel.subSystems(ba) = {'ExtraCell'};
    core_ss = [core_ss; {'ExtraCell'}];
    fprintf('Generic ExtraCellular Metabolites have been added to the core subsystems\n')
elseif strcmp(AddExtracellularSubsystem,'no')
    fprintf('Did not add The ExtraCellular Metabolites as a core subsystem\n')
else
    error('Wrong option!')
end


%% For the process of lumping, we can use the entire genomescale model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GSM_ForLumping = GEMmodel;

if strcmp(ImposeThermodynamics, 'no')
    % Remove unnecessary thermo fields (if they exist). If not, functions
    % like the extractSubnetwork might have an error.
    GSM_ForLumping = removeThermoFields(GSM_ForLumping);
end

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_1.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

%% Connecting Subsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we only want to perform the lumping on a predefined subsystem, or set
% of subsystems (i.e. we do not want to connect these subsystem, and/or
% include more reactions) then we can directly jump to the lumping part.
% The above implies that D=0, as we do not wish to expand the

% To find the adjacency matrices and perform all the graph search to obtain
% the desired reaction network we need a genomescale model without biomass,
% and without the metabolite pairs:
GSM_ForAdjMat = GEMmodel;
% Find and remove the biomass reactions
[~, BMIds] = ismember(Biomass_rxns, GSM_ForAdjMat.rxns);

GSM_ForAdjMat   = removeRxns(GSM_ForAdjMat, GSM_ForAdjMat.rxns(BMIds));

% We directly remove cofactors and small ionrganic metabolites from the S matrix
Smatrix = removeMetPairsFromS(GSM_ForAdjMat, met_pairs_to_remove, InorgMetSEEDIDs);
GSM_ForAdjMat.S = Smatrix;

f = filesep;

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_2.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

if L~=0 % if we want to connect subsystems
    fprintf('Calculating paths that connect subsystems\n')
    
    % To obtain paths that connect subsystems, create the adjacency matrix and associated matrices
    [~, DirAdjMatWn, DirAdjMatC] = adjacentFromStoich_150218(GSM_ForAdjMat);
    
    [L_DirAdjMatWn, L_DirAdjMatC, L_DirAdjMatMets] = allpaths(L, DirAdjMatWn, DirAdjMatC);
    kk=1;
    ind=[];
    for i=2:L
        for j=1:size(L_DirAdjMatC,1)
            for k=1:size(L_DirAdjMatC,2)
                rxns=L_DirAdjMatC(j,k,i);
                mets=L_DirAdjMatMets(j,k,i);
                if ~(isempty(cell2mat(rxns)))
                    rxns=cell2mat(rxns);
                    mets=cell2mat(mets);
                    rxns_new=rxns;
                    mets_new=mets;
                    for m=1:size(rxns,2)
                        rxns_small=rxns(:,m);
                        rxns_small_un=unique(rxns_small);
                        if length(rxns_small_un)~=length(rxns_small)
                            ind(kk,1)=m;
                            kk=kk+1;
                        end
                    end
                    if ~isempty(ind)
                        rxns_new(:, ind)=[];
                        mets_new(:, ind)=[];
                    end
                    L_DirAdjMatC{j,k,i}=rxns_new;
                    L_DirAdjMatMets{j,k,i}=mets_new;
                    
                    if length(ind)==size(rxns,2)
                        L_DirAdjMatWn(j,k,i)=0;
                    end
                    ind=[];
                    kk=1;
                end
            end
        end
    end
    eval(['structAllpaths_',Organism,'_L',num2str(D),'.ListForInorganicMets = ListForInorganicMets;'])
    eval(['structAllpaths_',Organism,'_L',num2str(D),'.ListForCofactorPairs = ListForCofactorPairs;'])
    eval(['save ' connectingpaths_folder f 'structAllpaths_',Organism,'_L',num2str(D),'.mat structAllpaths_',Organism,'_L',num2str(D),' DirAdjMatWn L_DirAdjMatWn L_DirAdjMatC L_DirAdjMatMets'])
    clear DirAdjMatWn DirAdjMatC %these two are included in countD and connectD
    
    
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_3.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
    % to select the reactions that join subsystems, extract the information.
    [selectedPaths, selectedPathsExternal, connecting_reactions, rxns_ss, mets_ss, otherReactions] = ...
        extractAdjData(GSM_ForAdjMat,core_ss,L_DirAdjMatWn, L_DirAdjMatC, L_DirAdjMatMets,D,startFromMin,...
                       OnlyConnectExclusiveMets,ConnectIntracellularSubsystems,ApplyShortestDistanceOfSubsystems,ThrowErrorOnDViolation);
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_4.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
else %% No Connections
    rxns_ss = cell(length(core_ss),1);
    % For each of the subsystems:
    for i=1:length(core_ss)
        % - find the reactions of this subsystem, and their indices
        [~, rxns_ss{i}]=ismember(GSM_ForAdjMat.subSystems, core_ss{i});
        rxns_ss{i}=find(rxns_ss{i});
    end
    connecting_reactions = [];
end


all_core_rxns = unique([connecting_reactions;vertcat(rxns_ss{:})]);

submodel = extractSubNetwork(GEMmodel, GSM_ForAdjMat.rxns(all_core_rxns));
% We go through all GS-reactions, and check if there exists a reactions
% that involves ONLY core metabolites. If yes, we keep it, and add it as a
% core (i.e. tightening)
k = 1;
for i=1:length(GEMmodel.rxns)
    mets = GEMmodel.mets(find(GEMmodel.S(:, i)));
    [~, ba] = ismember(mets, submodel.mets);
    if isempty(find(ba==0))
        core_rxns(k,1) = i;
        k=k+1;
    end
end
clear k;


if strcmp(performLUMPGEM, 'yes')
    
    otherReactions = setdiff(1:length(GSM_ForLumping.rxns),core_rxns);
    
    % Transform indices to GEM model (GSM_ForLumping) indexing
    [~, otherReactionsGSMForLump_idx] = ismember(GSM_ForLumping.rxns(otherReactions), GSM_ForLumping.rxns);
    
    %% Prevent thermodynamic relaxation of ETC reactions and/or ATP-synthase
    ATPsynth_RxnNames = [];
    ATPsynth_RxnIds = [];
    % In the case of thermodynamic feasibility when we add the thermodynaic
    % constraints, there might exist reactions in the model for which we would
    % prefer not to allow relaxation of their corresopnding DGo bounds. Such
    % reactions are
    % - ATP synthases, e.g.:
    % 'adp_c + pi_c + 4 h_p <=> atp_c + 3 h_c + h2o_c'
    % adp    'cpd00008'
    % pi     'cpd00009'
    % h      'cpd00067'
    % atp    'cpd00002'
    % h2o    'cpd00001'
    SEEDIDs_atpSynth = {'cpd00008','cpd00009','cpd00067','cpd00002','cpd00001'}';
    % Stoich matrix with binary entries
    Stoic01 = GSM_ForLumping.S~=0;
    for i = 1:size(Stoic01, 2)
        SEEDIDs_Rxn_i = unique(GSM_ForLumping.metSEEDID(Stoic01(:,i)));
        if isempty(setxor(SEEDIDs_Rxn_i, SEEDIDs_atpSynth))
            ATPsynth_RxnIds = [ATPsynth_RxnIds i];
        end
    end
    ATPsynth_RxnNames = GSM_ForLumping.rxns(ATPsynth_RxnIds);
    % - and reactions of the ETC sybsystem
    ETC_RxnNames = GSM_ForLumping.rxns(find_cell('ETC_Rxns', GSM_ForLumping.subSystems));
    RxnNames_PrevThermRelax =  [ETC_RxnNames; ATPsynth_RxnNames];
    
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_5.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
    if strcmp(performREDGEMX,'yes')
        fprintf('Connecting the metabolites from the extracellular medium to the core\n')
        % connect extracellular subsystem to core
        [ConnectExtrCell_rxns_all, ConnectExtrCell_id_all, sol_all] = ...
            redGEMX(rxns_ss, GSM_ForLumping, GSM_ForAdjMat, OriginalGEM, ...
                    GEMmodel, Organism, GEMname, NumOfConnections, CplexParameters, DB_AlbertyUpdate, ImposeThermodynamics,output_PATH);
        
        otherReactionsGSMForLump_idx = setdiff(otherReactionsGSMForLump_idx,unique(ConnectExtrCell_id_all));
        
        % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
        [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
        eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_5a.mat;'])   %
        % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
        
    elseif strcmp(performREDGEMX,'no')
        fprintf('The metabolites from the extracellular medium were not connected to the core\n')
    else
        error('Wrong option!')
    end
    
    
    [activeRxns, LumpedRxnFormulas, bbbNames, DPsAll, IdNCNTNER, relaxedDGoVarsValues_ForEveryLumpedRxn] = ...
        addBIOMASS(GSM_ForLumping, otherReactionsGSMForLump_idx, DB_AlbertyUpdate, BBBsToExclude, AerobicAnaerobic, ...
                   Organism, AlignTransportsUsingMatFile, TimeLimitForSolver, NumOfLumped, CplexParameters, ...
                   GEMname, RxnNames_PrevThermRelax, Biomass_rxns, ATPsynth_RxnNames, addGAM, PercentOfmuMaxForLumping, ...
                   ImposeThermodynamics, output_PATH);
    
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_6.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
    [~, indCORE] = ismember(1:length(GSM_ForLumping.rxns), IdNCNTNER);
    indCORE = find(indCORE==0);
    [RedModel, rxns ,info_LMDP] = ...
        add_lumped_reactions(GSM_ForLumping, LumpedRxnFormulas, bbbNames, DB_AlbertyUpdate, ...
                             GSM_ForLumping.rxns(indCORE), ImposeThermodynamics, CplexParameters, Organism, GEMname, activeRxns, output_PATH);
    
    if ~isempty(rxns)
        error('Additional reactions are required for growth!!')
    end
    
    if strcmp(performREDGEMX, 'yes')
        RedModel.info_Econnect = [OriginalGEM.ExtracellularMedium_connect ConnectExtrCell_rxns_all'];
    end
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_7.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
    %add ATPm reaction. With the tightening ATPM is probably already incuded
    if isempty(find_cell({'ATPM'},RedModel.rxns))
        RedModel = addATPm(RedModel, GSM_ForLumping);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Extract the gene information from the genome-scale model %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some of the genes might be unused now, because their corresponding
    % reactions are not in the reduced model
    
    % - Get the index if the Lumped reactions in the model
    idLUMPEDRxnsInRedGEM = find( ~cellfun(@(x) isempty(x), regexpi(RedModel.rxns ,'LMPD_')) );
    % - Get the index of the reduced model reactions that belong to the GEM. The indices
    % should correspond to all reactions until just before the first lumped reaction:
    idGEMRxnsInRedGEM = 1 : 1 : (idLUMPEDRxnsInRedGEM(1) - 1);
    
    % Checkpoint! If any reaction that is not lumped does not belong to the
    % GEM, display an error:
    if ~isempty(setdiff(RedModel.rxns(idGEMRxnsInRedGEM), GEMmodel_orig.rxns))
        setdiff(RedModel.rxns(idGEMRxnsInRedGEM), GEMmodel_orig.rxns)
        error('This reaction(s) is not part of the original GEM!!')
    end
    
    RedModelRightAfterLumping = RedModel;
    % This is the reduced model right after redGEM (without post-processing)
    RedModel.RedModelRightAfterLumping = RedModelRightAfterLumping;
    %%% UP TO HERE IS THE REDUCED MODEL GENERATION. AFTER THIS POINT IS JUST
    %%% POST PROCESSING
    
elseif strcmp(performLUMPGEM, 'no')
    fprintf('We do not perform lumping')
    % Set empty values to outputs
    activeRxns = [];
    LumpedRxnFormulas = [];
    bbbNames = [];
    DPsAll = [];
    rxns = [];
else
    error('Wrong option')
end


% Store the core reduced model stoichiometry in the output (without the lumped reactions and without the transports)
RedModel.core_model = extractSubNetwork(GEMmodel, GEMmodel.rxns(core_rxns));

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_8.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ->->->->->->->  P O S T    P R O C E S S I N G ->->->->->-> %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Blocked reactions (in FBA and TFA)

if strcmp(performPostProcessing, 'yes')
    
    % Run fluxVariability analysis, and remove all those reactions that are
    % unable to carry any flux, while the model is still able to produce biomass:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set growth bounds of reduced model to original GEM growth bounds:
    origGEM_growth_LB = OriginalGEM.lb(find_cell(RedModel.rxns(find(RedModel.c)), OriginalGEM.rxns));
    origGEM_growth_UB = OriginalGEM.ub(find_cell(RedModel.rxns(find(RedModel.c)), OriginalGEM.rxns));
    if origGEM_growth_UB>100
        origGEM_growth_UB = 100;
    end
    RedModel.lb(find(RedModel.c)) = origGEM_growth_LB;
    RedModel.ub(find(RedModel.c)) = origGEM_growth_UB;
    
    % We set a MINIMAL BOUND FOR GROWTH of 10% max growth to avoid deleting reactions
    % preventing growth !!!!!!
    temp_FBASol = solveFBAmodelCplex(RedModel,scalPar, feasTol, emphPar);
    RedModel.lb(find(RedModel.c))=temp_FBASol.f*0.1;
    
    FBA_minmax = runMinMax(RedModel,[],[],scalPar, feasTol, emphPar);
    FVA_BlockedRxnIds = find(abs(FBA_minmax(:,1))<feasTol & abs(FBA_minmax(:,2))<feasTol);
    FVA_BlockedRxnNames = RedModel.rxns(FVA_BlockedRxnIds);
    warning(['Using FVA we find that ',num2str(size(RedModel.rxns(FVA_BlockedRxnIds),1)),' reactions appear to be blocked'])
    warning('We remove these reactions!')
    rxns_to_remove = RedModel.rxns(FVA_BlockedRxnIds);
    sol_obj = solveFBAmodelCplex(RedModel,scalPar, feasTol, emphPar);
    k = 1;
    problematic_zero_minmax = {};
    % ATTENTION!! We need to finalize the removeRxns and removeRxnsThermo. Now they are
    % mixed and used in different contexts without being sure which fields to
    % include or not. In the meantime I just complete the following two fields
    % (rules & rxnDeltaGerr) to avoid the error:
    RedModel.rules = [RedModel.rules;repmat({''},size(RedModel.rxns,1)-size(RedModel.rules,1),1)];
    if strcmp(ImposeThermodynamics, 'yes')
        RedModel.rxnDeltaGRerr = [RedModel.rxnDeltaGRerr;repmat(10^7, size(RedModel.rxns,1)-size(RedModel.rxnDeltaGRerr,1) ,1)];
    end
    for i = 1:length(RedModel.rxns(FVA_BlockedRxnIds))
        RedModel_orig = RedModel;
        RedModel = removeRxns(RedModel,rxns_to_remove(i));
        sol = solveFBAmodelCplex(RedModel,scalPar, feasTol, emphPar);
        if sol.f < roundsd(sol_obj.f, 4, 'floor')
            RedModel = RedModel_orig;
            problematic_zero_minmax(k,1) = rxns_to_remove(i);
            k = k+1;
        end
    end
    clear k
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_9.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
    %% Convert the reduced model to TFBA
    if strcmp(ImposeThermodynamics, 'yes')
        
        RedFBASol = solveFBAmodelCplex(RedModel,scalPar, feasTol, emphPar);
        RedFBASol = roundsd(RedFBASol.f, 4, 'floor');
        
        RedModel = prepModelforTFA(RedModel, DB_AlbertyUpdate, GEMmodel.CompartmentData,  false, false);
        [RedModel, relaxedDGoVarsValues_inPP1a] = convToTFA(RedModel, DB_AlbertyUpdate, [], 'DGo', RxnNames_PrevThermRelax, 0.1*RedFBASol);
        
        RedTFASol = solveTFAmodelCplex(RedModel,[],[], mipTolInt, emphPar, feasTol, scalPar);
        RedTFASol = roundsd(RedTFASol.val, 4, 'floor');
        
        GEM_FBASol = solveFBAmodelCplex(GEMmodel,scalPar, feasTol, emphPar);
        GEM_FBASol = roundsd(GEM_FBASol.f, 4, 'floor');
        
        if abs(RedTFASol - GEM_FBASol)/GEM_FBASol < 0.05
            fprintf('The reduced model with the thermodynamic constraints is in very good agreement with the original GEM:\n')
            fprintf('- GEM-FBA: %s\n', num2str(GEM_FBASol))
            fprintf('- redGEM-FBA: %s\n', num2str(RedFBASol))
            fprintf('- redGEM-TFA: %s\n', num2str(RedTFASol))
        end
        
        % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
        [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
        eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_10.mat;'])   %
        % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
        
        % IN THE GENERATED-OUTPUT MODEL: Aligning the transport reactions that transport the same metabolite
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        checkgrowth = 1;
        RedModel = AlignTransportHelperFun(RedModel,AlignTransportsUsingMatFile, checkgrowth, CplexParameters, Biomass_rxns, ATPsynth_RxnNames);
        
        % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
        [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
        eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_11.mat;'])   %
        % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
        
        
        % Run Thermodynamic-fluxVariability analysis, and remove all those
        % reactions that are unable to carry any flux, while the model is still
        % able to produce biomass:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First we need to add the net-flux variables. They are
        RedModel = addNetFluxVariables(RedModel);
        
        temp_TFASol = solveTFAmodelCplex(RedModel,[],[],mipTolInt,emphPar,feasTol,scalPar);
        RedModel.var_lb(find(RedModel.f))=temp_TFASol.val*0.7;
        
        NFids = getAllVar(RedModel,{'NF'});
        TFVA_BlockedRxnNames = [];
        Tminmax = runTMinMax(RedModel,RedModel.varNames(NFids),30,mipTolInt,emphPar,feasTol,scalPar);
        TFVA_BlockedRxnIds = find(abs(Tminmax(:,1))<feasTol & abs(Tminmax(:,2))<feasTol);
        if ~isempty(TFVA_BlockedRxnIds)
            warning(['Using TFVA we find that ',num2str(size(RedModel.rxns(TFVA_BlockedRxnIds),1)),' reactions appear to be blocked, so we remove these reactions!'])
            TFVA_BlockedRxnNames = cellfun(@(x) x(4:end), RedModel.varNames(NFids(TFVA_BlockedRxnIds)), 'UniformOutput', false);
            k = 1;
            problematic_zero_minmax = {};
            for i = 1:length(RedModel.rxns(TFVA_BlockedRxnIds))
                RedModel_orig = RedModel;
                RedModel = removeRxns(RedModel,TFVA_BlockedRxnNames(i));
                sol = solveFBAmodelCplex(RedModel,scalPar, feasTol, emphPar);
                if sol.f < roundsd(temp_TFASol.val, 4, 'floor')
                    RedModel = RedModel_orig;
                    problematic_zero_minmax(k,1) = TFVA_BlockedRxnNames(i);
                    k = k+1;
                end
            end
        end
        clear k
        
        % Now convert again to tFA
        RedModel = prepModelforTFA(RedModel, DB_AlbertyUpdate, GEMmodel.CompartmentData, false, false);
        verboseFlag = false;
        [RedModel, relaxedDGoVarsValues_inPP1b] = convToTFA(RedModel, DB_AlbertyUpdate, [], 'DGo', [], roundsd(sol_obj.f, 2, 'floor'), [], [], verboseFlag);
        
        % Add net-flux variables
        RedModel = addNetFluxVariables(RedModel);
        
        % DGo that had to be relaxed for thermodynamic feasibility
        % - for thermodynamic feasibility of each lumped
        RedModel.relaxedDGoVarsValues_ForEveryLumpedRxn = relaxedDGoVarsValues_ForEveryLumpedRxn;
        RedModel.relaxedDGoVarsValues_ForAddingAllLumpedRxns = [relaxedDGoVarsValues_inPP1a; relaxedDGoVarsValues_inPP1b];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Keep track of the reactions that were removed from the reduced model
        % because they had zero FVA & tFVA fluxes:
        RedModel.PostProcessing.TFVA_BlockedRxnNames = TFVA_BlockedRxnNames;
    end
    RedModel.PostProcessing.FVA_BlockedRxnNames = FVA_BlockedRxnNames;
    
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_12.mat;'])   %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
elseif strcmp(performPostProcessing, 'no')
    fprintf('No post-processing!\n')
else
    error('Wrong post-processing option!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OOOOO    U    U  TTTTT  PPPPP  U    U  TTTTT
%  O     O   U    U    T    P    P U    U    T
%  O     O   U    U    T    P    P U    U    T
%  O     O   U    U    T    PPPPP  U    U    T
%  O     O   U    U    T    P      U    U    T
%   OOOOO     UUUU     T    P       UUUU     T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ->->->->->->  PREPARE OUTPUT ->->->->-> %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add date to the generation of the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dateStr, timeStr] = getDateTimeStrings(date,clock);

% Create a model name based on the date and time of the generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_name = [RedModelName,'_',dateStr,'_',timeStr];

%% Appending necessary informative fields to the model
% Keep the original GEM from which the reduced model has been derived
RedModel.OriginalGEM = OriginalGEM;

% Correct the genes field in the reduction
RedModel = genes_redGEM(RedModel,GEMmodel_orig);

% Store in the model some of the properties of the reduction
RedModel.ReductionProperties = paramRedGEM;

% Store in the model the core/initial subsystems that are selected for the
% reduction
RedModel.initialSubSystems = core_ss;

% Add a description field. ATTENTION: this is very important for the ORACLE
% workflow later on!!
RedModel.description = [RedModelName,'_Date_',dateStr,'_',timeStr,'_Generated_by_redGEMv1'];

% workflow later on!!
RedModel.solverPath = CPLEX_PATH;

% Add a field with the name of the user that generated this model
user_str = regexp(pwd,'[\\/]Users[\\/](\w+)[\\/]','tokens');
user_str = user_str{1}{1};
RedModel.GeneratedByUser = user_str;

%% Saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How much time it took for all these?
rG_RT = toc;

% Save the runtime of the generation
RedModel.ReductionRuntime = rG_RT;

% Save the model in the output folder
eval([model_name,' = RedModel;']);
eval(['save ',output_PATH,'/UserOutputs/Models/',Organism,'/',model_name,'.mat ', model_name,';']);


end
