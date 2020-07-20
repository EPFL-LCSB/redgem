function [activeRxns, LumpedRxnFormulas, bbbNames, DPsAll, IdNCNTNER, relaxedDGoVarsValues_ForEveryLumpedRxn] = ...
        addBIOMASS(GSM_ForLumping, otherReactionsGSMForLump_idx, DB_AlbertyUpdate, BBBsToExclude, AerobicAnaerobic,...
                   Organism, AlignTransportsUsingMatFile, TimeLimitForSolver, NumOfLumped, CplexParameters, ...
                   GEMname, RxnNames_PrevThermRelax,biomassRxnNames,ATPsynth_RxnNames, addGAM, PercentOfmuMaxForLumping, ...
                   ImposeThermodynamics, output_PATH)
% This function
%
% INPUTS
% - GSM_ForLumping: Genome scale model used for the lumping
% - otherReactionsGSMForLump_idx: Indices of non-core reactions. These indices
%   are according to GSM_ForLumping indexing.
% - InorgMetSEEDIDs:
% - met_pairs_to_remove:
% - DB_AlbertyUpdate:
% - AerobicAnaerobic:
%
% OUTPUTS
% - rxns_total: these are the reactions that participate in each lumped
%   reaction. For each bbb, is a set of reactions that are required to
%   constitute a lumped reaction.
% - activeRxns: Is a cell, where the first column corresponds to the
%   indices of the reactions that participate in each lumprd reaction, and
%   the second column to the actual reaction formulas of these reactions.
% - LumpedRxnFormulas: The reaction formulas of every lumped reaction
% - bbbNames: Names of the bbb that is produced by each lumped
% - DPsAll: Solutions of the system corresponding to distinct DPs. It is a
%   cell with number of cells equal to the number of biomass building
%   blocks. For each bbb there are so many DPs as defined by NumDPperBBB.


% This is the maximum theoretical yield under the particular media.
sol_obj = solveFBAmodelCplex(GSM_ForLumping);

% We put a lower bound to the muMax of biomass for the Lumping
muMax = roundsd((PercentOfmuMaxForLumping/100)*sol_obj.f, 5, 'floor');

% muMax = sol_obj.f; %0.05;%CHANGED FOR HUMAN!!!!! THIS SHOULD BE AN INPUT IN REDGEM
if muMax == 0
    muMax = 0.1;
end
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

if strcmp(NumOfLumped,'OnePerBBB')
    %flag: for each lumped reaction, how many alternative lumped reactions should the solver look for?
    NumAltLUMPedPerBBB = 1;
    SizeLUMPNetwork = 0.5;
elseif strcmp(NumOfLumped,'Smin')
    NumAltLUMPedPerBBB = 10^4;
    SizeLUMPNetwork = 0.5;
elseif strcmp(NumOfLumped,'Sminp1')
    NumAltLUMPedPerBBB = 10^4;
    SizeLUMPNetwork = 1.5;
elseif strcmp(NumOfLumped,'Sminp2')
    NumAltLUMPedPerBBB = 10^4;
    SizeLUMPNetwork = 2.5;
elseif strcmp(NumOfLumped,'Sminp3')
    NumAltLUMPedPerBBB = 10^4;
    SizeLUMPNetwork = 3.5;
else
    error('Wrong option!')
end



DPsAll = {};

% For easier operations we just rename the other reactions to indices of
% Non-core reactions (IdNCR)
IdNCR = otherReactionsGSMForLump_idx;
IdCR = setdiff(1:size(GSM_ForLumping.rxns,1),IdNCR);
% make sure biomass isnt included in the non-core list (which is used for lumping)
[ab, ba] = ismember(find(GSM_ForLumping.c), IdNCR);
if ab==1
    IdNCR(ba)=[];
end

% create outward drains for substrates of biomass rxn, and inward
% forproducts (why not the other way around?)
bbb_metnames = GSM_ForLumping.mets(find(GSM_ForLumping.S(:,find(GSM_ForLumping.c))<0));%bbb metnames are all the bbb

LumpCstModel = GSM_ForLumping;
LumpCstModel = checkifDrain(LumpCstModel);

for i = 1:size(bbb_metnames,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bbb_i_DM_RxnName = ['DM_' regexprep(bbb_metnames{i},'-','_')];
    stoichCoeffOfDrainedBBB = -1;
    revFlag    = true;
    lowerBound = 0;
    upperBound = 60;
    objCoeff   = 0;
    subSystem  = 'Demand';
    [LumpCstModel, DMrxnIDexists ] = addReaction(LumpCstModel, bbb_i_DM_RxnName, bbb_metnames(i), stoichCoeffOfDrainedBBB, revFlag, lowerBound, upperBound, objCoeff, subSystem);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(DMrxnIDexists)
        % if this demand reaction already exists, change the  name to the predefined one, so that
        % we can call it correctly below when we perform he check of lumping for each bbb
        fprintf('Renaming reaction: %s  to  %s\n',LumpCstModel.rxns{DMrxnIDexists},bbb_i_DM_RxnName)
        LumpCstModel.rxns{DMrxnIDexists} = bbb_i_DM_RxnName;
        % Change also the bounds to [0 60] to be sure that the bbb can be excreted
        LumpCstModel.lb(DMrxnIDexists) = 0;
        LumpCstModel.ub(DMrxnIDexists) = 60;
    else
        if strcmp(ImposeThermodynamics, 'yes')
            LumpCstModel.rxnMapResult{end+1} = 'drain flux';
        end
        LumpCstModel.isDrain(end+1) = 1;
        LumpCstModel.isTrans(end+1) = 0;
    end
end

% - LumpCstModel is based on GSM_ForLumping, and has the following:
%   (i) Adds demand-reactions (DM_) for all the biomass building blocks
%   (ii) Includes the DM_GAM_c reaction for feasibility
%   (iii) Adds thermodynamic constraints for the core system, using
%   convToTFA_BBB (generates thermodynamic variables ONLY for the core system)
%   (iv) sets the upper bound of all fluxes based on the maximum allowable
%   in the system, i.e.: 6* glucose uptake (i.e. 6*10=60 mmol/gDWhr)
%   (v) includes or not constraints to align all transports for same met in
%   same direction
%   (vi) if it is AerobicAnaerobic (AerobicAnaerobic=1), allow oxygen
%   (vii) create new binary variable that we will use to minimize the # of
%   reactions in the lumped reaction
%   (viii) now create constraint that is 1*FU + 1*BU + BFUSE < 1
%   interpretation: if BFUSE is on, FU and BU are off. We will maximize #
%   of BFUSE to get the fewest active reactions.

% This part added the hydrolysis reaction, but we have to make it automatic
% for other genome scale models: Growth associated maintenance (GAM)
if strcmp(addGAM, 'yes')
    [GAMequation, ppi_equation] = extractGAMcoefficientsFromBiomass(LumpCstModel);
    LumpCstModel = addReaction(LumpCstModel, 'DM_GAM_c', GAMequation);
    LumpCstModel.lb(length(LumpCstModel.rxns)) = 0;
    LumpCstModel.ub(length(LumpCstModel.rxns)) = 0;
    GAM_c = 'GAM_c';
    GAM_StCoeff = -1;
else
    GAM_c = [];
    GAM_StCoeff = [];
end

if isfield(LumpCstModel,'A')
    solTFA = solveFBAmodelCplex(LumpCstModel);
    minObjSolVal = roundsd(0.9*solTFA.f, 1, 'floor');
else
    minObjSolVal = sol_obj.f;% 1e-3;
end

if strcmp(ImposeThermodynamics, 'yes')
    LumpCstModel = prepModelforTFA(LumpCstModel, DB_AlbertyUpdate, LumpCstModel.CompartmentData, false, false, false);
    verboseFlag = false;
    [LumpCstModel, relaxedDGoVarsValues] = convToTFA(LumpCstModel, DB_AlbertyUpdate, [], 'DGo', RxnNames_PrevThermRelax, minObjSolVal, [], [], verboseFlag);
    LumpCstModel.relaxedDGoVarsValues = relaxedDGoVarsValues;
elseif strcmp(ImposeThermodynamics, 'no')
    fprintf('In the current case we do not impose thermodynamic constraints\n')
    LumpCstModel = convFBA2MILP(LumpCstModel, DB_AlbertyUpdate);  
else
    error('Wrong option for ImposeThermodynamics')
end

% Below we would like to investigate the possibility of the model to produce each of the
% bbbs without any production of biomass. Therefore we set max growth to 0:
LumpCstModel.var_ub(find(LumpCstModel.f)) = 0;
% It could be however that the lb of the bimoass is greater than zero. In
% this case we would have to set it to zero for this "exercise". Otherwise
% we would have inconsistent bound definition:
if LumpCstModel.var_lb(find(LumpCstModel.f)) > 0
    LumpCstModel.var_lb(find(LumpCstModel.f)) = 0;
end
LumpCstModel.lb(find(LumpCstModel.c))=0;
LumpCstModel.ub(find(LumpCstModel.c))=0;

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > >  > > > > > > >  > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_1.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < 

%% Aligning the transport reactions that transport the same metabolite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
checkgrowth = 0;
LumpCstModel = AlignTransportHelperFun(LumpCstModel, AlignTransportsUsingMatFile, checkgrowth, CplexParameters, biomassRxnNames, ATPsynth_RxnNames);

% Mt = LumpCstModel;
% Mt.var_ub(Mt.f==1) = 10;
% sol_obj_Mt = solveTFAmodelCplex(Mt);
% 
% if sol_obj_Mt.val < sol_obj.f
%     % This is the maximum theoretical yield under the particular media.
%     sol_obj = solveTFAmodelCplex(Mt);
%     % We put a lower bound to the muMax of biomass for the Lumping
%     muMax = roundsd((PercentOfmuMaxForLumping/100)*sol_obj_Mt.val, 5, 'floor');
% end


%% Aerobic/Anaerobic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set it from R_DM_o2_e, because a model doesnt have to have net fluxes.

for i = 1:length(LumpCstModel.rxns)
    mets_rxns = LumpCstModel.mets(find(LumpCstModel.S(:,i)));
    if ismember('o2_e', mets_rxns)
        if length(mets_rxns) == 1
            oxygen_exchange = LumpCstModel.rxns(i);
            [~, id_oxy] = ismember(strcat('R_', oxygen_exchange), LumpCstModel.varNames);
            id_oxyFBA=i;
            break;
        end
    end
end
if strcmp(AerobicAnaerobic,'aerobic')
    LumpCstModel.var_ub(id_oxy) = 20;
    LumpCstModel.lb(id_oxyFBA) = -20;
elseif strcmp(AerobicAnaerobic,'anaerobic')
    LumpCstModel.var_ub(id_oxy) = 0;
    LumpCstModel.lb(id_oxyFBA) = 0;
else
    error('Wrong option')
end

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > >  > > > > > > >  > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_2.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < 


% Now keep all reactions in non-core metabolism that are not transports
LumpCstModel.f = zeros(length(LumpCstModel.varNames),1);

% Extract a subnetwork with only the core-reactions
mod_core = extractSubNetwork(LumpCstModel, LumpCstModel.rxns(IdCR));
mets_core = [mod_core.mets; bbb_metnames];
LumpCstModel = isTrans_GEMS_MoreFields(LumpCstModel, mets_core);
is_core_trans = LumpCstModel.isCoreTrans(IdNCR);
% Keep only those Non-core reactions, that are also not transport reactions
IdNCNTR   = IdNCR(is_core_trans==0);

%keep in ind only those that arent exchange reactions 

EX_rxns = extract_drains(LumpCstModel);
for i=1:length(LumpCstModel.rxns(IdNCNTR))
    if ismember(LumpCstModel.rxns(IdNCNTR(i)), EX_rxns)
        exch_pattern(i) = 0;
    else
        exch_pattern(i) = 1;
    end
end

% These are the indices of the LumpCstModel that are NON-core, NOT
% transport, and also NOT exchange reactions
IdNCNTNER = IdNCNTR(find(exch_pattern));

% Get variable indeces for F and R of the IdNCNTNER
forward = getAllVar(LumpCstModel, {'F'});
forward = forward(IdNCNTNER);
backward = getAllVar(LumpCstModel, {'R'});
backward = backward(IdNCNTNER);

%create new binary variable that we will use to minimize the # of reactions in
%the lumped reaction
for i=1:length(IdNCNTNER)
    LumpCstModel.varNames(length(LumpCstModel.varNames)+1) = strcat('BFUSE_', LumpCstModel.rxns(IdNCNTNER(i)));
    LumpCstModel.var_ub(length(LumpCstModel.varNames)) = 1;
    LumpCstModel.var_lb(length(LumpCstModel.varNames)) = 0;
    LumpCstModel.vartypes(length(LumpCstModel.varNames)) = {'B'};
    LumpCstModel.f(length(LumpCstModel.varNames)) = 1; %this line makes objective function maximize BFUSE vars
end

% Find the indidces of the FUSE and BUSE variables:
ind_fu = getAllVar(LumpCstModel, {'FU'});
ind_bu = getAllVar(LumpCstModel, {'BU'});

ind_fu = ind_fu(IdNCNTNER);
ind_bu = ind_bu(IdNCNTNER);

% Add constraints BFMER5_ : FUSE + BUSE + BFUSE < 1
ind_bfuse=getAllVar(LumpCstModel, {'BFUSE'});
[num_constr,~] = size(LumpCstModel.A);
if length(forward)==length(backward)
    for i=1:length(forward)
        LumpCstModel.rhs(num_constr+i,1) =  1;
        LumpCstModel.constraintNames{num_constr+i,1} = strcat('BFMER5_',num2str(i));
        LumpCstModel.constraintType{num_constr+i,1} = '<';
        LumpCstModel.A(num_constr+i,ind_bfuse(i)) = 1;
        LumpCstModel.A(num_constr+i, ind_fu(i)) = 1;
        LumpCstModel.A(num_constr+i, ind_bu(i)) = 1;
    end
else
    error('length(forward) not equal to length(backward)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude water and bbbs that appear in more than one compartments in the
% biomass equation
bbb_not_to_lump = [{'h2o_c'} BBBsToExclude];
% Get indeces of forward and backward drains of the BBBs
% Create forward/reverse drain reactions for all biomass building blocks:
% adding the hydrolysis as a bbb_element that needs to be also balanced
% And also for the additional hydrolysis reaction

bbb_drains_forward  = strcat('F_DM_', [bbb_metnames; GAM_c]);
bbb_drains_backward = strcat('R_DM_', [bbb_metnames; GAM_c]);
bbb_drains_forward  = strrep(bbb_drains_forward,  '-', '_');
bbb_drains_backward = strrep(bbb_drains_backward, '-', '_');

[~, ind_bbb_forward]  = ismember(bbb_drains_forward,  LumpCstModel.varNames);
[~, ind_bbb_backward] = ismember(bbb_drains_backward, LumpCstModel.varNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LumpCstModel.S = full(LumpCstModel.S);
[~, bbb] = ismember(bbb_metnames, LumpCstModel.mets);
bbb = bbb(find(bbb));
stoich_bbb = LumpCstModel.S(bbb, find(LumpCstModel.c));
bbb_metnames = [bbb_metnames; GAM_c];
stoich_bbb = [stoich_bbb; GAM_StCoeff];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time limit for the solver
TimeLimit = 10*60; % 10 mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin lumping for each BBB

allBFUSEind = getAllVar(LumpCstModel, {'BFUSE'});
rxnsNamesTotal = {};
rxnFormulasTotal = {};

% We would like to avoid that any bbb (except atp, h2o, and GAM, that are
% not really bbbs) will be overproduced and hence accumulated to system,
% leading to infeasibilities.
% - Exclude the GAM related components of the bbbs
% not_to_add = {'atp_c' 'h2o_c' 'GAM_c' 'gly_c' 'cys_L_c'};
not_to_add = [{'atp_c' 'h2o_c'} GAM_c]; % NOTE:this is the default list of
%mets not_to_add, though in the human models sometimes we need to include
%others.
[~, bnta] = ismember(bbb_metnames, not_to_add);

ind_bta=find(bnta==0);

LumpCstModel.var_ub(ind_bbb_forward(ind_bta)) =  roundsd(-muMax*stoich_bbb(ind_bta), 1, 'ceil');
LumpCstModel.var_ub(ind_bbb_forward(find(bnta))) =  0;
bbb_drains=strcat('DM_', bbb_metnames);
bbb_drains=strrep(bbb_drains, '-', '_');
[~, bmetttt]=ismember(bbb_drains, LumpCstModel.rxns);
LumpCstModel.ub(bmetttt)=LumpCstModel.var_ub(ind_bbb_forward);

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_3.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < 

% LumpCstModel.b_for_lumping = roundsd(-muMax * stoich_bbb, 1, 'floor');

for i = 1:length(bbb_metnames)%parfor i = 1:length(bbb_metnames)
    bbb_metnames{i}
    
    if ~ismember(bbb_metnames{i}, bbb_not_to_lump)
        
        dmodel = LumpCstModel;
        dmodel.var_ub(ind_bbb_backward(i)) = 0;
        dmodel.var_ub(ind_bbb_forward(i))  = 10;
        dmodel.var_lb(ind_bbb_forward(i)) = roundsd(-muMax * stoich_bbb(i), 5, 'floor'); %negative because the coef is neg
        
        if strcmp(TimeLimitForSolver,'yes')
            sol = solveTFAmodelCplex(dmodel, TimeLimit, [], mipTolInt, emphPar, feasTol, scalPar, []);
        elseif strcmp(TimeLimitForSolver,'no')
            sol = solveTFAmodelCplex(dmodel, []       , [], mipTolInt, emphPar, feasTol, scalPar, []);
        else
            error('Wrong option!')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Remember here we're maximizing BFUSE vars -> minimizing non-core metab rxns
        % lumps that can't solve are saved in error and nosol
        if isempty(sol.x)
            rxnNames    = {'NA'};
            rxnFormulas = {'NA'};
        elseif isnan(sol.x) % don't we need one entry here for the if??
            rxnNames    = {'Crashed'};
            rxnFormulas = {'Crashed'};
        else
            bfuse = sol.x(allBFUSEind);
            rxnsBFUSE = find(bfuse<1e-9);%index of every bfuse var (reaction) that participates in the current lump
            if isempty(rxnsBFUSE) % lump is empty
                rxnNames    = {'ProducedByCore'};
                rxnFormulas = {'ProducedByCore'};
            else % there are reactions in the lump
                
                if NumAltLUMPedPerBBB>1
                    if strcmp(TimeLimitForSolver,'yes')
                        [DPs, ~, ~] = findDP_YIELD_TimeLimit(dmodel, NumAltLUMPedPerBBB, sol, allBFUSEind, TimeLimit, SizeLUMPNetwork, CplexParameters);
                    elseif strcmp(TimeLimitForSolver,'no')
                        [DPs, ~, ~] = findDP_YIELD_TimeLimit(dmodel, NumAltLUMPedPerBBB, sol, allBFUSEind, []       , SizeLUMPNetwork, CplexParameters);
                    else
                        error('Wrong option!')
                    end
                    DPsAll{i} = DPs;
                else
                    DPsAll{i} = sol.x;
                end
                
                rxnNames = dmodel.varNames(allBFUSEind(rxnsBFUSE));% for example rxns={'BFUSE_PDH' 'BFUSE_PPCK'}
                rxnNames = strrep(rxnNames, 'BFUSE_', '');   % for example rxns={'PDH'       'PPCK'}
                [~, ba] = ismember(rxnNames, dmodel.rxns); %ba is the reaction indeces in model dmodel
                rxnFormulas = printRxnFormula(dmodel, dmodel.rxns(ba), false); %rxnns is the formulas of the reactions in model dmodel
            end
        end
        rxnsNamesTotal{i}   = rxnNames;
        rxnFormulasTotal{i} = rxnFormulas;
    end
end

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_4.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <


%%%%%%%%%%%%%%%% Printout For Troubleshooting %%%%%%%%%%%%%%%%
summaryOfLUMPS = rxnFormulasTotal';
id_NA      = find(cellfun(@(x) isequal(x,{'NA'}), summaryOfLUMPS));
id_Crashed = find(cellfun(@(x) isequal(x,{'Crashed'}), summaryOfLUMPS));
id_ProducedByCore = find(cellfun(@(x) isequal(x,{'ProducedByCore'}),summaryOfLUMPS));

summaryOfLUMPS(id_NA) = repmat({'NA'},size(id_NA,1),1);
summaryOfLUMPS(id_Crashed) = repmat({'ProducedByCore'},size(id_Crashed,1),1);
summaryOfLUMPS(id_ProducedByCore) = repmat({'ProducedByCore'},size(id_ProducedByCore,1),1);

fprintf('Summary table about reaction sets for the LUMPED of each bbb :\n')
[bbb_metnames summaryOfLUMPS]

fprintf('It appears that the following mets are produced by the core, therefore there is no associated set of reactions to generate a lumped reaction:\n')
BBBsProducedByTheCore = bbb_metnames(id_ProducedByCore)

fprintf('It appears that the following mets cannot be produced by the core, therefore there is no associated set of reactions to generate a lumped reaction:\n')
BBBsNotProduced = bbb_metnames(id_NA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LumpedRxnFormulas = {};
bbbNames = [];
activeRxns = {};
relaxedDGoVarsValues_ForEveryLumpedRxn = {};
KeepTrackOfTransToBeCore = {};

for i = 1:size(DPsAll,2)%parfor i = 1:size(DPsAll,2)
    i
    if ~isempty(DPsAll{i})
        DPs_i = DPsAll{i};
        for j = 1:size(DPs_i, 2)
            % This function tests the thermodynamic feasibility of the
            % pathways obtained, and gives back only the feasible ones as
            % output. It also makes the coefficients of the LUMPED
            % reactions integers
            DPs_ij = DPs_i(:,j);
            [bbbName_i, id_ActRxnNames_i, ActiveRxnsFormulas, LumpedRxnFormulas_i, relaxedDGoVarsValues_ForEveryLumpedRxn_i, KeepTrackOfTransToBeCore_i] = ...
                testForIntegerAndThermo(LumpCstModel, GSM_ForLumping, DB_AlbertyUpdate, allBFUSEind, DPs_ij, muMax, IdNCNTNER, i, ...
                otherReactionsGSMForLump_idx, CplexParameters, bbb_metnames, RxnNames_PrevThermRelax, addGAM, ImposeThermodynamics);
            KeepTrackOfTransToBeCore = [KeepTrackOfTransToBeCore; KeepTrackOfTransToBeCore_i];
            bbbNames = [bbbNames; bbbName_i];
            relaxedDGoVarsValues_ForEveryLumpedRxn = [relaxedDGoVarsValues_ForEveryLumpedRxn {relaxedDGoVarsValues_ForEveryLumpedRxn_i}];
            activeRxns = [activeRxns; [{id_ActRxnNames_i} {ActiveRxnsFormulas}]];
            LumpedRxnFormulas = [LumpedRxnFormulas; LumpedRxnFormulas_i];
        end
    end
end

if ~isequal(unique(KeepTrackOfTransToBeCore), {'NA'})
KeepTrackOfTransToBeCore_unique=setdiff(KeepTrackOfTransToBeCore, {'NA'});
[~, ind_nc_tpt_to_be_core]=ismember(KeepTrackOfTransToBeCore_unique, GSM_ForLumping.rxns);
IdNCNTNER=setdiff(IdNCNTNER, ind_nc_tpt_to_be_core);
end

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_5.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

% Remove the Lumped reactions that were not feasible
idToRemove = find(strcmp(LumpedRxnFormulas,'NA'));
LumpedRxnFormulas(idToRemove) = [];
bbbNames(idToRemove) = [];
activeRxns(idToRemove,:) = [];
relaxedDGoVarsValues_ForEveryLumpedRxn(idToRemove)=[];

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_6.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

end
