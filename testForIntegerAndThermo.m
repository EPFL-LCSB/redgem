function [bbbName_i, id_ActRxnNames_i, ActiveRxnsFormulas_i, LumpedRxnFormulas_i, relaxedDGoVarsValues_i, KeepTrackOfTransToBeCore] = ...
    testForIntegerAndThermo(LumpCstModel, GSM_ForLumping, DB_AlbertyUpdate, allBFUSEind, DPs_ij, muMax, IdNCNTNER, i, ...
    otherReactionsGSMForLump_idx, CplexParameters, bbb_metnames, RxnNames_PrevThermRelax, addGAM, ImposeThermodynamics)
% INPUT
% - LumpCstModel:
% - GSM_ForLumping:
% - DB_AlbertyUpdate:
% - allBFUSEind: Indices of binary variable that are used to minimize the
%   number of reactions in the lumped reaction.
% - DPs_i:
% - muMax:
% - IdNCNTNER:
% - i:
% - otherReactionsGSMForLump_idx: Indices of non-core reactions. These indices
%   are according to GSM_ForLumping indexing.
% - CplexParameters
%
% OUTPUT
% - LumpedRxnFormulas_i:
% - id_ActRxnNames_i:
% - ActiveRxnsFormulas_i:
% - relaxedDGoVarsValues_i:
 
changeCobraSolver('cplex_direct')
 

flag=0;
% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);
 
% Get the non-core reactions that are active for this (j-th) solution of the current (i-th) BBB
% - Find the id of the binary variables that are active
idActBVars_in_allBFUSEind = find(DPs_ij(allBFUSEind)<0.5);
% - Get the names of the binary variables that are active
ActBVarNames = LumpCstModel.varNames(allBFUSEind(idActBVars_in_allBFUSEind));
% - Get the names of the corresponding reactions that are active
ActRxnNames = strrep(ActBVarNames, 'BFUSE_', '');
% - Get the index of the active reactions in LumpCstModel
[~, id_ActRxnNames_i] = ismember(ActRxnNames, LumpCstModel.rxns);

% Find the reactions of the lumped-subnetwork that are exchange transports
% (e <-> something)

LumpedSubnetworkModel = extractSubNetwork(LumpCstModel, ActRxnNames);
LumpedSubnetworkModelWithexchangeTrans = isTrans_GEMS_MoreFields(LumpedSubnetworkModel, LumpedSubnetworkModel.mets);

% If any of the active reactions is a drain make its index zero
drains = extract_drains(GSM_ForLumping);
for k = 1:length(ActRxnNames(find(id_ActRxnNames_i)))
    if ismember(ActRxnNames{k}, drains) 
        id_ActRxnNames_i(k) = 0;
    end
end
 
% Initialize the model
kmodel = LumpCstModel;
 
 
% Extract the metabolites of the core model. To do this we find the complement of the non-core
% reactions in the kmodel (LumpCstModel).
% - Find the reactions in the kmodel (LumpCstModel) that exist in the
%   non-core reactions of the GSM_ForLumping.
[~, id_nonCoreInKmodel] = ismember(1:length(kmodel.rxns), otherReactionsGSMForLump_idx);
% - Extract the core model
kmodelCore = extractSubNetwork(kmodel, kmodel.rxns(find(id_nonCoreInKmodel==0)));
% We add to the metabolites of the core-model all the bbbs
mets_core = [kmodelCore.mets; bbb_metnames];
% Add information about transports and drains to the kmodel
kmodel = isTrans_GEMS_MoreFields(kmodel, mets_core);
kmodel = checkifDrain(kmodel);
 
% Create a mini model that includes only:
% (1) The active reactions of this (j-th) solution of the current (i-th) BBB
% (2) Drain reactions for each BBB
% (3) core, transports or exchange reactions
%
% Get the reactions in kmodel that are the complement of IdNCNTNER, so either
% core, transport, or exchange reactions
[~, IdNCNTNER_inKmodel] = ismember(kmodel.rxns(1:length(LumpCstModel.rxns)), kmodel.rxns(IdNCNTNER));
rxns_coreTransOrExch = kmodel.rxns(find(IdNCNTNER_inKmodel==0));
% Create the abovementioned mini model, which consists of (1), (2) & (3)
mini_model = extractSubNetwork(kmodel, [ActRxnNames; rxns_coreTransOrExch; strcat('DM_', bbb_metnames)]);
 
% Add the Growth Associated Maintenance (GAM) reaction to this mini model:
if strcmp(addGAM, 'yes')
    GAM_c = 'GAM_c';
    % Extract the Growth Associated Maintenance (GAM) reaction to add it
    % later to the mini model:
    [GAMequation, ppi_equation] = extractGAMcoefficientsFromBiomass(LumpCstModel);
    bbb_metnames_and_GAM_c = [bbb_metnames; GAM_c];
    mini_model = addReaction(mini_model, 'DM_GAM_c', GAMequation);
    mini_model.lb(length(mini_model.rxns)) = 0;
    if strcmp('GAM_c',bbb_metnames_and_GAM_c{i})
        mini_model.ub(length(mini_model.rxns)) = 100;
    else
        mini_model.ub(length(mini_model.rxns)) = 0;
    end
else
    bbb_metnames_and_GAM_c = bbb_metnames;
	GAM_c = [];
    GAMequation = [];
end

if strcmp(ImposeThermodynamics, 'yes')
    % Assure that the mini model will have thermodynamic information in the same way as the GSM_ForLumping.
    LumpCstModel.rxnThermo(1:length(GSM_ForLumping.rxns)) = GSM_ForLumping.rxnThermo;
    [~, idRxnsInLumpCstModel] = ismember(mini_model.rxns, LumpCstModel.rxns);
    mini_model.rxnThermo(find(idRxnsInLumpCstModel)) = LumpCstModel.rxnThermo(idRxnsInLumpCstModel(find(idRxnsInLumpCstModel)));
    mini_model.thermo_units = 'kcal/mol';
end
 
% The function above puts thermo only to the core. The one below puts thermo to the entire system.
mini_model.c = zeros(size(mini_model.c));
mini_model.c(find_cell(['DM_',regexprep(bbb_metnames_and_GAM_c{i},'-','_')], mini_model.rxns)) = 1;
mini_model.ub(find(mini_model.c)) = 10;
FBAsolMiniModel = solveFBAmodelCplex(mini_model);
FBAsolMiniModel = roundsd(0.99*FBAsolMiniModel.f, 4, 'floor');
if FBAsolMiniModel < 10^-6
    disp(['Why is the maximum attainable FBA-flux for the drain of the bbb %s so small (<10^-6)?!?!', bbb_metnames_and_GAM_c{i}])
end
 

if strcmp(ImposeThermodynamics, 'yes')
    verboseFlag = false;
    % When we converted the model to TFA (at the top of addBIOMASS.m) in order
    % to find the maximum yield, we had to relax the thermodynamics of the
    % following reactions:
    if ~isempty(LumpCstModel.relaxedDGoVarsValues)
        relaxedDGoOfRxns = regexprep(LumpCstModel.relaxedDGoVarsValues(:,1),'^DGo_','');
    else
        relaxedDGoOfRxns = [];
    end
    % Now we convert again the model to TFA, without adding thermodynamic
    % constraints for these reactions:
    [tmini_model, relaxedDGoVarsValues_i] = convToTFA(mini_model, DB_AlbertyUpdate, relaxedDGoOfRxns, 'NO', RxnNames_PrevThermRelax, FBAsolMiniModel, [], [], verboseFlag);
elseif strcmp(ImposeThermodynamics, 'no')
    fprintf('In the current case we do not impose thermodynamic constraints\n')
    tmini_model = convFBA2MILP(mini_model, DB_AlbertyUpdate);
    relaxedDGoVarsValues_i = [];
else
    error('Wrong option for ImposeThermodynamics')
end
 
bbb_drains_forward = strcat('F_DM_', [bbb_metnames; GAM_c]);
bbb_drains_forward = strrep(bbb_drains_forward, '-', '_');
[~, ind_bbb_forward] = ismember(bbb_drains_forward, tmini_model.varNames);
 
% Extract the stoichiometric coefficients for each bbb including as last
% the GAM_c. ATTENTION!! Restore the objective to the Growth reaction!!
tmini_model.c = zeros(size(tmini_model.c));
tmini_model = changeObjective(tmini_model, GSM_ForLumping.rxns(find(GSM_ForLumping.c == 1)));
tmini_model.S = full(tmini_model.S);
[~, id_bbbs] = ismember(bbb_metnames, tmini_model.mets);
id_bbbs = id_bbbs(find(id_bbbs));
stoich_bbb = tmini_model.S(id_bbbs, find(tmini_model.c));
stoich_bbb = [stoich_bbb; -1];
 
% Set the upper bound of the forward drain flux of the current bbb to the highest alowed value
tmini_model.var_ub(ind_bbb_forward(i)) = 10;
% Set the lower bound of the forward drain flux of the current bbb to a value equal to its stoichiometric coefficient (times a muMax, but this is 1 without problems so far).
value = muMax * (-1*stoich_bbb(i));
% If this value is very small (i.e. <10^-6) do not use that value, but 5*10^-5
if value > 10^-8
    tmini_model.var_lb(ind_bbb_forward(i)) = roundsd(muMax*(-1*stoich_bbb(i)), 6, 'floor');
else
    warning('The stoichiometric coefficient of the %s bbb is <10^-6!', bbb_metnames{i})
    tmini_model.var_lb(ind_bbb_forward(i)) = 5*10^-7;
end
 
% Here we form an optimization problem where we assume that the forward and
% backward fluxes (of the active reactions of this j-th solution of the
% current i-th BBB) are integers, and we try to minimize their sum.
% Therefore we set up the objective of the mini model as the forward and
% backward fluxes of these active reactions.

% Remove the exchange transport from the active reactions that compose the
% lumped
if any(LumpedSubnetworkModelWithexchangeTrans.isExchangeTrans == 1)
    [~, ind_nc_tpt] = ismember(LumpedSubnetworkModelWithexchangeTrans.rxns(find(LumpedSubnetworkModelWithexchangeTrans.isExchangeTrans == 1)), ActRxnNames);    
    ActRxnNames(ind_nc_tpt)=[];
    KeepTrackOfTransToBeCore = LumpedSubnetworkModelWithexchangeTrans.rxns(find(LumpedSubnetworkModelWithexchangeTrans.isExchangeTrans == 1));
else
    KeepTrackOfTransToBeCore = {'NA'};
end

tmini_model.f = zeros(length(tmini_model.f), 1);
forward  = strcat('F_', ActRxnNames);
backward = strcat('R_', ActRxnNames);
[~, id_FB] = ismember([forward; backward], tmini_model.varNames);
tmini_model.f(id_FB) = 1;
tmini_model.objtype = 1;

% minimization
[~, id_F] = ismember(forward,  tmini_model.varNames);
[~, id_B] = ismember(backward, tmini_model.varNames);
 
% Exchange the upper bounds of all forward and reverse fluxes (that
% satisfy: 10^-6 < (F_x or R_x) < 1001) with 10^7
indfr=getAllVar(tmini_model, {'F' 'R'});
tmini_model_orig=tmini_model;

% The hope is that the problem will be feasible with integer coefficients:
tmini_model.var_ub(indfr)=1000000*tmini_model.var_ub(indfr);
tmini_model.var_lb(indfr)=1000000*tmini_model.var_lb(indfr);

tmini_model.vartypes(id_FB) = {'I'};
sol = solveTFAmodelCplex(tmini_model, 300, [], mipTolInt, emphPar, feasTol, scalPar, []);
 
% However, sometimes this is not possible (e.g. reactions with rational coefficients involved)
% And then we switch back to continuous variables to get a solution.
if isempty(sol.x)
    disp('Integer solution was not found!! You can try increasing the 10^6 in line 188-189 to higher values...\n')
    tmini_model=tmini_model_orig;
    sol = solveTFAmodelCplex(tmini_model, 300, [], mipTolInt, emphPar, feasTol, scalPar, []);
    flag=1;
end
 
if isempty(sol.x)
    LumpedRxnFormulas_i{1, 1} = 'NA';
else
    % Get the net flux of each reaction
    NetFlux_ofActRxns = sol.x(id_F) - sol.x(id_B);
    
    tmini_model.S = full(tmini_model.S);
    stoich = [];
    
    for j = 1:length(id_F)
        if length(id_F)==1
            [~, id_ActRxn_j] = ismember(ActRxnNames, tmini_model.rxns);
        else
            [~, id_ActRxn_j] = ismember(ActRxnNames(j), tmini_model.rxns);
        end
        % Get the stoichiometry column for each of the active reactions j
        % and multiply it with the (hopefully integer) minimal net-flux of
        % that reaction
        sto = tmini_model.S(:, id_ActRxn_j).*NetFlux_ofActRxns(j); % CAN'T WE JUST REPLACE THIS m WITH j?
        % We append and collect these stoichiometries in a matrix
        stoich = [stoich sto];
    end
    
    % Sum the stoichiometries of all active reactions to create a
    % stoichiometry of a lumped reaction
    sumStoich = sum(stoich,2);
    % Check if there exist stoichiometric coefficients that are very low
    % (i.e. <10^-8) in this lumped stoichiometric reaction.
    id_LowStCoef = find(abs(sumStoich)<10^-9);
    % If there exist, remove them (normally for the integer-case this should not happen. If it happens we need to check cplex settings)
    sumStoich(id_LowStCoef) = 0;
    % Check if there exist any non-zero stoichiometric coefficients
    id_NonZeroStCoef = find(sumStoich);
    if isempty(id_NonZeroStCoef)
        LumpedRxnFormulas_i{1, 1} = 'NA';
    else
        % Normalize wrt the gcd of stoichiometric coefficient to avoid
        % large numbers
        sumStoich_cc=sumStoich(find(sumStoich));
        if flag==0
            gcd_stoich = gcd_vect(round(sumStoich_cc));
            sumStoich = round(sumStoich)/gcd_stoich;
        else
            gcd_stoich=min(abs(sumStoich_cc));
            sumStoich = sumStoich/gcd_stoich;
        end
        % Create a dummy model which will be the same as the tmini_model,
        % but will also include the Lumped reaction j:
        [dModel, rxnIDexists] = addReaction(tmini_model, ['LMPD_', bbb_metnames{i}], tmini_model.mets(find(sumStoich)), sumStoich(find(sumStoich)));
        % Save the RxnFormula whether the reaction stoichiometry already
        % exists in the model or not:
        if ~isempty(rxnIDexists)
            LumpedRxnFormulas_i(1,1) = printRxnFormula(dModel, dModel.rxns(rxnIDexists));
        else
            LumpedRxnFormulas_i(1,1) = printRxnFormula(dModel, dModel.rxns(length(dModel.rxns)));
        end
    end
end
 
[~, id_ActRxnNames_i] = ismember(ActRxnNames, LumpCstModel.rxns);
 
if ~isempty(id_ActRxnNames_i)
    ActiveRxnsFormulas_i = printRxnFormula(LumpCstModel, ActRxnNames);
else
    ActiveRxnsFormulas_i = [];
end
 
bbbName_i = bbb_metnames_and_GAM_c(i);
 
end
