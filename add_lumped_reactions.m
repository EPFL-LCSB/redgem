function  [model_red, rxns ,info_LMDP, missingMets] = ...
    add_lumped_reactions(GSM_ForLumping, LumpedRxnFormulas,bbbNames, DB_AlbertyUpdate, ...
    all_core_rxns, ImposeThermodynamics, CplexParameters, Organism, GEMname, activeRxns, output_PATH)
% This function incorporates all the lumped reactions into a reduced model
%
% INPUTS
% - GSM_ForLumping: Genomescale without biomass. If the option rmvComp
%   is 1, then the compartments are merged, and this model has also no
%   periplasm
% - LumpCstModel:
% - GSM_ForAdjMat: Genomescale without biomass, and without the
%   metabolite pairs. if the option rmvComp is 1, then the compartments
%   are merged, and this model has also no periplasm.
% - overall_am_full:
% - rxns_core_subsystems_full:
% - connecting_reactions:
% - DB_AlbertyUpdate: Alberty thermodynamic database
%
% OUTPUTS
%
% - model_red:
% - rxns:
% - info_LMDP: [{LMDPname } {SubNetworkReactions} {SbNetworkFormulas}]

% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                        %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_1.mat;'])   %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < 

% Make the formulas of all the lumped reactions uni-directional
LumpedRxnFormulas = strrep(LumpedRxnFormulas, '<=>', '->');
% Extract from the GSM_ForLumping the core reactions
model_red = extractSubNetwork(GSM_ForLumping, all_core_rxns);

% Lumped reactions naming and tracking:
% LMPD_iBBB_iBBBnumberOfLumps_sizeOfLump

str_allLMPDs={};
cell_SubNetworks_LMPD={};
BBBs=unique(bbbNames,'stable'); 
for i_BBB=1:numel(BBBs)  
    id_LMPD=find_cell(BBBs(i_BBB),bbbNames);
    for i_BBBlumps=1:numel(id_LMPD)
        lump_SubNetwork=GSM_ForLumping.rxns(activeRxns{id_LMPD(i_BBBlumps),1});
        str_LMPD=strcat('LMPD_',BBBs(i_BBB),'_',num2str(i_BBBlumps),'_',num2str(numel(lump_SubNetwork)));
        str_allLMPDs=[str_allLMPDs;str_LMPD];
        cell_SubNetworks_LMPD=[cell_SubNetworks_LMPD ; {lump_SubNetwork}];
    end 
end
% Keep this information for future redGEM models so we have information for
% analysis on lumped reactions
if ~ isempty(bbbNames)
    info_LMDP=[str_allLMPDs bbbNames LumpedRxnFormulas cell_SubNetworks_LMPD activeRxns(:,2)];
else
    info_LMDP = [];
end

for i=1:length(LumpedRxnFormulas)
    model_red = addReaction(model_red, str_allLMPDs{i}, LumpedRxnFormulas{i});
    model_red.lb(length(model_red.rxns)) = 0; % have to be produced
end

sol = solveFBAmodelCplex(model_red, scalPar, feasTol, emphPar);
if sol.f > 1e-4
    rxns = {};
    missingMets = {};
else
    [missingMets, presentMets] = biomassPrecursorCheck(model_red);
    % If the model is not able to produce biomass, then we use the entire
    % GS-model, together with the core, connecting, and exchange reactions,
    % as well as the lumped reactions to get the missing reactions
    
    % ttmodel is the whold GEM used for lumping
    ttmodel=GSM_ForLumping;
    % Here we are adding the lumped reactions to the whole genome scale
    % model:
    for i=1:length(LumpedRxnFormulas)
        ttmodel = addReaction(ttmodel, str_allLMPDs{i}, LumpedRxnFormulas{i});
        ttmodel.lb(length(ttmodel.rxns)) = 0; % have to be produced
    end
    % Find indices of core, connecting, and exchange reactions in the
    % Genome-scale PLUS the lumped reactions
    GSM_ForLumping_idDrains = find(sum(abs(full(GSM_ForLumping.S)),1)==1);
    GSM_ForLumping_idDrains = [GSM_ForLumping_idDrains find(GSM_ForLumping.c)];
    [~, ba] = ismember(GSM_ForLumping.rxns, all_core_rxns);
    % Find indices of non-core, non-connecting, non-exchange reactions in
    % the Genome-scale PLUS the lumped reactions
    Ind_OtherReacs = find(ba==0);
    % This converts the model to irreversible-form, and adds the TFA
    % fields for constraints. ATTENTION!! It does NOT convert to TFA in
    % the classical way (DG's,etc.)
    tttmodel = convToTFBA_BNICE(ttmodel, DB_AlbertyUpdate);
    [~, rxns, ~, ~] = addBIOMASSminRea(tttmodel, Ind_OtherReacs, DB_AlbertyUpdate, 0.9, ImposeThermodynamics, CplexParameters);
    model_red = extractSubNetwork(GSM_ForLumping, [rxns(:,1);all_core_rxns]);
    for i = 1:length(LumpedRxnFormulas)
        model_red = addReaction(model_red, str_allLMPDs{i}, LumpedRxnFormulas{i});
        model_red.lb(length(model_red.rxns))=0;
    end
end


model_red.info_LMDP=info_LMDP;

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                        %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_2.mat;'])   %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

end
