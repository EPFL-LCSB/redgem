function [rxns_all, id_all, sol_all] = redGEMX(rxns_ss, GSM_ForLumping, GSM_ForAdjMat, OriginalGEM, GEMmodel, Organism, GEMname, NumOfConnections, CplexParameters, ReactionDB, ImposeThermodynamics,output_PATH)
    
if strcmp(ImposeThermodynamics, 'yes')
    GSM_ForLumping = prepModelforTFA(GSM_ForLumping, ReactionDB, GSM_ForLumping.CompartmentData);
    GSM_ForLumping = convToTFA(GSM_ForLumping, ReactionDB);
elseif strcmp(ImposeThermodynamics, 'no')
    GSM_ForLumping = convFBA2MILP(GSM_ForLumping, ReactionDB);
else
    error('Wrong option for ImposeThermodynamics')
end

if strcmp(NumOfConnections,'OnePerMetE')
    %flag: for each lumped reaction, how many alternative lumped reactions should the solver look for?
    NumAltPerMetE = 1;
    SizeLUMPNetwork = 0.5;
elseif strcmp(NumOfConnections,'SminMetE')
    NumAltPerMetE = 10^4;
    SizeLUMPNetwork = 0.5;
end



[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

ExtraCellSubsystem = GSM_ForLumping.ExtracellularMedium;

model_for_extra = GSM_ForLumping;

model_for_extra.rev=ones(length(model_for_extra.rxns),1);
model_for_extra=createintegermodel(model_for_extra);
[~, bfor]=ismember(strcat('F_', ExtraCellSubsystem), model_for_extra.varNames);
bfor = bfor(find(bfor));
[~, bback]=ismember(strcat('R_', ExtraCellSubsystem), model_for_extra.varNames);
bback = bback(find(bback));
model_for_extra.var_ub(find(model_for_extra.f))=0;
d0_rxns=GSM_ForAdjMat.rxns(vertcat(rxns_ss{:}));
other_rxns=setdiff(GSM_ForLumping.rxns, d0_rxns);
model_for_extra.f=zeros(length(model_for_extra.f),1);

model_for_extra.var_lb(bfor) = 0;
a = model_for_extra.var_ub(bfor);
a(a>0 & a<10) = 100;
model_for_extra.var_ub(bfor) = a;
model_for_extra.var_lb(bback) = 0;
b = model_for_extra.var_ub(bback);
b(b>0 & b<10) = 100;
model_for_extra.var_ub(bback) = b;

%% Block the extracellular reactions
ExtracellRxns = {};
for i = 1:length(model_for_extra.rxns)
    rMets = model_for_extra.mets(find(model_for_extra.S(:,i)));
    [~,ii] = ismember(rMets,model_for_extra.mets);
    if isempty(setdiff(unique(model_for_extra.metCompSymbol(ii)),'e')) && length(rMets)~=1
        ExtracellRxns = [ExtracellRxns;model_for_extra.rxns(i)];
    end
end

[~,FeRxns] = ismember(strcat('F_', ExtracellRxns),model_for_extra.varNames);
[~,ReRxns] = ismember(strcat('R_', ExtracellRxns),model_for_extra.varNames);

model_for_extra.var_ub(FeRxns) = 0;
model_for_extra.var_ub(ReRxns) = 0;

%% Block the 'drains' for non-core mets
% 1. Block all drains
drains = extract_drains(model_for_extra);

[~,Fdrains] = ismember(strcat('F_', drains),model_for_extra.varNames);
[~,Rdrains] = ismember(strcat('R_', drains),model_for_extra.varNames);

model_for_extra.var_ub(Fdrains) = 0;
model_for_extra.var_ub(Rdrains) = 0;

% 2. Open drains for core mets
coreNetwork = extractSubNetwork(GEMmodel, GSM_ForAdjMat.rxns(vertcat(rxns_ss{:})));
coreMets = coreNetwork.mets;
[~,ind_coreMets] = ismember(coreMets,model_for_extra.mets);
for i=1:length(coreNetwork.mets)
    coreMets_noComp{i}=coreNetwork.mets{i}(1:end-2);
end
coreMets_noComp = unique(coreMets_noComp');
coreMets_e = strcat(coreMets_noComp,'_e');

[~,ind_coreMets_e] = ismember(coreMets_e,model_for_extra.mets);
ind_coreMets_e = ind_coreMets_e(find(ind_coreMets_e));

checkMets = [ind_coreMets;ind_coreMets_e];
checkMets = unique(checkMets);

coreMets_drains={};

for j = 1:length(checkMets)
    rxnsID = find(model_for_extra.S(checkMets(j), :));
    for k = 1:length(rxnsID)
        no_of_mets = length(find(model_for_extra.S(:,rxnsID(k))));
        if no_of_mets == 1
            coreMets_drains = [coreMets_drains; model_for_extra.rxns(rxnsID(k))];
        end
    end
end

[~,Fcore_drains] = ismember(strcat('F_', coreMets_drains),model_for_extra.varNames);
[~,Rcore_drains] = ismember(strcat('R_', coreMets_drains),model_for_extra.varNames);

model_for_extra.var_ub(Fcore_drains) = OriginalGEM.var_ub(Fcore_drains);
model_for_extra.var_ub(Rcore_drains) = OriginalGEM.var_ub(Rcore_drains);

% 3. Open drains for medium metabolites
model_for_extra.var_ub(bfor) = OriginalGEM.var_ub(bfor);
model_for_extra.var_ub(bback) = OriginalGEM.var_ub(bback);


%%
[~, fu_rxn]=ismember(strcat('FU_',other_rxns), model_for_extra.varNames); 
[~, bu_rxn]=ismember(strcat('BU_',other_rxns), model_for_extra.varNames); 
 
[~, num_vars] = size(model_for_extra.A);
 
for i=1:length(fu_rxn)
    model_for_extra.varNames(length(model_for_extra.varNames)+1) = strcat('BFUSE_', other_rxns(i));
    model_for_extra.var_ub(length(model_for_extra.varNames)) = 1;
    model_for_extra.var_lb(length(model_for_extra.varNames)) = 0;
    model_for_extra.vartypes(length(model_for_extra.varNames)) = {'B'};
    model_for_extra.f(length(model_for_extra.varNames)) = 1; 
end


ind_bfuse=getAllVar(model_for_extra, {'BFUSE'});
[num_constr,~] = size(model_for_extra.A);

% constraint FUse+BUse+BFUSE <= 1; If BFUSE = 1 then the corresponding
% reaction cannot carry flux.
for i=1:length(fu_rxn)
    model_for_extra.rhs(num_constr+i,1) =  1;
    model_for_extra.constraintNames{num_constr+i,1} = strcat('Constraint_BFUSE_',num2str(i));
    model_for_extra.constraintType{num_constr+i,1} = '<';
    model_for_extra.A(num_constr+i,ind_bfuse(i)) = 1;
    model_for_extra.A(num_constr+i, fu_rxn(i)) = 1;
    model_for_extra.A(num_constr+i, bu_rxn(i)) = 1;
end

id_all = [];
ExtraCellSubsystem_connect = GSM_ForLumping.ExtracellularMedium_connect;
[~, bfor]=ismember(strcat('F_', ExtraCellSubsystem_connect), model_for_extra.varNames);
[~, bback]=ismember(strcat('R_', ExtraCellSubsystem_connect), model_for_extra.varNames);
DPsAll = {};

for i=1:length(ExtraCellSubsystem_connect)
    [i ExtraCellSubsystem_connect(i)]
    model_for_extra_temp=model_for_extra;

    mm_F = runTMinMax(model_for_extra_temp,model_for_extra_temp.varNames(bfor(i)),300);
    mm_R = runTMinMax(model_for_extra_temp,model_for_extra_temp.varNames(bback(i)),300);

    if min(mm_F(2),mm_R(2)) ~= 0 && min(mm_F(2),mm_R(2)) < 3
        rhs_i = min(mm_F(2),mm_R(2));       
    elseif min(mm_F(2),mm_R(2)) == 0 && max(mm_F(2),mm_R(2)) ~= 0 && max(mm_F(2),mm_R(2)) < 3
        rhs_i = max(mm_F(2),mm_R(2));
    elseif min(mm_F(2),mm_R(2)) > 3 || max(mm_F(2),mm_R(2)) > 3
        rhs_i = 1e-3;
    end 
     
    [~, f_rxn]=ismember(strcat('F_',other_rxns), model_for_extra_temp.varNames); 
    [~, r_rxn]=ismember(strcat('R_',other_rxns), model_for_extra_temp.varNames); 

    [num_constr,~] = size(model_for_extra_temp.A);

    % constraint F+R+BFUSE >= rhs_i; If BFUSE=0 then the corresponding reaction
    % has to carry a basal flux (to avoid suboptimal networks).
    for j = 1:length(f_rxn)
        model_for_extra_temp.rhs(num_constr+j,1) =  1e-6;
        model_for_extra_temp.constraintNames{num_constr+j,1} = strcat('BFMER_',num2str(j));
        model_for_extra_temp.constraintType{num_constr+j,1} = '>';
        model_for_extra_temp.A(num_constr+j,ind_bfuse(j)) = 1;
        model_for_extra_temp.A(num_constr+j, f_rxn(j)) = 1;
        model_for_extra_temp.A(num_constr+j, r_rxn(j)) = 1;
    end 
 
    % constraint F(ex_met)+R(ex_met) >= rhs_i; force flux of specific
    % metabolite
    [num_constr,num_vars] = size(model_for_extra_temp.A);
    model_for_extra_temp.rhs(num_constr+1,1) =  0.95*rhs_i;
    model_for_extra_temp.constraintNames{num_constr+1,1} = strcat('Path_',ExtraCellSubsystem_connect{i});   
    model_for_extra_temp.constraintType{num_constr+1,1} = '>';
    model_for_extra_temp.A(num_constr+1,bfor(i)) = 1;
    model_for_extra_temp.A(num_constr+1,bback(i)) = 1;

    sol = solveTFAmodelCplex(model_for_extra_temp, 300, [], mipTolInt, emphPar, feasTol, scalPar, []);
    
    if sol.val==0
        rxns_all{i,1} = 0;
        sol_all{i,1}  = 0;
    elseif isempty(sol.x)
        rxns_all{i,1} = 0;
        sol_all{i,1}  = 0;
    else
        sol_all{i,1} = sol.x;
        if NumAltPerMetE > 1
            [DPs, ~, ~] = findDP_YIELD_TimeLimit(model_for_extra_temp, NumAltPerMetE, sol, ind_bfuse, 300, SizeLUMPNetwork, CplexParameters);
            DPsAll{i} = DPs;
            save('RedGEMX_DPs.mat','DPsAll','sol_all')
        else
            DPsAll{i} = sol.x;
        end      
    end
end

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > 
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                                  %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_1.mat;'])  %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

for i = 1:length(DPsAll)
    if ~isempty(DPsAll{i})
        for j = 1:size(DPsAll{i},2) 
                this_dp = DPsAll{i}(:,j);
                id_rxns = find(this_dp(ind_bfuse)<0.1);
                rxnNames = other_rxns(id_rxns);
                rxns_all{i}{j} = rxnNames;
        end
    else
            rxns_all{i} = {};
    end
end

id_all_dps_min = [];
for i = 1:length(rxns_all)
    size_subntwk = [];
    for j = 1:length(rxns_all{i})
        size_subntwk = [size_subntwk;length(rxns_all{i}{j})];
    end
    min_size = min(size_subntwk);
    jj = find(size_subntwk==min_size);
    
    for k = 1:length(jj)
        rxns_all_min{i}{k} = rxns_all{i}{jj(k)};
        [~,id_rxns_inGSM] = ismember(rxns_all{i}{jj(k)},GSM_ForLumping.rxns);
        id_all_dps_min = [id_all_dps_min;id_rxns_inGSM];
    end   
end
rxns_all = rxns_all_min;
id_all = id_all_dps_min;

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                                  %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_2.mat;'])  %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
% 
end
