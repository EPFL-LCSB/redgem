function [selectedPaths, selectedPathsExternal, connecting_reactions, rxns_ss, mets_ss, otherReactions, ThrowErrorOnDViolation] = extractAdjData(model,core_ss,L_DirAdjMatWn,L_DirAdjMatC,L_DirAdjMatMets,D,startFromMin,OnlyConnectExclusiveMets,ConnectIntraSubsystem,ApplyShortestDistanceOfSubsystems,ThrowErrorOnDViolation)
% This function
%
% INPUTS
% - model: metabolic model of COBRAtoolbox format. (in the redGEM case,
%   this is the model that was used for the calculation of the Adj Mat)
% - core_ss: selected core subsystems (defined in redGEM_15*)
% - L_DirAdjMatWn: This is like DirAdjMatWn, but it refers to connections of
%   degree L.
% - L_DirAdjMatC: This is like DirAdjMatC, but it refers to connections of
%   degree L.
% - L_DirAdjMatMets: This is like DirAdjMatC, but it explicitly lists all
%   the participating metabolites for each connection.
% - D: (Defined in redGEM_150112)
%   Given that the algorithm has searched paths of up to length L, D
%   specifies how many lengths to include in the new model. if not
%   specified, this is equal to L.
% - startFromMin: (Defined in redGEM_150112)
%   default =0. if 1, indicates that the D's that are
%   included in the model start at the shortest distance between subsystem
%   pairs. if zero, the model is built with the first L=D levels regardless
%   of the presence of paths.
% - OnlyConnectExclusiveMets: if this is 1, than only metabolites that are
%   exlusive to each subsystem will be connected. if 0, even if a
%   metabolite is shared, it will still be connected. default = 1
% - ConnectIntraSubsystem: if this is 1, than metabolites that are
%   within a subsystem will be connected to eachother. if 0, they will not
%   be connected. default = 0;
%   ----------------------------------------------------------------------
%   Example: L=10 D=4 startFromMin=1
%   the adjacency matrix will be calculated up to length 10, but only the
%   first 4 lengths that connect subsystems will be included. If subsystem
%   A and B have a L=1 connection, then the paths of length 1:4 will be
%   included in the model. However, the shortest distance between
%   subsystems A and C is 3 so all paths of length 3:7 will be included.
%   ----------------------------------------------------------------------
% - stats: (Defined in redGEM_150112)
%   options to plot statistics on the connecting pathways. default=0
%
% OUTPUTS
% - selectedPaths: This is a structure containing the information about the path
%   from subsystem to subsystem, with fields:
%   .sub2subNumP: Number of pathways from source subsystem to target
%     subsystem
%   .sub2subNumR: Total number of reactions from source subsystem to
%     target subsystem
%   .joiningReacs: The actual reactions
%   .joiningPaths: The actual pathways
%   .joiningMets: The actual metabolites
%   .ModelReacs: All reactions per level (D)
%   .ModelMets: All metabolites per level (D)
% - connecting_reactions:
% - rxns_ss:
% - mets_ss:
% - otherReactions: Indexes of all the non-D1 (in general: non-core)
%   connections (rest-stuff/unused)


% The following will always exist
% %% flags
% if ~exist('startFromMin', 'var') || isempty(startFromMin)
%     startFromMin = 'yes';% flag. if 1, indicates that we should start counting levels starting at the shortest distance between subsystem pairs. if zero, it returns first D levels
% end
%
% if ~exist('D', 'var') || isempty(D)
%     D=1;% flag. it tells the algorithm to extract info from the first D levels from the input matrices
% end
%
% if ~exist('OnlyConnectExclusiveMets', 'var') || isempty(OnlyConnectExclusiveMets)
%     OnlyConnectExclusiveMets = 1;
% end
%
% if ~exist('ConnectIntraSubsystem', 'var') || isempty(ConnectIntraSubsystem)
%     ConnectIntraSubsystem = 0;
% end



%% Find the metabolites for each starting subsystem
% Initialize cell structures
rxns_ss = cell(length(core_ss),1);
mets_ss = cell(length(core_ss),1);
% For each of the subsystems:
for i=1:length(core_ss)
    % - find the reactions of this subsystem, and their indices
    [~, rxns_ss{i}]=ismember(model.subSystems, core_ss{i});
    rxns_ss{i}=find(rxns_ss{i});
    % - find the metabolite indices of those metabolites involved in the
    %   reactions of each subsystems
    mets_ss{i}=find(any(model.S(:,rxns_ss{i}),2));
end

% if connecting the extracellular mets, do this before defining pairs!!!!!
extraMets=[];
place=0;
if any(strcmp('ExtraCell',core_ss)) %If i chose to do extracell mets, save them in 'extraMets' and deal with them separately
    place = find(strcmp('ExtraCell',core_ss));
    placeCentral = find_cell(setdiff(core_ss,'ExtraCell'),core_ss);
    extraMets = mets_ss(place);
else
    placeCentral = find(find_cell(core_ss,core_ss));
end
CentralMets = {unique(vertcat(mets_ss{placeCentral}))};


oldmets_ss = mets_ss;   
oldcores_ss = core_ss; 
if place>0
mets_ss(place) = [];
core_ss(place)= [];
end


% Make one metabolite list as source, and the other as target:
% - assign numbers to each subsystem, and form all pairs
pairs = combnk(1:length(mets_ss),2);


% - initialize cell structures for source/target metabolite groups of each
%   subsystem.
source_met=cell(size(pairs,1),1);
target_met=cell(size(pairs,1),1);

% For each subsystem pair
if strcmp(OnlyConnectExclusiveMets, 'yes')
    for i=1:size(pairs,1)
        % - Find the source metabolites: mets that are in SS1 but not in SS2
        source_met{i}=setdiff(mets_ss{pairs(i,1)}, mets_ss{pairs(i,2)})';
        % - Find the target metabolites: mets that are in SS2 but not in SS1
        target_met{i}=setdiff(mets_ss{pairs(i,2)}, mets_ss{pairs(i,1)})';
    end
else
    for i=1:size(pairs,1)
        % - source metabolites: mets that are in SS1
        source_met{i}=mets_ss{pairs(i,1)};
        % - target metabolites: mets that are in SS2
        target_met{i}=mets_ss{pairs(i,2)};
    end
end

%if we want to connect the metabolites within a subsystem to eachother
if strcmp(ConnectIntraSubsystem, 'yes')
    %extend the pairs
    pairs = [pairs; repmat((1:length(core_ss))',1,2)];
    
    for i=1:length(core_ss)
        % - the source metabolites are the mets that are in SS1
        source_met{end+1}=mets_ss{i};
        % - we are connecting them to mets that are also in SS1
        target_met{end+1}=mets_ss{i};
    end
end


%% get unique reactions and metabolites in the new model
selectedPaths = pathsNum(L_DirAdjMatWn, L_DirAdjMatC, L_DirAdjMatMets, pairs, source_met, target_met, startFromMin, D, ApplyShortestDistanceOfSubsystems,ThrowErrorOnDViolation);

if ~isempty(extraMets)
    selectedPaths.ModelMets
    selectedPathsExternal = pathsNumExtra(L_DirAdjMatWn, L_DirAdjMatC, L_DirAdjMatMets, extraMets, 1, CentralMets, ApplyShortestDistanceOfSubsystems);
else
    selectedPathsExternal=[];
end

%% final selection

allstraight = cellfun(@(c)c(:),  selectedPaths.joiningReacs(:,:,1:D)  , 'UniformOutput',false);
% the connecting reactions between the subsystems should also include the
% reactions that connect the "extracellular sybsystem" with the rest
% subsystems.
if isfield(selectedPathsExternal,'ModelMets')
    connecting_reactions= unique([vertcat( allstraight{:,:,:}) ; selectedPathsExternal.ModelReacs{1}]);
else
    connecting_reactions= unique(vertcat(allstraight{:,:,:}));
end
% everything other than: 1)the core-system-reactions, and 2) the connecting reactions:
withconnectors=setdiff(setdiff(1:length(model.rxns),connecting_reactions),vertcat(rxns_ss{:}));

% dont forget to add the transports
if isfield(selectedPathsExternal,'ModelMets')
    otherReactions=setdiff(withconnectors,selectedPathsExternal.ModelReacs{1});%THIS WILL CAUSE AN ERROR IF WE DO MORE THAN ONE D FOR THE EXTERNAL CONNECTING PATHS
else
    otherReactions=withconnectors;
end



