function selectedPaths = pathsNumExtra(L_DirAdjMatWn, L_DirAdjMatC, L_DirAdjMatMets,  extraCellMets, D, CentralMets, ApplyShortestDistanceOfSubsystems)
% This function finds the shortest distance from subsystem to subsystem.
% ATTENTION: THIS IS DIRECTION SPECIFIC !!!
% 
% INPUTS
% - L_DirAdjMatWn: (Calculated in allpaths)
%   This is like DirAdjMatWn, but it refers to connections of degree L.
% - L_DirAdjMatC: (Calculated in allpaths)
%   This is like DirAdjMatC, but it refers to connections of degree L.
% - L_DirAdjMatMets: (Calculated in allpaths)
%   This is like DirAdjMatC, but it explicitly lists all the participating
%   metabolites for each connection.
% - pairs:
% - source_mets:
% - target_mets:
% - shortestPlusL:
% - D: (Defined in redGEM_150112)
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

% direction = 'bothways'; % 'bothways' or 'eachdirection' here we specify if the shortest distance is for both ways (a to b and b to a) or for each direction
direction = ApplyShortestDistanceOfSubsystems;

shortestPlusL=true

if ~exist('D', 'var') | isempty(D)
    D=1;% flag. it returns first D levels
end


source_mets=[  repmat({CentralMets{1}'},length(extraCellMets{1}),1)  ]
target_mets =  num2cell(extraCellMets{1})
pairs = [repmat(1,length(extraCellMets{1}),1) (2:(length(target_mets))+1)']





numSS=length(unique(pairs));
sub2subNumP=zeros(numSS,numSS,D);

if shortestPlusL
    shortestLsub2subNumP=zeros(numSS,numSS);%find shortest L between SS
    
    switch direction
        case {'eachdirection'}
            for i = 1:size(pairs,1)
                d=1;
                while ~any(L_DirAdjMatWn(source_mets{i},target_mets{i},d))
                    d=d+1;
                end
                shortestLsub2subNumP(pairs(i,1),pairs(i,2))=d;
            end
            
            for i = 1:size(pairs,1)
                d=1;
                while ~any(L_DirAdjMatWn(target_mets{i},source_mets{i},d))%switch direction!
                    d=d+1;
                end
                shortestLsub2subNumP(pairs(i,2),pairs(i,1))=d;     %switch direction!
            end
        case {'bothways'}
            for i = 1:size(pairs,1)
                d=1;
                while all(~any(L_DirAdjMatWn(target_mets{i},source_mets{i},d))) & all(~any(L_DirAdjMatWn(source_mets{i},target_mets{i},d)))
                    d=d+1;
                end
                shortestLsub2subNumP(pairs(i,1),pairs(i,2))=d;
                shortestLsub2subNumP(pairs(i,2),pairs(i,1))=d;     %switch direction!
            end
    end
else
    shortestLsub2subNumP=ones(numSS,numSS);%find shortest L between SS
end
clear direction d i
%% get the reactions and pathways joining the subsystems
% this section simply condenses all the met to met pathways(reactions) and represents
% them as subsystem to subsystem

joiningReacs=cell(numSS,numSS,D);
joiningMets=cell(numSS,numSS,D);
joiningPaths=cell(numSS,numSS,D);
ModelPaths=cell(D,1);
ModelReacs=cell(D,1);
ModelMets=cell(D,1);

for el=1:D
    for i=1:size(pairs,1)
        d=el-1+shortestLsub2subNumP(pairs(i,1),pairs(i,2));
        a=L_DirAdjMatC(source_mets{i}, target_mets{i},d);
        emptyCells = cellfun('isempty',a);
        a(emptyCells) = [];
        uniqueReacs=unique([a{:}],'stable');    
        joiningReacs{pairs(i,1),pairs(i,2),el}=uniqueReacs;
        joiningPaths{pairs(i,1),pairs(i,2),el}=[a{:}];
        sub2subNumP(pairs(i,1),pairs(i,2),el)=size([a{:}],2);
        sub2subNumR(pairs(i,1),pairs(i,2),el)=length(uniqueReacs);
        ModelReacs{el}=unique(   [ModelReacs{el}; uniqueReacs(:) ],'stable');   
        %collect metabolites
        m=L_DirAdjMatMets(source_mets{i}, target_mets{i},d);
        m(emptyCells) = [];
        uniqueMets=unique([m{:}],'stable');
        joiningMets{pairs(i,1),pairs(i,2),el}=uniqueMets;
        ModelMets{el}=unique(   [ModelMets{el}; uniqueMets(:) ],'stable');

        %other direction
        d=el-1+shortestLsub2subNumP(pairs(i,2),pairs(i,1));  %switch!
        a=L_DirAdjMatC(target_mets{i},source_mets{i},d);        %switch!
        emptyCells = cellfun('isempty',a);
        a(emptyCells) = [];
        uniqueReacs=unique([a{:}],'stable');
        joiningReacs{pairs(i,2),pairs(i,1),el}=uniqueReacs;  %switch!
        joiningPaths{pairs(i,2),pairs(i,1),el}=[a{:}];           %switch!
        sub2subNumP(pairs(i,2),pairs(i,1),el)=size([a{:}],2);
        sub2subNumR(pairs(i,2),pairs(i,1),el)=length(uniqueReacs);
        ModelReacs{el}=unique(   [ModelReacs{el}; uniqueReacs(:) ],'stable');
        %collect metabolites  
        m=L_DirAdjMatMets(target_mets{i}, source_mets{i},d);
        m(emptyCells) = [];
        uniqueMets=unique([m{:}],'stable');
        joiningMets{pairs(i,2),pairs(i,1),el}=uniqueMets;
        ModelMets{el}=unique(   [ModelMets{el}; uniqueMets(:) ],'stable');
    end
end

selectedPaths.sub2subNumP	= sub2subNumP;
selectedPaths.sub2subNumR	= sub2subNumR;
selectedPaths.joiningReacs	= joiningReacs;
selectedPaths.joiningPaths	= joiningPaths;
selectedPaths.joiningMets	= joiningMets;
selectedPaths.ModelReacs	= ModelReacs;
selectedPaths.ModelMets	= ModelMets;

clear a emptyCells d i

