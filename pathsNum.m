function selectedPaths = pathsNum(L_DirAdjMatWn, L_DirAdjMatC, ...
    L_DirAdjMatMets, pairs, source_mets, target_mets, startFromMin, D, ...
    ApplyShortestDistanceOfSubsystems, ThrowErrorOnDViolation)
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
% - startFromMin:
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

numSS=length(unique(pairs));
sub2subNumP=zeros(numSS,numSS,D);

% Adding d_max for the cases where subsystems are not connected even after
% L iterations (it happens), causing the while loop to go look for indices
% larger than L in L_DirAdjMatWn (which is a tensor with 1:L as 3rd
% dimension index)
d_max = size(L_DirAdjMatWn,3);

if strcmp(startFromMin,'yes')
    shortestLsub2subNumP=zeros(numSS,numSS);%find shortest L between SS
    
    switch direction
        case {'eachdirection'}
            for k = 1:size(pairs,1)
                d=1;   
                while d < d_max &&...
                    all(~any(L_DirAdjMatWn(source_mets{k},target_mets{k},d)))
                    d=d+1;
                end
                if d <= d_max
                    shortestLsub2subNumP(pairs(k,1),pairs(k,2))=d;
                end
            end

            for k = 1:size(pairs,1)
                d=1;   
                while d < d_max &&...
                    all(~any(L_DirAdjMatWn(target_mets{k},source_mets{k},d)))%switch direction!
                    d=d+1;
                end
                if d <= d_max
                    shortestLsub2subNumP(pairs(k,2),pairs(k,1))=d;     %switch direction!
                end
            end
        case {'bothways'}
            for k = 1:size(pairs,1)
                d=1;   
                while d <= d_max &&...
                        (all(~any(L_DirAdjMatWn(target_mets{k},source_mets{k},d))) &&...
                        all(~any(L_DirAdjMatWn(source_mets{k},target_mets{k},d))))
                    d=d+1;
                end
                if d <= d_max
                    shortestLsub2subNumP(pairs(k,1),pairs(k,2))=d;
                    shortestLsub2subNumP(pairs(k,2),pairs(k,1))=d;     %switch direction!
                end
            end
    end
elseif strcmp(startFromMin,'no')
    shortestLsub2subNumP=ones(numSS,numSS);%find shortest L between SS
else 
    error('Wrong option for startFromMin')
end


% Try to find any subsystems that are not connected by L
[source_ix, target_ix] = find(shortestLsub2subNumP == 0);
if ~isempty(source_ix) || ~isempty(target_ix)
    if strcmp(ThrowErrorOnDViolation,'error')
        error('L value exceeded in matching subsystems :')
        disp(numSS(source_ix));
        fprintf('and\n');
        disp(numSS(target_ix));
    else
        warning('L value exceeded in matching subsystems :')
        disp(source_ix);
        fprintf('and\n');
        disp(target_ix);
    end
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
    for k=1:size(pairs,1)
        % Added a if condition, in the case that subsystems are not linked,
        % d will stay at 0
        if shortestLsub2subNumP(pairs(k,1),pairs(k,2)) > 0
            d=el-1+shortestLsub2subNumP(pairs(k,1),pairs(k,2));
            a=L_DirAdjMatC(source_mets{k}, target_mets{k},d);
            emptyCells = cellfun('isempty',a);
            a(emptyCells) = [];
            uniqueReacs=unique([a{:}],'stable');    
            joiningReacs{pairs(k,1),pairs(k,2),el}=uniqueReacs;
            joiningPaths{pairs(k,1),pairs(k,2),el}=[a{:}];
            sub2subNumP(pairs(k,1),pairs(k,2),el)=size([a{:}],2);
            sub2subNumR(pairs(k,1),pairs(k,2),el)=length(uniqueReacs);
            ModelReacs{el}=unique(   [ModelReacs{el}; uniqueReacs(:) ],'stable');   
            %collect metabolites
            m=L_DirAdjMatMets(source_mets{k}, target_mets{k},d);
            m(emptyCells) = [];
            uniqueMets=unique([m{:}],'stable');
            joiningMets{pairs(k,1),pairs(k,2),el}=uniqueMets;
            ModelMets{el}=unique(   [ModelMets{el}; uniqueMets(:) ],'stable');
        end
        if shortestLsub2subNumP(pairs(k,2),pairs(k,1)) > 0
            %other direction
            d=el-1+shortestLsub2subNumP(pairs(k,2),pairs(k,1));  %switch!
            a=L_DirAdjMatC(target_mets{k},source_mets{k},d);        %switch!
            emptyCells = cellfun('isempty',a);
            a(emptyCells) = [];
            uniqueReacs=unique([a{:}],'stable');
            joiningReacs{pairs(k,2),pairs(k,1),el}=uniqueReacs;  %switch!
            joiningPaths{pairs(k,2),pairs(k,1),el}=[a{:}];           %switch!
            sub2subNumP(pairs(k,2),pairs(k,1),el)=size([a{:}],2);
            sub2subNumR(pairs(k,2),pairs(k,1),el)=length(uniqueReacs);
            ModelReacs{el}=unique(   [ModelReacs{el}; uniqueReacs(:) ],'stable');
            %collect metabolites  
            m=L_DirAdjMatMets(target_mets{k}, source_mets{k},d);
            m(emptyCells) = [];
            uniqueMets=unique([m{:}],'stable');
            joiningMets{pairs(k,2),pairs(k,1),el}=uniqueMets;
            ModelMets{el}=unique(   [ModelMets{el}; uniqueMets(:) ],'stable');
        end
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

