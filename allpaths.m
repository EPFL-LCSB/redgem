function [L_DirAdjMatWn, L_DirAdjMatC, L_DirAdjMatMets] = allpaths(L, DirAdjMatWn, DirAdjMatC)
% This function generates all pathways of length L
% 
% INPUTS
% - L: (defined in redGEM)
%   is an adjustable parameter indicating how far from each subsystem is
%   the algorithm allowed to search in order to find connections among
%   the subsystems. More precicely, it indicates how distant, in terms of
%   reactions, can the metabolites of any two core sybsystems be from any
%   potential interconnecting metabolite.
% - DirAdjMatWn : (calculated in adjacentFromStoich)
%   This is like DirAdjMatW1, however the weight (W) is the number
%   of reactions (n) that connect the two metabolites. This makes it useful
%   for calculating the total number of connections.
% - DirAdjMatC : (calculated in adjacentFromStoich)
%   This is like DirAdjMatWn, but is a cell array instead of a
%   matrix. Each cell contains the separate reactions between two
%   metabolites. This makes it useful for the distinctive enumeration of
%   all paths.
% 
% OUTPUTS
% - L_DirAdjMatWn: This is like DirAdjMatWn, but it refers to connections of
%   degree L.
% - L_DirAdjMatC: This is like DirAdjMatC, but it refers to connections of
%   degree L.
% - L_DirAdjMatMets: This is like DirAdjMatC, but it explicitly lists all
%   the participating metabolites for each connection.

tic

[r1, c1]=find(DirAdjMatWn);
place=cell(size(DirAdjMatC,1),1);

L_DirAdjMatMets=cell( [size(DirAdjMatC) L]);
for i=1:size(DirAdjMatWn,1)  
    place{i}=unique(c1(ismember(r1,i)));
    for j=place{i}'
        L_DirAdjMatMets{i,j,1}=repmat([i;j],[1,length(DirAdjMatC{i,j})]);  
    end
end
clear i j

L_DirAdjMatWn=zeros( [size(DirAdjMatWn) L]);
L_DirAdjMatWn(:,:,1)=DirAdjMatWn;

L_DirAdjMatC=cell( [size(DirAdjMatWn) L]);
L_DirAdjMatC(:,:,1)=DirAdjMatC;


for i = 2:L  
    if i==2
        rm = r1; cm= c1;
    else
        [rm, cm]=find(L_DirAdjMatWn(:,:,i-1));
    end
    for firstMet = 1:size(DirAdjMatC,1)
        for currentMet=unique(cm(ismember(rm,firstMet)))'

            prevMets=L_DirAdjMatMets{firstMet,currentMet,i-1};
            prevReacs=L_DirAdjMatC{firstMet,currentMet,i-1};
            rep2=size(prevReacs,2);% # of first half connections

            nextMets=place{currentMet};
            nextMets=setdiff(nextMets,union(prevMets, firstMet));%this prevents the path from repeating nodes
            if ~isempty(nextMets)

                for nextMet=nextMets'  
                secondReacs=L_DirAdjMatC{currentMet,nextMet,1}; %second half of the transition taken from L_DirAdjMatWn 1
                rep1=size(secondReacs,2);% # of second half connections

                L_DirAdjMatC{firstMet,nextMet,i}=[L_DirAdjMatC{firstMet,nextMet,i} [repmat(prevReacs,[1 rep1]) ; kron(secondReacs, ones(1,rep2))]];
                L_DirAdjMatMets{firstMet,nextMet,i}=[L_DirAdjMatMets{firstMet,nextMet,i} [repmat(prevMets,[1 rep1]) ; repmat(nextMet, [1,rep1*rep2] )]];
                L_DirAdjMatWn(firstMet,nextMet,i)=L_DirAdjMatWn(firstMet,nextMet,i)+size(repmat(prevMets,[1 rep1]),2);
                end

            end
        end
    end
    toc
end

end

