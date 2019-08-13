function [DirAdjMatW1, DirAdjMatWn, DirAdjMatC] = adjacentFromStoich_150218(model)
% This function creates a series of directed adjacency matrices (DirAdjMat). 
% ATTENTION: Directionality is defined based on lb and ub information!!!
% 
% INPUT 
% - model: metabolic model of COBRAtoolbox format.
% 
% OUTPUT
% - DirAdjMatW1 : a directed adjacency matrix of weights = 1 (W1)
% - DirAdjMatWn : This is like DirAdjMatW1, however the weight (W) is the
%   number of reactions (n) that connect the two metabolites. This makes it
%   useful for calculating the total number of connections.
% - DirAdjMatC : This is like DirAdjMatWn, but is a cell array instead of a
%   matrix. Each cell contains the distinct reactions between two
%   metabolites. This makes it useful for the distinctive enumeration of
%   all paths.
% 
% Note: at somepoint, we can integrate element balance to not allow unrelated
% metabolites to lead to each other. because of directionality, it is 
% not the same thing to go from A to B as it is from B to A.
%
% Modified: Tuure and GeorgiosF on 18th Feb 2015 modified the definition of
% the erroneous to exclude the cases of  [0 0] bounds. They are now
% automatically adjusted to [-100 100].

[num_mets,num_rxns]=size(model.S);
DirAdjMatW1 = zeros(num_mets,num_mets);%set memory

% Some bounds might be completely erroneous/wrong, if for example they have
% the lower bound of the flux being higher than the upper bound:
erroneous_bounds = find(model.lb>model.ub);

if sum(erroneous_bounds)>0
    warning('Bounds error: lb>ub!!! Automatically adjusted to -100 +100')
    model.lb(erroneous_bounds)=-100;
    model.ub(erroneous_bounds)= 100;
end

zero_zero_bounds = find((model.lb==0 & model.ub==0)==1);

if sum(zero_zero_bounds)>0
    warning('Bounds error: lb=ub=0!!! Automatically adjusted to -100 +100')
    model.lb(zero_zero_bounds)=-100;
    model.ub(zero_zero_bounds)= 100;
end

% Classify reactions as reversible, scrictly pos, or strictly neg
irrevneg = model.ub<=0; % strictly negative
irrevpos = model.lb>=0; % strictly positive

    
low = model.lb<0; 
up  = model.ub>0; 
reversible = (low & up); %reversible

% Now this should be zero!!! We fixed it above
erroneous_bounds = (irrevneg & irrevpos);


% Double check for consistency
if (sum(reversible) + sum(irrevneg) + sum(irrevpos))~=(num_rxns) && sum(erroneous_bounds)>0
    error('bounds arent correctly adjusted')
end

clear low up erroneous_bounds

DirAdjMatC=cell(num_mets, num_mets); %to save all reactions

%positive flow
posS=model.S;
posS(:,irrevneg)=0;

%% test case
%posS=[-1 -1 -1;1 1 0;1 0 0;1 0 0;0 1 0; 0 0 1];
%%

for i=1:num_mets % set the corresponding coefs for each met-met combination     
    reactions=find(posS(i,:)<0);    
    if ~isempty(reactions)
        [row col]=find(posS(:,reactions)>0);  
        
        
            ux = unique(row); %count the num of reactants
            if length(ux) == 1, DirAdjMatWn = length(row);
            else DirAdjMatWn = hist(row,ux);
            end
        DirAdjMatW1(i,ux)=DirAdjMatWn;
        
        for c=1:length(col) %add to connections array
        DirAdjMatC(i,row(c)) = {[DirAdjMatC{i,row(c)} reactions(col(c))]};
        end
    end
end

%negative flow
negS=model.S;
negS(:,irrevpos)=0;

for i=1:num_mets % set the corresponding coefs for each met-met combination
    reactions=find(negS(i,:)>0);
    if ~isempty(reactions)
    [row col]=find(negS(:,reactions)<0);  
        if~isempty(row)
        ux = unique(row);%count the num of reactions
            if length(ux) == 1, DirAdjMatWn = length(row);
            else DirAdjMatWn = hist(row,ux);
            end

        DirAdjMatW1(i,ux)=DirAdjMatW1(i,ux)+DirAdjMatWn;    
        end

        for c=1:length(col)%add to connections array
        DirAdjMatC(i,row(c)) = {[DirAdjMatC{i,row(c)} reactions(col(c))]};
        end
    end
end

DirAdjMatWn=DirAdjMatW1;
DirAdjMatW1(find(DirAdjMatW1))=1;

