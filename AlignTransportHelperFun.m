function RedModel = AlignTransportHelperFun(RedModel, AlignTransportsUsingMatFile, checkgrowth, CplexParameters,biomassRxnNames,ATPsynth_RxnNames)

if strcmp(AlignTransportsUsingMatFile,'yesautomatic')
    % Find the parallel-transport reactions, and the metabolites that these transport
    fprintf('Align parallel transports using a function that finds the parallel transports')
    % Find if there exist any LUMPED reactions or other reactions that we would
    % like to exclude from the search for the identification of
    % transported metabolites
    idLMPD = find(~cellfun(@isempty,regexpi(RedModel.rxns,'^LMPD')));
    idBiomass = find(RedModel.c);
    idRxnsToExclude = union(idLMPD,idBiomass);
    if ~isempty(idRxnsToExclude)
        ModelToFindTrans = removeRxns(RedModel,RedModel.rxns(idRxnsToExclude),false,false);
        % Find those unused metabolites and keep the indices
        UnusedMets = ModelToFindTrans.mets(any(sum(abs(ModelToFindTrans.S),2) == 0,2));
        if (~isempty(UnusedMets))
            ModelToFindTrans = removeMetabolites(ModelToFindTrans, UnusedMets, false);
        end
    end
    [AllTransports, TransportNoCouples, CoupledTransports, ImportantTransports, directions, TransportGroups] = identifyTransportRxns(ModelToFindTrans,biomassRxnNames,ATPsynth_RxnNames);
    rxns = [ImportantTransports(:,1) strcat(ImportantTransports(:,2), '_', ImportantTransports(:,3))];
    directions = ImportantTransports(:,4);
    directions = cell2mat(directions);
    
    RedModel = removefutiles(RedModel, rxns, directions, checkgrowth, CplexParameters);
elseif strcmp(AlignTransportsUsingMatFile,'yesusingmatfile')
    fprintf('Align parallel transports using a saved .mat file that is supposed to contain this information\n')
    RedModel = removefutilegs_Apr2016(RedModel);
elseif strcmp(AlignTransportsUsingMatFile,'no')
    fprintf('Proceed without additional directionality constraints')
else
    error('Wrong option!')
end


end