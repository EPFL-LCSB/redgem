function [GAMequation, ppi_equation] = extractGAMcoefficientsFromBiomass(model)

idBiomass = find(model.c);

% '54.12 atp_c + 48.7529 h2o_c <=> 53.95 adp_c + 53.95 h_c + 53.9459 pi_c + 0.74983 ppi_c

GAMmets = {'atp_c'
    'h2o_c'
    'adp_c'
    'h_c'
    'pi_c'};

ppi = {'ppi_c'};

FullS = full(model.S);

if any(abs(FullS(find_cell(GAMmets,model.mets),idBiomass))==0)
    error('Some of the metabolites have zero coefficient in the biomass!!')
end
GAMMetCoeffs = cellfun(@(x) num2str(x),num2cell(abs(FullS(find_cell(GAMmets,model.mets),idBiomass))),'UniformOutput',false);


ppiCoeff = cellfun(@(x) num2str(x),num2cell(abs(FullS(find_cell(ppi,model.mets),idBiomass))),'UniformOutput',false);

GAMequation = [[GAMMetCoeffs repmat({' '},size(GAMmets),1)  GAMmets]';...
                [repmat({' '},1,size(GAMmets,1)-1),{''} ];
                [{'+'},{'<=>'},repmat({'+'},1,2),{''}];
                repmat({' '},1,size(GAMmets,1)-1),{''}];
GAMequation = cell2mat(GAMequation(:)');

ppi_equation = cell2mat([' <=> ',ppiCoeff,' ',ppi]);



end