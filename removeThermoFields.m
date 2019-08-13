function model = removeThermoFields(model)

if (isfield(model,'metDeltaGFstd'))
    model = rmfield(model, 'metDeltaGFstd');
end
if (isfield(model,'metDeltaGFerr'))
    model = rmfield(model, 'metDeltaGFerr');
end
if (isfield(model,'metCharge'))
    model = rmfield(model, 'metCharge');
end
if (isfield(model,'metMass'))
    model = rmfield(model, 'metMass');
end
if (isfield(model,'metDeltaGFtr'))
    model = rmfield(model, 'metDeltaGFtr');
end
if (isfield(model,'struct_cues'))
    model = rmfield(model, 'struct_cues');
end
if (isfield(model,'rxnThermo'))
    model = rmfield(model, 'rxnThermo');
end
if (isfield(model,'rxnDeltaGR'))
    model = rmfield(model, 'rxnDeltaGR');
end
if (isfield(model,'rxnDeltaGRerr'))
    model = rmfield(model, 'rxnDeltaGRerr');
end
if (isfield(model,'rxnComp'))
    model = rmfield(model, 'rxnComp');
end
if (isfield(model,'rxnMapResult'))
    model = rmfield(model, 'rxnMapResult');
end
if (isfield(model,'isCoreTrans'))
    model = rmfield(model, 'isCoreTrans');
end

end