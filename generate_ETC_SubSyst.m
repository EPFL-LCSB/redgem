function model = generate_ETC_SubSyst(model, OxPhosSubsystem)

all_mets={'cpd15560'	'cpd15561';
'cpd15499'	'cpd15500';
'cpd15352'	'cpd15353';
'cpd00109'	'cpd00110';
'cpd01351'  'cpd11665'
'cpd00986'  'cpd00097'};

for i=1:length(model.rxns)
    mets_match={};
    metSEED=model.metSEEDID(find(model.S(:,i)));
    for j=1:size(all_mets,1)
        [~, ba]=ismember(all_mets(j,:), metSEED);
        if length(find(ba))==2
            mets_match=[mets_match all_mets(j,:)];
        end
    end
    if exist('mets_match', 'var')
        
        if length(mets_match) > 1
            model.isETC(i,1)=1;
        else
            model.isETC(i,1)=0;
        end
    else
        model.isETC(i,1)=0;
    end
end


%Update the subsystems vector in the model
model.subSystems = MakeCellArraySearchable(model.subSystems,'None');
model.subSystems(model.isETC==1) = {'ETC_Rxns'};

[~, boxidative]=ismember(model.subSystems,OxPhosSubsystem);
model.subSystems(find(boxidative))={'ETC_Rxns'};

end
