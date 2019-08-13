function [RedGEMOpts,paramNames] = redGEMQuestionnaire(varargin)
% This function takes care of the parameter set-up, according to which
% redGEM will function.

% Initialize a structure with all the parameters set to empty. In this way,
% if they are not defined in the input (varargin) questions will be triggered.

RedGEMOpts =  struct (...
    'D',                                   [],... %(1)% 1,2,...
    'L',                                   [],... %(2)% 1,2,...
    'startFromMin',                        [],... %(3)% yes, no
    'viewStats',                           [],... %(4)% yes, no
    'Organism',                            [],... %(5)% ecoli
    'GEMname',                             [],... %(6)% iMM904, iJO1366, etc.
    'ListForInorganicMets',                [],... %(7)% curated, automatic
    'ListForCofactorPairs',                [],... %(8)% curated, automatic
    'SelectedSubsystems',                  [],... %(9)% default (default is defined in the organism&GEM case-file), customly defined in a cell, e.g. {{'x';'y';'z'}}
    'AddETCAsSubsystem',                   [],... %(10)% yes, no
    'AddExtracellularSubsystem',           [],... %(11)% % no, automatic, default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'DM_x';'DM_y';'DM_z'}}
    'ConnectExtracellularToCore',          [],... %(12)%yes, no                    
    'NumOfConnections',                    [],... %(13)% OnePerMetE ,SminMetE
    'AerobicAnaerobic',                    [],... %(14)% aerobic, anaerobic
    'AlignTransportsUsingMatFile',         [],... %(15)% yesusingmatfile, yesusingReducedmatfile,yesautomatic,no
    'TimeLimitForSolver',                  [],... %(16)% yes, no
    'RemovePeriplasm',                     [],... %(17)% yes, no
    'OnlyConnectExclusiveMets',            [],... %(18)% yes, no
    'ConnectIntracellularSubsystems',      [],... %(19)% yes, no
    'ApplyShortestDistanceOfSubsystems',   [],... %(20)% bothways,eachdirection'
    'RedModelName',                        [],... %(21)% choose a name
    'ThrowErrorOnDViolation',              [],... %(22)% error/continue If any two subsystems cannont connect to the level D
    'performLUMPING',                      [],... %(23)% do we want to perform lumping or not?
    'PreventBBBuptake',                    [],... %(24)% yes/no: if yes, do not allow uptake through any bbb drains.
    'ZeroZeroGEMbounds',                   [],... %(25)% Original, DefineCustom, OpenTo100
    'FluxUnits',                           [],... %(26)% mmol, mumol, other
    'NumOfLumped',                         [],... %(27)% OnePerBBB, Smin, Sminp1, Sminp2, Sminp3
    'performPostProcessing',               [],... %(28)% no, PP_forMCA, PP_removeBlockedRxns
    'CplexParameters',                     [],... %(29)% CplexDefault, LCSBDefault, redGEM_m7, redGEM_m8
    'GITpath',                             [],... %(30)% Provide your GIT-path
    'CPLEX_PATH',                          []...  %(31)% Provide your CPLEX path
    );

%     'Tighten',                             [],... % yes/no: add reactions that connect metabolites without adding mass balances

% Create a string containing all the parameter names. This will be used in
% the main function to set values to all the parameter names individually.
paramNames = fieldnames(RedGEMOpts);

% If the input is not empty
if ~strcmp(varargin,'NoInput')
    % Extract the parameter names that were set in the input
    varargin = varargin{1};
    fieldnamesVarargin = fieldnames(varargin);
    
    % Check: if the parameter names are consistent in the entire parameter
    % vector. If not, then trigger an error
    if size(intersect(fieldnamesVarargin,paramNames),1)~=size(fieldnamesVarargin,1)
        error('parameter input names are not contained in the entire parameter vector')
    end
    
    % Overwrite the empty parameter values with those that have been defined by
    % the input (varargin)
    for i=1:size(fieldnamesVarargin,1)
        eval(['RedGEMOpts.',fieldnamesVarargin{i},' = getfield(varargin,''',fieldnamesVarargin{i},''')'])
    end
end


% For the parameters that have not been defined by the user in the input,
% ask the user to define them using a questionnaire:

%%%%%%%%% 1: 'D' %%%%%%%%%%%%%%%%%%%%
prompt{1}   = 'Please specify the D (0,1,2,...)';
choices{1}  = 'GE0Int'; % integer greater or equal to zero
%%%%%%%%% 2: 'L' %%%%%%%%%%%%%%%%%%%%
prompt{2}   = 'Please specify the L (0,1,2,...)';
choices{2}  = 'GE0Int'; % integer greater or equal to zero
%%%%%%%%% 3: 'startFromMin' %%%%%%%%%%%%%%%%%%%%
prompt{3}   = 'Please specify: Do you want to start from the min? \n- yes \n- no';
choices{3}  = {'yes', 'no'};
%%%%%%%%% 4: 'viewStats' %%%%%%%%%%%%%%%%%%%%
prompt{4}   = 'Would you like to perform statistic analysis? \n- yes \n- no';
choices{4}  = {'yes', 'no'};
%%%%%%%%% 5: 'Organism' %%%%%%%%%%%%%%%%%%%%
prompt{5}   = 'The GEM of which organism do you want to reduce?\n- ecoli \n- yeast \n- human \n- plasmodium';
choices{5}  = {'ecoli','yeast','human','plasmodium'};
%%%%%%%%% 6: 'GEMname' %%%%%%%%%%%%%%%%%%%%
prompt{6}   = 'Which GEM of the selected organism do you want to reduce?\n- iJO1366 \n- iAF1260 \n- iMM904 \n- yeast7';
choices{6}  = {'iJO1366','iAF1260','iMM904','yeast7'};
%%%%%%%%% 7: 'ListForInorganicMets' %%%%%%%%%%%%%%%%%%%%
prompt{7}   = 'Which list of inorganic metabolites do you want to use?\n- automatic: automatically calculated list\n- curated: manually curated list from matfile';
choices{7}  = {'automatic', 'curated'};
%%%%%%%%% 8: 'ListForCofactorPairs' %%%%%%%%%%%%%%%%%%%%
prompt{8}   = 'Which list of cofactor pairs do you want to remove?\n- automatic: automatically calculated list \n- curated: manually curated list from matfile';
choices{8}  = {'automatic', 'curated'};
%%%%%%%%% 9: 'SelectedSubsystems' %%%%%%%%%%%%%%%%%%%%
prompt{9}   = 'Would you like to define the subsystems or to be set to a default set specified in the organism&GEM case-file, or customly defined in a cell e.g. {{''x'';''y'';''z''}}';
choices{9}  = {'default', 'custom' };
%%%%%%%%% 10: 'AddETCAsSubsystem' %%%%%%%%%%%%%%%%%%%%
prompt{10}   = 'Would you like to add ETC to the core subsystems of the GEM?\n- yes \n- no';
choices{10}  = {'yes', 'no'};
%%%%%%%%% 11: 'AddExtracellularSubsystem' %%%%%%%%%%%%%%%%%%%%
prompt{11}   = 'Would you like to add an extracellular subsystem to the core subsystems? If yes, would you like to list the reactions in a cell?\n- no \n- default \n custom : {{''DM_x'',''DM_y'',''DM_z'',''etc...''}}';
choices{11}  = {'no','default','custom','automatic'};
%%%%%%%%% 12: 'AddExtracellularSubsystem' %%%%%%%%%%%%%%%%%%%%
prompt{12}   = 'Would you like to connect the medium to the core?\n- yes \n- no';
choices{12}  = {'yes', 'no'};
%%%%%%%%% 13: 'AerobicAnaerobic' %%%%%%%%%%%%%%%%%%%%
prompt{13}  = 'How many [e] connections should be generated?\n- One per [e]-met \n- Entire Smin network';
choices{13} = {'OnePerMetE' ,'SminMetE'};
%%%%%%%%% 14: 'AerobicAnaerobic' %%%%%%%%%%%%%%%%%%%%
prompt{14}  = 'Is the condition aerobic or anaerobic?\n- aerobic \n- anaerobic';
choices{14} = {'aerobic', 'anaerobic'};
%%%%%%%%% 15: 'AlignTransportsUsingMatFile' %%%%%%%%%%%%%%%%%%%%
prompt{15}  = 'Do you want to align the transport the same metabolites?\n- yes : find these reactions and add constraints to align them\n- no : proceed without additional constraints';
choices{15} = {'yesusingmatfile','yesusingReducedmatfile','yesautomatic','no'};
%%%%%%%%% 16: 'TimeLimitForSolver' %%%%%%%%%%%%%%%%%%%%
prompt{16}  = 'Do you want to use a time-limit for the solver?\n- yes\n- no';
choices{16} = {'yes', 'no'};
%%%%%%%%% 17: 'RemovePeriplasm' %%%%%%%%%%%%%%%%%%%%
prompt{17}  = 'Would you like to remove the periplasmic compartment from the reduced model?\n- yes \n- no';
choices{17} = {'yes', 'no'};
%%%%%%%%% 18: 'OnlyConnectExclusiveMets' %%%%%%%%%%%%%%%%%%%%
prompt{18}  = 'Do you want to connect only the exclusive metabolites among subsystems?\n- yes \n- no';
choices{18} = {'yes', 'no'};
%%%%%%%%% 19: 'ConnectIntracellularSubsystems' %%%%%%%%%%%%%%%%%%%%
prompt{19}  = 'Do you want to connect the metabolites within the subsystem with each other?\n- yes \n- no';
choices{19} = {'yes', 'no'};
%%%%%%%%% 20: 'ApplyShortestDistanceOfSubsystems' %%%%%%%%%%%%%%%%%%%%
prompt{20}  = 'Do you want the shortest distance between subsystems to be applied for both ways (a to b and b to a) or for each direction?\n- bothways \n- eachdirection';
choices{20} = {'bothways','eachdirection'};
%%%%%%%%% 21: 'RedModelName' %%%%%%%%%%
prompt{21}  = 'Type a desired name for the reduced model?';
choices{21} = {'StringComment'};
%%%%%%%%% 22: 'ThrowErrorOnDViolation' %%%%%%%%%%%%%%%%%%%%
prompt{22}  = 'Subsystem distance greater than specified D. Want to throw an error or continue up to feasible D?\n- error \n- continue';
choices{22} = {'error', 'continue'};
%%%%%%%%% 23: 'performLUMPING' %%%%%%%%%%%%%%%%%%%%
prompt{23}  = 'Do you want to perform lumping?\n- yes \n- no';
choices{23} = {'yes', 'no'};
%%%%%%%%% 24: 'PreventBBBuptake' %%%%%%%%%%%%%%%%%%%%
prompt{24}  = 'Do you want to prevent uptake through the drains that go to bbbs directly?\n- yes \n- no';
choices{24} = {'yes', 'no'};
%%%%%%%%% 25: 'ZeroZeroGEMbounds' %%%%%%%%%%%%%%%%%%%%
prompt{25}  = 'For the zero-zero GEM bounds, would you like to proceed with these bounds, open these bounds to +-100, or define custom?\n- Original \n- OpenTo100 \n- DefineCustom';
choices{25} = {'Original', 'OpenTo100', 'DefineCustom'};
%%%%%%%%% 26: 'FluxUnits' %%%%%%%%%%%%%%%%%%%%
prompt{26}  = 'What are the units of the fluxes?\n- mmol/(gDW*h) \n- mu-mmol/(gDW*h)';
choices{26} = {'mmol', 'mumol'};
%%%%%%%%% 27: 'NumOfLumped' %%%%%%%%%%%%%%%%%%%%
prompt{27}  = 'How many lumped reactions should be generated?\n- One per bbb \n- Entire Smin network \n- Entire Smin+1 network \n- Entire Smin+2 network';
choices{27} = {'OnePerBBB', 'Smin', 'Sminp1', 'Sminp2', 'Sminp3'};
%%%%%%%%% 28: 'CplexParameters' %%%%%%%%%%%%%%%%%%%%
prompt{28}  = 'Would you like to perform post-processing after the generation of the redGEM? Do you want just to  to prepare the model for MCA analysis? \n- no \n- PP_forMCA \n- PP_removeBlockedRxns';
choices{28} = {'no', 'PP_forMCA', 'PP_removeBlockedRxns'};
%%%%%%%%% 29: 'PostProcessingForMCA' %%%%%%%%%%%%%%%%%%%%
prompt{29}  = 'What parameters would you like to have for the cplex-solver?\n- Default cplex (e.g. tol 1e-6)? \n- Default LCSB (e.g. tol 1e-9)?';
choices{29} = {'CplexDefault', 'LCSBDefault', 'redGEM_m7', 'redGEM_m8'};
%%%%%%%%% 30: 'GITpath' %%%%%%%%%%%%%%%%%%%%
prompt{30}  = 'Provide the path of your GIT-folder?';
choices{30} = {'StringComment'};
%%%%%%%%% 31: 'CPLEX_PATH' %%%%%%%%%%%%%%%%%%%%
prompt{31}  = 'Provide the path of your CPLEX @Cplex folder?';
choices{31} = {'StringComment'};

% Add a short string asking to press enter after the input choice:
prompt = cellfun(@(x) [x,'\nchoose option and press enter: '],prompt,'UniformOutput',false)';

% Check if the size of the paramNames is the same as the choices and
% prompts:
if size(paramNames,1)~= size(prompt,1)
    error('Inconsistent number of questions and parameters!')
end


for i=1:size(prompt,1)
    % If the choice is a zero or positive integer:
    if strcmp(choices{i},'GE0Int')
        if isempty(eval(['RedGEMOpts.',paramNames{i}]))
            fprintf('==========================================\n')
            while ~(isposintscalar(eval(['RedGEMOpts.',paramNames{i}])) || issame(eval(['RedGEMOpts.',paramNames{i}]),0))
                eval(['RedGEMOpts.',paramNames{i},' = input(''',prompt{i},''');']);
            end
            fprintf('==========================================\n')
        end
        % If the choice is a string
    elseif strcmp(choices{i},'StringComment')
        if strcmp(paramNames{i},'CPLEX_PATH') && strcmp(prompt{i},'')
            CPLEX_PATH = fullfile(RedGEMOpts.GITpath,'CPLEX_Studio1251');
        end
        if isempty(eval(['RedGEMOpts.',paramNames{i}]))
            fprintf('==========================================\n')
            while isempty(['RedGEMOpts.',paramNames{i}])
                eval(['RedGEMOpts.',paramNames{i},' = input(''',prompt{i},''',''s'');']);
                % Allow for uppercase letters in the name
            end
            fprintf('==========================================\n')
        end
    else
        if isempty(eval(['RedGEMOpts.',paramNames{i}]))
            fprintf('==========================================\n')
            while isempty(intersect(eval(['RedGEMOpts.',paramNames{i}]),choices{i}))
                eval(['RedGEMOpts.',paramNames{i},' = input(''',prompt{i},''',''s'');']);
                eval(['RedGEMOpts.',paramNames{i},' = lower(RedGEMOpts.',paramNames{i},');']);
            end
            fprintf('==========================================\n')
        end
    end
end


end
