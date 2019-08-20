function [RedGEMOpts,paramNames] = redGEMQuestionnaire(varargin)
% This function takes care of the parameter set-up, according to which
% redGEM will function.

% Initialize a structure with all the parameters set to empty. In this way,
% if they are not defined in the input (varargin) questions will be triggered.

RedGEMOpts =  struct (...
    'Organism',                            [],... %(1)% human, ecoli, putida
    'GEMname',                             [],... %(2)% name of the GEM used (as it is saved in GEMs folder)
    'RedModelName',                        [],... %(3)% choose a name for the reduction
    'SelectedSubsystems',                  [],... %(4)% default (default is defined in the organism&GEM case-file), customly defined in a cell, e.g. {{'x';'y';'z'}}
    'AddETCAsSubsystem',                   [],... %(5)% yes, no
    'AddExtracellularSubsystem',           [],... %(6)% % no, automatic, default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'DM_x';'DM_y';'DM_z'}}
    'AerobicAnaerobic',                    [],... %(7)% aerobic, anaerobic
    'ListForInorganicMets',                [],... %(8)% curated, automatic
    'ListForCofactorPairs',                [],... %(9)% curated, automatic
    'ZeroZeroGEMbounds',                   [],... %(10)% Original, DefineCustom, OpenTo100
    'FluxUnits',                           [],... %(11)% mmol, mumol, other
    'L',                                   [],... %(12)% 1,2,...
    'D',                                   [],... %(13)% 1,2,...
    'startFromMin',                        [],... %(14)% yes, no
    'ThrowErrorOnDViolation',              [],... %(15)% error/continue If any two subsystems cannont connect to the level D
    'OnlyConnectExclusiveMets',            [],... %(16)% yes, no
    'ConnectIntracellularSubsystems',      [],... %(17)% yes, no
    'ApplyShortestDistanceOfSubsystems',   [],... %(18)% bothways,eachdirection'
    'performREDGEMX',                      [],... %(19)%yes, no                    
    'NumOfConnections',                    [],... %(20)% OnePerMetE ,SminMetE
    'performLUMPGEM',                      [],... %(21)% do we want to perform lumping or not?
    'PreventBBBuptake',                    [],... %(22)% yes/no: if yes, do not allow uptake through any bbb drains.
    'NumOfLumped',                         [],... %(23)% OnePerBBB, Smin, Sminp1, Sminp2, Sminp3
    'AlignTransportsUsingMatFile',         [],... %(24)% yesusingmatfile, yesusingReducedmatfile,yesautomatic,no
    'performPostProcessing',               [],... %(25)% no, PP_forMCA, PP_removeBlockedRxns
    'TimeLimitForSolver',                  [],... %(26)% yes, no
    'CplexParameters',                     [],... %(27)% CplexDefault, LCSBDefault, redGEM_m7, redGEM_m8
    'CPLEX_PATH',                          [],... %(28)% Provide your CPLEX path
    'TFA_PATH',                            [],... %(29)% Provide your GIT-path
    'thermo_data_PATH',                    []...  %(30)% Provide the path to the thermodynamic data (including the file name)
    );


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

%%%%%%%%% 1: 'Organism' %%%%%%%%%%%%%%%%%%%%
prompt{1}   = 'The GEM of which organism do you want to reduce?\n- ecoli \n- yeast \n- human \n- plasmodium';
choices{1}  = {'ecoli','yeast','human','plasmodium'};
%%%%%%%%% 2: 'GEMname' %%%%%%%%%%%%%%%%%%%%
prompt{2}   = 'Which GEM of the selected organism do you want to reduce?\n- iJO1366 \n- iAF1260 \n- iMM904 \n- yeast7';
choices{2}  = {'iJO1366','iAF1260','iMM904','yeast7'};
%%%%%%%%% 3: 'RedModelName' %%%%%%%%%%
prompt{3}  = 'Type a desired name for the reduced model?';
choices{3} = {'StringComment'};
%%%%%%%%% 4: 'SelectedSubsystems' %%%%%%%%%%%%%%%%%%%%
prompt{4}   = 'Would you like to define the subsystems or to be set to a default set specified in the organism&GEM case-file, or customly defined in a cell e.g. {{''x'';''y'';''z''}}';
choices{4}  = {'default', 'custom' };
%%%%%%%%% 5: 'AddETCAsSubsystem' %%%%%%%%%%%%%%%%%%%%
prompt{5}   = 'Would you like to add ETC to the core subsystems of the GEM?\n- yes \n- no';
choices{5}  = {'yes', 'no'};
%%%%%%%%% 6: 'AddExtracellularSubsystem' %%%%%%%%%%%%%%%%%%%%
prompt{6}  = 'Would you like to add an extracellular subsystem to the core subsystems? If yes, would you like to list the reactions in a cell?\n- no \n- default \n custom : {{''DM_x'',''DM_y'',''DM_z'',''etc...''}}';
choices{6} = {'no','default','custom','automatic'};
%%%%%%%%% 7: 'AerobicAnaerobic' %%%%%%%%%%%%%%%%%%%%
prompt{7}  = 'Is the condition aerobic or anaerobic?\n- aerobic \n- anaerobic';
choices{7} = {'aerobic', 'anaerobic'};
%%%%%%%%% 8: 'ListForInorganicMets' %%%%%%%%%%%%%%%%%%%%
prompt{8}   = 'Which list of inorganic metabolites do you want to use?\n- automatic: automatically calculated list\n- curated: manually curated list from matfile';
choices{8}  = {'automatic', 'curated'};
%%%%%%%%% 9: 'ListForCofactorPairs' %%%%%%%%%%%%%%%%%%%%
prompt{9}   = 'Which list of cofactor pairs do you want to remove?\n- automatic: automatically calculated list \n- curated: manually curated list from matfile';
choices{9}  = {'automatic', 'curated'};
%%%%%%%%% 10: 'ZeroZeroGEMbounds' %%%%%%%%%%%%%%%%%%%%
prompt{10}  = 'For the zero-zero GEM bounds, would you like to proceed with these bounds, open these bounds to +-100, or define custom?\n- Original \n- OpenTo100 \n- DefineCustom';
choices{10} = {'Original', 'OpenTo100', 'DefineCustom'};
%%%%%%%%% 11: 'FluxUnits' %%%%%%%%%%%%%%%%%%%%
prompt{11}  = 'What are the units of the fluxes?\n- mmol/(gDW*h) \n- mu-mmol/(gDW*h)';
choices{11} = {'mmol', 'mumol'};
%%%%%%%%% 12: 'L' %%%%%%%%%%%%%%%%%%%%
prompt{12}   = 'Please specify the L (0,1,2,...)';
choices{12}  = 'GE0Int'; % integer greater or equal to zero
%%%%%%%%% 13: 'D' %%%%%%%%%%%%%%%%%%%%
prompt{13}   = 'Please specify the D (0,1,2,...)';
choices{13}  = 'GE0Int'; % integer greater or equal to zero
%%%%%%%%% 14: 'startFromMin' %%%%%%%%%%%%%%%%%%%%
prompt{14}   = 'Please specify: Do you want to start from the min? \n- yes \n- no';
choices{14}  = {'yes', 'no'};
%%%%%%%%% 15: 'ThrowErrorOnDViolation' %%%%%%%%%%%%%%%%%%%%
prompt{15}  = 'Subsystem distance greater than specified D. Want to throw an error or continue up to feasible D?\n- error \n- continue';
choices{15} = {'error', 'continue'};
%%%%%%%%% 16: 'OnlyConnectExclusiveMets' %%%%%%%%%%%%%%%%%%%%
prompt{16}  = 'Do you want to connect only the exclusive metabolites among subsystems?\n- yes \n- no';
choices{16} = {'yes', 'no'};
%%%%%%%%% 17: 'ConnectIntracellularSubsystems' %%%%%%%%%%%%%%%%%%%%
prompt{17}  = 'Do you want to connect the metabolites within the subsystem with each other?\n- yes \n- no';
choices{17} = {'yes', 'no'};
%%%%%%%%% 18: 'ApplyShortestDistanceOfSubsystems' %%%%%%%%%%%%%%%%%%%%
prompt{18}  = 'Do you want the shortest distance between subsystems to be applied for both ways (a to b and b to a) or for each direction?\n- bothways \n- eachdirection';
choices{18} = {'bothways','eachdirection'};
%%%%%%%%% 19: 'performREDGEMX' %%%%%%%%%%%%%%%%%%%%
prompt{19}  = 'Would you like to connect the medium to the core?\n- yes \n- no';
choices{19} = {'yes', 'no'};
%%%%%%%%% 20: 'NumOfConnections' %%%%%%%%%%%%%%%%%%%%
prompt{20}  = 'How many connections should be generated \n- OnePerMetE \n- SminMetE';
choices{20} = {'OnePerMetE', 'SminMetE'};
%%%%%%%%% 21: 'performLUMPGEM' %%%%%%%%%%%%%%%%%%%%
prompt{21}  = 'Do you want to perform lumping?\n- yes \n- no';
choices{21} = {'yes', 'no'};
%%%%%%%%% 22: 'PreventBBBuptake' %%%%%%%%%%%%%%%%%%%%
prompt{22}  = 'Do you want to prevent uptake through the drains that go to bbbs directly?\n- yes \n- no';
choices{22} = {'yes', 'no'};
%%%%%%%%% 23: 'NumOfLumped' %%%%%%%%%%%%%%%%%%%%
prompt{23}  = 'How many lumped reactions should be generated?\n- One per bbb \n- Entire Smin network \n- Entire Smin+1 network \n- Entire Smin+2 network';
choices{23} = {'OnePerBBB', 'Smin', 'Sminp1', 'Sminp2', 'Sminp3'};
%%%%%%%%% 24: 'AlignTransportsUsingMatFile' %%%%%%%%%%%%%%%%%%%%
prompt{24}  = 'Do you want to align the transport the same metabolites?\n- yes : find these reactions and add constraints to align them\n- no : proceed without additional constraints';
choices{24} = {'yesusingmatfile','yesusingReducedmatfile','yesautomatic','no'};
%%%%%%%%% 25: 'PostProcessingForMCA' %%%%%%%%%%%%%%%%%%%%
prompt{25}  = 'Would you like to perform post-processing after the generation of the redGEM? Do you want just to  to prepare the model for MCA analysis? \n- no \n- PP_forMCA \n- PP_removeBlockedRxns';
choices{25} = {'no', 'PP_forMCA', 'PP_removeBlockedRxns'};
%%%%%%%%% 26: 'TimeLimitForSolver' %%%%%%%%%%%%%%%%%%%%
prompt{26}  = 'Do you want to use a time-limit for the solver?\n- yes\n- no';
choices{26} = {'yes', 'no'};
%%%%%%%%% 27: 'CplexParameters' %%%%%%%%%%%%%%%%%%%%
prompt{27}  = 'What parameters would you like to have for the cplex-solver?\n- Default cplex (e.g. tol 1e-6)? \n- Default LCSB (e.g. tol 1e-9)?';
choices{27} = {'CplexDefault', 'LCSBDefault', 'redGEM_m7', 'redGEM_m8'};
%%%%%%%%% 28: 'TFA_PATH' %%%%%%%%%%%%%%%%%%%%
prompt{28}  = 'Provide the path of your TFA-folder?';
choices{28} = {'StringComment'};
%%%%%%%%% 29: 'CPLEX_PATH' %%%%%%%%%%%%%%%%%%%%
prompt{29}  = 'Provide the path of your CPLEX @Cplex folder?';
choices{29} = {'StringComment'};
%%%%%%%%% 30: 'thermo_data_PATH' %%%%%%%%%%%%%%%%%%%%
prompt{30}  = 'Provide the path of the thermodynamic data (including the file name)?';
choices{30} = {'StringComment'};

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
