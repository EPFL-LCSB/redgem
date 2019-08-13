function getEcoli_D1(inputUserInitials)
restoredefaultpath

if ~exist('inputUserInitials','var')||isempty(inputUserInitials)
    prompt  = 'Please type your initials (select from the list) and press enter. Which one are you? \n-BZ \n-GF \n-MA \n-MM \n ... ';
    choices  = {'GF', 'BL', 'MA','MM'};
    currentInitials = [];
    while isempty(currentInitials)||~ismember(currentInitials, choices)
        currentInitials = input(prompt,'s');
    end
else
    currentInitials = inputUserInitials;
end

switch currentInitials
    case 'GF'
        redGEMlumpGEMpath = '/Users/georgiosfengos/Dropbox/SharedFolders/EPFL-LCSB/Projects/RedAndLumpGEM';
        OUTPUTpath = '/Users/georgiosfengos/TEMP/redGEMlumpGEMOUTPUT';
        CPLEX_PATH = '/Users/georgiosfengos/Applications/IBM/ILOG/CPLEX_Studio1271';
    case 'BL'
        redGEMlumpGEMpath = '/Users/Beatriz/Dropbox/RedAndLumpGEM';
        OUTPUTpath = '/Users/Beatriz/Documents/redGEMOUTPUT';
        CPLEX_PATH = '/Users/Beatriz/Applications/IBM/ILOG/CPLEX_Studio1271';
    case 'MA'
        redGEMlumpGEMpath = '';
        OUTPUTpath = '';
        CPLEX_PATH = '';
    case 'MM'
        redGEMlumpGEMpath = '';
        OUTPUTpath = '';
        CPLEX_PATH = '';
    case 'LM'
        redGEMlumpGEMpath = '';
        OUTPUTpath = '';
        CPLEX_PATH = '';
    otherwise
        error('Unknown user initials!!')
end


paramEcoli    = struct ('L',                                  6                                      ,... % 1,2,...
                        'D',                                  1                                      ,... % 1,2,...
                        'startFromMin',                       'no'                                   ,... % yes, no
                        'viewStats',                          'no'                                   ,... % yes, no
                        'Organism',                           'ecoli'                                ,... % ecoli, putida, putida butanol
                        'GEMname',                            'iJO1366'                              ,... % iMM904, iJO1366, etc.
                        'ListForInorganicMets',               'automatic'                            ,... % curated, automatic
                        'ListForCofactorPairs',               'curated'                              ,... % curated, automatic
                        'SelectedSubsystems',                 {{'Citric Acid Cycle';
                                                                'Pentose Phosphate Pathway';
                                                                'Glycolysis/Gluconeogenesis';
                                                                'Pyruvate Metabolism';
                                                                'Glyoxylate Metabolism'}}            ,... % default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'x';'y';'z'}}
                        'AddETCAsSubsystem',                  'yes'                                  ,... % yes, no
                        'AddExtracellularSubsystem',          {{'DM_succ_e';
                                                                'DM_ac_e';
                                                                'DM_etoh_e';
                                                                'DM_glyc_e';
                                                                'DM_lac-D_e';
                                                                'DM_akg_e';
                                                                'DM_for_e';
                                                                'DM_pyr_e';
                                                                'DM_fum_e';
                                                                'DM_co2_e';
                                                                'DM_mal-L_e'}}                       ,... % no, default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'DM_x';'DM_y';'DM_z'}}
                        'AerobicAnaerobic',                   'aerobic'                              ,... % aerobic, anaerobic
                        'AlignTransportsUsingMatFile',        'yesautomatic'                         ,... % yesusingmatfile, yesusingReducedmatfile, yesautomatic,no
                        'TimeLimitForSolver',                 'yes'                                  ,... % yes, no
                        'RemovePeriplasm',                    'no'                                   ,... % yes, no
                        'OnlyConnectExclusiveMets',           'no'                                   ,... % yes, no
                        'ConnectIntracellularSubsystems',     'yes'                                  ,... % yes, no
                        'ApplyShortestDistanceOfSubsystems',  'bothways'                             ,... % bothways,eachdirection'
                        'RedModelName',                       'EcoliD1Smin'                          ,... % choose a name
                        'ThrowErrorOnDViolation',             'error'                                ,... % error/continue If any two subsystems cannont connect to the level D                                                            
                        'performLUMPING',                     'yes'                                  ,... % do we want to perform lumping or not?
                        'PercentOfmuMaxForLumping',            100                                   ,... % Please specify the percentage of muMax that we should impose for the lumped reactions (100, 90, ..., 10, 0)
                        'PreventBBBuptake',                   'no'                                   ,... % yes/no: if yes, do not allow uptake through any bbb drains.                        
                        'ZeroZeroGEMbounds',                  'Original'                             ,... % Original, DefineCustom, OpenTo100      
                        'ImposeThermodynamics',               'no'                                   ,... % would you loke to impose thermodynamic cnstraints? yes, no
                        'NumOfLumped',                        'Smin'                                 ,... % OnePerBBB, Smin, Sminp1, Sminp2, Sminp3
                        'CplexParameters',                    'LCSBDefault'                          ,... % CplexDefault, LCSBDefault
                        'performPostProcessing',              'PP_removeBlockedRxns'                 ,... % No, PP_removeBlockedRxns
                        'redGEMlumpGEMpath',                  redGEMlumpGEMpath                      ,... % Provide the full path of your redGEMlumpGEMpath folder
                        'OUTPUTpath',                         OUTPUTpath                             ,... % Provide the full output path
                        'CPLEX_PATH',                         CPLEX_PATH);
cd(paramEcoli.redGEMlumpGEMpath);
redGEM(paramEcoli);
 
end
