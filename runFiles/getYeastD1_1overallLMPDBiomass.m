function getYeastD1_1overallLMPDBiomass

currentmfilenamepath = mfilename('fullpath');

id_slash = regexpi(currentmfilenamepath, '/');
cd(currentmfilenamepath(1:id_slash(end)));

cplex_PATH                 = '/Users/georgiosfengos/Applications/IBM/ILOG/CPLEX_Studio1271';
mattfa_PATH                = '/Users/georgiosfengos/GIT_Folders/mattfa';
ThermodynamicDatabase_PATH = '/Users/georgiosfengos/GIT_Folders/mattfa/thermoDatabases/thermo_data.mat';
output_PATH                = '/Users/georgiosfengos/TEMP/redGEMlumpGEMOUTPUT';

cd('./../')  %check where it is saved
paramRedGEM    = struct (...%model parameters
                        'Organism',                           'yeast'                                ,... % human, ecoli, putida
                        'GEMname',                            'iMM904bmDrain'                        ,... % file name of the GEM used (as it is saved in GEMs folder)
                        'RedModelName',                       'YeastBMDrainD1Smin1overallLmpdBM'     ,... % choose a name for the reduction
                        'SelectedSubsystems',                 {{'Citric Acid Cycle';
                                                                'Pentose Phosphate Pathway';
                                                                'Glycolysis/Gluconeogenesis';
                                                                'Pyruvate Metabolism';
                                                                'Glyoxylate Metabolism'}}            ,... % default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'x';'y';'z'}}
                        'AddETCAsSubsystem',                  'yes'                                  ,... % yes, no
                        'AddExtracellularSubsystem',          {{'EX_succ(e)';
                                                                'EX_ac(e)';
                                                                'EX_etoh(e)';
                                                                'EX_glyc(e)';
                                                                'EX_lac-D(e)';
                                                                'EX_akg(e)';
                                                                'EX_for(e)';
                                                                'EX_pyr(e)';
                                                                'EX_co2(e)';
                                                                'EX_mal-L(e)'}}                       ,... % no, default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'DM_x';'DM_y';'DM_z'}}
                        'AerobicAnaerobic',                   'aerobic'                              ,... % aerobic, anaerobic
                        'ListForInorganicMets',               'automatic'                            ,... % curated, automatic
                        'ListForCofactorPairs',               'curated'                              ,... % curated, automatic
                        'ZeroZeroGEMbounds',                  'Original'                             ,... % Original, DefineCustom, OpenTo100      
                        'FluxUnits',                          'mmol'                                 ,... % mmol, mumol, other      
                        ...%redGEM parameters
                        'L',                                  4                                      ,... % 1,2,...
                        'D',                                  1                                      ,... % 1,2,...
                        'startFromMin',                       'no'                                   ,... % yes, no
                        'ThrowErrorOnDViolation',             'error'                                ,... % error/continue If any two subsystems cannont connect to the level D                           
                        'OnlyConnectExclusiveMets',           'no'                                   ,... % yes, no
                        'ConnectIntracellularSubsystems',     'yes'                                  ,... % yes, no
                        'ApplyShortestDistanceOfSubsystems',  'bothways'                             ,... % bothways,eachdirection'
                        ...%redGEMX parameters 
                        'performREDGEMX',                     'no'                                   ,... % yes, no
                        'NumOfConnections',                   'OnePerMetE'                           ,... % OnePerMetE ,SminMetE
                        ...%lumpGEM parameters
                        'performLUMPGEM',                     'yes'                                  ,... % yes, no (do we want to perform lumping or not?)   
                        'PercentOfmuMaxForLumping',            100                                   ,... % Please specify the percentage of muMax that we should impose for the lumped reactions (100, 90, ..., 10, 0)
                        'addGAM',                             'no'                                   ,... % Would you like to extract and add to the reduced model a growth associated maintenance (GAM) reaction? yes, no 
                        'PreventBBBuptake',                   'no'                                   ,... % yes/no: if yes, do not allow uptake through any bbb drains.     
                        'NumOfLumped',                        'Smin'                                 ,... % OnePerBBB, Smin, Sminp1, Sminp2, Sminp3
                        'AlignTransportsUsingMatFile',        'yesautomatic'                         ,... % yesusingmatfile, yesusingReducedmatfile, yesautomatic,no
                        'ImposeThermodynamics',               'yes'                                  ,... % would you like to impose thermodynamic cnstraints? yes, no
                        ...%postprocessing parameters
                        'performPostProcessing',              'yes'                                  ,... % yes, no
                        ...%solver parameters
                        'TimeLimitForSolver',                 'yes'                                  ,... % yes, no
                        'CplexParameters',                    'LCSBDefault'                          ,... % CplexDefault, LCSBDefault
                        ....%paths to folders
                        'CPLEX_PATH',                         cplex_PATH                             ,... %provide path to CPLEX         
                        'TFA_PATH',                           mattfa_PATH                            ,... %provide path to matTFA
                        'thermo_data_PATH',                   ThermodynamicDatabase_PATH             ,... % Provide the path to the thermodynamic data (including the file name)
                        'output_PATH',                        output_PATH                            );   % Provide the full output path
redGEM(paramRedGEM);

end
