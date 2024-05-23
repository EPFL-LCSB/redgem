function getEcoli_D1
cd('./../')  %check where it is saved
paramEcoli    = struct (...%model parameters
                        'Organism',                           'ecoli'                                ,... % human, ecoli, putida
                        'GEMname',                            'iJO1366'                              ,... % name of the GEM used (as it is saved in GEMs folder)
                        'RedModelName',                       'EcoliD1Smin         '                 ,... % choose a name for the reduction
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
                        'ListForInorganicMets',               'automatic'                            ,... % curated, automatic
                        'ListForCofactorPairs',               'curated'                              ,... % curated, automatic
                        'ZeroZeroGEMbounds',                  'Original'                             ,... % Original, DefineCustom, OpenTo100      
                        'case_filename',                      'case_ecoli_iJO1366'                           ,... % name of the matlab function for the specific organism      
                        ...%redGEM parameters
                        'L',                                  6                                      ,... % 1,2,...
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
                        'addGAM',                             'yes'                                  ,... % Would you like to extract and add to the reduced model a growth associated maintenance (GAM) reaction? yes, no 
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
                        'CPLEX_PATH',                         'PATH/TO/CPLEX'                        ,... %provide path to CPLEX         
                        'TFA_PATH',                           'PATH/TO/TFA'                          ,... %provide path to matTFA
                        'thermo_data_PATH',                   'PATH/TO/THERMODATA/thermo_data.mat'   ,... % Provide the path to the thermodynamic data (including the file name)
                        'output_PATH',                        'OUPUT/PATH');                              % provide output path where the files will be saved
redGEM(paramEcoli);

end
