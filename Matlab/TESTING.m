    clear;
    %   stage_final     = 1:    outputs the clonal evolution
    %   stage_final     = 2:    also outputs the CN profiles of a sample
    %   stage_final     = 3:    also outputs the phylogeny of the sample
%-----------------------ADD PATH FOR MODEL VARIABLE FILES FROM VIGNETTES
    current_folder  = pwd;
    idcs            = strfind(current_folder,'/');
    mother_folder   = current_folder(1:idcs(end)-1);
    R_folder        = [mother_folder '/vignettes'];
    path(path,R_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEUTRAL FALLOPIAN TUBES - ONE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MODEL           = 'FALLOPIAN-TUBES-NEUTRAL';
    stage_final     = 3;

    package_output  = SIMULATOR_FULL_PROGRAM_one_simulation(MODEL,stage_final);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEUTRAL FALLOPIAN TUBES - MANY SIMULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MODEL           = 'FALLOPIAN-TUBES-NEUTRAL';
    stage_final     = 1;
    N_simulations   = 100;

    SIMULATOR_FULL_PROGRAM_many_simulations(MODEL,stage_final,N_simulations);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SELECTION OF BULK HGSOC - ONE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model                               = 'HGSOC-BULK';
%   Find how many drivers are in the library
    global driver_library
    SIMULATOR_VARIABLES_for_simulation(model)
    N_drivers                           = size(driver_library,1);
%   Set up individual parameters for fitting
    vec_selection_rates                 = 0.05*ones(1,N_drivers);

    rate_driver                         = 1e-16;

    prob_CN_whole_genome_duplication    = 1e-6;
    prob_CN_missegregation              = 1e-6;
    prob_CN_chrom_arm_missegregation    = 1e-6;
    prob_CN_focal_amplification         = 1e-6;
    prob_CN_focal_deletion              = 1e-6;
    prob_CN_cnloh_interstitial          = 1e-6;
    prob_CN_cnloh_terminal              = 1e-6;

    prob_CN_focal_amplification_length  = 0.1;
    prob_CN_focal_deletion_length       = 0.1;
    prob_CN_cnloh_interstitial_length   = 0.1;
    prob_CN_cnloh_terminal_length       = 0.1;
%   Package the parameters into one set
    parameter_set{1}                    = vec_selection_rates;
    parameter_set{2}                    = rate_driver;
    parameter_set{3}                    = prob_CN_whole_genome_duplication;
    parameter_set{4}                    = prob_CN_missegregation;
    parameter_set{5}                    = prob_CN_chrom_arm_missegregation;
    parameter_set{6}                    = prob_CN_focal_amplification;
    parameter_set{7}                    = prob_CN_focal_deletion;
    parameter_set{8}                    = prob_CN_cnloh_interstitial;
    parameter_set{9}                    = prob_CN_cnloh_terminal;
    parameter_set{10}                   = prob_CN_focal_amplification_length;
    parameter_set{11}                   = prob_CN_focal_deletion_length;
    parameter_set{12}                   = prob_CN_cnloh_interstitial_length;
    parameter_set{13}                   = prob_CN_cnloh_terminal_length;
%   Run the ABC model for this parameter set
    ABC_one_simulation(model,parameter_set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SELECTION OF BULK HGSOC - MANY SIMULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model                               = 'HGSOC-BULK';
%   Find how many drivers are in the library
    global driver_library
    SIMULATOR_VARIABLES_for_simulation(model)
    N_drivers                           = size(driver_library,1);
%   Set up individual parameters for fitting
    vec_selection_rates                 = 0.05*ones(1,N_drivers);

    rate_driver                         = 1e-16;

    prob_CN_whole_genome_duplication    = 1e-6;
    prob_CN_missegregation              = 1e-6;
    prob_CN_chrom_arm_missegregation    = 1e-6;
    prob_CN_focal_amplification         = 1e-6;
    prob_CN_focal_deletion              = 1e-6;
    prob_CN_cnloh_interstitial          = 1e-6;
    prob_CN_cnloh_terminal              = 1e-6;

    prob_CN_focal_amplification_length  = 0.1;
    prob_CN_focal_deletion_length       = 0.1;
    prob_CN_cnloh_interstitial_length   = 0.1;
    prob_CN_cnloh_terminal_length       = 0.1;
%   Package the parameters into one set
    parameter_set{1}                    = vec_selection_rates;
    parameter_set{2}                    = rate_driver;
    parameter_set{3}                    = prob_CN_whole_genome_duplication;
    parameter_set{4}                    = prob_CN_missegregation;
    parameter_set{5}                    = prob_CN_chrom_arm_missegregation;
    parameter_set{6}                    = prob_CN_focal_amplification;
    parameter_set{7}                    = prob_CN_focal_deletion;
    parameter_set{8}                    = prob_CN_cnloh_interstitial;
    parameter_set{9}                    = prob_CN_cnloh_terminal;
    parameter_set{10}                   = prob_CN_focal_amplification_length;
    parameter_set{11}                   = prob_CN_focal_deletion_length;
    parameter_set{12}                   = prob_CN_cnloh_interstitial_length;
    parameter_set{13}                   = prob_CN_cnloh_terminal_length;
%   Run the ABC model for this parameter set
    N_simulations                       = 100;

    ABC_many_simulations(model,parameter_set,N_simulations)
