    clear;
%-----------------------ADD PATH FOR MODEL VARIABLE FILES FROM VIGNETTES
    current_folder  = pwd;
    idcs            = strfind(current_folder,'/');
    mother_folder   = current_folder(1:idcs(end)-1);
    R_folder        = [mother_folder '/vignettes'];
    path(path,R_folder);
%----------------------------TEST THE COPY-NUMBER ALLELE TRACKING SYSTEM
    model           = 'FALLOPIAN-TUBES-NEUTRAL';
    stage_final     = 2;
    N_simulations   = 1;
tic
    TEST_SIMULATIONS_FOR_FIRST_MANUSCRIPT(model,stage_final,N_simulations);
toc




%=========================================SIMULATIONS FOR 1ST MANUSCRIPT
function TEST_SIMULATIONS_FOR_FIRST_MANUSCRIPT(model,stage_final,N_simulations)
% %-------------------------Create folder to store simulation output files
%     folder_name                 = [model '_MATLAB'];
%     if ~isfolder(folder_name)
% %       If the folder with same name doesn't already exist, then create it
%         mkdir(folder_name)
%     end
%-----------------------------------------------------Create simulations
    for i_simulation=1:N_simulations
%-------Create one simulation
        fprintf('Simulation %d/%d\n',i_simulation,N_simulations);
%       Set up the variables for simulations
        SIMULATOR_VARIABLES_for_simulation(model);
%       Create one simulation
        package_output                      = SIMULATOR_FULL_PROGRAM_one_simulation(model,stage_final);
        if stage_final==1
            package_clonal_evolution        = package_output{1};
        end
    end

end
