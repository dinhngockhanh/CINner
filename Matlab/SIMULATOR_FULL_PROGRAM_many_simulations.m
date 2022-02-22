%========================CREATE MANY SIMULATIONS WITH DEFAULT PARAMETERS
function  SIMULATOR_FULL_PROGRAM_many_simulations(model,stage_final,N_simulations)
    for i_simulation=1:N_simulations
        disp(i_simulation)
        SIMULATOR_FULL_PROGRAM_one_simulation(model,stage_final);
    end
end
