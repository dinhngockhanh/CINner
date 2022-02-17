%==========================CREATE MANY SIMULATIONS WITH INPUT PARAMETERS
function  ABC_many_simulations(model,parameter_set,N_simulations)
    for i_simulation=1:N_simulations
        disp(i_simulation)
        ABC_one_simulation(model,parameter_set)
    end
end
