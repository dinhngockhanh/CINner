%=====================================CLEAN VECTOR OF CLONAL POPULATIONS
function SIMULATOR_FULL_PHASE_1_clonal_population_cleaning()
    global clonal_population_current clonal_population_next clonal_ID_current
    vec_delete                              = find(clonal_population_next==0);
    clonal_population_current(vec_delete)   = [];
    clonal_population_next(vec_delete)      = [];
    clonal_ID_current(vec_delete)           = [];
end
