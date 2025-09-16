function [data_table, ct_sim] = add_benchmark_to_table(ct_sim,sim_res_folder, ...
    data_headers, data_table, slope, id_study, Lsim, msim)
%add_benchmark_to_table reads results of simulation and associated
%experiments and add it to a table with outcomes


[R, benchmark] = load_sim_file(sim_res_folder);
g=9.81;

% add to data table -- frequency
data_table(ct_sim,strcmp(data_headers,'sim_stride_frequency')) = ...
    R.spatiotemp.stride_freq;
if ~isempty(benchmark.stride_frequency)
    data_table(ct_sim,strcmp(data_headers,'exp_stride_frequency')) = ...
        benchmark.stride_frequency*(sqrt(g/Lsim));
end

% add to data table -- metabolic power
t = R.time.mesh_GC;
dt = t(end)-t(1);
Pmetab = R.metabolics.Bhargava2004.Edot_gait;
metab_work  = trapz(t(1:end-1)',Pmetab);
P_mean = sum(metab_work)./dt;
data_table(ct_sim,strcmp(data_headers,'sim_metabolic_power')) = P_mean;
if ~isempty(benchmark.Pmetab_mean)
    data_table(ct_sim,strcmp(data_headers,'exp_metabolic_power')) = ...
        benchmark.Pmetab_mean * (msim*sqrt(Lsim)*g^1.5);
end

% study information
data_table(ct_sim,strcmp(data_headers,'id_study')) = id_study;
data_table(ct_sim,strcmp(data_headers,'speed')) = R.S.misc.forward_velocity;
data_table(ct_sim,strcmp(data_headers,'slope')) = slope;
ct_sim = ct_sim +1;
end