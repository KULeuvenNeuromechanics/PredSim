function [R, benchmark] = load_sim_file(sim_res_folder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mat_files = dir(fullfile(sim_res_folder,'*.mat'));
if length(mat_files) ~= 1
    disp('warning mutiple mat files in folder')
    disp(sim_res_folder);
    disp(['assumes that file ' mat_files(1).name, ...
        'contains the simulation results'])
end
sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
load(sim_res_file,'R','benchmark');

end