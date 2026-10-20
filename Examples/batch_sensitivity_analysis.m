% --------------------------------------------------------------------------
% batch_sensitivity_analysis
%   Example on running a parameter sweep for sensitivity analysis.
%   Requires MATLAB's parallel computing toolbox.
%   In this example, the parameter of interest is the number of mesh
%   intervals.
%
% Original author: Lars D'Hondt
% Original date: 16/September/2024
% --------------------------------------------------------------------------

clear
close all
clc

[pathExDir,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathExDir);
[pathRepoFolder,~,~] = fileparts(pathRepo);

addpath(fullfile(pathRepo,'DefaultSettings'))
addpath(pathRepo)

%% Initialize S

[S] = initializeSettings('gait1018');

%% Settings

% name of the subject
S.subject.name = 'gait1018';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathExDir,'ExampleResults','BatchMesh'); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Run simulations as batch jobs, such that multiple simulations can run at
% the same time.
S.solver.run_as_batch_job = true;


% loop over different values for the parameter of interest
for N_meshes = [40,50,60,75,100]

    % assign ione value to the appropriate setting
    S.solver.N_meshes = N_meshes;

    % run a simulation for this value
    [savename] = runPredSim(S, osim_path);

end



