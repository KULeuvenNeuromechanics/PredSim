% example on external force at pelvis
%--------------------------------------

% simple example of a simulation with an external force at the pelvis
% this example is mainly used to learn how to work with the orthosis class
% the final goal is to also include an option to optimize controls of this
% force

%% path information
clear
close all
clc

[pathExDir,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathExDir);
[pathRepoFolder,~,~] = fileparts(pathRepo);

addpath(fullfile(pathRepo,'DefaultSettings'))
addpath(pathRepo)

%% Initialize S

[S] = initializeSettings('Falisse_et_al_2022');

%% Settings

% name of the subject
S.subject.name = 'Falisse_et_al_2022_Ftorso';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathExDir,'ExampleResults','TorsoForce_v2');  

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
% S.solver.IG_selection_gaitCyclePercent = 100;
S.solver.IG_selection = 'quasi-random';

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
S.OpenSimADOptions.verbose_mode = true;

S.bounds.activation_all_muscles.lower  = 0.01;

% run as a batch job
S.solver.run_as_batch_job = 1;
S.solver.N_threads      = 2;
S.solver.N_meshes       = 50;
S.solver.par_cluster_name = 'Cores3';

S.solver.CasADi_path = 'C:\Users\Maarten\Documents\Software\downloads\casadi_363';

%% Add external force at pelvis

% % proportional to muscle activity
% S.solver.run_as_batch_job = 0;
fext.function_name = 'Fext_torso_bumpm';
% fext.extForce = [100, 0, 0]'; % 100 N in x-direction
fext.r_origin = [0, 0.2, 0]';
% S.orthosis.settings{1} = fext;
% [savename] = runPredSim(S, osim_path);

%% Run predictive simulations

S.misc.save_folder  = fullfile(pathExDir,'ExampleResults',...
    'TorsoForce_Bumpm3D_v2');
S.misc.forward_velocity   = 1.25;
F_torso = 0:10:120;
fext.function_name = 'Fext_torso_bumpm';
fext.r_origin = [0, 0.2, 0]';
for i = 1:length(F_torso)
    fext.extForce = [F_torso(i), 0, 0]';
    S.orthosis.settings{1} = fext;
    [savename] = runPredSim(S, osim_path);
end
