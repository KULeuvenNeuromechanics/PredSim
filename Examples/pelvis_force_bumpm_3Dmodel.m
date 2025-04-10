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
S.subject.name = 'Falisse_et_al_2022';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathExDir,'ExampleResults','PelvisForce_Bumpm3D');  

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
S.OpenSimADOptions.verbose_mode = true;

%% Add external force at pelvis

% proportional to muscle activity
fext.function_name = 'Fext_pelvis_bumpm';
fext.extForce = [100, 0, 0]'; % 100 N in x-direction
fext.r_origin = [0, 0, 0]';

%% Run predictive simulations

F_pelvis = 0:10:100;
for i = 1:length(F_pelvis)
    fext.extForce = [F_pelvis(i), 0, 0]';
    S.orthosis.settings{1} = fext;
    [savename] = runPredSim(S, osim_path);
end
