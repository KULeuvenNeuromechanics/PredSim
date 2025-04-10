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

[S] = initializeSettings('gait1018');

%% Settings

% name of the subject
S.subject.name = 'gait1018';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathExDir,'ExampleResults','PelvisForce_Bumpm2D');  

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
S.OpenSimADOptions.verbose_mode = true;

%% Add external force at pelvis

% gait cycle independent force
%fext.function_name = 'Fext_pelvis';

% gait cycle dependent force
%fext.function_name = 'Fext_pelvis_dt';

% run for a range with forces

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





% 
% %% plot results
% R = load(fullfile(S.misc.save_folder,[savename '.mat']));
% 
% % plot results
% result_paths{1} = fullfile(S.misc.save_folder, 'gait1018_v1.mat');
% result_paths{2} = fullfile(S.misc.save_folder, 'gait1018_v2.mat');
% result_paths{3} = fullfile(S.misc.save_folder, 'gait1018_v3.mat');
% result_paths{4} = fullfile(S.misc.save_folder, 'gait1018_v4.mat');
% result_paths{5} = fullfile(S.misc.save_folder, 'gait1018_v5.mat');
% result_paths{6} = fullfile(S.misc.save_folder, 'gait1018_v6.mat');
% 
% legend_names = {'damper','motor','sin','emg-sol','emg-tibant','defaultsim'};
% figure_settings(1).name = 'all_angles';
% figure_settings(1).dofs = {'all_coords'};
% figure_settings(1).variables = {'Qs'};
% figure_settings(1).savepath = [];
% figure_settings(1).filetype = {};
% 
% figure_settings(2).name = 'all_activations';
% figure_settings(2).dofs = {'muscles_r'};
% figure_settings(2).variables = {'a'};
% figure_settings(2).savepath = [];
% figure_settings(2).filetype = {};
% 
% % call plotting function
% plot_figures(result_paths, legend_names, figure_settings);
% 
% % plot metabolic energy consumption
% figure();
% for i =1:length(result_paths)
%     R = load(result_paths{i});
%     bar(i, R.R.metabolics.Bhargava2004.COT); hold on;
%     legend(legend_names)
% end