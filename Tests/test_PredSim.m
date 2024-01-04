
clear
close all
clc

%%
% name of the reference model used for testing
reference_name = 'Falisse_et_al_2022';

% path to folder with program to create dll files from opensim model (will
% be downloaded automatically if it is not there)
PathCpp2Dll_Exe = 'C:\GBW_MyPrograms\Osim2Dll_exe';

% Do you want to run the simulation as a batch job (parallel computing toolbox)
run_as_batch_job = 0;

% Path to folder with casadi
CasADi_path = 'C:\GBW_MyPrograms\casadi_3_5_5';

% Value to be passed to S.Cpp2Dll.compiler. Most users can leave this empty.
Cpp2Dll_compiler = [];

%% Configure settings for simulation - Do not change!

[pathTestDir,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathTestDir);

cd(pathRepo);

pathDefaultSettings = fullfile(pathRepo,'DefaultSettings');
addpath(pathDefaultSettings)

% Initialize S
[S] = initializeSettings(reference_name);
S.misc.main_path = pathRepo;

addpath(fullfile(S.misc.main_path,'VariousFunctions'))

% Store result of test simulation in Debug folder
if ~isfolder(fullfile(S.misc.main_path,'Debug'))
    mkdir(fullfile(S.misc.main_path,'Debug'));
end
if ~isfolder(fullfile(S.misc.main_path,'Debug','TestResults'))
    mkdir(fullfile(S.misc.main_path,'Debug','TestResults'));
end
S.subject.save_folder  = fullfile(S.misc.main_path,'Debug','TestResults'); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.subject.IG_selection_gaitCyclePercent = 100;
% S.subject.IG_selection = 'quasi-random';

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% path to folder with program to create dll files from opensim model (will
% be downloaded automatically if it is not there)
S.Cpp2Dll.PathCpp2Dll_Exe = PathCpp2Dll_Exe;

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = run_as_batch_job;

S.solver.CasADi_path = CasADi_path;

if exist('Cpp2Dll_compiler','var') && ~isempty(Cpp2Dll_compiler)
    S.Cpp2Dll.compiler = Cpp2Dll_compiler;
end


result_filename = [reference_name '_test_' datestr(datetime,30)];
S.post_process.result_filename = result_filename;

S.Cpp2Dll.PathCpp2Dll_Exe = InstallOsim2Dll_Exe(S.Cpp2Dll.PathCpp2Dll_Exe);
S.Cpp2Dll.verbose_mode = 0; % 0 for no outputs from cmake


% Start simulation
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    run_pred_sim(S,osim_path);
end

cd(pathTestDir);

%% compare simulation results

R_ref = load(fullfile(pathTestDir,'ReferenceResults',reference_name,[reference_name,'_v1.mat']),'R');
R_ref = R_ref.R;
R_test = load(fullfile(S.misc.main_path,'Debug','TestResults',[result_filename,'.mat']),'R');
R_test = R_test.R;

[max_diff] = maxdiff(R_ref, R_test, {'kinematics.Qs_rad','kinematics.Qdots_rad','kinematics.Qddots_rad',...
    'muscles.a','muscles.da','muscles.FTtilde','muscles.dFTtilde',...
    'torque_actuators.a','torque_actuators.e','objective.absoluteValues'});


%% compare solver printout (logfile)

%%
function [max_diff] = maxdiff(R_ref, R_test, fields)
    max_diff = nan(length(fields),1);
    for i=1:length(fields)
        idx_dot = strfind(fields{i},'.');
        fieldname1 = fields{i}(1:idx_dot-1);
        fieldname2 = fields{i}(idx_dot+1:end);
        max_diff(i) = max(abs( R_test.(fieldname1).(fieldname2) - R_ref.(fieldname1).(fieldname2) ), [], 'all');
    end
end