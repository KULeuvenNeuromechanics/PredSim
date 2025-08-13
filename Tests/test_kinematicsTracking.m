% --------------------------------------------------------------------------
% test_PredSim
%   Run tests to ensure kinematics tracking works correctly.
%
%   Reference kinematics data taken from Van Criekinge et al., 2023
%       Van Criekinge, T., Saeys, W., Truijen, S. et al. A full-body motion 
%       capture gait dataset of 138 able-bodied adults across the life span 
%       and 50 stroke survivors. Sci Data 10, 852 (2023). 
%       https://doi.org/10.1038/s41597-023-02767-y
%
% Original author: Menthy Denayer
% Original date: 08 June 2025
% --------------------------------------------------------------------------

clear
close all
clc
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

cd(pathRepoFolder);

%% Parametersp
model_name = 'gait1018';

%% Initialize S
addpath(fullfile(pathRepoFolder,'DefaultSettings'))

[S] = initializeSettings(model_name);

%% Settings

% CasADi path
S.solver.CasADi_path = "";

% name of the subject
S.subject.name = model_name;

% tracking data
S.subject.TrackSim = true;
S.subject.TrackingFile = fullfile(S.misc.main_path,'Tests','vanCriekinge_referenceCoordinates.mot'); % mot file used for tracking 
% S.subject.TrackingJoints = {'knee_angle_r', 'knee_angle_l'}; % list of joints to track
S.subject.TrackingJoints = 'all';
S.weights.kinematicsTracking = 10;                  
S.misc.gaitmotion_type = 'FullGaitCycle';                                   % required for selected mot file

% set weights to zero to test out tracking
S.weights.E = 0;
S.weights.E_exp = 0;
S.weights.q_dotdot = 0;
S.weights.e_torqAct = 0;
S.weights.pass_torq = 0;
S.weights.pass_torq_includes_damping = 0;
S.weights.a = 0;
S.weights.a_exp = 0;
S.weights.slack_ctrl = 0;

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',[S.subject.name '_tracking']); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepoFolder,'Subjects',S.subject.name,[S.subject.name '.osim']);

% choose speed
S.misc.forward_velocity = 1.2;

% solver settings
S.solver.run_as_batch_job = false;

%% Run predictive simulations
[savename] = runPredSim(S, osim_path);

%% Plot Tracking Results
tracking_result = load(fullfile(S.misc.save_folder,savename));

Qref_tot = getIK(S.subject.TrackingFile,tracking_result.model_info);        % load .mot file with tracking data
Qref_time = Qref_tot.time;                                                  % extract tracking data time
Nmeshes = tracking_result.R.S.solver.N_meshes;                              % define number of meshes

% check joints to track
if(strcmp(S.subject.TrackingJoints,'all'))       
    % if tracking all joints
    desir_coo_names = string(fieldnames(tracking_result.model_info.ExtFunIO.coordi));    
else
    % if tracking only selected joints
    desir_coo_names = string(S.subject.TrackingJoints);                     
end

% create reference data
NtrackJoints = length(desir_coo_names);                                     % number of joints to track
Ndata = length(Qref_time);                                                  % size of the experimental data
Qref = zeros(Ndata,NtrackJoints);                                           % matrix to store tracking data  
for jointIdx = 1:NtrackJoints
    Qref(:,jointIdx) = Qref_tot.(desir_coo_names(jointIdx));                % fill matrix with data
end

% resample to be ( Nmeshes ) x ( number of joints to track )
Qrefsync = interp1(linspace(1,Nmeshes,Ndata),Qref,linspace(1,Nmeshes,Nmeshes),'spline','extrap');

isLinearCoo = contains(desir_coo_names,'tx') | ...                          % indices with non-rotational coordinates
    contains(desir_coo_names,'ty');
Qrefsync(:,~isLinearCoo) = Qrefsync(:,~isLinearCoo)*180/pi;                 % correct dimensions

% extract predicted kinematics
Qpred = tracking_result.R.kinematics.Qs;

isPelvisTx = contains(desir_coo_names,'tx');                                % index for pelvis translation
Qrefsync(:,~isPelvisTx) = circshift(Qrefsync(:,~isPelvisTx), ...            % shift predicted data to match tracking data
    -tracking_result.R.ground_reaction.idx_GC(end));

% plot results
ylabelList = repelem("coordinate value [ deg ]", NtrackJoints, 1);
ylabelList(isLinearCoo) = "coordinate value [ m ]";

for i = 1:NtrackJoints
    figure
    hold on
    grid on
    plot(linspace(1,Nmeshes,Nmeshes), Qrefsync(:,i), 'k--', 'LineWidth', 1)
    plot(linspace(1,Nmeshes,Nmeshes), Qpred(:,i), 'r', 'LineWidth', 1)
    title(replace(desir_coo_names(i),"_"," "))
    ylabel(ylabelList(i), "FontWeight", "bold")
    xlabel("Gait Cycle [ % ]", "FontWeight", "bold")
    legend(["Tracking Data", "Predicted Data"],"Location","best")
    hold off
end