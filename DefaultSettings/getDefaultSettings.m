function [S] = getDefaultSettings(S,model_info)
% --------------------------------------------------------------------------
%getDefaultSettings 
%   This functions sets default settings when the user didn't specify the
%   setting in main.m.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - IO -
%   * structure with all the model information
% 
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: FULL NAME
% Original date: DD/MM/YYYY
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 30/11/2021
% --------------------------------------------------------------------------

%% bounds

% minimal muscle activation
if ~isfield(S.bounds.a,'lower')
    S.bounds.a.lower = 0;
end

% minimal distance between orginins calcanei
if ~isfield(S.bounds.calcn_dist,'lower')
    S.bounds.calcn_dist.lower = 0.09;
end

% minimal distance between origins toes
if ~isfield(S.bounds.toes_dist,'lower')
    S.bounds.toes_dist.lower = 0.10;
end

% minimal distance between origins tibiae
if ~isfield(S.bounds.tibia_dist,'lower')
    S.bounds.tibia_dist.lower = 0.11;
end

% upper bound on left step length
if ~isfield(S.bounds.SLL,'upper')
    S.bounds.SLL.upper = [];
end

% upper bound on right step length
if ~isfield(S.bounds.SLR,'upper')
    S.bounds.SLR.upper = [];
end

% lower bound on distance travelled
if ~isfield(S.bounds.dist_trav,'lower')
    S.bounds.dist_trav.lower = [];
end

% upper bound on final time
if ~isfield(S.bounds.t_final,'upper')
    S.bounds.t_final.upper = [];
end

% lower bound on final time
if ~isfield(S.bounds.t_final,'lower')
    S.bounds.t_final.lower = [];
end


%% Initial Guess (IG)

% input is a string; "quasi-random" or the path to a .mot file
% check whether one of the two is inputted, otherwise give an error

if ~isfield(S.IG,'selection')
    error('A. Please specify what you want to use as an initial guess. Either choose "quasi-random" or specify an input .mot file.')
else
    [~,NAME,EXT] = fileparts(S.IG.selection);
    if EXT == ".mot" && isfile(S.IG.selection)
        disp(['Using ',char(S.IG.selection), ' as initial guess.'])
        
    elseif EXT == ".mot" && ~isfile(S.IG.selection)
        error('The motion file you specified does not exist.')
        
    elseif EXT == "" && NAME == "quasi-random"
         disp(['Using a ',char(S.IG.selection), ' guess as initial guess.'])
         
    elseif EXT == "" && NAME ~= "quasi-random"
        error('Please specify what you want to use as an initial guess. Either choose "quasi-random" or specify an input .mot file.')
    end
end


%% metabolicE

% hyperbolic tangent smoothing factor (used in metabolic cost)
if ~isfield(S.metabolicE,'tanh_b')
    S.metabolicE.tanh_b = 100;
end

% name of the metabolic energy model --> still to do!
if ~isfield(S.metabolicE,'model')
    S.metabolicE.model = ;
end

%% misc

% maximal contraction velocity identifier --> does still have to be in
% here?
if ~isfield(S.misc,'v_max_s')
    S.metabolicE.tanh_b = 100;
end

% type of gait simulation
if ~isfield(S.misc,'gaitmotion_type')
    S.misc.gaitmotion_type = 'HalfGaitCycle';
end

%% post_process

% boolean for making plots or not
if ~isfield(S.post_process,'make_plot')
    S.post_process.make_plot = [];
end

% name used for saving the results (choose custom or structurized savename)
if ~isfield(S.post_process,'savename')
    S.post_process.savename = 'structured';
end

%% solver

% solver algorithm used in the OCP
if ~isfield(S.solver,'linear_solver')
    S.solver.linear_solver = 'mumps';
end

% the power (10^-x) of the dual infeasibility for when the problem is 
% ‘solved’, higher number is more precise
if ~isfield(S.solver,'tol_ipopt')
    S.solver.ipopt = '4';
end

% the amount of iterations after which the solver stops
if ~isfield(S.solver,'max_iter')
    S.solver.max_iter = '10000';
end

% type of parellel computing
if ~isfield(S.solver,'parallel_mode')
    S.solver.parallel_mode = 'thread';
end

% number of threads in parallel mode
if ~isfield(S.solver,'N_threads')
    S.solver.N_threads = '4';
end

% number of mesh intervals
if ~isfield(S.solver,'N_meshes')
    S.solver.N_meshes = '50';
end

%% subject

% folder to store the subject specific results
if ~isfield(S.subject,'save_folder')
   S.subject.save_folder = []; 
end

% name of the subject, compare with opensim model
if ~isfield(S.subject,'name')
    error('Please provide a name for this subject. This name will be used to store the results. Specify the name in S.subject.name');
end

% mass of the subject
if ~isfield(S.subject,'mass')
   S.subject.mass = model_info.mass;
elseif S.subject.mass ~= model_info.mass
    message = sprintf(['The mass that you specified is different from the mass of the osim model. \n',...
        'The mass of the osim model: ', num2str(model_info.mass),'\n',...
        'The mass you specified: ', num2str(S.subject.mass)]);
    warning(message);
end

% height of the pelvis for the initial guess, compare with opensim model
if ~isfield(S.subject,'IG_pelvis_y')
   S.subject.IG_pelvis_y = model_info.pelvis_y;
elseif S.subject.IG_pelvis_y ~= model_info.pelvis_y
    message = sprintf(['The height that you specified for the pelvis is different from the height in the osim model. \n',...
        'The pelvis height in the osim model: ', num2str(model_info.pelvis_y),'\n',...
        'The pelvis height you specified: ', num2str(S.subject.IG_pelvis_y)]);
    warning(message);
end

% average velocity you want the model to have
if ~isfield(S.subject,'v_pelvis_x_trgt')
    S.subject.v_pelvis_x_trgt = 1.25;
end

% muscle strength, check with opensim model if muscles are present


% muscle stiffness, check with opensim model if muscles are present


% muscle symmetry, check with opensim model if muscles are present

% tendon stiffness, check with opensim model if muscles are present

% motion file for the initial guess

% type of mtp joint used in the model

% folder where you store all subject specific output

% muscle tendon properties


% will be added in the future:
% - spasticity
% - muscle coordination (co-contraction)


%% weights

% weight on metabolic energy rate

% exponent for the metabolic energy rate

% weight on joint accelerations

% weight on arm excitations

% weight on passive torques

% weight on activations

% weight on mtp excitations

% weight on slack controls



end

