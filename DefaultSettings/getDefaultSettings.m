function [S] = getDefaultSettings(S)
% --------------------------------------------------------------------------
%getDefaultSettings 
%   This functions checks if provided settings are valid and sets default 
%   settings when the user didn't specify them in main.m.
% 
% INPUT:
%   - S -
%   * setting structure S
% 
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: Bram Van Den Bosch
% Original date: 30/11/2021
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 21/01/2022
% --------------------------------------------------------------------------

%% bounds

% minimal muscle activation, a number between 0 and 1
if ~isfield(S.bounds.a,'lower')
    S.bounds.a.lower = 0;
end

% minimal distance between orginins calcanei, in meters
if ~isfield(S.bounds.calcn_dist,'lower')
    S.bounds.calcn_dist.lower = 0.09;
end

% minimal distance between origins toes, in meters
if ~isfield(S.bounds.toes_dist,'lower')
    S.bounds.toes_dist.lower = 0.10;
end

% minimal distance between origins tibiae, in meters
if ~isfield(S.bounds.tibia_dist,'lower')
    S.bounds.tibia_dist.lower = 0.11;
end

% upper bound on left step length, in meters
if ~isfield(S.bounds.SLL,'upper')
    S.bounds.SLL.upper = [];
end

% upper bound on right step length, in meters
if ~isfield(S.bounds.SLR,'upper')
    S.bounds.SLR.upper = [];
end

% lower bound on distance travelled, in meters
if ~isfield(S.bounds.dist_trav,'lower')
    S.bounds.dist_trav.lower = [];
end

% upper bound on final time, in seconds
if ~isfield(S.bounds.t_final,'upper')
    S.bounds.t_final.upper = [];
end

% lower bound on final time, in seconds
if ~isfield(S.bounds.t_final,'lower')
    S.bounds.t_final.lower = [];
end

%% metabolicE

% hyperbolic tangent smoothing factor (used in metabolic cost)
if ~isfield(S.metabolicE,'tanh_b')
    S.metabolicE.tanh_b = 100;
end

% name of the metabolic energy model used
if ~isfield(S.metabolicE,'model')
    S.metabolicE.model = 'Bhargava2004';
end

%% misc

% maximal contraction velocity identifier --TO CHECK--
if ~isfield(S.misc,'v_max_s')
    S.misc.v_max_s = 0;
end

% type of gait simulation
if ~isfield(S.misc,'gaitmotion_type')
    S.misc.gaitmotion_type = 'HalfGaitCycle';
end

% main path
if ~isfield(S.misc,'main_path')
    error('Please provide the path to the main folder in S.misc.main_path.');
elseif ~exist(S.misc.main_path)
    error('The main path you provided in S.misc.main_path does not exist.');
end

% type of equation to approximate musculo-skeletal geometry (moment arm and
% muscle-tendon lengths wrt. joint angle)
if ~isfield(S.misc,'msk_geom_eq')
    S.misc.msk_geom_eq = 'polynomials';
end

% minimal order of polynomial function
if ~isfield(S.misc.poly_order,'lower')
    S.misc.poly_order.lower = 3;
end

% maximal order of polynomial function
if ~isfield(S.misc.poly_order,'upper')
    S.misc.poly_order.upper = 9;
end

%% post_process

% boolean to plot post processing results
if ~isfield(S.post_process,'make_plot')
    S.post_process.make_plot = 0;
end

% name used for saving the resultfiles (choose custom or structurized savename)
if ~isfield(S.post_process,'savename')
    S.post_process.savename = 'structured';
end

%% solver

% solver algorithm used in the OCP
if ~isfield(S.solver,'linear_solver')
    S.solver.linear_solver = 'mumps';
end

% the power (10^-x) the dual infeasibility has to reach before the OCP can 
% be regarded as solved; a higher number gives a more precise answer
if ~isfield(S.solver,'tol_ipopt')
    S.solver.tol_ipopt = '4';
end

% maximal amount of itereations after wich the solver will stop
if ~isfield(S.solver,'max_iter')
    S.solver.max_iter = '10000';
end

% type of parallel computing
if ~isfield(S.solver,'parallel_mode')
    S.solver.parallel_mode = 'thread';
end

% ADD CHECK ST THIS IS ONLY  USED WHEN USING THREAD PARALLEL MODE?
% number of threads in parallel mode
if ~isfield(S.solver,'N_threads')
    S.solver.N_threads = '4';
end

% number of mesh intervals
if ~isfield(S.solver,'N_meshes')
    S.solver.N_meshes = '50';
end

%% subject

% name of the subject
if ~isfield(S.subject,'name')
    error('Please provide a name for this subject. This name will be used to store the results. Specify the name in S.subject.name.');
end

% folder path to store the results from the OCP in
if ~isfield(S.subject,'save_results')
    error('Please provide a name for this subject. This name will be used to store the results. Specify the name in S.subject.name.');
elseif ~exist(S.subject.save_results)
    mkdir(S.subject.save_results);
end

% folder path to store the intermediate subject specific results (muscle
% analysis etc.)
S.subject.save_folder = fullfile(S.misc.main_path, "Subjects", S.subject.name);
if ~exist(S.subject.save_folder)
    mkdir(S.subject.save_folder);
end

% mass of the subject, in kilograms
if ~isfield(S.subject,'mass')
   S.subject.mass = [];
end

% height of the pelvis for the initial guess, in meters
if ~isfield(S.subject,'IG_pelvis_y')
   S.subject.IG_pelvis_y = [];
end

% average velocity you want the model to have, in meters per second
if ~isfield(S.subject,'v_pelvis_x_trgt')
    S.subject.v_pelvis_x_trgt = 1.25;
end

% muscle strength
if ~isfield(S.subject,'muscle_strength')
    S.subject.muscle_strength = [];
end

% muscle stiffness
if ~isfield(S.subject,'muscle_stiff')
    S.subject.muscle_stiff = [];
end

% muscle symmetry
if ~isfield(S.subject,'muscle_sym')
    S.subject.muscle_sym = [];
end

% tendon stiffness
if ~isfield(S.subject,'tendon_stiff')
    S.subject.tendon_stiff = [];
end

% initial guess inputs
% input should be a string: "quasi-random" or the path to a .mot file
if ~isfield(S.subject,'IG_selection')
    error('Please specify what you want to use as an initial guess. Either choose "quasi-random" or specify the path of a .mot file in S.subject.IG_selection.')
else
    [~,NAME,EXT] = fileparts(S.subject.IG_selection);
    if EXT == ".mot" && isfile(S.subject.IG_selection)
        disp(['Using ',char(S.subject.IG_selection), ' as initial guess.'])
        
    elseif EXT == ".mot" && ~isfile(S.subject.IG_selection)
        error('The motion file path you specified does not exist. Check if the path exist and if you made a typo.')
        
    elseif EXT == "" && NAME == "quasi-random"
         disp(['Using a quasi-random guess as initial guess.'])
         
    elseif EXT == "" && NAME ~= "quasi-random"
        error('Please specify what you want to use as an initial guess. Either choose "quasi-random" or specify the path of a .mot file.')
    end
end

% initial guess bounds
if ~isfield(S.subject,'IG_bounds')
    error('Please provide a .mot file on which the IG bounds will be based. Specify the file in S.subject.IG_bounds.')
elseif ~isfile(S.subject.IG_bounds)
    error('The motion file you specified in S.subject.IG_bounds does not exist.')
end
disp([char(S.subject.IG_bounds), ' will be used to determine IG bounds.'])

% type of mtp joint used in the model
if ~isfield(S.subject,'mtp_type')
    S.subject.mtp_type = []; 
end

% muscle tendon properties
if ~isfield(S.subject,'MT_params')
    S.subject.MT_params = []; 
end

% muscle spasticity
if ~isfield(S.subject,'spasticity')
    S.subject.spasticity = []; 
end

% % muscle coordination
% if ~isfield(S.subject,'muscle_coordination')
%     S.subject.muscle_coordination = []; 
% end

%% weights

% weight on metabolic energy rate
if ~isfield(S.weights,'E')
    S.weights.E = 500; 
end

% exponent for the metabolic energy rate
if ~isfield(S.weights,'E_exp')
    S.weights.E_exp = 2; 
end

% weight on joint accelerations
if ~isfield(S.weights,'q_dotdot')
    S.weights.q_dotdot = 50000; 
end

% weight on arm excitations
if ~isfield(S.weights,'e_arm')
    S.weights.e_arm = 10^6; 
end

% weight on passive torques
if ~isfield(S.weights,'pass_torq')
    S.weights.pass_torq = 1000; 
end

% weight on muscle activations
if ~isfield(S.weights,'a')
    S.weights.a = 2000; 
end

% weight on mtp excitations
if ~isfield(S.weights,'e_mtp')
    S.weights.e_mtp = 10^6; 
end

% weight on slack controls
if ~isfield(S.weights,'slack_ctrl')
    S.weights.slack_ctrl = 0.001; 
end

end

