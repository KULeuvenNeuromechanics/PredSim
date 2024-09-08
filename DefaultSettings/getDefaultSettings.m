function [S] = getDefaultSettings(S,osim_path)
% --------------------------------------------------------------------------
% getDefaultSettings 
%   This functions sets default settings when the user didn't specify the
%   setting in main.m.
% 
% INPUT:
%   - S -
%   * setting structure S
%   
%   - osim_path -
%   * path to the osim model
%
% 
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: Bram Van Den Bosch
% Original date: 30/11/2021
% --------------------------------------------------------------------------

%% bounds
if ~isfield(S,'bounds')
    S.bounds = [];
end

% upper bound on left step length, in meters
if ~isfield(S.bounds,'SLL')
    S.bounds.SLL        = [];
end
if ~isfield(S.bounds.SLL,'upper')
    S.bounds.SLL.upper = [];
end

% upper bound on right step length, in meters
if ~isfield(S.bounds,'SLR')
    S.bounds.SLR        = [];
end
if ~isfield(S.bounds.SLR,'upper')
    S.bounds.SLR.upper = [];
end

% lower bound on distance travelled, in meters
if ~isfield(S.bounds.dist_trav,'lower')
    S.bounds.dist_trav.lower = [];
end

% bound on final time, in seconds
if ~isfield(S.bounds,'t_final')
    S.bounds.t_final    = [];
end
if ~isfield(S.bounds.t_final,'upper')
    S.bounds.t_final.upper = 2;
end
if ~isfield(S.bounds.t_final,'lower')
    S.bounds.t_final.lower = 0.1;
end

if ~isfield(S.bounds,'activation_all_muscles')
    S.bounds.activation_all_muscles = [];
end

% default coordinate bounds
if ~isfield(S.bounds,'default_coordinate_bounds')
    S.bounds.default_coordinate_bounds = 'Default_Coordinate_Bounds.csv';
end

% manually overwrite coordinate bounds
if ~isfield(S.bounds,'Qs')
    S.bounds.Qs = [];
end
if ~isfield(S.bounds,'Qdots')
    S.bounds.Qdots = [];
end
if ~isfield(S.bounds,'Qdotdots')
    S.bounds.Qdotdots = [];
end

% minimal muscle activation, a number between 0 and 1
if ~isfield(S.bounds.activation_all_muscles,'lower')
    S.bounds.activation_all_muscles.lower = 0.05;
end

% maximal muscle activation, a number between 0 and 1
if ~isfield(S.bounds.activation_all_muscles,'upper')
    S.bounds.activation_all_muscles.upper = 1;
end

% activation of selected muscles
if ~isfield(S.bounds,'activation_selected_muscles')
    S.bounds.activation_selected_muscles = [];
end

% when determining bounds based on range set in opensim model, qdot and
% qdotdot bounds are proportional to range of motion
if ~isfield(S.bounds,'Qdots_factor_RoM')
    S.bounds.Qdots_factor_RoM = 10;
end
if ~isfield(S.bounds,'Qdotdots_factor_RoM')
    S.bounds.Qdotdots_factor_RoM = 155;
end

% set pelvis ty bounds based on IG_pelvis_ty
if ~isfield(S.bounds,'factor_IG_pelvis_ty')
    S.bounds.factor_IG_pelvis_ty = [];
end
if ~isfield(S.bounds.factor_IG_pelvis_ty,'lower')
    S.bounds.factor_IG_pelvis_ty.lower = 0.5;
end
if ~isfield(S.bounds.factor_IG_pelvis_ty,'upper')
    S.bounds.factor_IG_pelvis_ty.upper = 1.2;
end

% constraints on distance between points
if ~isfield(S.bounds,'points')
    S.bounds.points = [];
end
if ~isfield(S.bounds,'distanceConstraints')
    S.bounds.distanceConstraints = [];
end

%% metabolicE
if ~isfield(S,'metabolicE')
    S.metabolicE = [];
end

% hyperbolic tangent smoothing factor (used in metabolic cost)
if ~isfield(S.metabolicE,'tanh_b')
    S.metabolicE.tanh_b = 10;
end

% name of the metabolic energy model used
if ~isfield(S.metabolicE,'model')
    S.metabolicE.model = 'Bhargava2004';
end

%% misc
if ~isfield(S,'misc')
    S.misc = [];
end

% subject folder to save intermediate data
S.misc.subject_path = fullfile(S.misc.main_path,'Subjects',S.subject.name);

% average velocity you want the model to have, in meters per second
if ~isfield(S.misc,'forward_velocity')
    S.misc.forward_velocity = 1.25;
end

% maximal contraction velocity identifier --TO CHECK--
if ~isfield(S.misc,'v_max_s')
    S.misc.v_max_s = 0;
end

% type of gait simulation
if ~isfield(S.misc,'gaitmotion_type')
    S.misc.gaitmotion_type = 'HalfGaitCycle';
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

% name to save musculoskeletal geometry CasADi function
[~, model_name, ~] = fileparts(osim_path);
model_name = char(strrep(model_name, ' ', '_'));
S.misc.msk_geom_name = [model_name '_f_lMT_vMT_dM'];
if strcmp(S.misc.msk_geom_eq,'polynomials') 
    S.misc.msk_geom_name = [S.misc.msk_geom_name '_poly_',...
        num2str(S.misc.poly_order.lower) '_' num2str(S.misc.poly_order.upper)];
end

% default coordinate bounds used to approximate musculoskeletal geometry
if ~isfield(S.misc,'default_msk_geom_bounds')
    S.misc.default_msk_geom_bounds = 'default_msk_geom_bounds.csv';
end

% manually overwrite coordinate bounds used to approximate musculoskeletal geometry
if ~isfield(S.misc,'msk_geom_bounds')
    S.misc.msk_geom_bounds = [];
end

% number of data points to sample musculoskeletal geometry
if ~isfield(S.misc,'msk_geom_n_samples')
    S.misc.msk_geom_n_samples = 5000;
end

% rmse threshold for muscle-tendon length approximation
if ~isfield(S.misc,'threshold_lMT_fit')
    S.misc.threshold_lMT_fit = 0.003;
end

% rmse threshold for muscle-tendon momentarm approximation
if ~isfield(S.misc,'threshold_dM_fit')
    S.misc.threshold_dM_fit = 0.003;
end

% visualize IG and bounds
if ~isfield(S.misc,'visualize_bounds')
    S.misc.visualize_bounds = 0;
end

% damping Coefficient in contraction dynamics
if ~isfield(S.misc,'dampingCoefficient') || isempty(S.misc.dampingCoefficient)
    S.misc.dampingCoefficient = 0.01;
end

% assume constant pennation angle instead of constant width
if ~isfield(S.misc,'constant_pennation_angle')
    S.misc.constant_pennation_angle = 0;
end

% default scale factors for variables and constraints
if ~isfield(S.misc,'default_scaling_NLP')
    S.misc.default_scaling_NLP = [];
end

% manually set scale factors
if ~isfield(S.misc,'scaling_Qs')
    S.misc.scaling_Qs = [];
end
if ~isfield(S.misc,'scaling_Qdots')
    S.misc.scaling_Qdots = [];
end
if ~isfield(S.misc,'scaling_Qdotdots')
    S.misc.scaling_Qdotdots = [];
end
if ~isfield(S.misc,'scaling_Moments')
    S.misc.scaling_Moments = [];
end

% folder path to store the subject specific results
if ~isfield(S.misc,'save_folder')
   error("Please provide a folder to store the results in. " + ...
       "Specify the folder path in S.misc.save_folder."); 
elseif ~isfolder(S.misc.save_folder)
    mkdir(S.misc.save_folder);
end

% name used for saving the resultfiles (choose custom or structurized savename)
if ~isfield(S.misc,'savename')
    S.misc.savename = 'structured';
end

% filename of the result
if ~isfield(S.misc,'result_filename')
    S.misc.result_filename = [];
end

%% post_process
if ~isfield(S,'post_process')
    S.post_process = [];
end

% rerun post-processing without solving OCP
if ~isfield(S.post_process,'rerun')
    S.post_process.rerun = 0;
end
if S.post_process.rerun && isempty(S.misc.result_filename)
    error(['Please provide the name of the result to post-process. ' ...
        '(S.post_process.result_filename)'])
end

% load w_opt and reconstruct R before rerunning the post-processing
% Advanced feature, for debugging only, you should not need this.
if ~isfield(S.post_process,'load_prev_opti_vars')
    S.post_process.load_prev_opti_vars = 0;
end
if S.post_process.load_prev_opti_vars && isempty(S.misc.result_filename)
    error(['Please provide the name of the result from which to load the ' ...
        'optimization variables. (S.post_process.result_filename)'])
end

%% solver
if ~isfield(S,'solver')
    S.solver = [];
end

% options for nlpsol
if ~isfield(S.solver,'nlpsol_options')
    S.solver.nlpsol_options = [];
end

% options for ipopt
if ~isfield(S.solver,'ipopt_options')
    S.solver.ipopt_options = [];
end

% solver algorithm used in the OCP
if ~isfield(S.solver,'linear_solver')
    S.solver.linear_solver = 'mumps';
end

% the power (10^-x) the error has to reach before the OCP can 
% be regarded as solved; a higher number gives a more precise answer, but
% requires more time
if ~isfield(S.solver,'tol_ipopt')
    S.solver.tol_ipopt = 4;
end

% constraint violation has to be below 10^(-x) at the solution
if ~isfield(S.solver,'constr_viol_tol_ipopt')
    S.solver.constr_viol_tol_ipopt = 6;
end

% maximal amount of itereations after wich the solver will stop
if ~isfield(S.solver,'max_iter')
    S.solver.max_iter = 10000;
end

% type of parallel computing
if ~isfield(S.solver,'parallel_mode')
    S.solver.parallel_mode = 'thread';
end

% ADD CHECK ST THIS IS ONLY  USED WHEN USING THREAD PARALLEL MODE?
% number of threads in parallel mode
if ~isfield(S.solver,'N_threads')
    S.solver.N_threads = 4;
end

% number of mesh intervals
if ~isfield(S.solver,'N_meshes')
    S.solver.N_meshes = 50;
    % by default, use double for full gait cycle
    if strcmp(S.misc.gaitmotion_type,'FullGaitCycle')
        S.solver.N_meshes = S.solver.N_meshes*2;
    end
end



% initial guess inputs
% input should be a string: "quasi-random" or the path to a .mot file
if ~isfield(S.solver,'IG_selection')
    error(['Please specify what you want to use as an initial guess. Either ' ...
        'choose "quasi-random" or specify the path of a .mot file in ' ...
        'S.solver.IG_selection.'])
else
    [~,NAME,EXT] = fileparts(S.solver.IG_selection);
    if EXT == ".mot" && isfile(S.solver.IG_selection)
        disp(['Using ',char(S.solver.IG_selection), ' as initial guess.'])
        if ~isfield(S.solver,'IG_selection_gaitCyclePercent')
            error(['Please specify what percent of gait cycle data is present ' ...
                'in the initial guess file in S.solver.IG_selection_gaitCyclePercent. ' ...
                'For example, use 50, 100, and 200 if the provided intial guess ' ...
                'file has half a gait cycle, full gait cycle or two full ' ...
                'gait cycles, respectively.'])
        end
        if ((strcmp(S.misc.gaitmotion_type,'FullGaitCycle')) ...
                && (S.solver.IG_selection_gaitCyclePercent < 100))
            error(['Cannot use an initial guess of an incomplete gait cycle ' ...
                'for predictive simulation of a full gait cycle. Please ' ...
                'adjust S.misc.gaitmotion_type or initial guess file.'])
        elseif ((strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')) ...
                && (S.solver.IG_selection_gaitCyclePercent < 50))
            error(['Cannot use an initial guess of less than a half gait ' ...
                'cycle for predictive simulation of a half gait cycle. Please ' ...
                'adjust S.misc.gaitmotion_type or initial guess file.'])
        end
        
    elseif EXT == ".mot" && ~isfile(S.solver.IG_selection)
        error(['The motion file path you specified does not exist. Check if ' ...
            'the path exists or if you made a typo.'])
        
    elseif EXT == "" && NAME == "quasi-random"
         disp('Using a quasi-random initial guess.')
         
    elseif EXT == "" && NAME ~= "quasi-random"
        error(['Please specify what you want to use as an initial guess. ' ...
            'Either choose "quasi-random" or specify the path of a .mot file.'])
    end
end

%% subject
if ~isfield(S,'subject')
    S.subject = [];
end

% name of the subject
if ~isfield(S.subject,'name')
    error(['Please provide a name for this subject. This name will be used ' ...
        'to store the results. Specify the name in S.subject.name.']);
end

% mass of the subject, in kilograms
if ~isfield(S.subject,'mass')
   S.subject.mass = [];
end

% height of the pelvis for the initial guess, in meters
if ~isfield(S.subject,'IG_pelvis_y')
   S.subject.IG_pelvis_y = [];
end

% muscle strength
if ~isfield(S.subject,'muscle_strength')
    S.subject.muscle_strength = [];
end

% muscle stiffness
if ~isfield(S.subject,'muscle_pass_stiff_shift')
    S.subject.muscle_pass_stiff_shift = [];
end
if ~isfield(S.subject,'muscle_pass_stiff_scale')
    S.subject.muscle_pass_stiff_scale = [];
end

% tendon stiffness
if ~isfield(S.subject,'tendon_stiff_scale')
    S.subject.tendon_stiff_scale = [];
end

% type of mtp joint used in the model
if ~isfield(S.subject,'mtp_type')
    S.subject.mtp_type = ''; 
end

% muscle tendon properties
if ~isfield(S.subject,'scale_MT_params')
    S.subject.scale_MT_params = []; 
end

% damping coefficient for all degrees of freedon
if ~isfield(S.subject,'damping_coefficient_all_dofs')
    S.subject.damping_coefficient_all_dofs = 0.1; 
end

% different damping coefficient for specific degrees of freedon
if ~isfield(S.subject,'set_damping_coefficient_selected_dofs')
    S.subject.set_damping_coefficient_selected_dofs = []; 
end

% stiffness coefficient for all degrees of freedon
if ~isfield(S.subject,'stiffness_coefficient_all_dofs')
    S.subject.stiffness_coefficient_all_dofs = 0; 
end

% different stiffness coefficient for specific degrees of freedon
if ~isfield(S.subject,'set_stiffness_coefficient_selected_dofs')
    S.subject.set_stiffness_coefficient_selected_dofs = []; 
end

% neutral position for stiffness coefficient for specific degrees of freedon
if ~isfield(S.subject,'set_stiffness_offset_selected_dofs')
    S.subject.set_stiffness_offset_selected_dofs = []; 
end

% default limit torque coefficient
if ~isfield(S.subject,'default_coord_lim_torq_coeff')
    S.subject.default_coord_lim_torq_coeff = 'default_coord_lim_torq_coeff.csv';
end

% scale factor for limit torque amplitude
if ~isfield(S.subject,'scale_default_coord_lim_torq')
    S.subject.scale_default_coord_lim_torq = [];
end

% limit torque coefficient for specific degrees of freedon
if ~isfield(S.subject,'set_limit_torque_coefficients_selected_dofs')
    S.subject.set_limit_torque_coefficients_selected_dofs = []; 
end

% ligament stiffness for all ligaments
if ~isfield(S.subject,'stiffness_all_ligaments')
    S.subject.stiffness_all_ligaments = 'ligamentGefen2002'; 
end

% ligament stiffness for selected ligaments
if ~isfield(S.subject,'set_stiffness_selected_ligaments')
    S.subject.set_stiffness_selected_ligaments =...
        {'PlantarFascia','plantarFasciaNatali2010'};
end

% joints that are considered base of a leg
if ~isfield(S.subject,'base_joints_legs')
    S.subject.base_joints_legs = 'hip';
end

% joints that are considered base of an arm
if ~isfield(S.subject,'base_joints_arms')
    S.subject.base_joints_arms = 'acromial';
end

%% weights
if ~isfield(S,'weights')
    S.weights = [];
end

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

% weight on actuator excitations
if ~isfield(S.weights,'e_torqAct')
    S.weights.e_torqAct = 10^6; 
end

% weight on passive torques
if ~isfield(S.weights,'pass_torq')
    S.weights.pass_torq = 1000; 
end

% damping can be included in the passive torques that go in the cost
% function
if ~isfield(S.weights,'pass_torq_includes_damping')
    S.weights.pass_torq_includes_damping = 0; 
end

% weight on muscle activations
if ~isfield(S.weights,'a')
    S.weights.a = 2000; 
end

% exponent for the muscle activations
if ~isfield(S.weights,'a_exp')
    S.weights.a_exp = 2; 
end

% weight on slack controls
if ~isfield(S.weights,'slack_ctrl')
    S.weights.slack_ctrl = 0.001; 
end


%% OpenSimADOptions

% settings for functions to convert .osim model to expression graph 
% file (.dll) to solve inverse dynamics
if ~isfield(S,'OpenSimADOptions')
    S.OpenSimADOptions = [];
end

% select compiler for cpp projects 
%   Visual studio 2015: 'Visual Studio 14 2015 Win64'
%   Visual studio 2017: 'Visual Studio 15 2017 Win64'
%   Visual studio 2017: 'Visual Studio 16 2019'
%   Visual studio 2017: 'Visual Studio 17 2022'
if ~isfield(S.OpenSimADOptions,'compiler')
    S.OpenSimADOptions.compiler = findVisualStudioInstallation;
end

% Input forces acting on bodies
if ~isfield(S.OpenSimADOptions,'input3DBodyForces')
    S.OpenSimADOptions.input3DBodyForces = [];
end

% Input forces acting on bodies
if ~isfield(S.OpenSimADOptions,'input3DBodyMoments')
    S.OpenSimADOptions.input3DBodyMoments = [];
end

% Export positions of points w.r.t. ground frame
if ~isfield(S.OpenSimADOptions,'export3DPositions')
    S.OpenSimADOptions.export3DPositions = [];
end

% Export velocities of points w.r.t. ground frame
if ~isfield(S.OpenSimADOptions,'export3DVelocities')
    S.OpenSimADOptions.export3DVelocities = [];
end

% If you want to choose the order of the joints and coordinate outputs
if ~isfield(S.OpenSimADOptions,'jointsOrder')
    S.OpenSimADOptions.jointsOrder = [];
end
if ~isfield(S.OpenSimADOptions,'coordinatesOrder')
    S.OpenSimADOptions.coordinatesOrder = [];
end

% Export total GRFs
if ~isfield(S.OpenSimADOptions,'exportGRFs')
    S.OpenSimADOptions.exportGRFs = true;
end

% Export separate GRFs of each contact element
if ~isfield(S.OpenSimADOptions,'exportSeparateGRFs')
    S.OpenSimADOptions.exportSeparateGRFs = true;
end

% Export GRMs.
if ~isfield(S.OpenSimADOptions,'exportGRMs')
    S.OpenSimADOptions.exportGRMs = true;
end

% Export contact sphere vertical deformation power.
if ~isfield(S.OpenSimADOptions,'exportContactPowers')
    S.OpenSimADOptions.exportContactPowers = true;
end

% Print outputs from compiler
if ~isfield(S.OpenSimADOptions,'verbose_mode')
    S.OpenSimADOptions.verbose_mode = false;
end 

% Test external function versus ID Tool
if ~isfield(S.OpenSimADOptions,'verify_ID')
    S.OpenSimADOptions.verify_ID = false;
end 


end