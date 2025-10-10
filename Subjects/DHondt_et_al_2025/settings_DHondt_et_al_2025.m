% --------------------------------------------------------------------------
% Settings for DHondt_et_al_2025 that deviate from the PredSim defaults
%
%   
%
% Original author: Lars D'Hondt
% Original date: 03/October/2025
% --------------------------------------------------------------------------

S.subject.name = 'DHondt_et_al_2025';

S.misc.forward_velocity = 1.33;

S.solver.N_meshes = 50;

S.solver.IG_selection = 'quasi-random';
S.subject.IG_pelvis_y = 0.9385;
S.subject.adapt_IG_pelvis_y = 0;

S.metabolicE.tanh_b = 100;


% S.misc.poly_order.upper = 5;

S.subject.base_joints_legs = {'hip'};
S.subject.base_joints_arms = {'sternoclavicular','scapulothoracic'};


% S.bounds.Qs = {'pelvis_ty',0.55,1.1, 'pelvis_tilt',-2.9302,nan};

% S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',1};
% 
% S.subject.muscle_pass_stiff_shift =...
%     {{'soleus','_gas','per_','tib_','_dig_','_hal_','FDB'},0.9};




S.subject.mtp_type = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',25};
S.subject.set_damping_coefficient_selected_dofs = {'mtp_angle',2};

S.subject.muscle_pass_stiff_shift =...
    {{'soleus','_gas','per_','tib_','_dig_','_hal_'},0.9};

S.subject.tendon_stiff_scale = {{'soleus','_gas'},0.5};


S.subject.default_coord_lim_torq_coeff =...
    'default_coord_lim_torq_coeff_incl_shoulder.csv';

% to prevent body segments from clipping into eachother
S.bounds.distanceConstraints(1).point1 = 'calcn_r';
S.bounds.distanceConstraints(1).point2 = 'calcn_l';
S.bounds.distanceConstraints(1).direction = 'xz';
S.bounds.distanceConstraints(1).lower_bound = 0.09;
S.bounds.distanceConstraints(1).upper_bound = 2;

S.bounds.distanceConstraints(2).point1 = 'hand_r';
S.bounds.distanceConstraints(2).point2 = 'femur_r';
S.bounds.distanceConstraints(2).direction = 'xz';
S.bounds.distanceConstraints(2).lower_bound = 0.18;
S.bounds.distanceConstraints(2).upper_bound = 2;

S.bounds.distanceConstraints(3).point1 = 'hand_l';
S.bounds.distanceConstraints(3).point2 = 'femur_l';
S.bounds.distanceConstraints(3).direction = 'xz';
S.bounds.distanceConstraints(3).lower_bound = 0.18;
S.bounds.distanceConstraints(3).upper_bound = 2;

S.bounds.distanceConstraints(4).point1 = 'tibia_r';
S.bounds.distanceConstraints(4).point2 = 'tibia_l';
S.bounds.distanceConstraints(4).direction = 'xz';
S.bounds.distanceConstraints(4).lower_bound = 0.11;
S.bounds.distanceConstraints(4).upper_bound = 2;

S.bounds.distanceConstraints(5).point1 = 'toes_r';
S.bounds.distanceConstraints(5).point2 = 'toes_l';
S.bounds.distanceConstraints(5).direction = 'xz';
S.bounds.distanceConstraints(5).lower_bound = 0.1;
S.bounds.distanceConstraints(5).upper_bound = 2;

% S.bounds.distanceConstraints(6).point1 = 'radius_r';
% S.bounds.distanceConstraints(6).point2 = 'torso';
% S.bounds.distanceConstraints(6).direction = 'xz';
% S.bounds.distanceConstraints(6).lower_bound = 0.2;
% S.bounds.distanceConstraints(6).upper_bound = 2;
% 
% S.bounds.distanceConstraints(7).point1 = 'radius_l';
% S.bounds.distanceConstraints(7).point2 = 'torso';
% S.bounds.distanceConstraints(7).direction = 'xz';
% S.bounds.distanceConstraints(7).lower_bound = 0.2;
% S.bounds.distanceConstraints(7).upper_bound = 2;
