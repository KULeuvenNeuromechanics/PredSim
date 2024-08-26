% --------------------------------------------------------------------------
% Settings for DHondt_et_al_2024_3seg that deviate from the PredSim defaults
%
%   L. D’Hondt, F. D. Groote, and M. Afschrift, “A dynamic foot model for 
%   predictive simulations of human gait reveals causal relations between 
%   foot structure and whole-body mechanics,” PLOS Computational Biology, 
%   vol. 20, no. 6, p. e1012219, Jun. 2024, doi: 10.1371/journal.pcbi.1012219.
%
% Original author: Lars D'Hondt
% Original date: 16/July/2024
% --------------------------------------------------------------------------

S.subject.name = 'DHondt_et_al_2024_3seg';

S.subject.v_pelvis_x_trgt   = 1.33;

S.solver.N_meshes = 100;

S.subject.IG_selection = 'quasi-random';
S.subject.IG_pelvis_y = 0.9385;
S.subject.adapt_IG_pelvis_y = 0;

S.metabolicE.tanh_b = 100;

S.bounds.Qs = {'pelvis_ty',0.55,1.1, 'pelvis_tilt',-2.9302,nan};

S.subject.mtp_type = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',25};
S.subject.set_damping_coefficient_selected_dofs = {'mtp_angle',2};

S.subject.muscle_pass_stiff_shift =...
    {{'soleus','_gas','per_','tib_','_dig_','_hal_','FDB'},0.9};
S.subject.tendon_stiff_scale = {{'soleus','_gas'},0.5};



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






