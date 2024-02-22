% --------------------------------------------------------------------------
% Settings for Falisse_et_al_2022 that deviate from the PredSim defaults
%
%   A. Falisse, M. Afschrift, and F. D. Groote, “Modeling toes contributes 
%   to realistic stance knee mechanics in three-dimensional predictive 
%   simulations of walking,” PLOS ONE, vol. 17, no. 1, p. e0256311, 
%   Jan. 2022, doi: 10.1371/journal.pone.0256311.
%
% Original author: Lars D'Hondt
% Original date: 11/December/2023
% --------------------------------------------------------------------------

S.subject.name = 'Falisse_et_al_2022';

S.subject.v_pelvis_x_trgt   = 1.33;

S.subject.mtp_type = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',25};
S.subject.set_damping_coefficient_selected_dofs = {'mtp_angle',2};

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

