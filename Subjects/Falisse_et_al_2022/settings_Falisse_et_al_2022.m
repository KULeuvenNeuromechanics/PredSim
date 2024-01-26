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

