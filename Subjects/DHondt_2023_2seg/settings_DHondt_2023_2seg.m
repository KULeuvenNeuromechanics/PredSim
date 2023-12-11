% --------------------------------------------------------------------------
% Settings for DHondt_2023_2seg that deviate from the PredSim defaults
%
%   L. D'Hondt, F. D. Groote, and M. Afschrift, "A dynamic foot model for 
%   predictive simulations of gait reveals causal relations between foot 
%   structure and whole body mechanics.” bioRxiv, p. 2023.03.22.533790, 
%   Mar. 24, 2023. doi: 10.1101/2023.03.22.533790.
%
% Original author: Lars D'Hondt
% Original date: 11/December/2023
% --------------------------------------------------------------------------

S.subject.name = 'DHondt_2023_2seg';

S.subject.v_pelvis_x_trgt   = 1.33;

S.solver.N_meshes = 100;

S.subject.IG_selection_gaitCyclePercent = 100;
S.subject.IG_selection = 'quasi-random';

S.metabolicE.tanh_b = 100;

S.bounds.coordinates = {'pelvis_ty',0.55,1.1, 'pelvis_tilt',-2.9302,nan};

S.subject.muscle_pass_stiff_shift = {{'soleus','_gas','per_','tib_','_dig_','_hal_'},0.9};
S.subject.tendon_stiff_scale = {{'soleus','_gas'},0.5};

S.subject.mtp_type = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',25};
S.subject.set_damping_coefficient_selected_dofs = {'mtp_angle',2};


