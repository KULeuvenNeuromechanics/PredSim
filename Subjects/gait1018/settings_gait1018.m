% --------------------------------------------------------------------------
% Settings for gait1018 (i.e. 2D model) that deviate from the PredSim defaults
%
% Original author: Lars D'Hondt
% Original date: 12/August/2024
% --------------------------------------------------------------------------

S.subject.name = 'gait1018';

% This model has no arms
S.subject.base_joints_arms = []; 

% Achilles tendon stiffness
S.subject.tendon_stiff_scale = {{'soleus','gastroc'},0.5};



