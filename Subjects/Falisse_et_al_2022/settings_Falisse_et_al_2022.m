function S = settings_Falisse_et_al_2022(S)
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
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 01/February/2024
% --------------------------------------------------------------------------

S.subject.name = 'Falisse_et_al_2022';

S.subject.v_pelvis_x_trgt   = 1.33;

S.subject.mtp_type = '2022paper';

switch S.personalization.jointparams
    case 'j1'
        S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',25};
        S.subject.set_damping_coefficient_selected_dofs = {'mtp_angle',2};
end

switch S.personalization.strength
    case 's1'
        S.subject.muscle_strength   = {...
            {'lat_gas_r' 'med_gas_r' 'soleus_r' },0.7};
    case 's2'
        S.subject.muscle_strength   = {...
            {'lat_gas_r' 'med_gas_r' 'soleus_r' },0.5};
end

switch S.personalization.MTparams
    case 'mt1'
        S.subject.scale_MT_params = {...
            {'bifemsh_r','bifemlh_r','grac_r','semimem_r','semiten_r','sar_r'},'lMo',0.815}; 
    case 'mt2'
        S.subject.scale_MT_params= {...
            {'soleus_l','soleus_r','med_gas_l','med_gas_r','lat_gas_l','lat_gas_r'},'FMo',1.2,...
            {'soleus_l','soleus_r','med_gas_l','med_gas_r','lat_gas_l','lat_gas_r'},'lMo',0.8,...
            {'soleus_l','soleus_r','med_gas_l','med_gas_r','lat_gas_l','lat_gas_r'},'lTs',1.05};
end

end