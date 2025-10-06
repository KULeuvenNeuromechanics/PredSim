function slow_twitch_fiber_ratio = getSlowTwitchRatios(muscleNames)
% --------------------------------------------------------------------------
% getSlowTwitchRatios
%   This function returns the ratio of slow twitch fibers in the muscles
%   The data come from Uchida et al. (2016).
%   We use 0.5 when no data were available.
%   
% INPUT:
%   - muscleNames -
%   * Cell array of muscle names
%
% OUTPUT:
%   - slow_twitch_fiber_ratio -
%   * Array with ratio of slow twitch fibers in each muscle. Default is 0.50
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

%% gait2392 muscles
% right muscles
pctst_data.glut_med1_r = 0.55;
pctst_data.glut_med2_r = 0.55;
pctst_data.glut_med3_r = 0.55;
pctst_data.glut_min1_r = 0.55;
pctst_data.glut_min2_r = 0.55;
pctst_data.glut_min3_r = 0.55;
pctst_data.semimem_r = 0.4925;
pctst_data.semiten_r = 0.425;    
pctst_data.bifemlh_r = 0.5425;
pctst_data.bifemsh_r = 0.529;
pctst_data.sar_r = 0.50;
pctst_data.add_mag1_r = 0.552;
pctst_data.add_mag2_r = 0.552;
pctst_data.add_mag3_r = 0.552;
pctst_data.tfl_r = 0.50;
pctst_data.pect_r = 0.50;
pctst_data.grac_r = 0.50;
pctst_data.glut_max1_r = 0.55;
pctst_data.glut_max2_r = 0.55;
pctst_data.glut_max3_r = 0.55;
pctst_data.iliacus_r = 0.50;
pctst_data.psoas_r = 0.50;
pctst_data.quad_fem_r = 0.50;
pctst_data.gem_r = 0.50;
pctst_data.peri_r = 0.50;
pctst_data.rect_fem_r = 0.3865;
pctst_data.vas_med_r = 0.503;
pctst_data.vas_int_r = 0.543;
pctst_data.vas_lat_r = 0.455;
pctst_data.med_gas_r = 0.566;
pctst_data.lat_gas_r = 0.507;
pctst_data.soleus_r = 0.803;
pctst_data.tib_post_r = 0.60;
pctst_data.flex_dig_r = 0.60;
pctst_data.flex_hal_r = 0.60;
pctst_data.tib_ant_r = 0.70;    
pctst_data.per_brev_r = 0.60;
pctst_data.per_long_r = 0.60;
pctst_data.per_tert_r = 0.75;
pctst_data.ext_dig_r = 0.75;
pctst_data.ext_hal_r = 0.75;    
pctst_data.ercspn_r = 0.60;
pctst_data.intobl_r = 0.56;
pctst_data.extobl_r = 0.58;    
pctst_data.add_long_r = 0.50;
pctst_data.add_brev_r = 0.50;

% left muscles
pctst_data.glut_med1_l = 0.55;
pctst_data.glut_med2_l = 0.55;
pctst_data.glut_med3_l = 0.55;
pctst_data.glut_min1_l = 0.55;
pctst_data.glut_min2_l = 0.55;
pctst_data.glut_min3_l = 0.55;
pctst_data.semimem_l = 0.4925;
pctst_data.semiten_l = 0.425;    
pctst_data.bifemlh_l = 0.5425;
pctst_data.bifemsh_l = 0.529;
pctst_data.sar_l = 0.50;
pctst_data.add_mag1_l = 0.552;
pctst_data.add_mag2_l = 0.552;
pctst_data.add_mag3_l = 0.552;
pctst_data.tfl_l = 0.50;
pctst_data.pect_l = 0.50;
pctst_data.grac_l = 0.50;
pctst_data.glut_max1_l = 0.55;
pctst_data.glut_max2_l = 0.55;
pctst_data.glut_max3_l = 0.55;
pctst_data.iliacus_l = 0.50;
pctst_data.psoas_l = 0.50;
pctst_data.quad_fem_l = 0.50;
pctst_data.gem_l = 0.50;
pctst_data.peri_l = 0.50;
pctst_data.rect_fem_l = 0.3865;
pctst_data.vas_med_l = 0.503;
pctst_data.vas_int_l = 0.543;
pctst_data.vas_lat_l = 0.455;
pctst_data.med_gas_l = 0.566;
pctst_data.lat_gas_l = 0.507;
pctst_data.soleus_l = 0.803;
pctst_data.tib_post_l = 0.60;
pctst_data.flex_dig_l = 0.60;
pctst_data.flex_hal_l = 0.60;
pctst_data.tib_ant_l = 0.70;    
pctst_data.per_brev_l = 0.60;
pctst_data.per_long_l = 0.60;
pctst_data.per_tert_l = 0.75;
pctst_data.ext_dig_l = 0.75;
pctst_data.ext_hal_l = 0.75;    
pctst_data.ercspn_l = 0.60;
pctst_data.intobl_l = 0.56;
pctst_data.extobl_l = 0.58;    
pctst_data.add_long_l = 0.50;
pctst_data.add_brev_l = 0.50;

%% gait1018 muscle groups
% Right muscle groups
pctst_data.hamstrings_r = 0.5425;
pctst_data.glut_max_r = 0.55;
pctst_data.iliopsoas_r = 0.50;
pctst_data.vasti_r = 0.543;
pctst_data.gastroc_r = 0.566;

%  reft muscle groups
pctst_data.hamstrings_l = 0.5425;
pctst_data.glut_max_l = 0.55;
pctst_data.iliopsoas_l = 0.50;
pctst_data.vasti_l = 0.543;
pctst_data.gastroc_l = 0.566;

%% shoulder muscles
% data from
% 1: 10.1002/ca.20349
% 2: https://www.sciencedirect.com/science/article/pii/0022510X73900233

% right side
pctst_data.BRA_r = 0.35; % cfr BIC and TRI
pctst_data.TrapeziusScapula_M_r = 0.537; %2
pctst_data.TrapeziusScapula_S_r = 0.537; %2
pctst_data.TrapeziusScapula_I_r = 0.537; %2
pctst_data.TrapeziusClavicle_S_r = 0.537; %2
% pctst_data.SerratusAnterior_I_r = % no data?
% pctst_data.SerratusAnterior_M_r = % no data?
% pctst_data.SerratusAnterior_S_r = % no data?
pctst_data.Rhomboideus_S_r = 0.446; %2
pctst_data.Rhomboideus_I_r = 0.446; %2
% pctst_data.LevatorScapulae_r = % no data?
pctst_data.Coracobrachialis_r = 0.29; %1
pctst_data.DeltoideusClavicle_A_r = 0.55; %1
pctst_data.DeltoideusScapula_P_r = 0.57; %1
pctst_data.DeltoideusScapula_M_r = 0.46; %1
pctst_data.LatissimusDorsi_S_r = 0.47; %1
pctst_data.LatissimusDorsi_M_r = 0.47; %1
pctst_data.LatissimusDorsi_I_r = 0.47; %1
pctst_data.PectoralisMajorClavicle_S_r = 0.36; %1
pctst_data.PectoralisMajorThorax_I_r = 0.36; %1
pctst_data.PectoralisMajorThorax_M_r = 0.36; %1
pctst_data.TeresMajor_r = 0.48; %1
pctst_data.Infraspinatus_I_r = 0.54; %1
pctst_data.Infraspinatus_S_r = 0.54; %1
pctst_data.PectoralisMinor_r = 0.36; %1
pctst_data.TeresMinor_r = 0.53; %1
pctst_data.Subscapularis_S_r = 0.42; %1
pctst_data.Subscapularis_M_r = 0.42; %1
pctst_data.Subscapularis_I_r = 0.42; %1
pctst_data.Supraspinatus_P_r = 0.51; %1
pctst_data.Supraspinatus_A_r = 0.51; %1
pctst_data.TRI_long_r = 0.43; %1
pctst_data.TRI_med_r = 0.33; %2
pctst_data.TRI_lat_r = 0.33; %2
pctst_data.BIC_long_r = 0.32; %1
pctst_data.BIC_brevis_r = 0.38; %1

% left side
pctst_data.BRA_l = 0.35; % cfr BIC and TRI
pctst_data.TrapeziusScapula_M_l = 0.537; %2
pctst_data.TrapeziusScapula_S_l = 0.537; %2
pctst_data.TrapeziusScapula_I_l = 0.537; %2
pctst_data.TrapeziusClavicle_S_l = 0.537; %2
% pctst_data.SerratusAnterior_I_l = % no data?
% pctst_data.SerratusAnterior_M_l = % no data?
% pctst_data.SerratusAnterior_S_l = % no data?
pctst_data.Rhomboideus_S_l = 0.446; %2
pctst_data.Rhomboideus_I_l = 0.446; %2
% pctst_data.LevatorScapulae_l = % no data?
pctst_data.Coracobrachialis_l = 0.29; %1
pctst_data.DeltoideusClavicle_A_l = 0.55; %1
pctst_data.DeltoideusScapula_P_l = 0.57; %1
pctst_data.DeltoideusScapula_M_l = 0.46; %1
pctst_data.LatissimusDorsi_S_l = 0.47; %1
pctst_data.LatissimusDorsi_M_l = 0.47; %1
pctst_data.LatissimusDorsi_I_l = 0.47; %1
pctst_data.PectoralisMajorClavicle_S_l = 0.36; %1
pctst_data.PectoralisMajorThorax_I_l = 0.36; %1
pctst_data.PectoralisMajorThorax_M_l = 0.36; %1
pctst_data.TeresMajor_l = 0.48; %1
pctst_data.Infraspinatus_I_l = 0.54; %1
pctst_data.Infraspinatus_S_l = 0.54; %1
pctst_data.PectoralisMinor_l = 0.36; %1
pctst_data.TeresMinor_l = 0.53; %1
pctst_data.Subscapularis_S_l = 0.42; %1
pctst_data.Subscapularis_M_l = 0.42; %1
pctst_data.Subscapularis_I_l = 0.42; %1
pctst_data.Supraspinatus_P_l = 0.51; %1
pctst_data.Supraspinatus_A_l = 0.51; %1
pctst_data.TRI_long_l = 0.43; %1
pctst_data.TRI_med_l = 0.33; %2
pctst_data.TRI_lat_l = 0.33; %2
pctst_data.BIC_long_l = 0.32; %1
pctst_data.BIC_brevis_l = 0.38; %1

%%
slow_twitch_fiber_ratio = zeros(length(muscleNames),1);
for i = 1:length(muscleNames)
    if isfield(pctst_data,muscleNames{i})
        slow_twitch_fiber_ratio(i,1) = pctst_data.(muscleNames{i});
    else
%         disp(['Default ratio of Slow twitch fibers for muscle ' muscleNames{i} ' selected (0.50)']);
        slow_twitch_fiber_ratio(i,1) = 0.5;
    end
end

end   
