% This function returns the percentage of slow twitch fibers in the muscles
% The data come from Uchida et al. (2016).
% We used 0.5 when no data were available.
%
% Author: Antoine Falisse
% Date: 12/19/2018
%     
function pctst = getSlowTwitchRatios_lr(muscleNames)

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
    
    pctst = zeros(length(muscleNames),1);
    for i = 1:length(muscleNames)
        pctst(i,1) = pctst_data.(muscleNames{i});
    end
    
end   
