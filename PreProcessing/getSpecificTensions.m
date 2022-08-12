function specific_tension = getSpecificTensions(muscleNames)
% --------------------------------------------------------------------------
% getSpecificTensions
%    This function returns the specific tensions in the muscles.
%    The data comes from Uchida et al. (2016).
%   
% INPUT:
%   - muscleNames -
%   * Cell array of muscle names
%
% OUTPUT:
%   - specific_tension -
%   * Array with specific tension of each muscle. Default is 0.70
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
%   Edit: add values for muscle groups
% Last edit by: Lars D'Hondt
% Last edit date: 12/Aug/2022
% --------------------------------------------------------------------------

% Right muscles
sigma_data.glut_med1_r = 0.74455;
sigma_data.glut_med2_r = 0.75395;
sigma_data.glut_med3_r = 0.75057;
sigma_data.glut_min1_r = 0.75;
sigma_data.glut_min2_r = 0.75;
sigma_data.glut_min3_r = 0.75116;
sigma_data.semimem_r = 0.62524;
sigma_data.semiten_r = 0.62121;    
sigma_data.bifemlh_r = 0.62222;
sigma_data.bifemsh_r = 1.00500;
sigma_data.sar_r = 0.74286;
sigma_data.add_mag1_r = 0.55217;
sigma_data.add_mag2_r = 0.55323;
sigma_data.add_mag3_r = 0.54831;
sigma_data.tfl_r = 0.75161;
sigma_data.pect_r = 0.76000;
sigma_data.grac_r = 0.73636;
sigma_data.glut_max1_r = 0.75395;
sigma_data.glut_max2_r = 0.74455;
sigma_data.glut_max3_r = 0.74595;
sigma_data.iliacus_r = 1.2477;
sigma_data.psoas_r = 1.5041;
sigma_data.quad_fem_r = 0.74706;
sigma_data.gem_r = 0.74545;
sigma_data.peri_r = 0.75254;
sigma_data.rect_fem_r = 0.74936;
sigma_data.vas_med_r = 0.49961;
sigma_data.vas_int_r = 0.55263;
sigma_data.vas_lat_r = 0.50027;
sigma_data.med_gas_r = 0.69865;
sigma_data.lat_gas_r = 0.69694;
sigma_data.soleus_r = 0.62703;
sigma_data.tib_post_r = 0.62520;
sigma_data.flex_dig_r = 0.5;
sigma_data.flex_hal_r = 0.50313;
sigma_data.tib_ant_r = 0.75417;    
sigma_data.per_brev_r = 0.62143;
sigma_data.per_long_r = 0.62450;
sigma_data.per_tert_r = 1;
sigma_data.ext_dig_r = 0.75294;
sigma_data.ext_hal_r = 0.73636;    
sigma_data.ercspn_r = 0.25;
sigma_data.intobl_r = 0.25;
sigma_data.extobl_r = 0.25;    
sigma_data.add_long_r = 0.74643;
sigma_data.add_brev_r = 0.75263;

% Left muscles
sigma_data.glut_med1_l = 0.74455;
sigma_data.glut_med2_l = 0.75395;
sigma_data.glut_med3_l = 0.75057;
sigma_data.glut_min1_l = 0.75;
sigma_data.glut_min2_l = 0.75;
sigma_data.glut_min3_l = 0.75116;
sigma_data.semimem_l = 0.62524;
sigma_data.semiten_l = 0.62121;    
sigma_data.bifemlh_l = 0.62222;
sigma_data.bifemsh_l = 1.00500;
sigma_data.sar_l = 0.74286;
sigma_data.add_mag1_l = 0.55217;
sigma_data.add_mag2_l = 0.55323;
sigma_data.add_mag3_l = 0.54831;
sigma_data.tfl_l = 0.75161;
sigma_data.pect_l = 0.76000;
sigma_data.grac_l = 0.73636;
sigma_data.glut_max1_l = 0.75395;
sigma_data.glut_max2_l = 0.74455;
sigma_data.glut_max3_l = 0.74595;
sigma_data.iliacus_l = 1.2477;
sigma_data.psoas_l = 1.5041;
sigma_data.quad_fem_l = 0.74706;
sigma_data.gem_l = 0.74545;
sigma_data.peri_l = 0.75254;
sigma_data.rect_fem_l = 0.74936;
sigma_data.vas_med_l = 0.49961;
sigma_data.vas_int_l = 0.55263;
sigma_data.vas_lat_l = 0.50027;
sigma_data.med_gas_l = 0.69865;
sigma_data.lat_gas_l = 0.69694;
sigma_data.soleus_l = 0.62703;
sigma_data.tib_post_l = 0.62520;
sigma_data.flex_dig_l = 0.5;
sigma_data.flex_hal_l = 0.50313;
sigma_data.tib_ant_l = 0.75417;    
sigma_data.per_brev_l = 0.62143;
sigma_data.per_long_l = 0.62450;
sigma_data.per_tert_l = 1;
sigma_data.ext_dig_l = 0.75294;
sigma_data.ext_hal_l = 0.73636;    
sigma_data.ercspn_l = 0.25;
sigma_data.intobl_l = 0.25;
sigma_data.extobl_l = 0.25;    
sigma_data.add_long_l = 0.74643;
sigma_data.add_brev_l = 0.75263;

% Right muscle groups
sigma_data.hamstrings_r = 0.62222;
sigma_data.glut_max_r = 0.74455;
sigma_data.iliopsoas_r = 1.5041;
sigma_data.vasti_r = 0.55263;
sigma_data.gastroc_r = 0.69865;

%  Left muscle groups
sigma_data.hamstrings_l = 0.62222;
sigma_data.glut_max_l = 0.74455;
sigma_data.iliopsoas_l = 1.5041;
sigma_data.vasti_l = 0.55263;
sigma_data.gastroc_l = 0.69865;

% Get specific tension for selected muscles
specific_tension = zeros(length(muscleNames),1);
for i = 1:length(muscleNames)
    if isfield(sigma_data,muscleNames{i})
        specific_tension(i,1) = sigma_data.(muscleNames{i});
    else
        disp(['Default specific tension for muscle ' muscleNames{i} ' selected (0.70)']);
        specific_tension(i,1) = 0.7;
    end
end

end    
