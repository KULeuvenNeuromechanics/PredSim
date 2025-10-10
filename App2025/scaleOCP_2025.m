[S] = scaleOCP_2025(S,U,mass_genmodel,length_genmodel);
% Ellis Van Can October 2025
% % coordinate limit torques
%              S.subject.set_limit_torque_coefficients_selected_dofs = {{'lumbar_extension'},[-0.7644*S.subject.St,11.2154,1.2788*S.subject.St,-7.2704], [-0.3716,0.1068],...
%              {'lumbar_bending'}, [-1.7301*S.subject.St,8.4705,1.6990*S.subject.St,-8.5249],[-0.2257,0.2245], ...
%              {'lumbar_rotation'}, [-1.7301*S.subject.St,8.4705,1.6990*S.subject.St,-8.5249],[-0.2257,0.2245],...
%              {'hip_flexion_r','hip_flexion_l'},[-2.44*S.subject.St,5.05,1.51*S.subject.St,-21.88],[-0.6981,1.81],...
%              {'hip_adduction_r','hip_adduction_l'},[-0.03*S.subject.St,14.94,0.03*S.subject.St,-14.94],[-0.5,0.5],...
%              {'hip_rotation_r','hip_rotation_l'},[-0.03*S.subject.St,14.94,0.03*S.subject.St,-14.94],[-0.92,0.92],...
%              {'knee_angle_r','knee_angle_l'},[-6.09*S.subject.St,33.94,11.03*S.subject.St,-11.33],[-2.4,0.13],... 
%              {'ankle_angle_r','ankle_angle_l'},[-2.03*S.subject.St,38.11,0.18*S.subject.St,-12.12],[-0.74,0.52],...
%              {'subtalar_angle_r','subtalar_angle_l'},[-60.21*S.subject.St,16.32,60.21*S.subject.St,-16.32],[-0.65,0.65],...
%              {'mtp_angle_r','mtp_angle_l'},[-0.9*S.subject.St,14.87,0.18*S.subject.St,-70.08],[0, 65/180*pi]};
% 
% % cost function
% if  isfield(S.settings,'scaling')
%     switch S.settings.scaling
%         case 'sc1' % based on Dhondt20236_2seg model and length of scaled model and weight based on CE CP_AFO_3_T0
%             S.subject.St = 0.70349;
%             S.subject.Swt =2.0206194;
%             S.subject.Swe1 = 0.686435;
%             S.subject.Swe2 = 1.456801;
%     end
%     S.weights.E = 500*S.subject.Swe2;
%     S.weights.pass_torq = 1000*S.subject.Swt;
% end
% 
% % joint params
% if  isfield(S.settings,'jointparams')
%     switch S.settings.jointparams
%         case 'j1'
%             S.subject.damping_coefficient_all_dofs = 0.1*S.subject.St;
%             S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25*S.subject.St};
%             S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2*S.subject.St,...
%             {'arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_r','arm_flex_l','arm_add_l','arm_rot_l','elbow_flex_l'},0.5*S.subject.St};
%     end
% end
% 


% sf_legLength = (sf.low_leg + sf.upp_leg)/2;
% % sf_mass = U.Mass/75;
% sf_force = sf_mass^(2/3);

% 
% actuator torques (~ size^3)
S.subject.scale_actuator_torque = {
    'lumbar_extension',sf_force*(sqrt(sf.torso*sf.shoulder)*sf.torso*sf.upp_leg),... % ~ torso and pelvis
    {'arm_flex_r','arm_flex_l'},sf_force*(sqrt(sf.torso*sf.shoulder)*sf.shoulder^2) % ~ torso and shoulder
    };
% 
% S.subject.set_damping_coefficient_selected_dofs = {
% %     'lumbar_extension',2*sf_force*(sqrt(sf.torso*sf.shoulder)*sf.torso*sf.upp_leg),...
%     {'arm_flex_r','arm_flex_l'},0.2*sf_force*(sqrt(sf.torso*sf.shoulder)*sf.shoulder^2),...
%     };
% 
% pass_torq_arm = 10*sf_force*(sqrt(sf.torso*sf.shoulder)*sf.shoulder^2);
% 
% S.subject.set_limit_torque_coefficients_selected_dofs = ...
%     {{'arm_flex_r','arm_flex_l'},[-pass_torq_arm; 22; pass_torq_arm; -22], [-1, 1]};
% 
% S.subject.scale_default_coord_lim_torq = sf_mass*sf_legLength;
% 
% S.subject.tendon_stiff_scale = {{'soleus_l','soleus_r','gastroc_r','gastroc_l'},0.7};
% 
% % gravity 
% import org.opensim.modeling.*
% 
% % % S.weights
% S.weights.E         = 0.05;
% % S.weights.E_exp     = ;
% S.weights.q_dotdot  = 1;
% S.weights.e_arm     = 10;
% S.weights.pass_torq = 0;
% S.weights.a         = 1;
% S.weights.slack_ctrl = 0.001;
% 
% 

% 
% S.subject.IG_selection = 'quasi-random';
S.subject.IG_selection = [pathRepo, '\OCP\IK_Guess_Full_GC.mot'];
S.subject.IG_selection_gaitCyclePercent = 100;

% S.subject.adapt_IG_pelvis_y =1;

S.subject.adapt_IG_pelvis_y = 1;
S.subject.IG_pelvis_y = U.IG_pelvis_y;
% 
% 
S.bounds.Qs = {
    'pelvis_tx',0,4*sf_legLength,...
    {'arm_flex_r','arm_flex_l'},-60,60
    };

S.bounds.Qdots = {
    'pelvis_tx',0.01,max(S.subject.v_pelvis_x_trgt)*2
    };

if max(S.subject.v_pelvis_x_trgt) > 4
    S.bounds.Qdotdots = {'pelvis_tilt',-2000,2000};%v=4.5;
end

S.misc.scaling_Moments = {'all',1.2*sf_mass*sf_legLength};
S.subject.damping_coefficient_all_dofs = 0.1*sf_mass*sf_legLength;