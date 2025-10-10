function [S] = scale4PredSim_2025(S,U,mass_genmodel,length_genmodel)
% Ellis Van Can October 2025
%%% scalefactors
S.subject.St = (U.Height*U.Mass)/(length_genmodel*mass_genmodel);
S.subject.Swt = 1/S.subject.St^2;
S.subject.Swe1 = (U.Mass/mass_genmodel)*(U.Height/length_genmodel)^2;
S.subject.Swe2 = 1/S.subject.Swe1;
sf_mass =   (U.Mass/mass_genmodel);
sf_force = sf_mass^(2/3);
sf_Length =(U.Height/length_genmodel);

% % cost function weights
% S.weights.E = 500*S.subject.Swe2;
% S.weights.pass_torq = 1000*S.subject.Swt;
% % S.weights
S.weights.E         = 0.05;
% S.weights.E_exp     = ;
S.weights.q_dotdot  = 1;
S.weights.e_arm     = 10;
S.weights.pass_torq = 0;
S.weights.a         = 1;
S.weights.slack_ctrl = 0.001;


% coordinate limit torques
             S.subject.set_limit_torque_coefficients_selected_dofs = {{'lumbar_extension'},[-0.7644*S.subject.St,11.2154,1.2788*S.subject.St,-7.2704], [-0.3716,0.1068],...
             {'hip_flexion_r','hip_flexion_l'},[-2.44*S.subject.St,5.05,1.51*S.subject.St,-21.88],[-0.6981,1.81],...
             {'knee_angle_r','knee_angle_l'},[-6.09*S.subject.St,33.94,11.03*S.subject.St,-11.33],[-2.4,0.13],... 
             {'ankle_angle_r','ankle_angle_l'},[-2.03*S.subject.St,38.11,0.18*S.subject.St,-12.12],[-0.74,0.52]};

% joint params
S.subject.damping_coefficient_all_dofs = 0.1*S.subject.St;
S.subject.scale_actuator_torque = {...
    {'lumbar_extension'},sf_force*(sqrt(sf_Length*sf_Length)*sf_Length*sf_Length)};
S.subject.tendon_stiff_scale = {{'soleus_l','soleus_r','gastroc_r','gastroc_l'},0.7};

S.misc.scaling_Moments = {'all',1.2*sf_mass*sf_Length};
end