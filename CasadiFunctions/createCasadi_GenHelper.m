function [f_casadi] = createCasadi_GenHelper(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_GenHelper
%   Function to create general Casadi functions.
%
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_casadi -
%   * Struct that contains all casadi functions.
%
% Original authors: Dhruv Gupta, Lars D'Hondt, Tom Buurke
% Original date: 01/12/2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

N_muscles = model_info.muscle_info.NMuscle;
N_arms_dof = model_info.ExtFunIO.jointi.nq.arms;
N_noarms_dof = model_info.ExtFunIO.jointi.nq.noArms;
N_torq_act = model_info.ExtFunIO.jointi.nq.torqAct;
N_pass_dof = model_info.ExtFunIO.jointi.nq.limTorq;

%% Normalized sum of squared values
if N_arms_dof > 0
    e_temp_arms_dof = SX.sym('e_temp_arms_dof',N_arms_dof);
    J_temp_arms_dof = 0;
    for i=1:length(e_temp_arms_dof)
        J_temp_arms_dof = J_temp_arms_dof + e_temp_arms_dof(i).^2;
    end
    J_temp_arms_dof = J_temp_arms_dof/N_arms_dof;
    f_casadi.J_arms_dof = Function('f_J_arms_dof',{e_temp_arms_dof},{J_temp_arms_dof});
end

if N_torq_act > 0
    e_temp_torq_act = SX.sym('e_temp_torq_act',N_torq_act);
    J_temp_torq_act = 0;
    for i=1:length(e_temp_torq_act)
        J_temp_torq_act = J_temp_torq_act + e_temp_torq_act(i).^2;
    end
    J_temp_torq_act = J_temp_torq_act/N_torq_act;
    f_casadi.J_torq_act = Function('f_J_torq_act',{e_temp_torq_act},{J_temp_torq_act});
end

e_temp_noarms_dof = SX.sym('e_temp_noarms_dof',N_noarms_dof);
J_temp_noarms_dof = 0;
for i=1:length(e_temp_noarms_dof)
    J_temp_noarms_dof = J_temp_noarms_dof + e_temp_noarms_dof(i).^2;
end
J_temp_noarms_dof = J_temp_noarms_dof/N_noarms_dof;
f_casadi.J_not_arms_dof = Function('f_J_not_arms_dof',{e_temp_noarms_dof},{J_temp_noarms_dof});

e_temp_pass_dof = SX.sym('e_temp_pass_dof',N_pass_dof);
J_temp_pass_dof = 0;
for i=1:length(e_temp_pass_dof)
    J_temp_pass_dof = J_temp_pass_dof + e_temp_pass_dof(i).^2;
end
J_temp_pass_dof = J_temp_pass_dof/N_pass_dof;
f_casadi.J_lim_torq = Function('f_J_lim_torq',{e_temp_pass_dof},{J_temp_pass_dof});

% Function for 2 elements
etemp2 = SX.sym('etemp2',2);
Jtemp2 = 0;
for i=1:length(etemp2)
    Jtemp2 = Jtemp2 + etemp2(i).^2;
end
Jtemp2 = Jtemp2/2;
f_casadi.J_2 = Function('f_J_2',{etemp2},{Jtemp2});

% Function for all muscles
e_temp_N_muscles = SX.sym('e_temp_N_muscles',N_muscles);
J_temp_N_muscles = 0;
for i=1:length(e_temp_N_muscles)
    J_temp_N_muscles = J_temp_N_muscles + e_temp_N_muscles(i).^2;
end
J_temp_N_muscles = J_temp_N_muscles/N_muscles;
f_casadi.J_muscles = Function('f_J_muscles',{e_temp_N_muscles},{J_temp_N_muscles});

%% Sum of squared values (non-normalized)
% Function for for distance between 2 points (in a certain plane)
e_temp_2 = SX.sym('e_temp_2',2);
J_temp_2 = 0;
for i=1:length(e_temp_2)
    J_temp_2 = J_temp_2 + e_temp_2(i).^2;
end
f_casadi.J_nn_2 = Function('f_J_nn_2',{e_temp_2},{J_temp_2});

%% Sum of squared values (non-normalized)
% Function for for distance between 2 points
e_temp_3 = SX.sym('e_temp_3',3);
J_temp_3 = 0;
for i=1:length(e_temp_3)
    J_temp_3 = J_temp_3 + e_temp_3(i).^2;
end
f_casadi.J_nn_3 = Function('f_J_nn_3',{e_temp_3},{J_temp_3});

%% Normalized sum of values to a certain power
% Function for number of muscles elements
e_temp_N_muscles_exp  = SX.sym('e_temp_N_muscles_exp',N_muscles);
expo        = SX.sym('exp',1);
J_temp_N_muscles_exp = 0;
for i=1:length(e_temp_N_muscles_exp)
    J_temp_N_muscles_exp = J_temp_N_muscles_exp + e_temp_N_muscles_exp(i).^expo;
end
J_temp_N_muscles_exp = J_temp_N_muscles_exp/N_muscles;
f_casadi.J_muscles_exp = Function('f_J_N_muscles_exp',{e_temp_N_muscles_exp,expo},{J_temp_N_muscles_exp});

%% Sum of products
% Function for number of muscles crossing a joint
sumCross = sum(model_info.muscle_info.muscle_spanning_joint_info);
N_musc_cross = setdiff(unique(sumCross),0);
for i = 1:length(N_musc_cross)
    ma_temp_musc_cross = SX.sym('ma_temp_musc_cross',N_musc_cross(i));
    ft_temp_musc_cross = SX.sym('ft_temp_musc_cross',N_musc_cross(i));
    J_sp_temp_musc_cross = 0;
    for j=1:length(ma_temp_musc_cross)
        J_sp_temp_musc_cross = J_sp_temp_musc_cross + ma_temp_musc_cross(j,1)*ft_temp_musc_cross(j,1);
    end
    f_casadi.(['musc_cross_' num2str(N_musc_cross(i))]) = Function(['musc_cross_' num2str(N_musc_cross(i))],{ma_temp_musc_cross,ft_temp_musc_cross},{J_sp_temp_musc_cross});
end
end