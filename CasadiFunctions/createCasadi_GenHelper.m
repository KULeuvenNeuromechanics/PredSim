function [f_casadi] = createCasadi_GenHelper(model_info,N_musc_cross)
%createCasadi_GenHelper
%   Function to create general Casadi functions.
%
% INPUT:
%   model_info
%   Struct with model info
%
%   N_musc_cross
%   Array with number of muscles crossing specific joinst, e.g. [27 13 12]
%   This will be part of model_info later, see issue on github.
%
% OUTPUT:
%   f_casadi
%   Struct that contains all casadi functions.
%
% Original authors: Dhruv Gupta, Lars D'Hondt, Tom Buurke
% Original date: 01/12/2021

import casadi.*

N_muscles = length(model_info.muscle_info.params.names);
N_non_musc_dof = []; %Zero elements in [# muscles crossing each joint]
N_musc_dof = []; %Non-zero elements in [# muscles crossing each joint]

%% Normalized sum of squared values
% Function for non muscle driven dofs ***(previously f_J8)***
e_temp_non_musc_dof = SX.sym('e_temp_non_musc_dof',N_non_musc_dof);
J_temp_non_musc_dof = 0;
for i=1:length(e_temp_non_musc_dof)
    J_temp_non_musc_dof = J_temp_non_musc_dof + e_temp_non_musc_dof(i).^2;
end
J_temp_non_musc_dof = J_temp_non_musc_dof/N_non_musc_dof;
f_casadi.J_non_musc_dof = Function('f_J_non_musc_dof',{e_temp_non_musc_dof},{J_temp_non_musc_dof});

% Function for muscle driven dofs ***(previously f_J23)***
e_temp_musc_dof = SX.sym('e_temp_musc_dof',N_musc_dof);
J_temp_musc_dof = 0;
for i=1:length(e_temp_musc_dof)
    J_temp_musc_dof = J_temp_musc_dof + e_temp_musc_dof(i).^2;
end
J_temp_musc_dof = J_temp_musc_dof/N_musc_dof;
f_casadi.J_musc_dof = Function('f_J_musc_dof',{e_temp_musc_dof},{J_temp_musc_dof});

% Function for all muscles
e_temp_N_muscles = SX.sym('e_temp_N_muscles',N_muscles);
J_temp_N_muscles = 0;
for i=1:length(e_temp_N_muscles)
    J_temp_N_muscles = J_temp_N_muscles + e_temp_N_muscles(i).^2;
end
J_temp_N_muscles = J_temp_N_muscles/N_muscles;
f_casadi.J_N_muscles = Function('f_J_N_muscles',{e_temp_N_muscles},{J_temp_N_muscles});

% Function for 2 elements
etemp2 = SX.sym('etemp2',2);
Jtemp2 = 0;
for i=1:length(etemp2)
    Jtemp2 = Jtemp2 + etemp2(i).^2;
end
Jtemp2 = Jtemp2/2;
f_casadi.J_2 = Function('f_J_2',{etemp2},{Jtemp2});

%% Sum of squared values (non-normalized)
% Function for for distance between 2 points
e_temp_2 = SX.sym('e_temp_2',2);
J_temp_2 = 0;
for i=1:length(e_temp_2)
    J_temp_2 = J_temp_2 + e_temp_2(i).^2;
end
f_casadi.J_nn_2 = Function('f_J_nn_2',{e_temp_2},{J_temp_2});

%% Normalized sum of values to a certain power
% Function for 92 elements
e_temp_N_muscles_exp  = SX.sym('e_temp_N_muscles_exp',N_muscles);
expo        = SX.sym('exp',1);
J_temp_N_muscles_exp = 0;
for i=1:length(e_temp_N_muscles_exp)
    J_temp_N_muscles_exp = J_temp_N_muscles_exp + e_temp_N_muscles_exp(i).^expo;
end
J_temp_N_muscles_exp = J_temp_N_muscles_exp/N_muscles;
f_casadi.J_N_muscles_exp = Function('f_J_N_muscles_exp',{e_temp_N_muscles_exp,expo},{J_temp_N_muscles_exp});

%% Sum of products
% Function for number of muscles crossing a joint
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