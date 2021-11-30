function [f_ArmActivationDynamics,f_TrunkActivationDynamics,...
    f_MtpActivationDynamics] = createCasadi_ActDynam()
%%createCasadi_ActDynam 
% Function to create Casadi functions for activation dynamics
% 
% INPUT:

% 
% OUTPUT:
%   f_ArmActivationDynamics
%   * description of this output
%
%   f_TrunkActivationDynamics
%   * description of this output
%
%   f_MtpActivationDynamics
%   * description of this output
% 
% Original author: Tom Buurke, KU Leuven
% Original date: 30/11/2021

%% Arm activation dynamics
e_a = SX.sym('e_a',nq.arms); % arm excitations
a_a = SX.sym('a_a',nq.arms); % arm activations
dadt = ArmActivationDynamics(e_a,a_a);
f_ArmActivationDynamics = ...
    Function('f_ArmActivationDynamics',{e_a,a_a},{dadt},...
    {'e','a'},{'dadt'});

%% Lumbar activation dynamics
e_lumb = SX.sym('e_a',nq.trunk); % Lumbar excitations
a_lumb = SX.sym('a_a',nq.trunk); % Lumbar activations
dlumbdt = ArmActivationDynamics(e_lumb,a_lumb);
f_TrunkActivationDynamics = ...
    Function('f_TrunkActivationDynamics',{e_lumb,a_lumb},{dlumbdt},...
    {'e','a'},{'dadt'});

%% Mtp activation dynamics
e_mtp = SX.sym('e_mtp',nq.mtp); % mtp excitations
a_mtp = SX.sym('a_mtp',nq.mtp); % mtp activations
dmtpdt = ArmActivationDynamics(e_mtp,a_mtp);
f_MtpActivationDynamics = ...
    Function('f_MtpActivationDynamics',{e_mtp,a_mtp},{dmtpdt},...
    {'e','a'},{'dadt'});
end

