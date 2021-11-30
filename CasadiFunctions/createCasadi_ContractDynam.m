function [f_forceEquilibrium_FtildeState_all_tendon,f_FiberLength_TendonForce_tendon...
    ,f_FiberVelocity_TendonForce_tendon,f_lT_vT] = createCasadi_ContractDynam(MainPath,NMuscle,IO)
%% createCasadi_ContractDynam.m
%Function to create Casadi functions for muscle contraction dynamics.
%
%INPUT
% - Mainpath
% - NMuscle
% - IO
%
%OUTPUT
%  - Casadi functions
%
% Authors: Ines Vandekerckhove & Tom Buurke, KU Leuven
% Date: 30-11-2021 

import casadi.*

%% Muscle contraction dynamics
% Function for Hill-equilibrium
FTtilde     = SX.sym('FTtilde',NMuscle); % Normalized tendon forces
a           = SX.sym('a',NMuscle); % Muscle activations
dFTtilde    = SX.sym('dFTtilde',NMuscle); % Time derivative tendon forces
lMT         = SX.sym('lMT',NMuscle); % Muscle-tendon lengths
vMT         = SX.sym('vMT',NMuscle); % Muscle-tendon velocities
tension_SX  = SX.sym('tension',NMuscle); % Tensions
% atendon_SX  = SX.sym('atendon',NMuscle); % Tendon stiffness
% shift_SX    = SX.sym('shift',NMuscle); % shift curve
Hilldiff    = SX(NMuscle,1); % Hill-equilibrium
FT          = SX(NMuscle,1); % Tendon forces
Fce         = SX(NMuscle,1); % Contractile element forces
Fiso        = SX(NMuscle,1); % Normalized forces from force-length curve
vMmax       = SX(NMuscle,1); % Maximum contraction velocities
massM       = SX(NMuscle,1); % Muscle mass
Fpass       = SX(NMuscle,1); % Passive element forces
% Parameters of force-length-velocity curves
load(fullfile(MainPath,'MuscleModel','Fvparam.mat'),'Fvparam');
load(fullfile(MainPath,'MuscleModel','Fpparam.mat'),'Fpparam');
load(fullfile(MainPath,'MuscleModel','Faparam.mat'),'Faparam');
% Parameters of force-length-velocity curves
for m = 1:NMuscle
    [Hilldiff(m),FT(m),Fce(m),Fpass(m),Fiso(m),vMmax(m),massM(m)] = ...
        ForceEquilibrium_FtildeState_all_tendon(a(m),FTtilde(m),...
        dFTtilde(m),lMT(m),vMT(m),IO.muscle.params.params(:,m),Fvparam,Fpparam,...
        Faparam,tension_SX(m),aTendon(m),shift(m),MuscMoAsmp);
end
f_forceEquilibrium_FtildeState_all_tendon = ...
    Function('f_forceEquilibrium_FtildeState_all_tendon',{a,FTtilde,...
    dFTtilde,lMT,vMT,tension_SX},{Hilldiff,FT,Fce,Fpass,Fiso,vMmax,massM},...
    {'a','FTtilde','dFTtilde','lMT','vMT','tension_SX'},...
    {'Hilldiff','FT','Fce','Fpass','Fiso','vMmax','massM'});

% Function to get (normalized) muscle fiber lengths
lM      = SX(NMuscle,1);
lMtilde = SX(NMuscle,1);
lT      = SX(NMuscle,1);
for m = 1:NMuscle
    [lM(m),lMtilde(m),lT(m)] = FiberLength_TendonForce_tendon(FTtilde(m),...
        IO.muscle.params.params(:,m),lMT(m),aTendon(m),shift(m),MuscMoAsmp);
end
f_FiberLength_TendonForce_tendon = Function(...
    'f_FiberLength_Ftilde_tendon',{FTtilde,lMT},{lM,lMtilde},...
    {'FTtilde','lMT'},{'lM','lMtilde'});

% Function to get (normalized) muscle fiber velocities
vM      = SX(NMuscle,1);
vMtilde = SX(NMuscle,1);
vT      = SX(NMuscle,1);
for m = 1:NMuscle
    [vM(m),vMtilde(m),vT(m)] = FiberVelocity_TendonForce_tendon(FTtilde(m),...
        dFTtilde(m),IO.muscle.params.params(:,m),lMT(m),vMT(m),aTendon(m),shift(m),MuscMoAsmp);
end
f_FiberVelocity_TendonForce_tendon = Function(...
    'f_FiberVelocity_Ftilde_tendon',{FTtilde,dFTtilde,lMT,vMT},...
    {vM,vMtilde},{'FTtilde','dFTtilde','lMT','vMT'},{'vM','vMtilde'});

f_lT_vT = Function('f_lT_vT',{FTtilde,dFTtilde,lMT,vMT},...
    {lT,vT},{'FTtilde','dFTtilde','lMT','vMT'},{'lT','vT'});
end