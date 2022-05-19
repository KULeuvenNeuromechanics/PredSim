function [R] = PostProcess_metabolic_energy(S,model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_metabolic_energy
%   This function computes the metabolic energy expenditure for the predicted
%   gait, according to different metabolic models.
%   Implemented models are:
%       * Bhargava et al. (2004)
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
%
%   - R -
%   * struct with simulation results
%
% OUTPUT:
%   - R -
%   * struct with simulation results
% 
% Original author: Lars D'Hondt
% Original date: 19/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

N = size(R.kinematics.Qs,1);
NMuscle = model_info.muscle_info.NMuscle;

%% Bhargava et al. (2004)
% Get metabolic energy rate 
pctsts = struct_array_to_double_array(model_info.muscle_info.parameters,'slow_twitch_fiber_ratio');
MuscleMass = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_mass');

R.metabolics.Bhargava2004.Edot_gait = zeros(N,NMuscle);
R.metabolics.Bhargava2004.Adot = zeros(N,NMuscle);
R.metabolics.Bhargava2004.Mdot = zeros(N,NMuscle);
R.metabolics.Bhargava2004.Sdot = zeros(N,NMuscle);
R.metabolics.Bhargava2004.Wdot = zeros(N,NMuscle);
R.metabolics.Bhargava2004.Edot_incl_basal = zeros(N,NMuscle);

for i=1:N
    [Edot_tot_i,Adot_i,Mdot_i,Sdot_i,Wdot_i,Edot_b_i] = f_casadi.getMetabolicEnergySmooth2004all(...
            R.muscles.a(i,:)',R.muscles.a(i,:)',R.muscles.lMtilde(i,:)',R.muscles.vM(i,:)',...
            R.muscles.Fce(i,:)',R.muscles.Fpass(i,:)',MuscleMass',pctsts,R.muscles.Fiso(i,:)',...
            S.subject.mass,S.metabolicE.tanh_b);

    R.metabolics.Bhargava2004.Edot_gait(i,:) = full(Edot_tot_i)';
    R.metabolics.Bhargava2004.Adot(i,:) = full(Adot_i)';
    R.metabolics.Bhargava2004.Mdot(i,:) = full(Mdot_i)';
    R.metabolics.Bhargava2004.Sdot(i,:) = full(Sdot_i)';
    R.metabolics.Bhargava2004.Wdot(i,:) = full(Wdot_i)';
    R.metabolics.Bhargava2004.Edot_incl_basal(i,:) = full(Edot_b_i)';

end



%%









