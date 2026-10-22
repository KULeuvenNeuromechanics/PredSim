function [R] = PostProcess_metabolic_energy(model_info,f_casadi,R)
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
% Adapted by: Tim van der Zee
% Adapted on: 29/Sep/2026
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

N = size(R.kinematics.Qs,1);
NMuscle = model_info.muscle_info.NMuscle;

all_models = {'Bhargava2004', 'Bhargava2004_simple','Bhargava2004_original','Bhargava2004_revised', 'Umberger2003', 'Umberger2010', 'Uchida2016', 'Margaria'};

%% Get metabolic energy rate 
pctsts = struct_array_to_double_array(model_info.muscle_info.parameters,'slow_twitch_fiber_ratio');
MuscleMass = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_mass');

for j = 1:length(all_models)

    S.metabolicE.model = all_models{j};
    
    [f_casadi.fgetMetabolicEnergy] = createCasadi_E_Metab(S,model_info);
    
    R.metabolics.(all_models{j}).Edot_gait = zeros(N,NMuscle);
    R.metabolics.(all_models{j}).Adot = zeros(N,NMuscle);
    R.metabolics.(all_models{j}).Mdot = zeros(N,NMuscle);
    R.metabolics.(all_models{j}).Sdot = zeros(N,NMuscle);
    R.metabolics.(all_models{j}).Wdot = zeros(N,NMuscle);
    R.metabolics.(all_models{j}).Edot_incl_basal = zeros(N,1);

    for i=1:N
        [Edot_tot_i,Adot_i,Mdot_i,Sdot_i,Wdot_i,Edot_b_i] = f_casadi.getMetabolicEnergy(...
                R.muscles.a(i,:)',R.muscles.a(i,:)',R.muscles.lMtilde(i,:)',R.muscles.vM(i,:)',...
                R.muscles.Fce(i,:)',R.muscles.Fpass(i,:)',MuscleMass',pctsts,R.muscles.Fiso(i,:)',...
                R.misc.body_mass,R.S.metabolicE.tanh_b);

        R.metabolics.(all_models{j}).Edot_gait(i,:) = full(Edot_tot_i)';
        R.metabolics.(all_models{j}).Adot(i,:) = full(Adot_i)';
        R.metabolics.(all_models{j}).Mdot(i,:) = full(Mdot_i)';
        R.metabolics.(all_models{j}).Sdot(i,:) = full(Sdot_i)';
        R.metabolics.(all_models{j}).Wdot(i,:) = full(Wdot_i)';
        R.metabolics.(all_models{j}).Edot_incl_basal(i) = full(Edot_b_i)';

    end

    % cost of transport
    E_sum_GC = trapz(R.time.mesh_GC(1:end-1),R.metabolics.(all_models{j}).Edot_incl_basal);
    R.metabolics.(all_models{j}).COT = E_sum_GC/R.misc.body_mass/R.spatiotemp.dist_trav;

end

