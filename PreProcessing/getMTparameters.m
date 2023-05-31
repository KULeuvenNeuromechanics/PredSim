function [FMo, lMo, lTs, alphao, vMmax, fiber_damping, muscle_pass_stiff_scale,...
    tendon_stiffness] = getMTparameters(modelPath,muscleNames)
% --------------------------------------------------------------------------
% getMTparameters
%   This function returns the muscle-tendon parameters of the muscles
%   specified in muscleNames from the model in modelPath.
%   
% INPUT:
%   - modelPath -
%   * path to the OpenSim model file (.osim)
% 
%   - muscleNames -
%   * muscles for which we want to read the muscle-tendon parameters
%
% OUTPUT:
%   - FMo -
%   * maximum isometric force (N)
%
%   - lMo -
%   * optimal fiber length (m)
%
%   - lTs -
%   * tendon slack length (m)
%
%   - alphao -
%   * pennation angle at optimal fiber length (rad)
%
%   - vMmax -
%   * maximum contraltion velocity (m/s)
%
%   - fiber_damping -
%   * damping coefficient of muscle fiber (-). Set to NaN if the muscle is
%   not a "DeGrooteFregly2016Muscle".
%
%   - muscle_pass_stiff_scale -
%   * scale factor for normalised passive force-length of muscle fiber (-). 
%   For "DeGrooteFregly2016Muscle", this is calculated based on 
%   "passive_fiber_strain_at_one_norm_force", for other muscle models the
%   output is 1;
%
%   - tendon_stiffness -
%   * stiffness parameter for normalised force-length of tendon (-). 
%   For "DeGrooteFregly2016Muscle", this is calculated based on 
%   "tendon_strain_at_one_norm_force", for other muscle models the
%   output is 35;
%
%
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: Lars D'Hondt
% Last edit date: 30/May/2023
% --------------------------------------------------------------------------

import org.opensim.modeling.*

model = Model(modelPath);

NMuscles = length(muscleNames);
FMo = zeros(1,NMuscles);
lMo = zeros(1,NMuscles);
lTs = zeros(1,NMuscles);
alphao = zeros(1,NMuscles);
vMmax = zeros(1,NMuscles);
fiber_damping = nan(1,NMuscles);
muscle_pass_stiff_scale = ones(1,NMuscles);
tendon_stiffness = 35*ones(1,NMuscles);

muscles = model.getMuscles();

for i = 1:NMuscles
   muscle = muscles.get(muscleNames{i});
   FMo(i) = muscle.getMaxIsometricForce();
   lMo(i) = muscle.getOptimalFiberLength();
   lTs(i) = muscle.getTendonSlackLength();
   alphao(i) = muscle.getPennationAngleAtOptimalFiberLength();
   vMmax(i) = muscle.getMaxContractionVelocity()*lMo(i);

    if strcmp(muscle.getConcreteClassName(),'DeGrooteFregly2016Muscle')
        muscle = DeGrooteFregly2016Muscle.safeDownCast(muscle);
        
        fiber_damping(i) = muscle.get_fiber_damping();

        eMtilde_at_FMo = muscle.get_passive_fiber_strain_at_one_norm_force();
        muscle_pass_stiff_scale_i = 0.6/eMtilde_at_FMo;
        % avoid changes due to rounding errors in osim file
        if abs(muscle_pass_stiff_scale_i-1) < 0.01
            muscle_pass_stiff_scale_i = 1;
        end
        muscle_pass_stiff_scale(i) = muscle_pass_stiff_scale_i;
        
        eTtilde_at_FMo = muscle.get_tendon_strain_at_one_norm_force();
        tendon_stiffness_i = calc_tendon_stiffness_from_strain(eTtilde_at_FMo);
        % avoid changes due to rounding errors in osim file
        if abs(tendon_stiffness_i-35) < 0.01
            tendon_stiffness_i = 35;
        end
        tendon_stiffness(i) = tendon_stiffness_i;

    end
end

end % end of function