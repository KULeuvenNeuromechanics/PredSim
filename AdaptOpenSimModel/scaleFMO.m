function [model] = scaleFMO(muscle,model,mass_genmodel,mass_subject)
% --------------------------------------------------------------------------
%scaleFMO
%   This function scales maximal isometric force based on the mass of the
%   subject and the mass of gen model
% 
% INPUT:
%   -muscle-
%   * char name of muscle in which you want to change fMo
%
%   -model-
%   * org.opensim.modeling.Model
%
%   -statemodel1-
%   * org.opensim.modeling.State of model1
%
%   -mass_genmodel-
%   * mass of the generic model, if using gait2392 this is 75.1646 kg
%
%   -mass_subject
%   * mass of the subject, in kilograms
% 
% OUTPUT:
%   -model-
%   * updated org.opensim.modeling.Model
% 
% Original author: Bram Van Den Bosch
% Original date: 23/05/2022
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 23/05/2022
% --------------------------------------------------------------------------

% clear all
% 
% muscle = 'rect_fem_r';
% mass_genmodel = 75.1646;
% mass_subject = 30.15;

import org.opensim.modeling.*
% rootdir_model = 'C:\Users\u0138016\OneDrive - KU Leuven\SimCP_2\Subjects\CP3\T0\Model';
% name_model = 'mri_temp_v2.osim';
% filename = fullfile(rootdir_model,name_model);
% model = Model(filename);
% statemodel = model.initSystem;

%scale factor muscle-tendon parameters, based on van der Krogt et al.,
%(2016), doi: 10.1186/s12984-016-0170-5
scaleFactor = (mass_subject/mass_genmodel)^(2/3);

%define the muscle in the modified model for changing
currentMuscle = model.getMuscles.get(muscle);
% currentMuscle.dump

%define the new tendon slack length by multiplying current 
% muscle-tendon parameter by the scale factor
currentMuscle.setMaxIsometricForce(currentMuscle.getMaxIsometricForce*scaleFactor);
% currentMuscle.dump

model.updForceSet;

% %added
% currentMuscle.delete

end