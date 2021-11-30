function [S] = getDefaultSettings(S)
% --------------------------------------------------------------------------
%getDefaultSettings 
%   This functions sets default settings when the user didn't specify the
%   setting in main.m.
% 
% INPUT:
%   - S -
%   * setting structure S
% 
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: FULL NAME
% Original date: DD/MM/YYYY
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 30/11/2021
% --------------------------------------------------------------------------

%% bounds

% lower bounds on activations
if ~isfield(S,'bouds.a.lower')
    S.bouds.a.lower    = 0;        
end

% minimal distance between orginins calcanei

% minimal distance between origins toes

% minimal distance between origins tibiae

% upper bound on left step length

% upper bound on right step length

% lower bound on distance travelled

% upper bound on final time

% lower bound on final time


%% Initial Guess (IG) --> ask what the conclusion was about this

% name of the IG file

% folder path where the IG file is stored


%% metabolicE

% hyperbolic tangent smoothing factor (used in metabolic cost)

% name of the metabolic energy model


%% misc

% maximal contraction velocity identifier

% type of gait simulation

%% post_process

% boolean for making plots or not

% name used for saving the results (choose custom or structurized savename)


%% solver

% solver algorithm used in the OCP

% the power (10^-x) of the dual infeasibility for when the problem is 
% ‘solved’, higher number is more precise

% the amount of iterations after which the solver stops

% type of parellel computing

% number of threads in parallel mode

% number of mesh intervals


%% subject

% name of the subject, compare with opensim model

% mass of the subject, compare with opensim model

% height of the pelvis, compare with opensim model

% average velocity you want the model to have

% muscle strength, check with opensim model if muscles are present

% muscle stiffness, check with opensim model if muscles are present

% muscle symmetry, check with opensim model if muscles are present

% tendon stiffness, check with opensim model if muscles are present

% motion file for the initial guess

% type of mtp joint used in the model

% folder where you store all subject specific output

% muscle tendon properties


% will be added in the future:
% - spasticity
% - muscle coordination (co-contraction)


%% weights

% weight on metabolic energy rate

% exponent for the metabolic energy rate

% weight on joint accelerations

% weight on arm excitations

% weight on passive torques

% weight on activations

% weight on mtp excitations

% weight on slack controls



end

