function [force] = plantarFasciaNatali2010(cross_section_area,slack_length,PF_length)
% --------------------------------------------------------------------------
% plantarFasciaNatali2010
%   This function contains the force-length characteristic of the plantar fascia, 
%   based on A. N. Natali, P. G. Pavan, and C. Stecco, "A constitutive model for 
%   the mechanical characterization of the plantar fascia,” Connective Tissue 
%   Research, vol. 51, no. 5, pp. 337–346, 2010, doi: 10.3109/03008200903389127.
%
% 
% INPUT:
%   - cross_section_area -
%   * cross section area of the plantar fascia, in mm^2
%
%   - slack_length -
%   * plantar fascia length at zero force, in m
%
%   - PF_length -
%   * plantar fascia length, in m
%
%
% OUTPUT:
%   - force -
%   * tensile force, in N
%
% 
% Original author: Lars D'Hondt
% Original date: 4/April/2023
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


elongation = PF_length - slack_length;
lambda = PF_length/slack_length;

% Poisson ratio
nu = 0.4;
A = cross_section_area*lambda.^(-nu*2); % actual cross-section

% model parameters
mu = 14.449; % (MPa)
k = 254.02; % (MPa)
alpha = 10.397; % (-)

sigma = mu*(lambda.^2 - 1./lambda) + k/(2*alpha) *(exp(alpha*(lambda.^2-1))-1).*lambda.^2; % Cauchy stress
F_PF = sigma.*A;

% Set compressive forces to zero, smoothed
force = F_PF.*(tanh(elongation*4e3-1.1)+1)/2;
    

end