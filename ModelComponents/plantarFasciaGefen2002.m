function [force] = plantarFasciaGefen2002(cross_section_area,slack_length,PF_length)
% --------------------------------------------------------------------------
% plantarFasciaGefen2002
%   This function contains the force-length characteristic of the plantar fascia, 
%   based on A. Gefen, “Stress analysis of the standing foot following 
%   surgical plantar fascia release,” Journal of biomechanics, vol. 35, 
%   no. 5, pp. 629–637, 2002.
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
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


elongation = PF_length - slack_length;
lambda = PF_length/slack_length;

% Poisson ratio
nu = 0.4;
A = cross_section_area*lambda.^(-nu*2); % actual cross-section

% model parameters
a1 = -488737.9;
a2 = 2648898.5;
a3 = -5736967.6;
a4 = 6206986.7;
a5 = -3354935.1;
a6 = 724755.5;

sigma = a1*lambda.^5 + a2*lambda.^4 + a3*lambda.^3 + a4*lambda.^2 + a5*lambda + a6 -0.100; % stress correction term, to make F=0 for l=ls
F_PF = sigma.*A; 

% Set compressive forces to zero, smoothed
force = F_PF.*(tanh(elongation*3e3-1)+1)/2;


    

end