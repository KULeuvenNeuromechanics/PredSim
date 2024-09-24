function [force] = ligamentGefen2002(cross_section_area,slack_length,lig_length,varargin)
% --------------------------------------------------------------------------
% ligamentGefen2002
%   This function contains the force-length characteristic of a ligament,
%   based on  A. Gefen, “Stress analysis of the standing foot following 
%   surgical plantar fascia release,” Journal of biomechanics, vol. 35, 
%   no. 5, pp. 629–637, 2002.
% 
%   We approximated the polynomial stress-strain curve with an exponential
%   expression, to prevent a drop in force when extrapolating. See
%   supplementary material of D’Hondt, L., De Groote, F., & Afschrift, M. 
%   2024 A dynamic foot model for predictive simulations of human gait 
%   reveals causal relations between foot structure and whole-body mechanics. 
%   PLOS Computational Biology, 20(6), e1012219. 
%   https://doi.org/10.1371/journal.pcbi.1012219.
%
% 
% INPUT:
%   - cross_section_area -
%   * cross section area of the ligament, in mm^2
%
%   - slack_length -
%   * ligament length at zero force, in m
%
%   - lig_length -
%   * ligament length, in m
%
%   - set_negative_force_to_zero (optional input) -
%   * Select how the negative forces should be treated
%       - 'exact' use if-statement (default if length is of type double)
%       - 'smooth' use smoothed if-statemen (default if length is not of
%       type double)
%       - 'none' leave them as is
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

%% Calculate force
% stretch ratio lambda
lambda = lig_length/slack_length;


coeff = [7.9884 8.2061 1.0018];
sigma = coeff(1)*(exp( coeff(2)*(lambda-coeff(3)) ) - 1);
force = sigma*cross_section_area;

%% Handle negative forces
% read preferred mode
set_negative_force_to_zero_options = {'exact','smooth','none'};
if length(varargin) >=1
    if sum(contains(set_negative_force_to_zero_options,varargin{1}))
        set_negative_force_to_zero = varargin{1};
    else
        error(["Invalid input for set_negative_force_to_zero. Please choose from 'exact','smooth','none'."])
    end
elseif isa(lig_length,'double')
    set_negative_force_to_zero = 'exact';
else
    set_negative_force_to_zero = 'smooth';
end

% apply
if strcmp(set_negative_force_to_zero,'exact')
    force(force<0) = 0;

elseif strcmp(set_negative_force_to_zero,'smooth')
    force = force.*smoothIf(sigma,0.1,0);

end





end