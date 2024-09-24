function [force] = ligamentForceLength_template(cross_section_area,slack_length,lig_length)
% --------------------------------------------------------------------------
% ligamentForceLength_template
%   This function serves as a template for defining your own ligament
%   force-length characteristic. To use this force-lenght, pass it as
%   argument to "S.subject.stiffness_all_ligaments" or
%   "S.subject.set_stiffness_selected_ligaments". 
%   Note that the inputs and outputs of this function should be preserved.
%
%   Your force-length should not include conditional statements or splines.
%   If you do not want compressive forces, you should handle this here.
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
%
% OUTPUT:
%   - force -
%   * tensile force, in N
%
% 
% Original author: (First name Last name)
% Original date: (Using "30/May/2022" format avoids confusion)
% --------------------------------------------------------------------------

Youngs_modulus = 500; % [MPa] (= [N/mm^2])

strain = (lig_length - slack_length)/slack_length; % [-]

stress = Youngs_modulus*strain; % [MPa]

F = cross_section_area*stress; % [N]
    
% make sure there is no force at negative elongation
% e.g.  strain < 0: force = 0
%       strain > 2%: force = F
force = F.*smoothIf(strain, 0.02, 0); % [N]


end