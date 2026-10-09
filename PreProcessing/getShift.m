function tendon_stiff_shift = getShift(tendon_stiff)
% --------------------------------------------------------------------------
% getShift
%   This script returns the value used to shift the tendon force-length curve
%   when changing the tendon stiffness. 
%   With the standard stiffness (35), the shift is 0. For a different
%   stiffness, the curve is shifted so that the normalized tendon force is
%   the same as with the standard stiffness when the normalized tendon length
%   is 1.
%   
% INPUT:
%   - tendon_stiff -
%   * tendon stiffness (default 35)
%
% OUTPUT:
%   - tendon_stiff_shift -
%   * value with which the tendon force-length curve is shifted to
%   compesate the new stiffness
% 
% Original author: Antoine Falisse
% Original date: 
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

kT35 = 35;
tendon_stiff_shift = 0;
lTtilde = 1;
fse = (exp(kT35.*(lTtilde - 0.995)))/5 - 0.25 + tendon_stiff_shift; 
fse_kt35 = fse;

fse = (exp(tendon_stiff.*(lTtilde - 0.995)))/5 - 0.25 + tendon_stiff_shift; 
fse_kt = fse;

tendon_stiff_shift = fse_kt35-fse_kt;
fse = (exp(tendon_stiff.*(lTtilde - 0.995)))/5 - 0.25 + tendon_stiff_shift;

% if sum(abs(fse - fse_kt35)>1e-12) ~= 0
%     shift = NaN;
%     warning('Error in shift tendon force-length curve');
% end

end

