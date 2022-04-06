function [S] = initializeSettings()
% --------------------------------------------------------------------------
%initializeSettings
%   This function creates the empty settings struct S up to the field 
%   above the field containing data. 
% 
% INPUT:
%   * no input arguments
% 
% OUTPUT:
%   - S -
%   * empty settings struct S
%
% Original author: Bram Van Den Bosch
% Original date: 01/12/2021
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 01/12/2021
% --------------------------------------------------------------------------

S = struct;

S.metabolicE   = [];
S.misc         = [];
S.post_process = [];
S.solver       = [];
S.subject      = [];
S.weights      = [];

% bounds have an .upper and .lower field
S.bounds.a          = [];
S.bounds.calcn_dist = [];
S.bounds.femur_hand_dist = [];
S.bounds.toes_dist  = [];
S.bounds.tibia_dist = [];
S.bounds.SLL        = [];
S.bounds.SLR        = [];
S.bounds.dist_trav  = [];
S.bounds.t_final    = [];

% polynomial order has .lower and .upper field
S.misc.poly_order = [];

end