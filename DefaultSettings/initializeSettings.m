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
% Last edit date: 13/April/2023
% --------------------------------------------------------------------------

S = struct;

S.metabolicE   = [];
S.misc         = [];
S.post_process = [];
S.solver       = [];
S.subject      = [];
S.weights      = [];
S.OpenSimADOptions  = [];

% bounds have an .upper and .lower field
S.bounds.a          = [];
S.bounds.SLL        = [];
S.bounds.SLR        = [];
S.bounds.dist_trav  = [];
S.bounds.t_final    = [];

% polynomial order has .lower and .upper field
S.misc.poly_order = [];

% initiate for warning
S.subject.adapt_IG_pelvis_y = 0;

end