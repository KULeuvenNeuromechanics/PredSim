function [w_opt,stats] = solve_NLPSOL(opti,optionssol)
% --------------------------------------------------------------------------
% solve_NLPSOL
%   (Explanation)
%   
% INPUT:
%   - opti -
%   * 
%
%   - optionssol -
%   * 
%
% OUTPUT:
%   - w_opt -
%   * 
%
%   - stats -
%   * 
% 
% Original author: 
% Original date: 
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*
% this will return the guess
X = opti.x;
guess=opti.debug.value(X,opti.initial);
F = opti.f;
Fguess=opti.debug.value(F,opti.initial);
G = opti.g;
Gguess=opti.debug.value(G,opti.initial);
% Sparsity pattern of the constraint jacobian (with respect to the variables)
% i.e. on which variables does each constraint depend
sp = sparse(DM(sparsity(jacobian(opti.g,opti.x)),1));

% Find the constraints that only depend on 1 variable
is_single = full(sum(sp,2)==1);
% Find constraints that contain (at least one) a linear dependency or no
% dependency
is_linear = ~which_depends(opti.g,opti.x,2,true)';
% Constraints that only depend linearly on 1 variable should be turned in
% to bounds. These are the ones we are looking for.
is_simple = is_single & is_linear;

% Find for each constraint that should become a bound which variable it is
% referring to
[col,~] = find(sp(is_simple,:)');

% Read out the values of the contraints from the opti instance
lbg = opti.lbg;
lbg = opti.value(lbg);
ubg = opti.ubg;
ubg = opti.value(ubg);

lb = lbg(find(is_simple));
ub = ubg(find(is_simple));

% % Detect  f2(p)x+f1(p)==0
% This is important if should you have scaled variables: x =
% 10*opti.variable(...) with a constraint -10 < x < 10. Because in the
% reformulation we read out the original variable and thus we need to scale
% the bounds appropriately.
g = opti.g;
gf = Function('gf',{opti.x,opti.p},{g(find(is_simple)),jtimes(g(find(is_simple)),opti.x,ones(opti.nx,1))});
[f1,f2] = gf(0,opti.p);
f1 = full(evalf(f1));
f2 = full(evalf(f2));

lb = (lbg(find(is_simple))-f1)./abs(f2);
ub = (ubg(find(is_simple))-f1)./abs(f2);

% Initialize the bound vector
lbx = -inf(opti.nx,1);
ubx = inf(opti.nx,1);
% Fill the bound vector at the correct indices with the correct values.
% For variable indices where the variables are unbounded we leave -inf, inf
% as bounds.
for i=1:numel(col)
    j = col(i);
    lbx(j) = max(lbx(j),lb(i));
    ubx(j) = min(ubx(j),ub(i));
end

lbx(col) = (lbg(find(is_simple))-f1)./abs(f2);
ubx(col) = (ubg(find(is_simple))-f1)./abs(f2);

% Generate a new constraint vector as now the variable constraints are
% turned into bounds
g = opti.g;
new_g = g(find(~is_simple));

% Generate the constraint value vector for the newly generated constraint
% vector
llb = lbg(find(~is_simple));
uub = ubg(find(~is_simple));

% Generate problem structure and intialize NLP solver
prob = struct('f', opti.f, 'x', opti.x, 'g',new_g );
solver = nlpsol('solver', 'ipopt',prob,optionssol);

%% Solve
sol = solver('x0',guess,'lbx',lbx,'ubx',ubx,'lbg',llb,'ubg',uub); 
w_opt = full(sol.x);
stats = solver.stats();
end
