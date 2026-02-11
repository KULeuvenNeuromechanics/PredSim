function [guess] = adaptInitialGuess(varargin)
% --------------------------------------------------------------------------
% adaptInitialGuess
%   Adapts the kinematics at the initial guess with the aim of improving
%   convergence of the main simulation.
%   
%   The position, velocity and acceleration guesses are adapted by adding
%   an offset to each value. To the vertical position of the floating base,
%   we add a second offset which is constant for all timepoints.
%
%   The values of these offsets are optimised to minimise a multi-objective
%   costfunction. The different terms are (the squared 2-norm of):
%       - residuals: inverse dynamics forces and moments on the floating
%       base (scaled by body weight)
%       - collocation: error on the relation between derivatives (scaled by
%       NLP scaling)
%       - Qs: offsets on the positions (scaled by NLP scaling) Note that the
%       constant offset on floating base vertical position is not included
%       in the cost.
%       - Qdots: offsets on the velocities (scaled by NLP scaling)
%       - Qdotdots: offsets on the accelerations (scaled by NLP scaling)
% 
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - scaling -
%   * scale factors for all optimisation variables
% 
%   - guess -
%   * initial guess values for all optimisation variables
%
%   - bounds -
%   * boundaries for all scaled optimisation variables
%
% OUTPUT:
%   - guess -
%   * initial guess values for all optimisation variables
% 
% Original author: Lars D'Hondt
% Original date: 8/May/2024
% --------------------------------------------------------------------------


% [advanced use] detect call with alternative inputs
if nargin==5
    S = varargin{1};
    model_info = varargin{2};
    scaling = varargin{3};
    guess = varargin{4};
    bounds = varargin{5};

else
    [S, model_info, scaling, guess, bounds] = inputHelper(varargin{:});

end

%%

fprintf("\tStarting optimisation to adapt initial guess...\n\n")
t0=tic;


import casadi.*


% load external function
F  = external('F',replace(fullfile(S.misc.subject_path,S.misc.external_function),'\','/'));


% indices of vertical GRFs
idx_GRFy = [];
GRF_labels = {};
for s=string(fieldnames(model_info.ExtFunIO.GRFs)')
%     if ~contains(s,'_total')
        idx_GRFy(end+1) = model_info.ExtFunIO.GRFs.(s)(2);
        GRF_labels{end+1} = char(s);
%     end
end

BW = model_info.mass*9.81;

if isempty(model_info.ExtFunIO.jointi.base_vertical)
    model_info.ExtFunIO.jointi.base_vertical = model_info.ExtFunIO.coordi.pelvis_ty;
end

%%
d = size(guess.Qs_col,1)/(size(guess.Qs,1)-1);
[~,C,~,~] = CollocationScheme(d,'radau');

h = guess.tf/(size(guess.Qs_col,1)/d);

x_SX = SX.sym('x',size(guess.Qs_col,2),3);
dx_SX = SX.sym('dx',size(guess.Qs_col,2),3);
xk_SX = SX.sym('x',size(guess.Qs_col,2),1);
dxk_SX = SX.sym('dx',size(guess.Qs_col,2),1);
ddx_SX = SX.sym('ddx',size(guess.Qs_col,2),3);

dxi_SX = [xk_SX,x_SX]*C(:,2:end);
err_dx = (h*dx_SX - dxi_SX)./scaling.Qs';

ddxi_SX = [dxk_SX,dx_SX]*C(:,2:end);
err_ddx = (h*ddx_SX - ddxi_SX)./scaling.Qdots';

f_coll = Function('f_coll',{x_SX,dx_SX,ddx_SX, xk_SX,dxk_SX},{[err_dx;err_ddx]});

%%

% get initial guess
Qs_0 = guess.Qs_col';
Qdots_0 = guess.Qdots_col';
Qdotdots_0 = guess.Qdotdots_col';

% define variables for offsets
opti = Opti();
offset_vertical = opti.variable(1,1);
offset_Qs = opti.variable(size(Qs_0,1),size(Qs_0,2));
offset_Qdots = opti.variable(size(Qdots_0,1),size(Qdots_0,2));
offset_Qdotdots = opti.variable(size(Qdotdots_0,1),size(Qdotdots_0,2));

% add offsets
Qs = Qs_0 + offset_Qs;
Qs(model_info.ExtFunIO.jointi.base_vertical,:) =...
    Qs(model_info.ExtFunIO.jointi.base_vertical,:) + offset_vertical;
Qdots = Qdots_0 + offset_Qdots;
Qdotdots = Qdotdots_0 + offset_Qdotdots;

% unscale
Qs_nsc = Qs.*(scaling.Qs'*ones(1,size(Qs,2)));
Qdots_nsc = Qdots.*(scaling.Qdots'*ones(1,size(Qdots,2)));
Qdotdots_nsc = Qdotdots.*(scaling.Qdotdots'*ones(1,size(Qdotdots,2)));

% contrain bounds
opti.subject_to(bounds.Qs.lower'*ones(1,size(Qs,2)) < Qs < ...
    bounds.Qs.upper'*ones(1,size(Qs,2)));
opti.subject_to(bounds.Qdots.lower'*ones(1,size(Qdots,2)) < Qdots < ...
    bounds.Qs.upper'*ones(1,size(Qdots,2)));
opti.subject_to(bounds.Qdotdots.lower'*ones(1,size(Qdotdots,2)) < Qdotdots < ...
    bounds.Qs.upper'*ones(1,size(Qdotdots,2)));

% Create input vector for external function
F_ext_input = MX(model_info.ExtFunIO.input.nInputs, size(Qs,2));
F_ext_input(model_info.ExtFunIO.input.Qs.all,:) = Qs_nsc;
F_ext_input(model_info.ExtFunIO.input.Qdots.all,:) = Qdots_nsc;
F_ext_input(model_info.ExtFunIO.input.Qdotdots.all,:) = Qdotdots_nsc;

% Evaluate external function
[F_ext_output] = F(F_ext_input);

GRF_y = F_ext_output(idx_GRFy,:)./BW;

% cost of residuals
residuals = F_ext_output(model_info.ExtFunIO.jointi.floating_base,:)./BW;
obj_res = sumsqr(residuals)/numel(residuals) *S.weights.adaptIG_residuals;

% cost of collocation error
Qsk_nsc = Qs_nsc(:,[1,d:d:end-d]);
Qdotsk_nsc = Qdots_nsc(:,[1,d:d:end-d]);
coll_err = f_coll(Qs_nsc,Qdots_nsc,Qdotdots_nsc,Qsk_nsc,Qdotsk_nsc);
obj_coll = sumsqr(coll_err)/numel(coll_err) *S.weights.adaptIG_collocation;


% cost of offsets
obj_offset_Qs = sumsqr(offset_Qs)/numel(offset_Qs) *S.weights.adaptIG_Qs;
obj_offset_Qdots = sumsqr(offset_Qdots)/numel(offset_Qdots) *S.weights.adaptIG_Qdots;
obj_offset_Qdotdots = sumsqr(offset_Qdotdots)/numel(offset_Qdotdots) *S.weights.adaptIG_Qdotdots;


obj = obj_res + obj_coll + obj_offset_Qs + obj_offset_Qdots + obj_offset_Qdotdots;
opti.minimize(obj)


optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 10^(-S.solver.adaptIG_tol_ipopt);
optionssol.ipopt.hessian_approximation = 'limited-memory';
optionssol.ipopt.nlp_scaling_method = 'none';
optionssol.ipopt.max_iter = S.solver.adaptIG_max_iter;

opti.solver('ipopt',optionssol);

sol = opti.solve_limited();

%%
res_0 = opti.value(residuals,opti.initial).*BW;
Qs_nsc_0 = opti.value(Qs_nsc,opti.initial);
Qdots_nsc_0 = opti.value(Qdots_nsc,opti.initial);
Qdotdots_nsc_0 = opti.value(Qdotdots_nsc,opti.initial);
GRF_0 = opti.value(GRF_y,opti.initial);

res_sol = sol.value(residuals).*BW;
Qs_nsc_sol = sol.value(Qs_nsc);
Qdots_nsc_sol = sol.value(Qdots_nsc);
Qdotdots_nsc_sol = sol.value(Qdotdots_nsc);
GRF_sol = sol.value(GRF_y);


obj_res_0 = opti.value(obj_res,opti.initial);
obj_coll_0 = opti.value(obj_coll,opti.initial);
obj_offset_Qs_0 = opti.value(obj_offset_Qs,opti.initial);
obj_offset_Qdots_0 = opti.value(obj_offset_Qdots,opti.initial);
obj_offset_Qdotdots_0 = opti.value(obj_offset_Qdotdots,opti.initial);

obj_res_sol = sol.value(obj_res);
obj_coll_sol = sol.value(obj_coll);
obj_offset_Qs_sol = sol.value(obj_offset_Qs);
obj_offset_Qdots_sol = sol.value(obj_offset_Qdots);
obj_offset_Qdotdots_sol = sol.value(obj_offset_Qdotdots);

%%

fprintf("\n\tAdapting initial guess done. (Time elapsed: %.2fs)\n\n",toc(t0))
fprintf("\t\tCost breakdown for optimising initial guess adaptation:\n\t\tTerm\t\t\t\t|\tOriginal IG\t\tAdapted IG\n")
fprintf("\t\t--------------------------------------------------------\n")
fprintf("\t\tResiduals\t\t\t|\t%0.3d\t\t%0.3d\n",obj_res_0,obj_res_sol)
fprintf("\t\tCollocation error\t|\t%0.3d\t\t%0.3d\n",obj_coll_0,obj_coll_sol)
fprintf("\t\tQs adjustment\t\t|\t%0.3d\t\t\t\t%0.3d\n",obj_offset_Qs_0,obj_offset_Qs_sol)
fprintf("\t\tQdots adjustment\t|\t%0.3d\t\t\t\t%0.3d\n",obj_offset_Qdots_0,obj_offset_Qdots_sol)
fprintf("\t\tQdotdots adjustment\t|\t%0.3d\t\t\t\t%0.3d\n\n\n",obj_offset_Qdotdots_0,obj_offset_Qdotdots_sol)



%%

figure
tiledlayout('flow')
for i=1:size(res_sol,1)
    nexttile
    hold on
    title(model_info.ExtFunIO.coord_names.all{model_info.ExtFunIO.jointi.floating_base(i)},'Interpreter','none')
    xlabel('collocation point')
    if any(model_info.ExtFunIO.jointi.rotations == model_info.ExtFunIO.jointi.floating_base(i))
        ylabel('Residual (Nm)')
    else
        ylabel('Residual (N)')
    end
    plot(res_0(i,:),'LineWidth',2)
    plot(res_sol(i,:),'LineWidth',2)
end

for i=1:size(GRF_sol,1)
    nexttile
    hold on
    title(GRF_labels{i},'Interpreter','none')
    xlabel('collocation point')
    ylabel('GRF (% BW)')
    plot(GRF_0(i,:)*100,'LineWidth',2)
    plot(GRF_sol(i,:)*100,'LineWidth',2)
end
lg=legend({'original initial guess','adapted initial guess'},'FontSize',10);
lg.Layout.Tile = i+size(res_sol,1)+1;

figure
Qs_0_deg = Qs_nsc_0;
Qs_0_deg(model_info.ExtFunIO.jointi.rotations,:) = Qs_0_deg(model_info.ExtFunIO.jointi.rotations,:)*180/pi;
Qs_sol_deg = Qs_nsc_sol;
Qs_sol_deg(model_info.ExtFunIO.jointi.rotations,:) = Qs_sol_deg(model_info.ExtFunIO.jointi.rotations,:)*180/pi;
tiledlayout('flow')
for i=1:size(Qs_sol_deg,1)
    nexttile
    hold on
    title(model_info.ExtFunIO.coord_names.all{i},'Interpreter','none')
    xlabel('collocation point')
    if any(model_info.ExtFunIO.jointi.rotations == i)
        ylabel('Angle (°)')
    else
        ylabel('Position (m)')
    end
    plot(Qs_0_deg(i,:),'LineWidth',2)
    plot(Qs_sol_deg(i,:),'LineWidth',2)
end
sgtitle('Qs')
lg=legend({'original initial guess','adapted initial guess'},'FontSize',10);
lg.Layout.Tile = i+1;

figure
tiledlayout('flow')
for i=1:size(Qdots_nsc_sol,1)
    nexttile
    hold on
    title(model_info.ExtFunIO.coord_names.all{i},'Interpreter','none')
    xlabel('collocation point')
    if any(model_info.ExtFunIO.jointi.rotations == i)
        ylabel('Velocity (rad/s)')
    else
        ylabel('Velocity (m/s)')
    end
    plot(Qdots_nsc_0(i,:),'LineWidth',2)
    plot(Qdots_nsc_sol(i,:),'LineWidth',2)
end
sgtitle('Qdots')
lg=legend({'original initial guess','adapted initial guess'},'FontSize',10);
lg.Layout.Tile = i+1;

figure
tiledlayout('flow')
for i=1:size(Qdotdots_nsc_sol,1)
    nexttile
    hold on
    title(model_info.ExtFunIO.coord_names.all{i},'Interpreter','none')
    xlabel('collocation point')
    if any(model_info.ExtFunIO.jointi.rotations == i)
        ylabel('Acceleration (rad/s^2)')
    else
        ylabel('Acceleration (m/s^2)')
    end
    plot(Qdotdots_nsc_0(i,:),'LineWidth',2)
    plot(Qdotdots_nsc_sol(i,:),'LineWidth',2)
end
sgtitle('Qdotdots')
lg=legend({'original initial guess','adapted initial guess'},'FontSize',10);
lg.Layout.Tile = i+1;

%% write adapted guess to output, but save original values

guess.Qs_og = guess.Qs;
guess.Qs = sol.value(Qs_nsc(:,[1,d:d:end]))';

guess.Qdots_og = guess.Qdots;
guess.Qdots = sol.value(Qdots_nsc(:,[1,d:d:end]))';

guess.Qs_col_og = guess.Qs_col;
guess.Qs_col = sol.value(Qs)';

guess.Qdots_col_og = guess.Qdots_col;
guess.Qdots_col = sol.value(Qdots)';

guess.Qdotdots_col_og = guess.Qdotdots_col;
guess.Qdotdots_col = sol.value(Qdotdots)';

guess.residuals_og = res_0;
guess.residuals = res_sol;


end % end of function

function [S, model_info, scaling, guess, bounds] = inputHelper(options)
% helper function to use alternative inputs when calling adaptInitialGuess
% from outside PredSim.
arguments
    options.prevRes char {mustBeFile} = [];
    options.S struct = []
    options.model_info struct = [];
    options.scaling struct = [];
    options.guess struct = [];
    options.bounds struct = [];
    options.weight_residuals (1,1) double {mustBeNonnegative} = 10;
    options.weight_collocation (1,1) double {mustBeNonnegative} = 1;
    options.weight_Qs (1,1) double {mustBeNonnegative} = 1;
    options.weight_Qdots (1,1) double {mustBeNonnegative} = 1e-3;
    options.weight_Qdotdots (1,1) double {mustBeNonnegative} = 1e-3;
    options.solver_max_iter (1,1) double {mustBeNonnegative} = 100;
    options.solver_tol_ipopt (1,1) double {mustBeNonnegative} = 2;
end

% initialise S
S.weights = [];
S.solver = [];

% load previous result
if ~isempty(options.prevRes)
    load(options.prevRes, 'setup','model_info','R');
    
    if exist('R','var')
        S = R.S;
    else
        load(options.prevRes, 'S');
    end
    guess = setup.guess;
    scaling = setup.scaling;
    bounds = setup.bounds
end

if ~isempty(options.S)
    S = options.S;
end
if ~isempty(options.model_info)
    model_info = options.model_info;
end
if ~isempty(options.scaling)
    scaling = options.scaling;
end
if ~isempty(options.guess)
    guess = options.guess;
end
if ~isempty(options.bounds)
    bounds = options.bounds;
end

% read weights and solver options
for s=string(fieldnames(options))
    % weights
    if contains(s,'weight_')
        s_w = replace(s,'weight_','');
        if ~isfield(S.weights,s_w)
            S.weights.(s_w) = options.(s);
        end
    end
    % solver
    if contains(s,'solver_')
        s_s = replace(s,'solver_','');
        if ~isfield(S.solver,s_s)
            S.solver.(s_s) = options.(s);
        end
    end
end

if ~exist('guess','var')
    N = R.S.solver.N_meshes;
    d = 3;
    ncoord = model_info.ExtFunIO.jointi.nq.all;

    guess.Qs = zeros(N+1,ncoord);
    guess.Qdots = zeros(N+1,ncoord);
    guess.Qdotdots = zeros(N+1,ncoord);

    guess.Qs_col = zeros(N*d,ncoord);
    guess.Qdots_col = zeros(N*d,ncoord);
    guess.Qdotdots_col = zeros(N*d,ncoord);

end

if ~exist('scaling','var')
    scaling.Qs = ones(1,ncoord);
    scaling.Qdots = ones(1,ncoord);
    scaling.Qdotdots = ones(1,ncoord);
end

if ~exist('bounds','var')
    bounds.Qs.lower = -inf(1,ncoord);
    bounds.Qs.upper = inf(1,ncoord);
    bounds.Qdots.lower = -inf(1,ncoord);
    bounds.Qdots.upper = inf(1,ncoord);
    bounds.Qdotdots.lower = -inf(1,ncoord);
    bounds.Qdotdots.upper = inf(1,ncoord);
end

end % end of inputHelper
