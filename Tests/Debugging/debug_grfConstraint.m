%--------------------------------------------------------------------------
% Script to test the ground reaction forces constraint
% Note: The reference result contains the stored GRFs at the collocation
% points for easy debugging, but these are not typically stored inside the 
% PredSim results structure.
%
% Author: Menthy Denayer
% Date: 18/November/2026
%
%--------------------------------------------------------------------------

clear all
clc
close all

%% Add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions'])

%% Load Reference Result
load(fullfile(pathRepo,'Tests','ReferenceResults','gait1018','gait1018_GRFconstraints.mat'));

%% Define Variables
S = R.S;                                                                    % settings structure
N = R.S.solver.N_meshes;                                                    % number of meshes
tau_root = [0, (4-sqrt(6))/10, (4+sqrt(6))/10, 1];                          % radau collocation points

tau_c = tau_root(2:end);
pct_window = ((0:N-1) + tau_c(:))/N; pct_window = pct_window(:);            % simulated window percent
pct_window = [0; pct_window(1:end-1)];

% GRFs at the collocation points 
GRFcoll_right = R.ground_reaction.rawGRF(1:2:end); 
GRFcoll_left = R.ground_reaction.rawGRF(2:2:end); 

% Shift to make sure first data point matches first mesh point of the 
% results. Otherwise, the final collocation point (t = 100%), which matches
% the first mesh point, receives the bound for the end of the gait cycle, 
% while it should receive the bound for t = 0 %.
GRFcoll_right_shift = [GRFcoll_left(end), GRFcoll_right(1:end-1)];          
GRFcoll_left_shift = [GRFcoll_right(end), GRFcoll_left(1:end-1)];

% Number of GRF constraints
Ncon = 0;
for c = 1:numel(S.bounds.GRFConstraints)
    con = S.bounds.GRFConstraints{c};
    Ncon = Ncon + numel(con.side);
end

%% Check Constraints
fig = figure;
set(gcf,"Units","centimeters")                                              % cm units for position
set(gcf,"Position",[0 0 8.89*3 8.89])                                       % IEEE 1-column: 8.89cm

t = tiledlayout(1,Ncon,"TileSpacing","compact");
t.InnerPosition = [0.04 0.08 0.95 0.75];

for c = 1:numel(S.bounds.GRFConstraints)
    con = S.bounds.GRFConstraints{c};
    conActive = find((pct_window >= con.gait_percent(1)) & (pct_window <= con.gait_percent(2)));

    % Change scale
    if(isfield(con, "unit"))
        if strcmp(con.unit,'BW')
            scale = model_info.mass*9.81; 
        elseif(strcmp(con.unit,'BM'))
            scale = model_info.mass;
        end
    else
        scale = 1;
    end

    % Choose side
    for s = 1:numel(con.side)
        if(strcmp(con.side(s), 'l'))
            GRFcoll = GRFcoll_left_shift;
            GRFcoll_con = GRFcoll_left_shift(conActive);                    % only when constraint is active
        elseif(strcmp(con.side(s),'r'))
            GRFcoll = GRFcoll_right_shift;
            GRFcoll_con = GRFcoll_right_shift(conActive);
        end
    
        % Create tiled figure
        nexttile
        hold on
        plot(pct_window, GRFcoll, "black", "LineWidth", 1.5)
        plot(pct_window(conActive), GRFcoll_con, "blue", "LineWidth", 1.5)
    
        % Apply lower bound constraints
        if isfield(con, "lower_bound")
            plotConstraint(con.lower_bound*scale, [], con.gait_percent, 1.5 * R.misc.body_weight)
            if(all(GRFcoll_con > con.lower_bound*scale - 10^(-R.S.solver.constr_viol_tol_ipopt)))
                fprintf('Constraint satisfied: Vertical GRF on %s side > %.2f N from %.1f%% to %.1f%% of simulated window\n', con.side(s), con.lower_bound * scale, con.gait_percent(1)*1e2, con.gait_percent(2)*1e2);
            else
                fprintf('Constraint violated: Vertical GRF on %s side < %.2f N from %.1f%% to %.1f%% of simulated window\n', con.side(s), con.lower_bound * scale, con.gait_percent(1)*1e2, con.gait_percent(2)*1e2);
            end
        end
        
        % Apply upper bound constraints
        if isfield(con, "upper_bound")
            plotConstraint([], con.upper_bound*scale, con.gait_percent, 1.5 * R.misc.body_weight)
            if(all(GRFcoll_con < con.upper_bound*scale + 10^(-R.S.solver.constr_viol_tol_ipopt)))
                fprintf('Constraint satisfied: Vertical GRF on %s side < %.2f N from %.1f%% to %.1f%% of simulated window\n', con.side(s), con.upper_bound * scale, con.gait_percent(1)*1e2, con.gait_percent(2)*1e2);
            else
                fprintf('Constraint violated: Vertical GRF on %s side > %.2f N from %.1f%% to %.1f%% of simulated window\n', con.side(s), con.upper_bound * scale, con.gait_percent(1)*1e2, con.gait_percent(2)*1e2);
            end
        end
    
        plot(R.time.mesh_GC(1:N)/R.time.mesh_GC(N+1), R.ground_reaction.(['GRF_' con.side(s)])(1:N,2), "Marker", "o", "MarkerFaceColor", "m", "LineStyle", "none", "MarkerSize", 3)
        xlim([0 1])
        ylim([-0.1 1.5 * R.misc.body_weight])
        xlabel("Simulated Window [%]", "FontWeight", "bold")
        if(c == 1)
            ylabel("Force [N]", "FontWeight", "bold")
        end
        title(['Vertical Force ' con.side(s)])
        hold off
    end
end

lg = legend(["Raw Force", "Active Constraint", "Bounds", "Mesh Points"], "Orientation", "horizontal", "Box", "off");
lg.Layout.Tile = "north";

% figure settings
set(findall(fig,'-property','FontSize'),'FontSize',8)                       % font size
set(0,"DefaultFigureColor","w")                                             % white background
set(0,"defaulttextinterpreter","tex")                                       % tex style font
set(0,"DefaultAxesFontName","Helvetica")                                    % times new roman font
set(gca,"Units","centimeters")                                              % cm units for position

%% Helper Functions
% Function to plot the constrained area
function plotConstraint(lower_bound, upper_bound, gait_percent, ylim)

    if(~isempty(lower_bound))
        fill(repelem(gait_percent,1,2), [0 lower_bound lower_bound 0], "red", "FaceAlpha", 0.3, "LineStyle", "none")
    end

    if(~isempty(upper_bound))
        fill(repelem(gait_percent,1,2), [ylim upper_bound upper_bound ylim], "red", "FaceAlpha", 0.3, "LineStyle", "none")
    end

end