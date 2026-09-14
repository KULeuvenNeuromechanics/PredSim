clear all
clc
close all

%% Test Joint Stiffness
% Script to plot the results of the joint stiffness computation. Plots are
% designed for the gait1018 model (i.e. 10 joints, 18 muscles). Results for
% other models (e.g. D'Hondt 2024 4 segments) can also be loaded, but (by
% default) plotted muscles will be limited to those from the gait1018
% model.
%
% Plots include:
%   o normalized derivatives of muscle-tendon-force curves
%   o partial derivative of muscle moment arms to joint angles over gait
%   o muscle stiffness over gait
%   o tendon stiffness over gait
%   o joint stiffness over gait
%   o contributions of muscles to joint stiffness over gait
%
% Joint Stiffness computation based on:
% Luis, I., Gutierrez-Farewik, E.M. & Afschrift, M. Why inverse simulations 
% overestimate optimal ankle exoskeleton assistance: the role of bi-articular 
% coordination and joint stiffness. J NeuroEngineering Rehabil 23, 224 (2026). 
% https://doi.org/10.1186/s12984-026-02106-3

%% Add Paths
[pathTests,~,~] = fileparts(mfilename('fullpath'));                         % path to the tests folder
[pathRepo,~,~] = fileparts(pathTests);                                      % path to the repository folder

addpath(fullfile(pathRepo,"CasadiFunctions/"))
addpath(fullfile(pathRepo,"PreProcessing/"))

%% Load Muscle-Force Parameters
load(fullfile(pathRepo,"CasadiFunctions/","Fpparam.mat"))
load(fullfile(pathRepo,"CasadiFunctions/","Ftparam.mat"))
load(fullfile(pathRepo,"CasadiFunctions/","Faparam.mat"))
load(fullfile(pathRepo,"CasadiFunctions/","Fvparam.mat"))

%% Define Plot Variables
fig_height = 8.89;                                                          % cm
fig_width = 8.89;                                                           % cm
lineWidth = 1.5;
colors = makeGroupColors(0.5, 10, 0.8, 0.3, 0.9);                           % plot colors

%% Choose Muscles
% Muscles from the gait1018 model, can be adapted to other muscles, but
% this might cause some plots to become less clear.
subMuscles = ["hamstrings", "bifemsh", "glut_max", "iliopsoas", ...
    "rect_fem", "vasti", "gastroc", "soleus", "tib_ant"];

%% Choose Results File
[res_file_name, res_file_dir] = uigetfile(".mat","Choose PredSim results file");

%% Load Data
res = load(fullfile(res_file_dir, res_file_name));

%% Define Variables
% headers
muscleNames = string(res.R.colheaders.muscles);
jointNames = string(res.R.colheaders.coordinates);

% muscle properties
tendon_stiff_scale = [res.model_info.muscle_info.parameters.tendon_stiff];
muscle_strength = [res.model_info.muscle_info.parameters.muscle_strength];
muscle_pass_stiff_shift = [res.model_info.muscle_info.parameters.muscle_pass_stiff_shift];
muscle_pass_stiff_scale = [res.model_info.muscle_info.parameters.muscle_pass_stiff_scale];
FMmax = [res.model_info.muscle_info.parameters.FMo];
lMopt = [res.model_info.muscle_info.parameters.lMo];
lTs = [res.model_info.muscle_info.parameters.lTs];
vMmax = [res.model_info.muscle_info.parameters.vMmax];

%% ========================================================================
%   Normalized partial derivatives of muscle-tendon force curves
% =========================================================================
%   These curves show the normalized muscle-tendon force curves and their
%   partial derivatives to the muscle fiber and tendon lengths, for
%   different activations and muscle contractile velocities, computed using
%   the function "ForceEquilibrium_dFtildeState_all_tendon" in PredSim.
% =========================================================================

% Define mesh   
Npoints = 100;                                                              % number of points
Niter = 10;                                                                 % number of cases per variable

lMtilde_dx = linspace(0, 2, Npoints);                                       % normalized muscle fibre length
vMtilde_dx = linspace(-1.5, 1.5, Niter);                                    % normalized muscle fibre velocities
vM_dx = vMtilde_dx .* vMmax(1);                                             % muscle fibre velocities
lTtilde_dx = linspace(0.99, 1.06, Npoints);                                 % normalized tendon lengths
lT_dx = lTtilde_dx .* lTs(1);                                               % tendon lengths
a_dx = linspace(0, 1, Niter);                                               % muscle activations

% save results
dFT = NaN(Npoints,Npoints);                                                 % tendon force partial derivative
dFM_dv = NaN(Npoints,Npoints);                                              % muscle force partial derivative (variable velocity)
dFM_da = NaN(Npoints,Npoints);                                              % muscle force partial derivative (variable activation)

FTtilde = NaN(Npoints,Npoints);                                             % normalized tendon force
FMtilde_dv = NaN(Npoints,Npoints);                                          % normalized muscle force (variable velocity)
FMtilde_da = NaN(Npoints,Npoints);                                          % normalized muscle force (variable activation)

% get MTU forces & partial derivatives
for v = 1:Niter
    for t = 1:Npoints
        % MTU forces (using function defined at the end of this code)
        [FTtilde(t,v), FMtilde_dv(t,v)] = f_muscle_mechanics(Fvparam, Fpparam, Faparam, Ftparam, lTtilde_dx(t), lMtilde_dx(t), vMtilde_dx(v), 1, 1,...
                1, 1, 1, res.R.S.misc.dampingCoefficient);

        [~, FMtilde_da(t,v)] = f_muscle_mechanics(Fvparam, Fpparam, Faparam, Ftparam, lTtilde_dx(t), lMtilde_dx(t), 0, a_dx(v), 1,...
            1, 1, 1, res.R.S.misc.dampingCoefficient);

        % Partial derivatives (using function implemented in PredSim)
        [dFT(t,v), dFM_dv(t,v)] = ...
            ForceEquilibrium_dFtildeState_all_tendon(1,lMtilde_dx(t),vM_dx(v),lT_dx(t),FMmax(1),lMopt(1),...
            lTs(1),vMmax(1),Ftparam,Fvparam,Fpparam,Faparam,1,1,1,35);

        [~, dFM_da(t,v)] = ...
            ForceEquilibrium_dFtildeState_all_tendon(a_dx(v),lMtilde_dx(t),0,lT_dx(t),FMmax(1),lMopt(1),...
            lTs(1),vMmax(1),Ftparam,Fvparam,Fpparam,Faparam,1,1,1,35);
    end
end

% normalize partial derivatives
dFTtilde = dFT / FMmax(1) * lTs(1);                                         % normalize tendon force
dFMtilde_da = dFM_da / FMmax(1) * lMopt(1);                                 % normalize muscle force
dFMtilde_dv = dFM_dv / FMmax(1) * lMopt(1);                                 % normalize muscle force

% create figure
fig = figure("Color","white");
set(gcf,"Units","centimeters")                                              % cm units for position
set(gcf,"Position",[0 0 fig_width*2 fig_height])                          
t = tiledlayout(2,3,"TileSpacing","tight");
t.InnerPosition = [0.08 0.1 0.87 0.85];

% muscle force (variable velocity)
nexttile
hold on
yline(0,"LineWidth",0.5,"Color","black")
for v = 1:10
    plot(lMtilde_dx, FMtilde_dv(:,v), "Color", colors(v,:))
end
xlim([0.3 1.8])
ylim([-0.2 4])
ylabel(["Normalized Muscle"; "Force [-]"], "FontWeight", "bold")
colormap(colors)
clim([-1.5 1.5]);
c = colorbar;
c.Label.String = "Normalized Muscle Velocity [-]";
c.Label.FontWeight = "bold";
hold off

% muscle force (variable activation)
nexttile
hold on
yline(0,"LineWidth",0.5,"Color","black")
for v = 1:10
    plot(lMtilde_dx, FMtilde_da(:,v), "Color", colors(v,:))
end
xlim([0.3 1.8])
ylim([-0.2 4])
ylabel(["Normalized Muscle"; "Force [-]"], "FontWeight", "bold")
colormap(colors)
clim([0, 1]);
c = colorbar;
c.Label.String = "Activation [-]";
c.Label.FontWeight = "bold";
hold off

% tendon force
nexttile
hold on
yline(0,"LineWidth",0.5,"Color","black")
plot(lTtilde_dx, FTtilde(:,1), "Color", colors(1,:))
ylabel(["Normalized Tendon"; "Force [-]"], "FontWeight", "bold")
ylim([-0.2 2])
hold off

% muscle force derivative (variable velocity)
nexttile
hold on
yline(0,"LineWidth",0.5,"Color","black")
for v = 1:10
    plot(lMtilde_dx, dFMtilde_dv(:,v), "Color", colors(v,:))
end
xlim([0.3 1.8])
xlabel("Normalized Muscle Length [-]", "FontWeight", "bold")
ylabel(["Normalized Muscle"; "Force Derivative [-]"], "FontWeight", "bold")
colormap(colors)
clim([-1.5 1.5]);
c = colorbar;
c.Label.String = "Normalized Muscle Velocity [-]";
c.Label.FontWeight = "bold";
hold off

% muscle force derivative (variable activation)
nexttile
hold on
yline(0,"LineWidth",0.5,"Color","black")
for v = 1:10
    plot(lMtilde_dx, dFMtilde_da(:,v), "Color", colors(v,:))
end
xlim([0.3 1.8])
xlabel("Normalized Muscle Length [-]", "FontWeight", "bold")
ylabel(["Normalized Muscle"; "Force Derivative [-]"], "FontWeight", "bold")
colormap(colors)
clim([0, 1]);
c = colorbar;
c.Label.String = "Activation [-]";
c.Label.FontWeight = "bold";
hold off

% tendon force derivative 
nexttile
hold on
yline(0,"LineWidth",0.5,"Color","black")
plot(lTtilde_dx, dFTtilde(:,1), "Color", colors(1,:))
xlabel("Normalized Tendon Length [-]", "FontWeight", "bold")
ylabel(["Normalized Tendon"; "Force Derivative [-]"], "FontWeight", "bold")
hold off

% figure settings
set(findall(fig,'-property','FontSize'),'FontSize',8)                       % font size
set(0,"DefaultFigureColor","w")                                             % white background
set(0,"defaulttextinterpreter","tex")                                       % tex style font
set(0,"DefaultAxesFontName","Helvetica")                                    % times new roman font

%% ========================================================================
%  Muscle moment arms & partial derivative
% =========================================================================
%   These curves show the muscle moment arms and their partial derivatives
%   to the relevant joint angles.
% =========================================================================

isMuscle = contains(muscleNames, subMuscles);
muscleNamesLim = muscleNames(isMuscle);

dM = res.R.muscles.dM(:,isMuscle,:);                                        % muscle moment arm
dMdq = res.R.joint_stiffness.drdtheta(:,isMuscle,:);                        % muscle moment arm partial derivative
MJmap = res.model_info.muscle_info.muscle_spanning_joint_info(isMuscle,:);  % muscles spanning joints info
isSpanning = MJmap > 0;                                                     % muscles spanning joints info
Ntiles = sum(MJmap,'all');                                                  % number of tiles in tiledlayout

% create figure
fig = figure("Color","white");
set(gcf,"Units","centimeters")                                              % cm units for position
set(gcf,"Position",[0 0 fig_width*3 fig_height*2])                         
t = tiledlayout(round(Ntiles/4),4,"TileSpacing","tight");
t.InnerPosition = [0.04 0.04 0.93 0.88];

% loop over muscles
for i = 1:length(muscleNamesLim)
    dMi = dM(:,i,isSpanning(i,:));
    dMdqi = dMdq(:,i,isSpanning(i,:));
    Nspanning = size(dMi,3);
    jointNamesi = jointNames(isSpanning(i,:));

    % loop over relevant joints
    for j = 1:Nspanning
        nexttile
        hold on
        plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100, dMi(:,j),"black","LineWidth",1.5);
        plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100, dMdqi(:,j),"red","LineWidth",1.5);
        title(strrep(muscleNamesLim(i),"_"," ") + "(" + strrep(jointNamesi(j),"_"," ") + ")")
        hold off
    end

end

lg = legend(["moment arm", "moment arm derivative"],"Box","off","Orientation","horizontal");
lg.Layout.Tile = "north";

% X-axis labels only on bottom row
xlabel(nexttile(Ntiles-3),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(Ntiles-2),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(Ntiles-1),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(Ntiles),  'Gait cycle [%]', 'FontWeight', 'bold');

% figure settings
set(findall(fig,'-property','FontSize'),'FontSize',8)                       % font size
set(0,"DefaultFigureColor","w")                                             % white background
set(0,"defaulttextinterpreter","tex")                                       % tex style font
set(0,"DefaultAxesFontName","Helvetica")                                    % times new roman font

%% ========================================================================
%  Muscle stiffness, force, activation, length and velocity over gait
% =========================================================================
%   These curves show the computed muscle stiffness, their force,
%   activation, length and velocity over gait. You can verify the muscle
%   stiffness qualitatively by comparing the muscle state with the plots of
%   muscle force and their partial derivatives.
% =========================================================================

% define variables
Fce  = res.R.muscles.Fce(:,isMuscle);                                       % muscle contractile force
Fpass = res.R.muscles.Fpass(:,isMuscle);                                    % muscle passive force
KM   = res.R.joint_stiffness.KM(:,isMuscle);                                % computed muscle stiffness
a    = res.R.muscles.a(:,isMuscle);                                         % muscle activations
lM   = res.R.muscles.lM(:,isMuscle);                                        % muscle fibre lengths
vM   = res.R.muscles.vM(:,isMuscle);                                        % muscle fibre velocities
FMtot = Fce + Fpass;                                                        % total muscle fibre force

FMmaxLim = FMmax(:,isMuscle);                                               % maximal isometric force for desired muscles
lMoptLim = lMopt(:,isMuscle);                                               % optimal fibre length for desired muscles
vMmaxLim = vMmax(:,isMuscle);                                               % maximal contractile velocity for desired muscles

% figure
figure("Color","white");
set(gcf,"Units","centimeters")                                              % cm units for position
set(gcf,"Position",[0 0 fig_width*3 fig_height*2])                          
t = tiledlayout(length(muscleNamesLim)/2,5,"TileSpacing","tight");
t.InnerPosition = [0.04 0.05 0.94 0.92];

% plot for right side muscles
for m = 1:length(muscleNamesLim)/2
    % muscle stiffness
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,KM(:,m)/FMmaxLim(m)*lMoptLim(m),'LineWidth',1.5,"Color",colors(1,:));
    xlim([0 100])
    ylim([0 1.5])
    ylabel(strrep(muscleNamesLim(m),"_"," "), 'FontWeight','bold');
    hold off

    if m == 1
        title('Muscle Stiffness K_M');
    end

    % force
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,FMtot(:,m)/FMmaxLim(m),'LineWidth',1.5,"Color",colors(end,:));
    xlim([0 100])
    ylim([0 1])
    hold off

    if m == 1
        title('Normalized Muscle Force');
    end

    % activation
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,a(:,m),'LineWidth',1.5,"Color",colors(end,:));
    xlim([0 100])
    ylim([0 1])
    hold off

    if m == 1
        title('Muscle Activation');
    end

    % length
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,lM(:,m)/lMoptLim(m),'LineWidth',1.5,"Color",colors(end,:));
    xlim([0 100])
    ylim([0.45 1.1])
    hold off

    if m == 1
        title('Normalized Muscle Fibre Length');
    end

    % velocity
    nexttile
    hold on
    yline(0,"LineWidth",0.5,"Color","black")
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,vM(:,m)/vMmaxLim(m),'LineWidth',1.5,"Color",colors(end,:));
    xlim([0 100])
    ylim([-0.3 0.3])
    hold off

    if m == 1
        title('Normalized Muscle Fibre Velocity');
    end

end

% X-axis labels only on bottom row
xlabel(nexttile(5*length(muscleNamesLim)/2-4),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(5*length(muscleNamesLim)/2-3),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(5*length(muscleNamesLim)/2-2),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(5*length(muscleNamesLim)/2-1),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(5*length(muscleNamesLim)/2),  'Gait cycle [%]', 'FontWeight', 'bold');

% figure settings
set(findall(fig,'-property','FontSize'),'FontSize',8)                       % font size
set(0,"DefaultFigureColor","w")                                             % white background
set(0,"defaulttextinterpreter","tex")                                       % tex style font
set(0,"DefaultAxesFontName","Helvetica")                                    % times new roman font

%% ========================================================================
%  Tendon force, stiffness, length over gait
% =========================================================================
%   These curves show the computed tendon stiffness, their force,
%   length and force vs. length over gait. You can verify the tendon
%   stiffness qualitatively by comparing the tendon state with the plots of
%   tendon force and their partial derivatives.
% =========================================================================

lTsLim      = lTs(:,isMuscle);                                              % tendon slack length for desired muscles
FTtilde     = res.R.muscles.FTtilde(:,isMuscle);                            % tendon force
lTtilde     = res.R.muscles.lT(:,isMuscle)./lTsLim;                         % normalized tendon length
KT          = res.R.joint_stiffness.KT(:,isMuscle)./FMmaxLim.*lTsLim;       % tendon stiffness

% figure
figure("Color","white");
set(gcf,"Units","centimeters")                                              % cm units for position
set(gcf,"Position",[0 0 fig_width*3 fig_height*2])                          
t = tiledlayout(length(muscleNamesLim)/2,4,"TileSpacing","tight");
t.InnerPosition = [0.04 0.05 0.94 0.92];

% plot for right side muscles
for m = 1:length(muscleNamesLim)/2
    % tendon stiffness
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,KT(:,m),'LineWidth',1.5,"Color",colors(1,:));
    xlim([0 100])
    ylim([0 50])
    ylabel(strrep(muscleNamesLim(m),"_"," "), 'FontWeight','bold');
    hold off

    if m == 1
        title('Tendon Stiffness K_T');
    end

    % force
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,FTtilde(:,m),'LineWidth',1.5,"Color",colors(end,:));
    xlim([0 100])
    ylim([0 1])
    hold off

    if m == 1
        title('Normalized Tendon Force');
    end

    % length
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100,lTtilde(:,m),'LineWidth',1.5,"Color",colors(end,:));
    xlim([0 100])
    ylim([1 1.1])
    hold off

    if m == 1
        title('Normalized Tendon Fibre Length');
    end

    % force/stiffness vs length
    nexttile
    hold on
    yyaxis left
    plot(lTtilde(:,m),FTtilde(:,m),'LineWidth',1.5,"Color",colors(end,:));
    ylim([0 1])

    yyaxis right
    plot(lTtilde(:,m),KT(:,m),'LineWidth',1.5,"Color",colors(1,:));
    ylim([0 30])
    xlim([min(lTtilde(:,m)), max(lTtilde(:,m))])
    hold off

    ax = gca;
    ax.YAxis(1).Color = colors(end,:);
    ax.YAxis(2).Color = colors(1,:);

    if m == 1
        title('Normalized Tendon Force/K_T');
    end

end

% X-axis labels only on bottom row
xlabel(nexttile(4*length(muscleNamesLim)/2-3),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(4*length(muscleNamesLim)/2-2),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(4*length(muscleNamesLim)/2-1),'Gait cycle [%]', 'FontWeight', 'bold');
xlabel(nexttile(4*length(muscleNamesLim)/2),  'Normalized Tendon Length', 'FontWeight', 'bold');

% figure settings
set(findall(fig,'-property','FontSize'),'FontSize',8)                       % font size
set(0,"DefaultFigureColor","w")                                             % white background
set(0,"defaulttextinterpreter","tex")                                       % tex style font
set(0,"DefaultAxesFontName","Helvetica")                                    % times new roman font

%% ========================================================================
%  Joint stiffness over gait
% =========================================================================
%   These curves show the computed joint stiffness over the gait cycle.
% =========================================================================

% choose joints for plotting joint stiffness
subJoints = ["hip_flexion_r", "knee_angle_r", "ankle_angle_r", "hip_flexion_l", "knee_angle_l", "ankle_angle_l"];
[~,isJoint] = ismember(subJoints,jointNames); isJoint = isJoint(isJoint > 0);

K_J = res.R.joint_stiffness.KJ;                                             % joint stiffness
K_J_lim = K_J(:,isJoint);                                                   % limited to desired joints
Nsubjoints = length(subJoints);

% figure
fig = figure("Color","white");
set(gcf,"Units","centimeters")                                              % cm units for position
set(gcf,"Position",[0 0 fig_width*2 fig_height])                            
t = tiledlayout(2,3,"TileSpacing","tight");
t.InnerPosition = [0.06 0.09 0.92 0.85];

for j = 1:Nsubjoints
    nexttile
    hold on
    plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100, K_J_lim(:,j)*pi/180,"black", "LineWidth", 1.5)
    title(strrep(subJoints(j),"_", " "))
    
    if(any(j == [1,4]))
        ylabel("Joint Stiffness [Nm/deg]", "FontWeight", "bold")
    end
    
    if(j > 3)
        xlabel("Gait Cycle [%]", "FontWeight", "bold")
    end
    xlim([0 100])
    ylim([0 2.5])
    hold off
end

% figure settings
set(findall(fig,'-property','FontSize'),'FontSize',8)                       % font size
set(0,"DefaultFigureColor","w")                                             % white background
set(0,"defaulttextinterpreter","tex")                                       % tex style font
set(0,"DefaultAxesFontName","Helvetica")                                    % times new roman font

%% ========================================================================
%  Joint stiffness muscle contributions 
% =========================================================================
%   These curves show the computed joint stiffness and the contribution
%   from the relevant muscles. Muscles with a contribution smaller than 10%
%   to the total joint stiffness are denoted as "other".
% =========================================================================

threshold = 10;                                                             % threshold for joint stiffness contribution (default 10)

% choose joints for plotting joint stiffness
subJoints = ["hip_flexion_r", "knee_angle_r", "ankle_angle_r", "hip_flexion_l", "knee_angle_l", "ankle_angle_l"];
[~,isJoint] = ismember(subJoints,jointNames); isJoint = isJoint(isJoint > 0);

K_M_J = res.R.joint_stiffness.KMJ;                                          % joint stiffness per muscle
K_J = res.R.joint_stiffness.KJ;                                             % joint stiffness
K_M_J_lim = K_M_J(:,:,isJoint);                                             % limited to desired joints
K_J_lim = K_J(:,isJoint);                                                   % limited to desired joints
Nsubjoints = length(subJoints);

% figure
fig = figure("Color","white");
set(gcf,"Units","centimeters")                                              % cm units for position
set(gcf,"Position",[0 0 fig_width*3 fig_height])                           
t = tiledlayout(2,Nsubjoints,"TileSpacing","tight");
t.InnerPosition = [0.04 0.13 0.95 0.80];

for j = 1:Nsubjoints
    isNonZero = ~all(K_M_J_lim(:,:,j) == 0, 1);                             % muscles contributing to joint stiffness
    nonZeroIdxs = find(isNonZero);                                          % indexes of relevant muscles

    if(~isempty(nonZeroIdxs))
        % compute area under curve of joint stiffness
        totalArea = trapz(res.R.time.mesh_GC(1:end-1), sum(K_M_J_lim(:,:,j),2)*pi/180);
        areaList = zeros(length(nonZeroIdxs), 1);
        for i = 1:length(nonZeroIdxs)
            areaList(i) = trapz(res.R.time.mesh_GC(1:end-1), K_M_J_lim(:,nonZeroIdxs(i),j)*pi/180);
        end

        % muscle contributions
        contributions = areaList/totalArea*100;
        nonZeroIdxsLim = nonZeroIdxs(contributions > threshold);                % save muscles that contribution > 10%
        contributionsLim = contributions(contributions > threshold);
        contributionsOther = sum(contributions(contributions <= threshold));    % classify contributions < 10% as other
    
        % sort to plot largest to smallest area
        [a,b] = sort(contributionsLim, 1, "descend");
        K_M_J_right_sorted = K_M_J_lim(:,nonZeroIdxsLim,j)*pi/180;
        K_M_J_right_sorted = K_M_J_right_sorted(:,b);
    
        colors = makeGroupColors(0.5, length(nonZeroIdxsLim), 0.8, 0.5, 0.9);
    
        Y = [];
        for i = 1:length(nonZeroIdxsLim)
            Y = [Y K_M_J_right_sorted(:,i)];
        end
        Y = [Y sum(K_M_J_lim(:, nonZeroIdxs(contributions<=threshold), j),2)*pi/180];

        % % plot stacked area curve
        nexttile
        hold on
        h = area(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100, Y);                           
        plot(res.R.time.mesh_GC(1:end-1)/res.R.time.mesh_GC(end-1)*100, K_J_lim(:,j)*pi/180, "black", "LineWidth", 1.5)
        
        for i = 1:length(nonZeroIdxsLim)
            h(i).FaceColor = colors(i,:);
        end
        h(end).FaceColor = [0.9 0.9 0.9];

        title(strrep(subJoints(j),"_", " "))
        xlim([0 100])
        ylim([0 2.5])
        xlabel("Gait Cycle [%]", "FontWeight", "bold")
        ylabel("Joint Stiffness [Nm/deg]", "FontWeight", "bold")
        hold off

        % make muscle contributions bar plot
        nexttile
        hold on
        bHandle = bar(1:length(nonZeroIdxsLim)+1, [contributionsLim(b); contributionsOther], "FaceColor", "flat");
        xticklabels([strrep(muscleNames(nonZeroIdxsLim(b)),"_"," ") "other"])
        xticks(1:length(nonZeroIdxsLim)+1)

        for i = 1:length(nonZeroIdxsLim)
            bHandle.CData(i,:) = colors(i,:);
        end
        bHandle.CData(end,:) = [0.9 0.9 0.9];
        ylabel("Contribution [%]", "FontWeight", "bold")
        ylim([-15 100])
        hold off
    end
end

% figure settings
set(findall(fig,'-property','FontSize'),'FontSize',8)                       % font size
set(0,"DefaultFigureColor","w")                                             % white background
set(0,"defaulttextinterpreter","tex")                                       % tex style font
set(0,"DefaultAxesFontName","Helvetica")                                    % times new roman font
set(gca,"Units","centimeters")                                              % cm units for position

%% Functions
% function to compute muscle-tendon forces
function [fse, Fmtilde] = f_muscle_mechanics(Fvparamf, Fpparamf, Faparamf,...
    Ftparam, lTtildef, lMtildef, vMtilde, a, tendon_stiff, strength, stiffness_shift, stiffness_scale, damping_coeff)

% tendon force-length characteristic
kt = Ftparam(1);
c1 = Ftparam(2);
c2 = Ftparam(3);
c3 = Ftparam(4);

tendon_stiff = kt*tendon_stiff;

fse = c1*exp((lTtildef - c2).*tendon_stiff) - c3 + getShift(tendon_stiff);

% Active muscle force-length characteristic
b11 = Faparamf(1);
b21 = Faparamf(2);
b31 = Faparamf(3);
b41 = Faparamf(4);
b12 = Faparamf(5);
b22 = Faparamf(6);
b32 = Faparamf(7);
b42 = Faparamf(8);
b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtildef-b23;
den3 = b33+b43*lMtildef;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);
num1 = lMtildef-b21;
den1 = b31+b41*lMtildef;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);
num2 = lMtildef-b22;
den2 = b32+b42*lMtildef;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);
FMltilde = FMtilde1+FMtilde2+FMtilde3;
% Fiso = strength.*FMltilde;

% Passive muscle force-length characteristic
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtildef - stiffness_shift) / (e0/stiffness_scale));
% Passive muscle force
Fpetilde = ((t5 - 0.10e1) - Fpparamf(1)) / Fpparamf(2);

% Muscle force-velocity characteristic
d1 = Fvparamf(1);
d2 = Fvparamf(2);
d3 = Fvparamf(3);
d4 = Fvparamf(4);

FMvtilde = d1 * log((d2*vMtilde + d3) + sqrt((d2*vMtilde + d3).^2 + 1)) + d4;

% Active muscle force
Fcetilde = strength.*a.*FMltilde.*FMvtilde + damping_coeff*vMtilde;

Fmtilde = Fcetilde + Fpetilde;

end

% function to define gradient of plot colors
% S = saturation, Vmin/Vmax control light → dark
function colors = makeGroupColors(hue, n, S, Vmin, Vmax)
    V = linspace(Vmin, Vmax, n);        
    colors = hsv2rgb([hue*ones(n,1), S*ones(n,1), V']);
end