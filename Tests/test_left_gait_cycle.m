%% Minimal test script for get_left_gait_cycle
clear; clc; close all

result_path = 'C:\GBW_MyPrograms\PredSimSHARED\Tests\ReferenceResults\DHondt_et_al_2024_4seg\DHondt_et_al_2024_4seg_paper.mat';

data = load(result_path, 'R', 'model_info');
R_right = data.R;
model_info = data.model_info;

R_left = get_left_gait_cycle(R_right, model_info);

%% Plots: right vs. left gait cycle -- kinematics (Qs)
coordNames = R_right.colheaders.coordinates;
n_coords = length(coordNames);
n_cols = 5;
n_rows = ceil(n_coords/n_cols);

figure('Name','Kinematics: right vs left gait cycle')
for j = 1:n_coords
    subplot(n_rows,n_cols,j)
    plot(R_right.kinematics.Qs(:,j),'b','LineWidth',1.5); hold on
    plot(R_left.kinematics.Qs(:,j),'r','LineWidth',1.5);
    title(coordNames{j},'Interpreter','none')
end
legend({'right','left'})
sgtitle('Kinematics: right (original) vs. left gait cycle','Interpreter','none')

%% Plots: right vs. left gait cycle -- kinetics (T_ID)
figure('Name','Kinetics: right vs left gait cycle')
for j = 1:n_coords
    subplot(n_rows,n_cols,j)
    plot(R_right.kinetics.T_ID(:,j),'b','LineWidth',1.5); hold on
    plot(R_left.kinetics.T_ID(:,j),'r','LineWidth',1.5);
    title(coordNames{j},'Interpreter','none')
end
legend({'right','left'})
sgtitle('Kinetics: right (original) vs. left gait cycle','Interpreter','none')

%% Plots: right vs. left gait cycle -- ground reaction forces (GRF)
n_grf_cols = size(R_right.ground_reaction.GRF_r, 2); % typically x/y/z

figure('Name','GRF: right vs left gait cycle')
for j = 1:n_grf_cols
    subplot(2,n_grf_cols,j)
    plot(R_right.ground_reaction.GRF_r(:,j),'b','LineWidth',1.5); hold on
    plot(R_left.ground_reaction.GRF_r(:,j),'r','LineWidth',1.5);
    title(sprintf('GRF\\_r col %d', j),'Interpreter','none')

    subplot(2,n_grf_cols,n_grf_cols+j)
    plot(R_right.ground_reaction.GRF_l(:,j),'b','LineWidth',1.5); hold on
    plot(R_left.ground_reaction.GRF_l(:,j),'r','LineWidth',1.5);
    title(sprintf('GRF\\_l col %d', j),'Interpreter','none')
end
legend({'right','left'})
sgtitle('Ground reaction forces: right (original) vs. left gait cycle','Interpreter','none')