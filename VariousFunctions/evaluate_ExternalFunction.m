
% this script performs a quasi-static evaluation of the external function
% for a given motion file. Bases on
% PredSim/Tests/compare_external_functions.m

% original author: Bram Van Den Bosch
% original date: 11/May/2023
%
% Last edited by: 
% Last edit: 

% TO DO: not rely on previous solution but create model_info and necessary
% S settings if they do not exist yet

clear
clc
addpath(genpath('C\:GBW_MyPrograms\casadi_3_5_5'))
import casadi.*

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/PreProcessing'])
addpath([pathRepo '/OCP'])

%% input
name_1 = 'Falisse_et_al_2022';
post_process.result_filename = 'Falisse_et_al_2022_v1';
F1  = external('F',fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '.dll'])); 

% load IO
load(fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '_IO.mat'])); 

% path to motion file to use
motion_file = fullfile(pathRepo,'OCP','IK_Guess_Full_GC.mot');

% outfile with model_info
S.subject.save_folder  = fullfile(pathRepo,'Tests',[name_1,'_Results']);

% load S and model_info
results_file = fullfile(S.subject.save_folder,[post_process.result_filename '.mat']);
load(results_file,'R','model_info');
S = R.S;

% coordinates
coord_names = fieldnames(IO.coordi);
n_coord = length(coord_names);

%% Get bounds and initial guess
[bounds,scaling] = getBounds(S,model_info);
Qs_IK = getIK(motion_file,model_info);

for i = 1:100
    Qs = Qs_IK.allfilt(i,2:end)';
    QsQdots1 = zeros(n_coord*2,1);
    QsQdots1(1:2:end) = Qs(:);
    Qsddots = zeros(n_coord,1);
    QsQsdotsQsddots = [QsQdots1; Qsddots];
    res1 = F1(QsQsdotsQsddots);
    evaluated(:,i) = full(res1); % the evaluated external function
end






