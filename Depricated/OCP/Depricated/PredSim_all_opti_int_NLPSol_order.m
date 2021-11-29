%%  Three-dimensional muscle-driven predictive simulations of human gaits
%
% Author: Antoine Falisse
% Date: 1/7/2019
%
clear all;
clc
close all;

%% User inputs
% This script can be run to solve the optimal control problems but also to
% analyze and process the results. The user does not need to re-run the
% optimal control problems to analyze the results. Therefore, the user can
% select some settings beforehand through the variable num_set. For
% example, when num_set(1) is set to 0, the script will not run the
% optimal control problem. Here is a brief description of num_set:
% num_set(1): set to 1 to solve problem
% num_set(2): set to 1 to analyze results
% num_set(3): set to 1 to load results
% num_set(4): set to 1 to save results
% num_set(5): set to 1 to visualize guess-bounds 
% num_set(6): set to 1 to write .mot file
% num_set(7): set to 1 to decompose cost
% Note that you should re-run the simulations to write out the .mot files
% and visualize the results in the OpenSim GUI.

num_set = [1,0,0,0,0,0,0]; % This configuration solves the problem
% num_set = [0,1,1,0,0,0,0]; % This configuration analyzes the results

% The variable settings in the following section will set some parameters 
% of the optimal control problems. Through the variable idx_ww, the user  
% can select which row of parameters will be used.
idx_ww = [22]; % Index row in matrix settings (1:198)

%% Settings
import casadi.*
subject = 'subject1';
parallelMode = 'thread';
NThreads = 8;

solveProblem    = num_set(1); % set to 1 to solve problem
analyseResults  = num_set(2); % set to 1 to analyze results
loadResults     = num_set(3); % set to 1 to load results
saveResults     = num_set(4); % set to 1 to save sens. results
checkBoundsIG   = num_set(5); % set to 1 to visualize guess-bounds 
writeIKmotion   = num_set(6); % set to 1 to write .mot file
decomposeCost   = num_set(7); % set to 1 to decompose cost

% settings describes the parameters used in the optimal control problems.
% settings(1): average speed: max 2 digits.
% settings(2): NLP error tolerance: 1*10^(-settings(2)).
% settings(3): number of mesh intervals
% settings(4): weight metabolic energy rate
% settings(5): weight joint accelerations
% settings(6): weight arm excitations
% settings(7): weight passive torques
% settings(8): weight muscle activations
% settings(9): power metabolic energy rate
% settings(10): initial guess identifier: 1 quasi-random, 2 data-informed
% settings(11): contact model identifier: 1 generic, 2 subject-specific
% settings(12): initial guess mode identifier for data-informed guesses
%   1 data from average walking motion, 2 data from average running motion 
%   3 partial solution from motion specified in settings(13)
%   4 full solution from motion specified in settings(13)
%   if 1-3 then settings(10) should be 2
%   if 4 then settings(10) does not matter
% settings(13): initial guess case identifier: row number in "settings" of
% the motion used as initial guess when settings(12)=3 || settings(12)=4
% settings(14): weakness hip muscles, 0 no weakness, 50 weakened by 50%,
% 75 weakened by 75%, 90 weakened by 90%,
% settings(15): maximal contraction velocity identifier
%   0 generic values, 1 values*2
% settings(16): weakness ankle plantarflexors, 0 no weakness, 50 weakened 
% by 50%, 75 weakened by 75%, 90 weakened by 90%,
% settings(17): metabolic energy model identifier: 0 Bhargava et al. (2004)
% 1 Umberger et al. (2003), 2 Umberger (2010), 3 Uchida et al. (2016), 
% 4 Umberger (2010) treating muscle lengthening heat rate as 
% Umberger et al. (2003), 5 Umberger (2010) treating negative mechanical 
% work as Umberger et al. (2003)
% settings(18): co-contraction identifier: 0 lower bound activation = 0.05,
% 1 lower bound activation = 0.1, 2 lower bound activation = 0.15,
% 3 lower bound activation = 0.2
settings = [    
    % A. Varying prescribed gait speed
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 1
    1.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 2
    1.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 3
    1.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 4
    0.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 5
    0.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 6
    0.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 7
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 8
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 9
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 10
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 11
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 12
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 13
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 14
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 15
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 16
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 17
    2.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 18
    2.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 19
    2.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 20
    2.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 21    
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 22
    1.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 23
    1.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 24
    1.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 25
    0.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 26
    0.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 27
    0.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 28
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 29
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 30
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 31
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 32
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 33
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 34
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 35
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 36
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 37
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 38
    2.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 39
    2.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 40
    2.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 41
    2.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 42
    % Data-informed (running) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 43
    1.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 44
    1.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 45
    1.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 46
    0.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 47
    0.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 48
    0.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 49
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 50
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 51
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 52
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 53
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 54
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 55
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 56
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 57
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 58
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 59
    2.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 60
    2.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 61
    2.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 62
    2.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0;    % 63
    % Data-informed (full solution at closest speed) initial guess
    1.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 1, 0, 0, 0, 0, 0;    % 64
    1.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 64, 0, 0, 0, 0, 0;   % 65
    1.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 65, 0, 0, 0, 0, 0;   % 66
    0.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 66, 0, 0, 0, 0, 0;   % 67
    0.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 67, 0, 0, 0, 0, 0;   % 68
    0.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 68, 0, 0, 0, 0, 0;   % 69
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 1, 0, 0, 0, 0, 0;    % 70
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 70, 0, 0, 0, 0, 0;   % 71
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 71, 0, 0, 0, 0, 0;   % 72
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 72, 0, 0, 0, 0, 0;   % 73
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 73, 0, 0, 0, 0, 0;   % 74
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 74, 0, 0, 0, 0, 0;   % 75
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 75, 0, 0, 0, 0, 0;   % 76
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 76, 0, 0, 0, 0, 0;   % 77
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 77, 0, 0, 0, 0, 0;   % 78
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 78, 0, 0, 0, 0, 0;   % 79
    2.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 79, 0, 0, 0, 0, 0;   % 80
    2.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 80, 0, 0, 0, 0, 0;   % 81
    2.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 81, 0, 0, 0, 0, 0;   % 82
    2.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 82, 0, 0, 0, 0, 0;   % 83 
    % Data-informed (partial solution at closest speed) initial guess
    1.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 1, 0, 0, 0, 0, 0;    % 84
    1.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 84, 0, 0, 0, 0, 0;   % 85
    1.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 85, 0, 0, 0, 0, 0;   % 86
    0.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 86, 0, 0, 0, 0, 0;   % 87
    0.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 87, 0, 0, 0, 0, 0;   % 88
    0.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 88, 0, 0, 0, 0, 0;   % 89
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 1, 0, 0, 0, 0, 0;    % 90
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 90, 0, 0, 0, 0, 0;   % 91
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 91, 0, 0, 0, 0, 0;   % 92
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 92, 0, 0, 0, 0, 0;   % 93
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 93, 0, 0, 0, 0, 0;   % 94
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 94, 0, 0, 0, 0, 0;   % 95
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 95, 0, 0, 0, 0, 0;   % 96
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 96, 0, 0, 0, 0, 0;   % 97
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 97, 0, 0, 0, 0, 0;   % 98
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 98, 0, 0, 0, 0, 0;   % 99
    2.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 99, 0, 0, 0, 0, 0;   % 100
    2.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 100, 0, 0, 0, 0, 0;  % 101
    2.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 101, 0, 0, 0, 0, 0;  % 102
    2.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 102, 0, 0, 0, 0, 0;  % 103    
    % B. Altering cost function
    % No metabolic energy rate term
    % Quasi-random initial guess
    1.33, 4, 50, 0, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;      % 104
    % Data-informed (walking) initial guess
    1.33, 4, 50, 0, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;      % 105
    % No muscle activity term
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;       % 106
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 0, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;       % 107
    % Not squaring metabolic energy rate term
    % Quasi-random initial guess
    1.33, 4, 50, 1000, 50000, 1000000, 1000, 2000, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0;   % 108
    % Data-informed (walking) initial guess
    1.33, 4, 50, 1000, 50000, 1000000, 1000, 2000, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0;   % 109
    % No passive torque term
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 0, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;       % 110
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 0, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;       % 111
    % Very low weight on joint acceleration term
    % Quasi-random initial guess
    1.33, 4, 50, 500, 1, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;        % 112
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 1, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;        % 113    
    % C. Subject-specific contact model
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0;    % 114
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0;    % 115      
    % D. More mesh points
    % Quasi-random initial guess
    1.33, 4, 100, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;   % 116
    % Data-informed (walking) initial guess
    1.33, 4, 100, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;   % 117    
    % E. Weak hip muscles
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 50, 0, 0, 0, 0;   % 118
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 50, 0, 0, 0, 0;   % 119
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 75, 0, 0, 0, 0;   % 120
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 75, 0, 0, 0, 0;   % 121
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 90, 0, 0, 0, 0;   % 122
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 90, 0, 0, 0, 0;   % 123         
    % F. Increased maximum muscle contraction velocities
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 124
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 125
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 126
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 127
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 128
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 129
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 130
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 131
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 132
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 133
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0;    % 134
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 135
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 136
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 137
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 138
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 139
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 140
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 141
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 142
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 143
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 144
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 1, 0, 0, 0;    % 145    
    % Data-informed (running) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 146
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 147
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 148
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 149
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 150
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 151
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 152
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 153
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 154
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 155
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 2, 0, 0, 1, 0, 0, 0;    % 156
    % Data-informed (full solution at closest speed) initial guess
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 124, 0, 1, 0, 0, 0;  % 157
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 157, 0, 1, 0, 0, 0;  % 158
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 158, 0, 1, 0, 0, 0;  % 159
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 159, 0, 1, 0, 0, 0;  % 160
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 160, 0, 1, 0, 0, 0;  % 161
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 161, 0, 1, 0, 0, 0;  % 162
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 162, 0, 1, 0, 0, 0;  % 163
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 163, 0, 1, 0, 0, 0;  % 164
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 164, 0, 1, 0, 0, 0;  % 165
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 4, 165, 0, 1, 0, 0, 0;  % 166
    % Data-informed (partial solution at closest speed) initial guess
    1.43, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 124, 0, 1, 0, 0, 0;  % 167
    1.53, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 167, 0, 1, 0, 0, 0;  % 168
    1.63, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 168, 0, 1, 0, 0, 0;  % 169
    1.73, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 169, 0, 1, 0, 0, 0;  % 170
    1.83, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 170, 0, 1, 0, 0, 0;  % 171
    1.93, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 171, 0, 1, 0, 0, 0;  % 172
    2.03, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 172, 0, 1, 0, 0, 0;  % 173
    2.13, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 173, 0, 1, 0, 0, 0;  % 174
    2.23, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 174, 0, 1, 0, 0, 0;  % 175
    2.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 3, 175, 0, 1, 0, 0, 0;  % 176    
    % G. Weakness ankle plantaflexors
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 50, 0, 0;   % 177   
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 75, 0, 0;   % 178   
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 90, 0, 0;   % 179
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 50, 0, 0;   % 180 
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 75, 0, 0;   % 181
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 90, 0, 0;   % 182    
    % H. Different metabolic energy models
    % Metabolic energy model from Umberger et al. (2003)
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 1, 0;    % 183
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 1, 0;    % 184
    % Metabolic energy model from Umberger (2010)
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 2, 0;    % 185
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 2, 0;    % 186
    % Metabolic energy model from Uchida et al. (2016)
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 3, 0;    % 187
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 3, 0;    % 188  
    % Metabolic energy model from Umberger (2010) treating muscle lengthening heat rate as Umberger et al. (2003)
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 4, 0;    % 189
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 4, 0;    % 190
    % Metabolic energy model from Umberger (2010) treating negative mechanical work as Umberger et al. (2003)
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 5, 0;    % 191
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 5, 0;    % 192    
    % I. Co-contraction
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 1;    % 193
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 1;    % 194
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 2;    % 195
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 2;    % 196    
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 3;    % 197
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 3;    % 198
];

%% Select settings
for www = 1:length(idx_ww)

%% Set parameters based on settings
ww = idx_ww(www);
% Variable parameters
v_tgt       = settings(ww,1);    % average speed
tol_ipopt   = settings(ww,2);    % tolerance (means 1e-(tol_ipopt))
N           = settings(ww,3);    % number of mesh intervals
W.E         = settings(ww,4);    % weight metabolic energy rate
W.Ak        = settings(ww,5);    % weight joint accelerations
W.ArmE      = settings(ww,6);    % weight arm excitations
W.passMom   = settings(ww,7);    % weight passive torques
W.A         = settings(ww,8);    % weight muscle activations
exp_E       = settings(ww,9);    % power metabolic energy
IGsel       = settings(ww,10);   % initial guess identifier
cm          = settings(ww,11);   % contact model identifier
IGm         = settings(ww,12);   % initial guess mode identifier
IGcase      = settings(ww,13);   % initial guess case identifier
h_weak      = settings(ww,14);   % weakness hip actuators
vMax_s      = settings(ww,15);   % maximal contraction velocity identifier
pf_weak     = settings(ww,16);   % weakness ankle plantaflexors
mE          = settings(ww,17);   % metabolic energy model identifier
coCont      = settings(ww,18);   % co-contraction identifier
% Fixed parameter
W.u = 0.001;
% The filename used to save the results depends on the settings 
v_tgt_id = round(v_tgt,2);
savename = ['_c',num2str(ww),'_v',num2str(v_tgt_id*100),...
    '_T',num2str(tol_ipopt),'_N',num2str(N),'_E',num2str(W.E),...
    '_Ak',num2str(W.Ak),'_AE',num2str(W.ArmE),'_P',num2str(W.passMom),...
    '_A',num2str(W.A),'_eE',num2str(exp_E),'_G',num2str(IGsel),...
    '_M',num2str(cm),'_Gm',num2str(IGm),...
    '_W',num2str(h_weak),'_vM',num2str(vMax_s),...
    '_pW',num2str(pf_weak),'_mE',num2str(mE),'_cc',num2str(coCont)];

% In some cases, the inital guess depends on results from other simulations
if IGm == 3 || IGm == 4    
ww_ig          = IGcase;               % Case identifier used as IG       
v_tgt_ig       = settings(ww_ig,1);    % average speed
tol_ipopt_ig   = settings(ww_ig,2);    % tolerance (means 1e-(tol_ipopt))
N_ig           = settings(ww_ig,3);    % number of mesh intervals
W.E_ig         = settings(ww_ig,4);    % weight metabolic energy
W.Ak_ig        = settings(ww_ig,5);    % weight joint accelerations
W.ArmE_ig      = settings(ww_ig,6);    % weight arm excitations
W.passMom_ig   = settings(ww_ig,7);    % weight passive torques
W.A_ig         = settings(ww_ig,8);    % weight muscle activations
exp_E_ig       = settings(ww_ig,9);    % power metabolic energy
IGsel_ig       = settings(ww_ig,10);   % initial guess identifier
cm_ig          = settings(ww_ig,11);   % contact model identifier
IGm_ig         = settings(ww_ig,12);   % initial guess mode identifier
% There is no need for index 13, since it is not in savename
h_weak_ig      = settings(ww_ig,14);   % weakness hip actuators
vMax_s_ig      = settings(ww_ig,15);   % maximal contraction velocity id
pf_weak_ig     = settings(ww_ig,16);   % weakness ankle plantarflexors
ME_ig          = settings(ww_ig,17);   % metabolic energy model identifier
coCont_ig      = settings(ww_ig,18);   % co-contraction identifier

% The filename used to load the results depends on the settings 
v_tgt_id_ig = round(v_tgt_ig,2);
savename_ig = ['_c',num2str(ww_ig),'_v',num2str(v_tgt_id_ig*100),...
    '_T',num2str(tol_ipopt_ig),'_N',num2str(N_ig),'_E',num2str(W.E_ig),...
    '_Ak',num2str(W.Ak_ig),'_AE',num2str(W.ArmE_ig),...
    '_P',num2str(W.passMom_ig),'_A',num2str(W.A_ig),...
    '_eE',num2str(exp_E_ig),'_G',num2str(IGsel_ig),'_M',num2str(cm_ig),...
    '_Gm',num2str(IGm_ig),'_W',num2str(h_weak_ig),'_vM',num2str(vMax_s_ig),...
    '_pW',num2str(pf_weak_ig),'_mE',num2str(ME_ig),'_cc',num2str(coCont_ig)];
end

%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathmain = pwd;
% We use different external functions, since we also want to access some 
% parameters of the model in a post-processing phase.
[pathRepo,~,~] = fileparts(pathmain);
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
% Loading external functions. 
cd(pathExternalFunctions);
setup.derivatives =  'AD'; % Algorithmic differentiation
if ispc    
    switch setup.derivatives
        case {'AD'}   
            if cm == 1
                F = external('F','PredSim.dll');   
                if analyseResults
                    F1 = external('F','PredSim_pp.dll');
                end
            elseif cm == 2
                F = external('F','PredSim_SSCM.dll');   
                if analyseResults
                    F1 = external('F','PredSim_SSCM_pp.dll');
                end
            end
    end
elseif ismac
    switch setup.derivatives
        case {'AD'}   
            if cm == 1
                F = external('F','PredSim.dylib');   
                if analyseResults
                    F1 = external('F','PredSim_pp.dylib');
                end
            elseif cm == 2
                F = external('F','PredSim_SSCM.dylib');   
                if analyseResults
                    F1 = external('F','PredSim_SSCM_pp.dylib');
                end
            end
    end
else
    disp('Platform not supported')
end
cd(pathmain);
% This is an example of how to call an external function with some
% numerical values.
% vec1 = -ones(87,1);
% res1 = full(F(vec1));
% res2 = full(F1(vec1));

%% Indices external function
% Indices of the elements in the external functions
% External function: F
% First, joint torques. 
jointi.pelvis.tilt  = 1; 
jointi.pelvis.list  = 2; 
jointi.pelvis.rot   = 3; 
jointi.pelvis.tx    = 4;
jointi.pelvis.ty    = 5;
jointi.pelvis.tz    = 6;
jointi.hip_flex.l   = 7;
jointi.hip_add.l    = 8;
jointi.hip_rot.l    = 9;
jointi.hip_flex.r   = 10;
jointi.hip_add.r    = 11;
jointi.hip_rot.r    = 12;
jointi.knee.l       = 13;
jointi.knee.r       = 14;
jointi.ankle.l      = 15;
jointi.ankle.r      = 16;
jointi.subt.l       = 17;
jointi.subt.r       = 18;
jointi.trunk.ext    = 19;
jointi.trunk.ben    = 20;
jointi.trunk.rot    = 21;
jointi.sh_flex.l    = 22;
jointi.sh_add.l     = 23;
jointi.sh_rot.l     = 24;
jointi.sh_flex.r    = 25;
jointi.sh_add.r     = 26;
jointi.sh_rot.r     = 27;
jointi.elb.l        = 28;
jointi.elb.r        = 29;
% Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all 
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
residuals_noarmsi   = jointi.pelvis.tilt:jointi.trunk.rot; % all but arms
roti                = [jointi.pelvis.tilt:jointi.pelvis.rot,...
    jointi.hip_flex.l:jointi.elb.r];
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all 
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.leg      = 9; % #joints needed for polynomials
% Second, origins bodies. 
% Calcaneus
calcOr.r    = 30:31;
calcOr.l    = 32:33;
calcOr.all  = [calcOr.r,calcOr.l];
NcalcOr     = length(calcOr.all);
% Femurs
femurOr.r   = 34:35;
femurOr.l   = 36:37;
femurOr.all = [femurOr.r,femurOr.l];
NfemurOr    = length(femurOr.all);
% Hands
handOr.r    = 38:39;
handOr.l    = 40:41;
handOr.all  = [handOr.r,handOr.l];
NhandOr     = length(handOr.all);
% Tibias
tibiaOr.r   = 42:43;
tibiaOr.l   = 44:45;
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
NtibiaOr    = length(tibiaOr.all);
% External function: F1 (post-processing purpose only)
% Ground reaction forces (GRFs)
GRFi.r      = 30:32;
GRFi.l      = 33:35;
GRFi.all    = [GRFi.r,GRFi.l];
NGRF        = length(GRFi.all);
% Origins calcaneus (3D)
calcOrall.r     = 36:38;
calcOrall.l     = 39:41;
calcOrall.all   = [calcOrall.r,calcOrall.l];
NcalcOrall      = length(calcOrall.all);

%% Model info
body_mass = 62;
body_weight = body_mass*9.81;

%% Collocation scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points. 
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Muscle-tendon parameters 
% Muscles from one leg and from the back
muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
        'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
        'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
        'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
        'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
        'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
        'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
        'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
        'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
        'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
        'intobl_l','extobl_l'};
% Muscle indices for later use
pathmusclemodel = [pathRepo,'/MuscleModel'];
addpath(genpath(pathmusclemodel)); 
% (1:end-3), since we do not want to count twice the back muscles
musi = MuscleIndices(muscleNames(1:end-3));
% Total number of muscles
NMuscle = length(muscleNames(1:end-3))*2;
% Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation 
% angles; Row 5: maximal contraction velocities
load([pathmusclemodel,'/MTparameters_',subject,'.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];
% Indices of the muscles actuating the different joints for later use
pathpolynomial = [pathRepo,'/Polynomials'];
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),...
    tl.muscle_spanning_joint_INFO(1:end-3,:));

% By default, the tendon stiffness is 35 and the shift is 0.
aTendon = 35*ones(NMuscle,1);
shift = zeros(NMuscle,1);

% Adjust the maximal isometric force of the hip actuators if needed.
if h_weak ~= 0
    MTparameters_m(1,[mai(1).mus.l,mai(1).mus.r]) = ...
        MTparameters_m(1,[mai(1).mus.l,mai(1).mus.r])-...
        h_weak/100*MTparameters_m(1,[mai(1).mus.l,mai(1).mus.r]);
end

% Adjust the maximal isometric force of the ankle plantarflexors if needed.
if pf_weak ~= 0
    idx_FD = find(strcmp(muscleNames,'flex_dig_r'));
    idx_FH = find(strcmp(muscleNames,'flex_hal_r'));
    idx_GL = find(strcmp(muscleNames,'lat_gas_r'));
    idx_GM = find(strcmp(muscleNames,'med_gas_r'));
    idx_PB = find(strcmp(muscleNames,'per_brev_r'));
    idx_PL = find(strcmp(muscleNames,'per_long_r'));
    idx_SO = find(strcmp(muscleNames,'soleus_r'));
    idx_TP = find(strcmp(muscleNames,'tib_post_r'));    
    idx_pf = [idx_FD,idx_FH,idx_GL,idx_GM,idx_PB,idx_PL,idx_SO,idx_TP];
    idx_pf_all = [idx_pf,idx_pf+NMuscle/2];        
    MTparameters_m(1,idx_pf_all) = MTparameters_m(1,idx_pf_all)-...
        pf_weak/100*MTparameters_m(1,idx_pf_all);
end

% Adjust the maximum contraction velocities if needed.
if vMax_s == 1
    % Maximum contraction velocities * 2
    MTparameters_m(end,:) = MTparameters_m(end,:).*2;
elseif vMax_s == 2
    % Maximum contraction velocities * 1.5
    MTparameters_m(end,:) = MTparameters_m(end,:).*1.5;
end

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
% (1:end-3), since we do not want to count twice the back muscles
tension = getSpecificTensions(muscleNames(1:end-3)); 
tensions = [tension;tension];
% (1:end-3), since we do not want to count twice the back muscles
pctst = getSlowTwitchRatios(muscleNames(1:end-3)); 
pctsts = [pctst;pctst];

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));
% We load some variables for the polynomial approximations
load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
load([pathpolynomial,'/MuscleInfo_',subject,'.mat']);
% For the polynomials, we want all independent muscles. So we do not need 
% the muscles from both legs, since we assume bilateral symmetry, but want
% all muscles from the back (indices 47:49).
musi_pol = [musi,47,48,49];
NMuscle_pol = NMuscle/2+3;
CasADiFunctions_all

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk
pathPassiveMoments = [pathRepo,'/PassiveMoments'];
addpath(genpath(pathPassiveMoments));
PassiveMomentsData

%% Experimental data
% We extract experimental data to set bounds and initial guesses if needed
pathData = [pathRepo,'/OpenSimModel/',subject];
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','lumbar_extension',...
    'lumbar_bending','lumbar_rotation','arm_flex_l','arm_add_l',...
    'arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_l',...
    'elbow_flex_r'};
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));
% Extract joint positions from average walking motion
motion_walk         = 'walking';
nametrial_walk.id   = ['average_',motion_walk,'_HGC']; 
nametrial_walk.IK   = ['IK_',nametrial_walk.id];
pathIK_walk         = [pathData,'/IK/',nametrial_walk.IK,'.mat'];
Qs_walk             = getIK(pathIK_walk,joints);
% Depending on the initial guess mode, we extract further experimental data
if IGm == 2
    % Extract joint positions from average running motion
    motion_run          = 'running';
    nametrial_run.id    = ['average_',motion_run,'_HGC'];  
    nametrial_run.IK    = ['IK_',nametrial_run.id];
    pathIK_run          = [pathData,'/IK/',nametrial_run.IK,'.mat'];
    Qs_run              = getIK(pathIK_run,joints);    
elseif IGm == 3
    % Extract joint positions from existing motion (previous results)
    p = mfilename('fullpath');
    [~,namescript,~] = fileparts(p);
    pathIK = [pathRepo,'/Results/',namescript,'/IK',savename_ig,'.mot'];
    Qs_ig = getIK(pathIK,joints);
    % When saving the results, we save a full gait cycle (2*N) so here we
    % only select 1:N to have half a gait cycle
    Qs_ig_sel.allfilt = Qs_ig.allfilt(1:N,:);
    Qs_ig_sel.time = Qs_ig.time(1:N,:);
    Qs_ig_sel.colheaders = Qs_ig.colheaders;
end    

%% Bounds
pathBounds = [pathRepo,'/Bounds'];
addpath(genpath(pathBounds));
[bounds,scaling] = getBounds_all(Qs_walk,NMuscle,nq,jointi,v_tgt);
% Simulate co-contraction by increasing the lower bound on muscle activations
if coCont == 1
    bounds.a.lower = 0.1*ones(1,NMuscle);
elseif coCont == 2
    bounds.a.lower = 0.15*ones(1,NMuscle);
elseif coCont == 3
    bounds.a.lower = 0.2*ones(1,NMuscle);
end

%% Initial guess
% The initial guess depends on the settings
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));
if IGsel == 1 % Quasi-random initial guess  
    guess = getGuess_QR_opti_int(N,nq,NMuscle,scaling,v_tgt,jointi,d);
elseif IGsel == 2 % Data-informed initial guess
    if IGm == 1 % Data from average walking motion
        time_IC = [Qs_walk.time(1),Qs_walk.time(end)];
        guess = getGuess_DI_opti_int(Qs_walk,nq,N,time_IC,NMuscle,jointi,...
            scaling,v_tgt,d);   
    elseif IGm == 2 % Data from average runing motion    
        time_IC = [Qs_run.time(1),Qs_run.time(end)];
        guess = getGuess_DI_opti_int(Qs_run,nq,N,time_IC,NMuscle,jointi,...
            scaling,v_tgt,d);
    elseif IGm == 3 % Data from selected motion
        time_IC = [Qs_ig_sel.time(1),Qs_ig_sel.time(end)];
        guess = ...
            getGuess_DI_t(Qs_ig_sel,nq,N,time_IC,NMuscle,jointi,scaling);    
    end 
end
% If co-contraction, the initial guess of muscles activations is increased
if coCont == 1
    guess.a = 0.15*ones(N,NMuscle);
elseif coCont == 2
    guess.a = 0.20*ones(N,NMuscle);
elseif coCont == 3
    guess.a = 0.25*ones(N,NMuscle);
end
% This allows visualizing the initial guess and the bounds
if checkBoundsIG
    pathPlots = [pathRepo,'/Plots'];
    addpath(genpath(pathPlots));
    plot_BoundsVSInitialGuess_all
end

%% Formulate the NLP
if solveProblem
    % Create opti instance
    opti = casadi.Opti();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static parameters
    % Final time
    tf = opti.variable();
    opti.subject_to(bounds.tf.lower < tf < bounds.tf.upper);
    opti.set_initial(tf, guess.tf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define states 
    % Muscle activations at mesh points
    a = opti.variable(NMuscle,N+1);
    opti.subject_to(bounds.a.lower'*ones(1,N+1) < a < ...
        bounds.a.upper'*ones(1,N+1));
    opti.set_initial(a, guess.a');
    % Muscle activations at collocation points
    a_col = opti.variable(NMuscle,d*N);
    opti.subject_to(bounds.a.lower'*ones(1,d*N) < a_col < ...
        bounds.a.upper'*ones(1,d*N));
    opti.set_initial(a_col, guess.a_col');    
    % Muscle-tendon forces at mesh points
    FTtilde = opti.variable(NMuscle,N+1);
    opti.subject_to(bounds.FTtilde.lower'*ones(1,N+1) < FTtilde < ...
        bounds.FTtilde.upper'*ones(1,N+1));
    opti.set_initial(FTtilde, guess.FTtilde');
    % Muscle-tendon forces at collocation points
    FTtilde_col = opti.variable(NMuscle,d*N);
    opti.subject_to(bounds.FTtilde.lower'*ones(1,d*N) < FTtilde_col < ...
        bounds.FTtilde.upper'*ones(1,d*N));
    opti.set_initial(FTtilde_col, guess.FTtilde_col');    
    % Qs at mesh points
    Qs = opti.variable(nq.all,N+1);
    % We want to constraint the pelvis_tx position at the first mesh point,
    % and avoid redundant bounds
    lboundsQsk = bounds.QsQdots.lower(1:2:end)'*ones(1,N+1);
    lboundsQsk(jointi.pelvis.tx,1) = ...
        bounds.QsQdots_0.lower(2*jointi.pelvis.tx-1);
    uboundsQsk = bounds.QsQdots.upper(1:2:end)'*ones(1,N+1);
    uboundsQsk(jointi.pelvis.tx,1) = ...
        bounds.QsQdots_0.upper(2*jointi.pelvis.tx-1);
    opti.subject_to(lboundsQsk < Qs < uboundsQsk);
    opti.set_initial(Qs, guess.QsQdots(:,1:2:end)');
    % Qs at collocation points
    Qs_col = opti.variable(nq.all,d*N);
    opti.subject_to(bounds.QsQdots.lower(1:2:end)'*ones(1,d*N) < Qs_col < ...
        bounds.QsQdots.upper(1:2:end)'*ones(1,d*N));
    opti.set_initial(Qs_col, guess.QsQdots_col(:,1:2:end)');   
    % Qdots at mesh points
    Qdots = opti.variable(nq.all,N+1);
    opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,N+1) < Qdots < ...
        bounds.QsQdots.upper(2:2:end)'*ones(1,N+1));
    opti.set_initial(Qdots, guess.QsQdots(:,2:2:end)');    
    % Qdots at collocation points
    Qdots_col = opti.variable(nq.all,d*N);
    opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,d*N) < Qdots_col < ...
        bounds.QsQdots.upper(2:2:end)'*ones(1,d*N));
    opti.set_initial(Qdots_col, guess.QsQdots_col(:,2:2:end)');
    % Arm activations at mesh points
    a_a = opti.variable(nq.arms,N+1);
    opti.subject_to(bounds.a_a.lower'*ones(1,N+1) < a_a < ...
        bounds.a_a.upper'*ones(1,N+1));
    opti.set_initial(a_a, guess.a_a');  
    % Arm activations at collocation points
    a_a_col = opti.variable(nq.arms,d*N);
    opti.subject_to(bounds.a_a.lower'*ones(1,d*N) < a_a_col < ...
        bounds.a_a.upper'*ones(1,d*N));
    opti.set_initial(a_a_col, guess.a_a_col');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define controls
    % Time derivative of muscle activations (states) at mesh points
    vA = opti.variable(NMuscle, N);
    opti.subject_to(bounds.vA.lower'*ones(1,N) < vA < ...
        bounds.vA.upper'*ones(1,N));
    opti.set_initial(vA, guess.vA');   
    % Arm excitations
    e_a = opti.variable(nq.arms, N);
    opti.subject_to(bounds.e_a.lower'*ones(1,N) < e_a < ...
        bounds.e_a.upper'*ones(1,N));
    opti.set_initial(e_a, guess.e_a');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define "slack" controls
    % Time derivative of muscle-tendon forces (states) at collocation points
    dFTtilde_col = opti.variable(NMuscle, d*N);
    opti.subject_to(bounds.dFTtilde.lower'*ones(1,d*N) < dFTtilde_col < ...
        bounds.dFTtilde.upper'*ones(1,d*N));
    opti.set_initial(dFTtilde_col, guess.dFTtilde_col');
    % Time derivative of Qdots (states) at collocation points
    A_col = opti.variable(nq.all, d*N);
    opti.subject_to(bounds.Qdotdots.lower'*ones(1,d*N) < A_col < ...
        bounds.Qdotdots.upper'*ones(1,d*N));
    opti.set_initial(A_col, guess.Qdotdots_col');          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unscale variables
    Qs_nsc = Qs.*(scaling.QsQdots(1:2:end)'*ones(1,N+1));
    Qs_col_nsc = Qs_col.*(scaling.QsQdots(1:2:end)'*ones(1,d*N));    
    Qdots_nsc = Qdots.*(scaling.QsQdots(2:2:end)'*ones(1,N+1));
    Qdots_col_nsc = Qdots_col.*(scaling.QsQdots(2:2:end)'*ones(1,d*N));    
    FTtilde_nsc = FTtilde.*(scaling.FTtilde'*ones(1,N+1));
    FTtilde_col_nsc = FTtilde_col.*(scaling.FTtilde'*ones(1,d*N));
    dFTtilde_col_nsc = dFTtilde_col.*scaling.dFTtilde;
    vA_nsc = vA.*scaling.vA;
    A_col_nsc = A_col.*(scaling.Qdotdots'*ones(1,d*N));    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Provide expression for the distance traveled
    dist_trav_tot = Qs_nsc((jointi.pelvis.tx),end) - ...
        Qs_nsc((jointi.pelvis.tx),1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Time step
    h = tf/N; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J = 0; % Initialize cost function
    % Loop over mesh points
    for k=1:N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Variables within current mesh interval
        % States      
        akj = [a(:,k), a_col(:,(k-1)*d+1:k*d)]; 
        FTtildekj_nsc = [FTtilde_nsc(:,k), FTtilde_col_nsc(:,(k-1)*d+1:k*d)];
        FTtildekj = [FTtilde(:,k), FTtilde_col(:,(k-1)*d+1:k*d)];        
        Qskj_nsc = [Qs_nsc(:,k), Qs_col_nsc(:,(k-1)*d+1:k*d)];
        Qskj = [Qs(:,k), Qs_col(:,(k-1)*d+1:k*d)];       
        Qdotskj_nsc = [Qdots_nsc(:,k), Qdots_col_nsc(:,(k-1)*d+1:k*d)];
        Qdotskj = [Qdots(:,k), Qdots_col(:,(k-1)*d+1:k*d)];
        a_akj = [a_a(:,k), a_a_col(:,(k-1)*d+1:k*d)];        
        QsQdotskj_nsc = MX(nq.all*2, d+1);
        QsQdotskj_nsc(1:2:end,:) = Qskj_nsc;
        QsQdotskj_nsc(2:2:end,:) = Qdotskj_nsc;        
        % Controls
        vAk = vA(:,k); vAk_nsc = vA_nsc(:,k);
        e_ak = e_a(:,k); 
        % "Slack" controls
        dFTtildej = dFTtilde_col(:,(k-1)*d+1:k*d); 
        dFTtildej_nsc = dFTtilde_col_nsc(:,(k-1)*d+1:k*d);
        Aj = A_col(:,(k-1)*d+1:k*d); Aj_nsc = A_col_nsc(:,(k-1)*d+1:k*d);                 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop over collocation points
        for j=1:d            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
            % Get muscle-tendon lengths, velocities, and moment arms
            % Left leg
            qinj_l = [Qskj_nsc(jointi.hip_flex.l,j+1),...
                Qskj_nsc(jointi.hip_add.l,j+1), ...
                Qskj_nsc(jointi.hip_rot.l,j+1), ...
                Qskj_nsc(jointi.knee.l,j+1), ...
                Qskj_nsc(jointi.ankle.l,j+1),...
                Qskj_nsc(jointi.subt.l,j+1),...
                Qskj_nsc(jointi.trunk.ext,j+1),...
                Qskj_nsc(jointi.trunk.ben,j+1),...
                Qskj_nsc(jointi.trunk.rot,j+1)];  
            qdotinj_l = [Qdotskj_nsc(jointi.hip_flex.l,j+1),...
                Qdotskj_nsc(jointi.hip_add.l,j+1),...
                Qdotskj_nsc(jointi.hip_rot.l,j+1),...
                Qdotskj_nsc(jointi.knee.l,j+1),...
                Qdotskj_nsc(jointi.ankle.l,j+1),...
                Qdotskj_nsc(jointi.subt.l,j+1),...
                Qdotskj_nsc(jointi.trunk.ext,j+1),...
                Qdotskj_nsc(jointi.trunk.ben,j+1),...
                Qdotskj_nsc(jointi.trunk.rot,j+1)];  
            [lMTj_l,vMTj_l,MAj_l] =  f_lMT_vMT_dM(qinj_l,qdotinj_l);    
            MAj.hip_flex.l   =  MAj_l(mai(1).mus.l',1);
            MAj.hip_add.l    =  MAj_l(mai(2).mus.l',2);
            MAj.hip_rot.l    =  MAj_l(mai(3).mus.l',3);
            MAj.knee.l       =  MAj_l(mai(4).mus.l',4);
            MAj.ankle.l      =  MAj_l(mai(5).mus.l',5);  
            MAj.subt.l       =  MAj_l(mai(6).mus.l',6); 
            % For the back muscles, we want left and right together: left
            % first, right second. In MuscleInfo, we first have the right
            % muscles (44:46) and then the left muscles (47:49). Since the back
            % muscles only depend on back dofs, we do not care if we extract
            % them "from the left or right leg" so here we just picked left.
            MAj.trunk_ext    =  MAj_l([47:49,mai(7).mus.l]',7);
            MAj.trunk_ben    =  MAj_l([47:49,mai(8).mus.l]',8);
            MAj.trunk_rot    =  MAj_l([47:49,mai(9).mus.l]',9);
            % Right leg
            qinj_r = [Qskj_nsc(jointi.hip_flex.r,j+1),...
                Qskj_nsc(jointi.hip_add.r,j+1),...
                Qskj_nsc(jointi.hip_rot.r,j+1),...
                Qskj_nsc(jointi.knee.r,j+1),...
                Qskj_nsc(jointi.ankle.r,j+1),...
                Qskj_nsc(jointi.subt.r,j+1),...
                Qskj_nsc(jointi.trunk.ext,j+1),...
                Qskj_nsc(jointi.trunk.ben,j+1),...
                Qskj_nsc(jointi.trunk.rot,j+1)];  
            qdotinj_r = [Qdotskj_nsc(jointi.hip_flex.r,j+1),...
                Qdotskj_nsc(jointi.hip_add.r,j+1),...
                Qdotskj_nsc(jointi.hip_rot.r,j+1),...
                Qdotskj_nsc(jointi.knee.r,j+1),...
                Qdotskj_nsc(jointi.ankle.r,j+1),...
                Qdotskj_nsc(jointi.subt.r,j+1),...
                Qdotskj_nsc(jointi.trunk.ext,j+1),...
                Qdotskj_nsc(jointi.trunk.ben,j+1),...
                Qdotskj_nsc(jointi.trunk.rot,j+1)];      
            [lMTj_r,vMTj_r,MAj_r] = f_lMT_vMT_dM(qinj_r,qdotinj_r);
            % Here we take the indices from left since the vector is 1:49
            MAj.hip_flex.r   =  MAj_r(mai(1).mus.l',1);
            MAj.hip_add.r    =  MAj_r(mai(2).mus.l',2);
            MAj.hip_rot.r    =  MAj_r(mai(3).mus.l',3);
            MAj.knee.r       =  MAj_r(mai(4).mus.l',4);
            MAj.ankle.r      =  MAj_r(mai(5).mus.l',5);
            MAj.subt.r       =  MAj_r(mai(6).mus.l',6);
            % Both legs
            % In MuscleInfo, we first have the right back muscles (44:46) and 
            % then the left back muscles (47:49). Here we re-organize so that
            % we have first the left muscles and then the right muscles.
            lMTj_lr = [lMTj_l([1:43,47:49],1);lMTj_r(1:46,1)];
            vMTj_lr = [vMTj_l([1:43,47:49],1);vMTj_r(1:46,1)];   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get muscle-tendon forces and derive Hill-equilibrium
            [Hilldiffj,FTj,Fcej,Fpassj,Fisoj,vMmaxj,massMj] = ...
                f_forceEquilibrium_FtildeState_all_tendon(akj(:,j+1),...
                    FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
                    lMTj_lr,vMTj_lr,tensions,aTendon,shift); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get metabolic energy rate if in the cost function   
            if W.E ~= 0    
                % Get muscle fiber lengths
                [~,lMtildej] = f_FiberLength_TendonForce_tendon(...
                    FTtildekj_nsc(:,j+1),lMTj_lr,aTendon,shift); 
                % Get muscle fiber velocities
                [vMj,~] = f_FiberVelocity_TendonForce_tendon(...
                    FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
                    lMTj_lr,vMTj_lr,aTendon,shift);
                % Get metabolic energy rate
                if mE == 0 % Bhargava et al. (2004)
                    [e_totj,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                        akj(:,j+1),akj(:,j+1),lMtildej,...
                        vMj,Fcej,Fpassj,massMj,pctsts,Fisoj,...
                        MTparameters_m(1,:)',body_mass,10);
                elseif mE == 1 % Umberger et al. (2003)
                    % vMtilde defined for this model as vM/lMopt
                    vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                    [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2003all(...
                        akj(:,j+1),akj(:,j+1),lMtildej,...
                        vMtildeUmbj,vMj,Fcej,massMj,...
                        pctsts,vMmaxj,Fisoj,body_mass,10);
                elseif mE == 2 % Umberger (2010)
                    % vMtilde defined for this model as vM/lMopt
                    vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                    [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2010all(...
                        akj(:,j+1),akj(:,j+1),lMtildej,...
                        vMtildeUmbj,vMj,Fcej,massMj,...
                        pctsts,vMmaxj,Fisoj,body_mass,10);  
                elseif mE == 3 % Uchida et al. (2016)
                    % vMtilde defined for this model as vM/lMopt
                    vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                    [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2016all(...
                        akj(:,j+1),akj(:,j+1),lMtildej,...
                        vMtildeUmbj,vMj,Fcej,massMj,...
                        pctsts,vMmaxj,Fisoj,body_mass,10);  
                elseif mE == 4 % Umberger (2010) treating muscle lengthening 
                    % heat rate as Umberger et al. (2003)
                    % vMtilde defined for this model as vM/lMopt
                    vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                    [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2010all_hl(...
                        akj(:,j+1),akj(:,j+1),lMtildej,...
                        vMtildeUmbj,vMj,Fcej,massMj,...
                        pctsts,vMmaxj,Fisoj,body_mass,10); 
                elseif mE == 5 % Umberger (2010) treating negative mechanical 
                    % work as Umberger et al. (2003)
                    % vMtilde defined for this model as vM/lMopt
                    vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                    [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2010all_neg(...
                        akj(:,j+1),akj(:,j+1),lMtildej,...
                        vMtildeUmbj,vMj,Fcej,massMj,...
                        pctsts,vMmaxj,Fisoj,body_mass,10); 
                end                
            end        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get passive joint torques
            Tau_passj.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
                theta.pass.hip.flex,Qskj_nsc(jointi.hip_flex.l,j+1),...
                Qdotskj_nsc(jointi.hip_flex.l,j+1));
            Tau_passj.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
                theta.pass.hip.flex,Qskj_nsc(jointi.hip_flex.r,j+1),...
                Qdotskj_nsc(jointi.hip_flex.r,j+1));
            Tau_passj.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
                theta.pass.hip.add,Qskj_nsc(jointi.hip_add.l,j+1),...
                Qdotskj_nsc(jointi.hip_add.l,j+1));
            Tau_passj.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
                theta.pass.hip.add,Qskj_nsc(jointi.hip_add.r,j+1),...
                Qdotskj_nsc(jointi.hip_add.r,j+1));
            Tau_passj.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
                theta.pass.hip.rot,Qskj_nsc(jointi.hip_rot.l,j+1),...
                Qdotskj_nsc(jointi.hip_rot.l,j+1));
            Tau_passj.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
                theta.pass.hip.rot,Qskj_nsc(jointi.hip_rot.r,j+1),...
                Qdotskj_nsc(jointi.hip_rot.r,j+1));
            Tau_passj.knee.l        = f_PassiveMoments(k_pass.knee,...
                theta.pass.knee,Qskj_nsc(jointi.knee.l,j+1),...
                Qdotskj_nsc(jointi.knee.l,j+1));
            Tau_passj.knee.r        = f_PassiveMoments(k_pass.knee,...
                theta.pass.knee,Qskj_nsc(jointi.knee.r,j+1),...
                Qdotskj_nsc(jointi.knee.r,j+1));
            Tau_passj.ankle.l       = f_PassiveMoments(k_pass.ankle,...
                theta.pass.ankle,Qskj_nsc(jointi.ankle.l,j+1),...
                Qdotskj_nsc(jointi.ankle.l,j+1));
            Tau_passj.ankle.r       = f_PassiveMoments(k_pass.ankle,...
                theta.pass.ankle,Qskj_nsc(jointi.ankle.r,j+1),...
                Qdotskj_nsc(jointi.ankle.r,j+1));        
            Tau_passj.subt.l       = f_PassiveMoments(k_pass.subt,...
                theta.pass.subt,Qskj_nsc(jointi.subt.l,j+1),...
                Qdotskj_nsc(jointi.subt.l,j+1));
            Tau_passj.subt.r       = f_PassiveMoments(k_pass.subt,...
                theta.pass.subt,Qskj_nsc(jointi.subt.r,j+1),...
                Qdotskj_nsc(jointi.subt.r,j+1));        
            Tau_passj.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
                theta.pass.trunk.ext,Qskj_nsc(jointi.trunk.ext,j+1),...
                Qdotskj_nsc(jointi.trunk.ext,j+1));
            Tau_passj.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
                theta.pass.trunk.ben,Qskj_nsc(jointi.trunk.ben,j+1),...
                Qdotskj_nsc(jointi.trunk.ben,j+1));
            Tau_passj.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
                theta.pass.trunk.rot,Qskj_nsc(jointi.trunk.rot,j+1),...
                Qdotskj_nsc(jointi.trunk.rot,j+1));        
            Tau_passj_all = [Tau_passj.hip.flex.l,Tau_passj.hip.flex.r,...
                Tau_passj.hip.add.l,Tau_passj.hip.add.r,...
                Tau_passj.hip.rot.l,Tau_passj.hip.rot.r,...
                Tau_passj.knee.l,Tau_passj.knee.r,Tau_passj.ankle.l,...
                Tau_passj.ankle.r,Tau_passj.subt.l,Tau_passj.subt.r,...
                Tau_passj.trunk.ext,Tau_passj.trunk.ben,...
                Tau_passj.trunk.rot]';
            % Expression for the state derivatives at the collocation points
            Qsp_nsc          = Qskj_nsc*C(:,j+1);
            Qdotsp_nsc       = Qdotskj_nsc*C(:,j+1);            
            FTtildep_nsc    = FTtildekj_nsc*C(:,j+1);
            ap              = akj*C(:,j+1);
            a_ap            = a_akj*C(:,j+1);
            % Append collocation equations
            % Dynamic constraints are scaled using the same scale
            % factors as the ones used to scale the states
            % Activation dynamics (implicit formulation)
            opti.subject_to((h*vAk_nsc - ap)./scaling.a == 0);
            % Contraction dynamics (implicit formulation)     
            opti.subject_to((h*dFTtildej_nsc(:,j) - FTtildep_nsc)./...
                scaling.FTtilde' == 0);
            % Skeleton dynamics (implicit formulation)               
            qdotj_nsc = Qdotskj_nsc(:,j+1); % velocity
            opti.subject_to((h*qdotj_nsc - Qsp_nsc)./...
                scaling.QsQdots(1:2:end)' == 0)
            opti.subject_to((h*Aj_nsc(:,j) - Qdotsp_nsc)./...
                scaling.QsQdots(2:2:end)' == 0)  
            % Arm activation dynamics (explicit formulation)  
            dadtj = f_ArmActivationDynamics(e_ak,a_akj(:,j+1)');
            opti.subject_to((h*dadtj - a_ap)./scaling.a_a == 0)   
            % Add contribution to the quadrature function
            if W.E == 0
            J = J + 1/(dist_trav_tot)*(...
                W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
                W.Ak*B(j+1)     *(f_J21(Aj(residuals_noarmsi,j)))*h + ... 
                W.passMom*B(j+1)*(f_J15(Tau_passj_all))*h + ...
                W.u*B(j+1)      *(f_J92(vAk))*h + ...
                W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...                
                W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);
            elseif W.A == 0
            J = J + 1/(dist_trav_tot)*(...
                W.E*B(j+1)      *(f_J92exp(e_totj,exp_E))/body_mass*h + ...                
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
                W.Ak*B(j+1)     *(f_J21(Aj(residuals_noarmsi,j)))*h + ... 
                W.passMom*B(j+1)*(f_J15(Tau_passj_all))*h + ...
                W.u*B(j+1)      *(f_J92(vAk))*h + ...
                W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...                
                W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h); 
            elseif W.passMom == 0
                J = J + 1/(dist_trav_tot)*(...
                W.E*B(j+1)      *(f_J92exp(e_totj,exp_E))/body_mass*h + ...                
                W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
                W.Ak*B(j+1)     *(f_J21(Aj(residuals_noarmsi,j)))*h + ... 
                W.u*B(j+1)      *(f_J92(vAk))*h + ...
                W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...                
                W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);
            else                
            J = J + 1/(dist_trav_tot)*(...
                W.E*B(j+1)      *(f_J92exp(e_totj,exp_E))/body_mass*h + ...
                W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
                W.Ak*B(j+1)     *(f_J21(Aj(residuals_noarmsi,j)))*h + ... 
                W.passMom*B(j+1)*(f_J15(Tau_passj_all))*h + ...
                W.u*B(j+1)      *(f_J92(vAk))*h + ...
                W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...                
                W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Call external function (run inverse dynamics)
            [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add path constraints
            % Null pelvis residuals
            opti.subject_to(Tj(ground_pelvisi,1) == 0);
            % Muscle-driven joint torques for the lower limbs and the trunk
            % Hip flexion, left
            Ft_hip_flex_l   = FTj(mai(1).mus.l',1);
            T_hip_flex_l    = f_T27(MAj.hip_flex.l,Ft_hip_flex_l);
            opti.subject_to(Tj(jointi.hip_flex.l,1)-(T_hip_flex_l + ...
                Tau_passj.hip.flex.l - ...
                0.1*Qdotskj_nsc(jointi.hip_flex.l,j+1)) == 0);
            % Hip flexion, right
            Ft_hip_flex_r   = FTj(mai(1).mus.r',1);
            T_hip_flex_r    = f_T27(MAj.hip_flex.r,Ft_hip_flex_r);
            opti.subject_to(Tj(jointi.hip_flex.r,1)-(T_hip_flex_r + ...
                Tau_passj.hip.flex.r - ...
                0.1*Qdotskj_nsc(jointi.hip_flex.r,j+1)) == 0);
            % Hip adduction, left
            Ft_hip_add_l    = FTj(mai(2).mus.l',1);
            T_hip_add_l     = f_T27(MAj.hip_add.l,Ft_hip_add_l);
            opti.subject_to(Tj(jointi.hip_add.l,1)-(T_hip_add_l + ...
                Tau_passj.hip.add.l - ...
                0.1*Qdotskj_nsc(jointi.hip_add.l,j+1)) == 0);
            % Hip adduction, right
            Ft_hip_add_r    = FTj(mai(2).mus.r',1);
            T_hip_add_r     = f_T27(MAj.hip_add.r,Ft_hip_add_r);
            opti.subject_to(Tj(jointi.hip_add.r,1)-(T_hip_add_r + ...
                Tau_passj.hip.add.r - ...
                0.1*Qdotskj_nsc(jointi.hip_add.r,j+1)) == 0);
            % Hip rotation, left
            Ft_hip_rot_l    = FTj(mai(3).mus.l',1);
            T_hip_rot_l     = f_T27(MAj.hip_rot.l,Ft_hip_rot_l);
            opti.subject_to(Tj(jointi.hip_rot.l,1)-(T_hip_rot_l + ...
                Tau_passj.hip.rot.l - ...
                0.1*Qdotskj_nsc(jointi.hip_rot.l,j+1)) == 0);
            % Hip rotation, right
            Ft_hip_rot_r    = FTj(mai(3).mus.r',1);
            T_hip_rot_r     = f_T27(MAj.hip_rot.r,Ft_hip_rot_r);
            opti.subject_to(Tj(jointi.hip_rot.r,1)-(T_hip_rot_r + ...
                Tau_passj.hip.rot.r - ...
                0.1*Qdotskj_nsc(jointi.hip_rot.r,j+1)) == 0);
            % Knee, left
            Ft_knee_l       = FTj(mai(4).mus.l',1);
            T_knee_l        = f_T13(MAj.knee.l,Ft_knee_l);
            opti.subject_to(Tj(jointi.knee.l,1)-(T_knee_l + ...
                Tau_passj.knee.l - 0.1*Qdotskj_nsc(jointi.knee.l,j+1)) == 0);   
            % Knee, right
            Ft_knee_r       = FTj(mai(4).mus.r',1);
            T_knee_r        = f_T13(MAj.knee.r,Ft_knee_r);
            opti.subject_to(Tj(jointi.knee.r,1)-(T_knee_r + ...
                Tau_passj.knee.r - 0.1*Qdotskj_nsc(jointi.knee.r,j+1)) == 0);
            % Ankle, left
            Ft_ankle_l      = FTj(mai(5).mus.l',1);
            T_ankle_l       = f_T12(MAj.ankle.l,Ft_ankle_l);
            opti.subject_to(Tj(jointi.ankle.l,1)-(T_ankle_l + ...
                Tau_passj.ankle.l - 0.1*Qdotskj_nsc(jointi.ankle.l,j+1)) == 0);
            % Ankle, right
            Ft_ankle_r      = FTj(mai(5).mus.r',1);
            T_ankle_r       = f_T12(MAj.ankle.r,Ft_ankle_r);
            opti.subject_to(Tj(jointi.ankle.r,1)-(T_ankle_r + ...
                Tau_passj.ankle.r - 0.1*Qdotskj_nsc(jointi.ankle.r,j+1)) == 0);
            % Subtalar, left
            Ft_subt_l       = FTj(mai(6).mus.l',1);
            T_subt_l        = f_T12(MAj.subt.l,Ft_subt_l);
            opti.subject_to(Tj(jointi.subt.l,1)-(T_subt_l + ...
                Tau_passj.subt.l - 0.1*Qdotskj_nsc(jointi.subt.l,j+1)) == 0);   
            % Subtalar, right
            Ft_subt_r       = FTj(mai(6).mus.r',1);
            T_subt_r        = f_T12(MAj.subt.r,Ft_subt_r);
            opti.subject_to(Tj(jointi.subt.r,1)-(T_subt_r + ...
                Tau_passj.subt.r - 0.1*Qdotskj_nsc(jointi.subt.r,j+1)) == 0);
            % Lumbar extension
            Ft_trunk_ext    = FTj([mai(7).mus.l,mai(7).mus.r]',1);
            T_trunk_ext     = f_T6(MAj.trunk_ext,Ft_trunk_ext);
            opti.subject_to(Tj(jointi.trunk.ext,1)-(T_trunk_ext + ...
                Tau_passj.trunk.ext - ...
                0.1*Qdotskj_nsc(jointi.trunk.ext,j+1)) == 0);
            % Lumbar bending
            Ft_trunk_ben    = FTj([mai(8).mus.l,mai(8).mus.r]',1);
            T_trunk_ben     = f_T6(MAj.trunk_ben,Ft_trunk_ben);
            opti.subject_to(Tj(jointi.trunk.ben,1)-(T_trunk_ben + ...
                Tau_passj.trunk.ben - ...
                0.1*Qdotskj_nsc(jointi.trunk.ben,j+1)) == 0);
            % Lumbar rotation
            Ft_trunk_rot    = FTj([mai(9).mus.l,mai(9).mus.r]',1);
            T_trunk_rot     = f_T6(MAj.trunk_rot,Ft_trunk_rot);
            opti.subject_to(Tj(jointi.trunk.rot,1)-(T_trunk_rot + ...
                Tau_passj.trunk.rot - ...
                0.1*Qdotskj_nsc(jointi.trunk.rot,j+1)) == 0);
            % Torque-driven joint torques for the arms
            % Arms
            opti.subject_to(Tj(armsi,1)/scaling.ArmTau - a_akj(:,j+1) + ...
                0.1/scaling.ArmTau*Qdotskj_nsc(armsi,j+1) == 0);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Activation dynamics (implicit formulation)            
            act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
            act2 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tact);
            opti.subject_to(act1 >= 0);
            opti.subject_to(act2 <= 1/tact);       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Contraction dynamics (implicit formulation)
            opti.subject_to(Hilldiffj == 0);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constraints to prevent parts of the skeleton to penetrate each
            % other.
            % Origins calcaneus (transv plane) at minimum 9 cm from each other.
            opti.subject_to(0.0081 < f_Jnn2(Tj(calcOr.r,1) - Tj(calcOr.l,1))<4);  
            % Constraint to prevent the arms to penetrate the skeleton       
            % Origins femurs and ipsilateral hands (transv plane) at minimum 
            % 18 cm from each other.
            opti.subject_to(0.0324 < f_Jnn2(Tj(femurOr.r,1)-Tj(handOr.r,1))<4);
            opti.subject_to(0.0324 < f_Jnn2(Tj(femurOr.l,1)-Tj(handOr.l,1))<4);
            % Origins tibia (transv plane) at minimum 11 cm from each other.   
            opti.subject_to(0.0121 < f_Jnn2(Tj(tibiaOr.r,1)-Tj(tibiaOr.l,1))<4);     
        end % End loop over collocation points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add equality constraints (next interval starts with end values of 
        % states from previous interval)
        opti.subject_to(a(:,k+1) == akj*D);
        opti.subject_to(FTtilde(:,k+1) == FTtildekj*D); % scaled
        opti.subject_to(Qs(:,k+1) == Qskj*D); % scaled
        opti.subject_to(Qdots(:,k+1) == Qdotskj*D); % scaled
        opti.subject_to(a_a(:,k+1) == a_akj*D);          
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Additional path constraints
    % Periodicity of the states
    % Qs and Qdots
    QsInvA = [jointi.pelvis.tilt,...
        jointi.pelvis.ty,...
        jointi.hip_flex.l:jointi.trunk.ext,...
        jointi.sh_flex.l:jointi.elb.r]';
    QsInvB = [jointi.pelvis.tilt,...
        jointi.pelvis.ty,...
        jointi.hip_flex.r:jointi.hip_rot.r,...
        jointi.hip_flex.l:jointi.hip_rot.l,...
        jointi.knee.r,...
        jointi.knee.l,...
        jointi.ankle.r,...
        jointi.ankle.l,...
        jointi.subt.r,...
        jointi.subt.l,...
        jointi.trunk.ext,...
        jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r,...
        jointi.elb.l]';   
    
    QdotsInvA = [jointi.pelvis.tilt,...
        jointi.pelvis.tx,jointi.pelvis.ty,...
        jointi.hip_flex.l:jointi.trunk.ext,...
        jointi.sh_flex.l:jointi.elb.r]';
    QdotsInvB = [jointi.pelvis.tilt,...
        jointi.pelvis.tx,jointi.pelvis.ty,...
        jointi.hip_flex.r:jointi.hip_rot.r,...
        jointi.hip_flex.l:jointi.hip_rot.l,...
        jointi.knee.r,...
        jointi.knee.l,...
        jointi.ankle.r,...
        jointi.ankle.l,...
        jointi.subt.r,...
        jointi.subt.l,...
        jointi.trunk.ext,...
        jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r,...
        jointi.elb.l]';     
    
    orderQsOpp = [jointi.pelvis.list:jointi.pelvis.list,...   
        jointi.pelvis.rot:jointi.pelvis.rot,...
        jointi.pelvis.tz:jointi.pelvis.tz,...
        jointi.trunk.ben:jointi.trunk.ben,...
        jointi.trunk.rot:jointi.trunk.rot];
    
    opti.subject_to(Qs(QsInvA,end) - Qs(QsInvB,1) == 0);
    opti.subject_to(Qdots(QdotsInvA,end) - Qdots(QdotsInvB,1) == 0);
    opti.subject_to(Qs(orderQsOpp,end) + Qs(orderQsOpp,1) == 0);  
    opti.subject_to(Qdots(orderQsOpp,end) + Qdots(orderQsOpp,1) == 0);  
    % Muscle activations
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
    opti.subject_to(a(:,end) - a(orderMusInv,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(:,end) - FTtilde(orderMusInv,1) == 0);
    % Arm activations
    orderArmInv = [jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r:jointi.elb.r,...
    jointi.elb.l:jointi.elb.l]-jointi.sh_flex.l+1;
    opti.subject_to(a_a(:,end) - a_a(orderArmInv,1) == 0);
    % Average speed
    vel_aver_tot = dist_trav_tot/tf; 
    opti.subject_to(vel_aver_tot - v_tgt == 0)  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    % Create NLP solver
    opti.minimize(J);
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy      = 'adaptive';
    options.ipopt.max_iter = 10000;
    options.ipopt.tol = 1*10^(-tol_ipopt);
    opti.solver('ipopt', options);  
    % Create and save diary
    p = mfilename('fullpath');
    [~,namescript,~] = fileparts(p);
    pathresults = [pathRepo,'/Results'];
    if ~(exist([pathresults,'/',namescript],'dir')==7)
        mkdir(pathresults,namescript);
    end
    if (exist([pathresults,'/',namescript,'/D',savename],'file')==2)
        delete ([pathresults,'/',namescript,'/D',savename])
    end 
    diary([pathresults,'/',namescript,'/D',savename]);  
    % Data-informed (full solution at closest speed) initial guess    
    if IGm == 4   
        disp('Not supported')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve problem
    % Opti does not use bounds on variables but constraints. This function
    % adjusts for that.    
    [w_opt,stats] = solve_NLPSOL(opti,options);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    diary off
    % Extract results 
    % Create setup
    setup.tolerance.ipopt = tol_ipopt;
    setup.bounds = bounds;
    setup.scaling = scaling;
    setup.guess = guess;
    % Save results and set    
    save([pathresults,'/',namescript,'/w',savename],'w_opt');
    save([pathresults,'/',namescript,'/stats',savename],'stats');
end

%% Analyze results
if analyseResults
    %% Load results
    if loadResults
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathRepo,'/Results'];
        load([pathresults,'/',namescript,'/w',savename]);
        load([pathresults,'/',namescript,'/stats',savename]);
    end 
    NParameters = 1;    
    tf_opt = w_opt(1:NParameters);
    starti = NParameters+1;
    a_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    a_col_opt = reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    FTtilde_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    FTtilde_col_opt =reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    Qs_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    Qs_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);
    Qdots_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    Qdots_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);    
    a_a_opt = reshape(w_opt(starti:starti+nq.arms*(N+1)-1),nq.arms,N+1)';
    starti = starti + nq.arms*(N+1);
    a_a_col_opt = reshape(w_opt(starti:starti+nq.arms*(d*N)-1),nq.arms,d*N)';
    starti = starti + nq.arms*(d*N);
    vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
    starti = starti + NMuscle*N;
    e_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
    starti = starti + nq.arms*N;   
    dFTtilde_col_opt=reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    qdotdot_col_opt =reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,(d*N))';
    starti = starti + nq.all*(d*N);
    % Combine results at mesh and collocation points
    a_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    a_mesh_col_opt(1:(d+1):end,:)= a_opt;
    FTtilde_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    FTtilde_mesh_col_opt(1:(d+1):end,:)= FTtilde_opt;
    Qs_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    Qs_mesh_col_opt(1:(d+1):end,:)= Qs_opt;
    Qdots_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    Qdots_mesh_col_opt(1:(d+1):end,:)= Qdots_opt;
    a_a_mesh_col_opt=zeros(N*(d+1)+1,nq.arms);
    a_a_mesh_col_opt(1:(d+1):end,:)= a_a_opt;
    for k=1:N
        rangei = k*(d+1)-(d-1):k*(d+1);
        rangebi = (k-1)*d+1:k*d;
        a_mesh_col_opt(rangei,:) = a_col_opt(rangebi,:);
        FTtilde_mesh_col_opt(rangei,:) = FTtilde_col_opt(rangebi,:);
        Qs_mesh_col_opt(rangei,:) = Qs_col_opt(rangebi,:);
        Qdots_mesh_col_opt(rangei,:) = Qdots_col_opt(rangebi,:);
        a_a_mesh_col_opt(rangei,:) = a_a_col_opt(rangebi,:);
    end
    
    %% Unscale results
    % States at mesh points
    % Qs (1:N-1)
    q_opt_unsc.rad = Qs_opt(1:end-1,:).*repmat(...
        scaling.Qs,size(Qs_opt(1:end-1,:),1),1); 
    % Convert in degrees
    q_opt_unsc.deg = q_opt_unsc.rad;
    q_opt_unsc.deg(:,roti) = q_opt_unsc.deg(:,roti).*180/pi;
    % Qs (1:N)
    q_opt_unsc_all.rad = Qs_opt.*repmat(scaling.Qs,size(Qs_opt,1),1); 
    % Convert in degrees
    q_opt_unsc_all.deg = q_opt_unsc_all.rad;
    q_opt_unsc_all.deg(:,roti) = q_opt_unsc_all.deg(:,roti).*180/pi;       
    % Qdots (1:N-1)
    qdot_opt_unsc.rad = Qdots_opt(1:end-1,:).*repmat(...
        scaling.Qdots,size(Qdots_opt(1:end-1,:),1),1);
    % Convert in degrees
    qdot_opt_unsc.deg = qdot_opt_unsc.rad;
    qdot_opt_unsc.deg(:,roti) = qdot_opt_unsc.deg(:,roti).*180/pi;
    % Qdots (1:N)
    qdot_opt_unsc_all.rad =Qdots_opt.*repmat(scaling.Qdots,size(Qdots_opt,1),1); 
    % Muscle activations (1:N-1)
    a_opt_unsc = a_opt(1:end-1,:).*repmat(...
        scaling.a,size(a_opt(1:end-1,:),1),size(a_opt,2));
    % Muscle-tendon forces (1:N-1)
    FTtilde_opt_unsc = FTtilde_opt(1:end-1,:).*repmat(...
        scaling.FTtilde,size(FTtilde_opt(1:end-1,:),1),1);
    % Arm activations (1:N-1)
    a_a_opt_unsc = a_a_opt(1:end-1,:).*repmat(...
        scaling.a_a,size(a_a_opt(1:end-1,:),1),size(a_a_opt,2));
    % Arm activations (1:N)
    a_a_opt_unsc_all = ...
        a_a_opt.*repmat(scaling.a_a,size(a_a_opt,1),size(a_a_opt,2));
    % Controls at mesh points
    % Time derivative of muscle activations (states)
    vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
    % Get muscle excitations from time derivative of muscle activations
    e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
    % Arm excitations
    e_a_opt_unsc = e_a_opt.*repmat(scaling.e_a,size(e_a_opt,1),...
        size(e_a_opt,2));
    % States at collocation points
    % Qs
    q_col_opt_unsc.rad = Qs_col_opt.*repmat(scaling.Qs,size(Qs_col_opt,1),1); 
    % Convert in degrees
    q_col_opt_unsc.deg = q_col_opt_unsc.rad;
    q_col_opt_unsc.deg(:,roti) = q_col_opt_unsc.deg(:,roti).*180/pi;
    % Qdots
    qdot_col_opt_unsc.rad = Qdots_col_opt.*repmat(...
        scaling.Qdots,size(Qdots_col_opt,1),1); 
    % Convert in degrees
    qdot_col_opt_unsc.deg = qdot_col_opt_unsc.rad;
    qdot_col_opt_unsc.deg(:,roti) = qdot_col_opt_unsc.deg(:,roti).*180/pi;
    % Muscle activations
    a_col_opt_unsc = a_col_opt.*repmat(...
        scaling.a,size(a_col_opt,1),size(a_col_opt,2));
    % Muscle-tendon forces
    FTtilde_col_opt_unsc = FTtilde_col_opt.*repmat(...
        scaling.FTtilde,size(FTtilde_col_opt,1),1);
    % Arm activations
    a_a_col_opt_unsc = a_a_col_opt.*repmat(...
        scaling.a_a,size(a_a_col_opt,1),size(a_a_col_opt,2));    
    % "Slack" controls at collocation points    
    % Time derivative of Qdots
    qdotdot_col_opt_unsc.rad = ...
        qdotdot_col_opt.*repmat(scaling.Qdotdots,size(qdotdot_col_opt,1),1);
    % Convert in degrees
    qdotdot_col_opt_unsc.deg = qdotdot_col_opt_unsc.rad;
    qdotdot_col_opt_unsc.deg(:,roti) = qdotdot_col_opt_unsc.deg(:,roti).*180/pi;    
    % Time derivative of muscle-tendon forces
    dFTtilde_col_opt_unsc = dFTtilde_col_opt.*repmat(...
        scaling.dFTtilde,size(dFTtilde_col_opt,1),size(dFTtilde_col_opt,2));    
    
    %% Time grid    
    % Mesh points
    tgrid = linspace(0,tf_opt,N+1);
    dtime = zeros(1,d+1);
    for i=1:4
        dtime(i)=tau_root(i)*(tf_opt/N);
    end
    % Mesh points and collocation points
    tgrid_ext = zeros(1,(d+1)*N+1);
    for i=1:N
        tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
    end
    tgrid_ext(end)=tf_opt;     
 
    %% Joint torques and ground reaction forces at mesh points (N-1), except #1
    Xk_Qs_Qdots_opt             = zeros(N,2*nq.all);
    Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc_all.rad(2:end,:);
    Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc_all.rad(2:end,:);
    Xk_Qdotdots_opt             = qdotdot_col_opt_unsc.rad(d:d:end,:);
    Foutk_opt = zeros(N,nq.all+NGRF+NcalcOrall);
    for i = 1:N
        [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
        Foutk_opt(i,:) = full(res);    
    end
    GRFk_opt = Foutk_opt(:,GRFi.all);
    % assertArmTmaxj should be 0
    assertArmTmaxk = max(max(abs(Foutk_opt(:,armsi) - ...
        (a_a_opt_unsc_all(2:end,:)*scaling.ArmTau - ...
        0.1*Xk_Qs_Qdots_opt(:,armsi*2))))); 
    if assertArmTmaxk > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing residual forces')
    end 

    %% Joint torques and ground reaction forces at collocation points
    Xj_Qs_Qdots_opt             = zeros(d*N,2*nq.all);
    Xj_Qs_Qdots_opt(:,1:2:end)  = q_col_opt_unsc.rad;
    Xj_Qs_Qdots_opt(:,2:2:end)  = qdot_col_opt_unsc.rad;
    Xj_Qdotdots_opt             = qdotdot_col_opt_unsc.rad;
    Foutj_opt = zeros(d*N,nq.all+NGRF+NcalcOrall);
    for i = 1:d*N
        [res] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)']);
        Foutj_opt(i,:) = full(res);    
    end
    GRFj_opt = Foutj_opt(:,GRFi.all);
    % assertArmTmaxj should be 0
    assertArmTmaxj = max(max(abs(Foutj_opt(:,armsi)-(a_a_col_opt_unsc*...
        scaling.ArmTau - 0.1*Xj_Qs_Qdots_opt(:,armsi*2))))); 
    if assertArmTmaxj > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing residual forces')
    end 

    %% Stride length and width  
    % For the stride length we also need the values at the end of the
    % interval so N+1 where states but not controls are defined
    Xk_Qs_Qdots_opt_all = zeros(N+1,2*size(q_opt_unsc_all.rad,2));
    Xk_Qs_Qdots_opt_all(:,1:2:end)  = q_opt_unsc_all.rad;
    Xk_Qs_Qdots_opt_all(:,2:2:end)  = qdot_opt_unsc_all.rad;
    % We just want to extract the positions of the calcaneus origins so we
    % do not really care about Qdotdot that we set to 0
    Xk_Qdotdots_opt_all = zeros(N+1,size(q_opt_unsc_all.rad,2));
    out_res_opt_all = zeros(N+1,nq.all+NGRF+NcalcOrall);
    for i = 1:N+1
        [res] = F1([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)']);
        out_res_opt_all(i,:) = full(res);    
    end
    % The stride length is the distance covered by the calcaneus origin
    % Right leg
    dist_r = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.r)-...
        out_res_opt_all(1,calcOrall.r)));
    % Left leg
    dist_l = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.l)-...
        out_res_opt_all(1,calcOrall.l)));
    % The total stride length is the sum of the right and left stride 
    % lengths after a half gait cycle, since we assume symmetry
    StrideLength_opt = full(dist_r + dist_l);    
    % The stride width is the medial distance between the calcaneus origins
    StepWidth_opt = full(abs(out_res_opt_all(:,calcOrall.r(3)) - ...
        out_res_opt_all(:,calcOrall.l(3))));
    stride_width_mean = mean(StepWidth_opt);
    stride_width_std = std(StepWidth_opt);
    
    %% Passive joint torques at optimal solution
    Tau_pass_opt_all = zeros(N,15);
    i = 1;
    for k = 1:N    
        for j=1:d
            Tau_pass_opt.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
               theta.pass.hip.flex,Xj_Qs_Qdots_opt(i,jointi.hip_flex.l*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.hip_flex.l*2));
            Tau_pass_opt.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
               theta.pass.hip.flex,Xj_Qs_Qdots_opt(i,jointi.hip_flex.r*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.hip_flex.r*2));
            Tau_pass_opt.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
               theta.pass.hip.add,Xj_Qs_Qdots_opt(i,jointi.hip_add.l*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.hip_add.l*2));
            Tau_pass_opt.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
               theta.pass.hip.add,Xj_Qs_Qdots_opt(i,jointi.hip_add.r*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.hip_add.r*2));
            Tau_pass_opt.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
               theta.pass.hip.rot,Xj_Qs_Qdots_opt(i,jointi.hip_rot.l*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.hip_rot.l*2));
            Tau_pass_opt.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
               theta.pass.hip.rot,Xj_Qs_Qdots_opt(i,jointi.hip_rot.r*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.hip_rot.r*2));
            Tau_pass_opt.knee.l        = f_PassiveMoments(k_pass.knee,...
               theta.pass.knee,Xj_Qs_Qdots_opt(i,jointi.knee.l*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.knee.l*2));
            Tau_pass_opt.knee.r        = f_PassiveMoments(k_pass.knee,...
               theta.pass.knee,Xj_Qs_Qdots_opt(i,jointi.knee.r*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.knee.r*2));
            Tau_pass_opt.ankle.l       = f_PassiveMoments(k_pass.ankle,...
               theta.pass.ankle,Xj_Qs_Qdots_opt(i,jointi.ankle.l*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.ankle.l*2));
            Tau_pass_opt.ankle.r       = f_PassiveMoments(k_pass.ankle,...
               theta.pass.ankle,Xj_Qs_Qdots_opt(i,jointi.ankle.r*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.ankle.r*2));
            Tau_pass_opt.subt.l       = f_PassiveMoments(k_pass.subt,...
               theta.pass.subt,Xj_Qs_Qdots_opt(i,jointi.subt.l*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.subt.l*2));
            Tau_pass_opt.subt.r       = f_PassiveMoments(k_pass.subt,...
               theta.pass.subt,Xj_Qs_Qdots_opt(i,jointi.subt.r*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.subt.r*2));
            Tau_pass_opt.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
               theta.pass.trunk.ext,Xj_Qs_Qdots_opt(i,jointi.trunk.ext*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.trunk.ext*2));
            Tau_pass_opt.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
               theta.pass.trunk.ben,Xj_Qs_Qdots_opt(i,jointi.trunk.ben*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.trunk.ben*2));
            Tau_pass_opt.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
               theta.pass.trunk.rot,Xj_Qs_Qdots_opt(i,jointi.trunk.rot*2-1),...
               Xj_Qs_Qdots_opt(i,jointi.trunk.rot*2));        
            Tau_pass_opt_all(i,:) = full([Tau_pass_opt.hip.flex.l,...
               Tau_pass_opt.hip.add.l,Tau_pass_opt.hip.rot.l,...             
               Tau_pass_opt.hip.flex.r,Tau_pass_opt.hip.add.r,...
               Tau_pass_opt.hip.rot.r,Tau_pass_opt.knee.l,...
               Tau_pass_opt.knee.r,Tau_pass_opt.ankle.l,...
               Tau_pass_opt.ankle.r,Tau_pass_opt.subt.l,...
               Tau_pass_opt.subt.r,Tau_pass_opt.trunk.ext,...
               Tau_pass_opt.trunk.ben,Tau_pass_opt.trunk.rot]);
           i = i + 1;
        end
    end    
    
    %% Assert average speed    
    dist_trav_opt = q_opt_unsc_all.rad(end,jointi.pelvis.tx) - ...
        q_opt_unsc_all.rad(1,jointi.pelvis.tx); % distance traveled
    time_elaps_opt = tf_opt; % time elapsed
    vel_aver_opt = dist_trav_opt/time_elaps_opt; 
    % assert_v_tg should be 0
    assert_v_tg = abs(vel_aver_opt-v_tgt);
    if assert_v_tg > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing average speed')
    end 
    
    %% Decompose optimal cost
if decomposeCost
    J_opt           = 0;
    E_cost          = 0;
    A_cost          = 0;
    Arm_cost        = 0;
    Qdotdot_cost    = 0;
    Pass_cost       = 0;
    GRF_cost        = 0;
    vA_cost         = 0;
    dFTtilde_cost   = 0;
    QdotdotArm_cost = 0;
    count           = 1;
    h_opt           = tf_opt/N;
    for k=1:N     
        for j=1:d 
            % Get muscle-tendon lengths, velocities, moment arms
            % Left leg
            qin_l_opt_all = [Xj_Qs_Qdots_opt(count,jointi.hip_flex.l*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.hip_add.l*2-1), ...
                Xj_Qs_Qdots_opt(count,jointi.hip_rot.l*2-1), ...
                Xj_Qs_Qdots_opt(count,jointi.knee.l*2-1), ...
                Xj_Qs_Qdots_opt(count,jointi.ankle.l*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.subt.l*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ext*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ben*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.rot*2-1)];  
            qdotin_l_opt_all = [Xj_Qs_Qdots_opt(count,jointi.hip_flex.l*2),...
                Xj_Qs_Qdots_opt(count,jointi.hip_add.l*2),...
                Xj_Qs_Qdots_opt(count,jointi.hip_rot.l*2),...
                Xj_Qs_Qdots_opt(count,jointi.knee.l*2),...
                Xj_Qs_Qdots_opt(count,jointi.ankle.l*2),...
                Xj_Qs_Qdots_opt(count,jointi.subt.l*2),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ext*2),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ben*2),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.rot*2)];  
            [lMTk_l_opt_all,vMTk_l_opt_all,~] = ...
                f_lMT_vMT_dM(qin_l_opt_all,qdotin_l_opt_all);    
            % Right leg
            qin_r_opt_all = [Xj_Qs_Qdots_opt(count,jointi.hip_flex.r*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.hip_add.r*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.hip_rot.r*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.knee.r*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.ankle.r*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.subt.r*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ext*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ben*2-1),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.rot*2-1)];  
            qdotin_r_opt_all = [Xj_Qs_Qdots_opt(count,jointi.hip_flex.r*2),...
                Xj_Qs_Qdots_opt(count,jointi.hip_add.r*2),...
                Xj_Qs_Qdots_opt(count,jointi.hip_rot.r*2),...
                Xj_Qs_Qdots_opt(count,jointi.knee.r*2),...
                Xj_Qs_Qdots_opt(count,jointi.ankle.r*2),...
                Xj_Qs_Qdots_opt(count,jointi.subt.r*2),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ext*2),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.ben*2),...
                Xj_Qs_Qdots_opt(count,jointi.trunk.rot*2)];      
            [lMTk_r_opt_all,vMTk_r_opt_all,~] = ...
                f_lMT_vMT_dM(qin_r_opt_all,qdotin_r_opt_all);
            % Both legs
            lMTk_lr_opt_all = ...
                [lMTk_l_opt_all([1:43,47:49],1);lMTk_r_opt_all(1:46,1)];
            vMTk_lr_opt_all = ...
                [vMTk_l_opt_all([1:43,47:49],1);vMTk_r_opt_all(1:46,1)];        
            % Metabolic energy rate
            [~,~,Fce_opt_all,Fpass_opt_all,Fiso_opt_all,vMmax_opt_all,...
                massM_opt_all] = f_forceEquilibrium_FtildeState_all_tendon(...
                    a_col_opt_unsc(count,:)',FTtilde_col_opt_unsc(count,:)',...
                    dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
                    full(vMTk_lr_opt_all),tensions,aTendon,shift);                  
            [~,lMtilde_opt_all] = f_FiberLength_TendonForce_tendon(...
                FTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),aTendon,...
                shift);                
            [vM_opt_all,~] = ...
                f_FiberVelocity_TendonForce_tendon(...
                FTtilde_col_opt_unsc(count,:)',...
                dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
                full(vMTk_lr_opt_all),aTendon,shift);   

            if mE == 0 % Bhargava et al. (2004)
                [e_tot_all,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    full(vM_opt_all),full(Fce_opt_all)',full(Fpass_opt_all)',...
                    full(massM_opt_all)',pctsts,full(Fiso_opt_all)',...
                    MTparameters_m(1,:)',body_mass,10); 
            elseif mE == 1 % Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2003all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);                
            elseif mE == 2 % Umberger (2010)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);            
            elseif mE == 3 % Uchida et al. (2016)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2016all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);      
            elseif mE == 4 % Umberger (2010) treating muscle lengthening 
                % heat rate as Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all_hl(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);  
            elseif mE == 5 % Umberger (2010) treating negative mechanical 
                % work as Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all_neg(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);  
            end
            e_tot_opt_all = full(e_tot_all)';
                            
            J_opt = J_opt + 1/(dist_trav_opt)*(...
                W.E*B(j+1)      *...
                    (f_J92exp(e_tot_opt_all,exp_E))/body_mass*h_opt + ...                   
                W.A*B(j+1)      *(f_J92(a_col_opt(count,:)))*h_opt +...      
                W.ArmE*B(j+1)   *(f_J8(e_a_opt(k,:)))*h_opt +...                    
                W.Ak*B(j+1)     *...
                    (f_J21(qdotdot_col_opt(count,residuals_noarmsi)))*h_opt +...                     
                W.passMom*B(j+1)*(f_J15(Tau_pass_opt_all(count,:)))*h_opt + ...                    
                W.u*B(j+1)      *(f_J92(vA_opt(k,:)))*h_opt + ...        
                W.u*B(j+1)      *(f_J92(dFTtilde_col_opt(count,:)))*h_opt + ...                           
                W.u*B(j+1)      *(f_J8(qdotdot_col_opt(count,armsi)))*h_opt);    
            
                E_cost = E_cost + W.E*B(j+1)*...
                    (f_J92exp(e_tot_opt_all,exp_E))/body_mass*h_opt;
                A_cost = A_cost + W.A*B(j+1)*...
                    (f_J92(a_col_opt(count,:)))*h_opt;
                Arm_cost = Arm_cost + W.ArmE*B(j+1)*...
                    (f_J8(e_a_opt(k,:)))*h_opt;
                Qdotdot_cost = Qdotdot_cost + W.Ak*B(j+1)*...
                    (f_J21(qdotdot_col_opt(count,residuals_noarmsi)))*h_opt;
                Pass_cost = Pass_cost + W.passMom*B(j+1)*...
                    (f_J15(Tau_pass_opt_all(count,:)))*h_opt;
                vA_cost = vA_cost + W.u*B(j+1)*...
                    (f_J92(vA_opt(k,:)))*h_opt;
                dFTtilde_cost = dFTtilde_cost + W.u*B(j+1)*...
                    (f_J92(dFTtilde_col_opt(count,:)))*h_opt;
                QdotdotArm_cost = QdotdotArm_cost + W.u*B(j+1)*...
                    (f_J8(qdotdot_col_opt(count,armsi)))*h_opt;                
                count = count + 1;  
        end
    end      
    J_optf = full(J_opt);     
    E_costf = full(E_cost);
    A_costf = full(A_cost);
    Arm_costf = full(Arm_cost);
    Qdotdot_costf = full(Qdotdot_cost);
    Pass_costf = full(Pass_cost);
    vA_costf = full(vA_cost);
    dFTtilde_costf = full(dFTtilde_cost);
    QdotdotArm_costf = full(QdotdotArm_cost);
    % assertCost should be 0 
    assertCost = J_optf - 1/(dist_trav_opt)*(E_costf+A_costf+Arm_costf+...
        Qdotdot_costf+Pass_costf+vA_costf+dFTtilde_costf+...
        QdotdotArm_costf);
    assertCost2 = stats.iterations.obj(end) - J_optf;
    if assertCost > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing optimal cost')
    end 
    if assertCost2 > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing optimal cost')
    end 
end 

    %% Reconstruct full gait cycle
    % We reconstruct the full gait cycle from the simulated half gait cycle
    % Identify heel strike
    threshold = 20; % there is foot-ground contact above the threshold
    % We adjust some thresholds manually
    if ww == 17
        threshold = 31;
    elseif ww == 133
        threshold = 21;    
    elseif ww == 159
        threshold = 42;
    elseif ww == 160
        threshold = 44;
    elseif ww == 161
        threshold = 57;
    elseif ww == 162
        threshold = 69;
    end
        
    if exist('HS1','var')
        clear HS1
    end
    % Check if heel strike is on the right side   
    phase_tran_tgridi = find(GRFk_opt(:,2)<threshold,1,'last');
    if ~isempty(phase_tran_tgridi)        
        if phase_tran_tgridi == N
            temp_idx = find(GRFk_opt(:,2)>threshold,1,'first');
            if ~isempty(temp_idx)
                if temp_idx-1 ~= 0 && ...
                        find(GRFk_opt(temp_idx-1,2)<threshold)
                    phase_tran_tgridi_t = temp_idx;             
                    IC1i = phase_tran_tgridi_t;
                    HS1 = 'r';
                end 
            else            
                IC1i = phase_tran_tgridi + 1; 
                HS1 = 'r';
            end
        else            
            IC1i = phase_tran_tgridi + 1; 
            HS1 = 'r';
        end        
    end
    if ~exist('HS1','var')
        % Check if heel strike is on the left side  
        phase_tran_tgridi = find(GRFk_opt(:,5)<threshold,1,'last');       
        if phase_tran_tgridi == N
            temp_idx = find(GRFk_opt(:,5)>threshold,1,'first');
            if ~isempty(temp_idx)  
                if temp_idx-1 ~= 0 && ...
                        find(GRFk_opt(temp_idx-1,5)<threshold)
                    phase_tran_tgridi_t = temp_idx;             
                    IC1i = phase_tran_tgridi_t;
                    HS1 = 'l';
                else
                    IC1i = phase_tran_tgridi + 1; 
                    HS1 = 'l';
                end 
            else
                IC1i = phase_tran_tgridi + 1; 
                HS1 = 'l';
            end
        else            
            IC1i = phase_tran_tgridi + 1; 
            HS1 = 'l';
        end        
    end
    
    % TODO: GRFk_opt is at mesh points starting from k=2, we should thus
    % add 1 to IC1i
    IC1i = IC1i + 1;
    
    if isempty(phase_tran_tgridi)
        disp('No heel strike detected, consider increasing the threshold');
        continue;
    end

    % Joint kinematics
    % Helper variables to reconstruct full gait cycle assuming symmetry
    QsSymA = [jointi.pelvis.tilt,jointi.pelvis.ty,...
        jointi.hip_flex.l:jointi.trunk.ext,...
        jointi.sh_flex.l:jointi.elb.r];
    QsSymB = [jointi.pelvis.tilt,jointi.pelvis.ty,...    
        jointi.hip_flex.r:jointi.hip_rot.r,...
        jointi.hip_flex.l:jointi.hip_rot.l,...
        jointi.knee.r,jointi.knee.l,...
        jointi.ankle.r,jointi.ankle.l,...
        jointi.subt.r,jointi.subt.l,...
        jointi.trunk.ext,...
        jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r,jointi.elb.l]; 
    QsOpp = [jointi.pelvis.list:jointi.pelvis.rot,jointi.pelvis.tz,...
        jointi.trunk.ben:jointi.trunk.rot];
    QsSymA_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
        jointi.pelvis.ty,...
        jointi.hip_flex.l:jointi.trunk.ext,...
        jointi.sh_flex.l:jointi.elb.r];
    QsSymB_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
        jointi.pelvis.ty,...    
        jointi.hip_flex.r:jointi.hip_rot.r,...
        jointi.hip_flex.l:jointi.hip_rot.l,...
        jointi.knee.r,jointi.knee.l,...
        jointi.ankle.r,jointi.ankle.l,...
        jointi.subt.r,jointi.subt.l,...
        jointi.trunk.ext,...
        jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r,jointi.elb.l]; 
    
    % Qs
    Qs_opt = zeros(N*2,size(q_opt_unsc.deg,2));
    Qs_opt(1:N-IC1i+1,:) = q_opt_unsc.deg(IC1i:end,:);   
    Qs_opt(N-IC1i+2:N-IC1i+1+N,QsSymA) = q_opt_unsc.deg(1:end,QsSymB);
    Qs_opt(N-IC1i+2:N-IC1i+1+N,QsOpp) = -q_opt_unsc.deg(1:end,QsOpp);    
    Qs_opt(N-IC1i+2:N-IC1i+1+N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:end,jointi.pelvis.tx) + ...
        q_opt_unsc_all.deg(end,jointi.pelvis.tx);        
    Qs_opt(N-IC1i+2+N:2*N,:) = q_opt_unsc.deg(1:IC1i-1,:);    
    Qs_opt(N-IC1i+2+N:2*N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:IC1i-1,jointi.pelvis.tx) + ...
        2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')         
        Qs_opt(:,QsSymA_ptx)  = Qs_opt(:,QsSymB_ptx);
        Qs_opt(:,QsOpp)       = -Qs_opt(:,QsOpp);    
    end  
    temp_Qs_opt_pelvis_tx = Qs_opt(1,jointi.pelvis.tx);
    Qs_opt(:,jointi.pelvis.tx) = Qs_opt(:,jointi.pelvis.tx)-...
        temp_Qs_opt_pelvis_tx;
    
    % Qdots
    Qdots_opt = zeros(N*2,size(Qs_opt,2));
    Qdots_opt(1:N-IC1i+1,:) = qdot_opt_unsc.deg(IC1i:end,:);
    Qdots_opt(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
        qdot_opt_unsc.deg(1:end,QsSymB_ptx);
    Qdots_opt(N-IC1i+2:N-IC1i+1+N,QsOpp) = ...
        -qdot_opt_unsc.deg(1:end,QsOpp);
    Qdots_opt(N-IC1i+2+N:2*N,:) = qdot_opt_unsc.deg(1:IC1i-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Qdots_opt(:,QsSymA_ptx) = Qdots_opt(:,QsSymB_ptx);
        Qdots_opt(:,QsOpp) = -Qdots_opt(:,QsOpp);    
    end

%     % Qdotdots
%     qdotdot_opt_GC = zeros(N*2,size(Qs_opt,2));
%     qdotdot_opt_GC(1:N-IC1i+1,:) = qdotdot_opt_unsc.deg(IC1i:end,:);
%     qdotdot_opt_GC(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
%         qdotdot_opt_unsc.deg(1:end,QsSymB_ptx);
%     qdotdot_opt_GC(N-IC1i+2:N-IC1i+1+N,QsOpp) = ...
%         -qdotdot_opt_unsc.deg(1:end,QsOpp);
%     qdotdot_opt_GC(N-IC1i+2+N:2*N,:) = qdotdot_opt_unsc.deg(1:IC1i-1,:);
%     % If the first heel strike was on the left foot then we invert so that
%     % we always start with the right foot, for analysis purpose
%     if strcmp(HS1,'l')
%         qdotdot_opt_GC(:,QsSymA_ptx) = qdotdot_opt_GC(:,QsSymB_ptx);
%         qdotdot_opt_GC(:,QsOpp) = -qdotdot_opt_GC(:,QsOpp);    
%     end

%     % Ground reaction forces
%     GRFs_opt = zeros(N*2,NGRF);
%     GRFs_opt(1:N-IC1i+1,:) = GRFj_opt(IC1i:end,1:6);
%     GRFs_opt(N-IC1i+2:N-IC1i+1+N,:) = GRFj_opt(1:end,[4:6,1:3]);
%     GRFs_opt(N-IC1i+2:N-IC1i+1+N,[3,6]) = ...
%         -GRFs_opt(N-IC1i+2:N-IC1i+1+N,[3,6]);
%     GRFs_opt(N-IC1i+2+N:2*N,:) = GRFj_opt(1:IC1i-1,1:6);
%     GRFs_opt = GRFs_opt./(body_weight/100);
%     % If the first heel strike was on the left foot then we invert so that
%     % we always start with the right foot, for analysis purpose
%     if strcmp(HS1,'l')
%         GRFs_opt(:,[4:6,1:3]) = GRFs_opt(:,:);
%         GRFs_opt(:,[3,6]) = -GRFs_opt(:,[3,6]);
%     end
    
%     % Joint torques
%     Ts_opt = zeros(N*2,size(Qs_opt,2));
%     Ts_opt(1:N-IC1i+1,1:nq.all) = ...
%         Foutj_opt(IC1i:end,1:nq.all)./body_mass;
%     Ts_opt(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
%         Foutj_opt(1:end,QsSymB_ptx)./body_mass;
%     Ts_opt(N-IC1i+2:N-IC1i+1+N,QsOpp) = ...
%         -Foutj_opt(1:end,QsOpp)./body_mass;
%     Ts_opt(N-IC1i+2+N:2*N,1:nq.all) = ...
%         Foutj_opt(1:IC1i-1,1:nq.all)./body_mass;
%     % If the first heel strike was on the left foot then we invert so that
%     % we always start with the right foot, for analysis purpose
%     if strcmp(HS1,'l')
%         Ts_opt(:,QsSymA_ptx) = Ts_opt(:,QsSymB_ptx);
%         Ts_opt(:,QsOpp) = -Ts_opt(:,QsOpp);    
%     end

    % Muscle-Tendon Forces
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
    FTtilde_opt_GC = zeros(N*2,NMuscle);
    FTtilde_opt_GC(1:N-IC1i+1,:) = FTtilde_opt_unsc(IC1i:end,:);
    FTtilde_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
        FTtilde_opt_unsc(1:end,orderMusInv);
    FTtilde_opt_GC(N-IC1i+2+N:2*N,:) = FTtilde_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        FTtilde_opt_GC(:,:) = FTtilde_opt_GC(:,orderMusInv);
    end

    % Muscle activations
    Acts_opt = zeros(N*2,NMuscle);
    Acts_opt(1:N-IC1i+1,:) = a_opt_unsc(IC1i:end,:);
    Acts_opt(N-IC1i+2:N-IC1i+1+N,:) = a_opt_unsc(1:end,orderMusInv);
    Acts_opt(N-IC1i+2+N:2*N,:) = a_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Acts_opt(:,:) = Acts_opt(:,orderMusInv);
    end

%     % Time derivative of muscle-tendon force
%     dFTtilde_opt_GC = zeros(N*2,NMuscle);
%     dFTtilde_opt_GC(1:N-IC1i+1,:) = dFTtilde_opt_unsc(IC1i:end,:);
%     dFTtilde_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
%         dFTtilde_opt_unsc(1:end,orderMusInv);
%     dFTtilde_opt_GC(N-IC1i+2+N:2*N,:) = dFTtilde_opt_unsc(1:IC1i-1,:);
%     % If the first heel strike was on the left foot then we invert so that
%     % we always start with the right foot, for analysis purpose
%     if strcmp(HS1,'l')
%         dFTtilde_opt_GC(:,:) = dFTtilde_opt_GC(:,orderMusInv);
%     end

    % Muscle excitations
    vA_opt_GC = zeros(N*2,NMuscle);
    vA_opt_GC(1:N-IC1i+1,:) = vA_opt_unsc(IC1i:end,:);
    vA_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = vA_opt_unsc(1:end,orderMusInv);
    vA_opt_GC(N-IC1i+2+N:2*N,:) = vA_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        vA_opt_GC(:,:) = vA_opt_GC(:,orderMusInv);
    end
    e_opt_GC = computeExcitationRaasch(Acts_opt,vA_opt_GC,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);   
    
    % Arm activations
    orderArmInv = [jointi.sh_flex.r:jointi.sh_rot.r,...
            jointi.sh_flex.l:jointi.sh_rot.l,...
            jointi.elb.r:jointi.elb.r,...
            jointi.elb.l:jointi.elb.l]-jointi.sh_flex.l+1; 
    a_a_opt_GC = zeros(N*2,nq.arms);
    a_a_opt_GC(1:N-IC1i+1,:) = a_a_opt_unsc(IC1i:end,:);
    a_a_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = a_a_opt_unsc(1:end,orderArmInv);
    a_a_opt_GC(N-IC1i+2+N:2*N,:) = a_a_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        a_a_opt_GC(:,:) = a_a_opt_GC(:,orderArmInv);
    end
    
    % Arm excitations
    e_a_opt_GC = zeros(N*2,nq.arms);
    e_a_opt_GC(1:N-IC1i+1,:) = e_a_opt_unsc(IC1i:end,:);
    e_a_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = e_a_opt_unsc(1:end,orderArmInv);
    e_a_opt_GC(N-IC1i+2+N:2*N,:) = e_a_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        e_a_opt_GC(:,:) = e_a_opt_GC(:,orderArmInv);
    end  

%     % Passive joint torques
%     Tau_pass_opt_inv = [4:6,1:3,8,7,10,9,12,11,13:15];
%     Tau_pass_opt_GC = zeros(N*2,15);
%     Tau_pass_opt_GC(1:N-IC1i+1,:) = Tau_pass_opt_all(IC1i:end,:); 
%     Tau_pass_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
%         Tau_pass_opt_all(1:end,Tau_pass_opt_inv);
%     Tau_pass_opt_GC(N-IC1i+2+N:2*N,:) = Tau_pass_opt_all(1:IC1i-1,:); 
%     % If the first heel strike was on the left foot then we invert so that
%     % we always start with the right foot, for analysis purpose
%     if strcmp(HS1,'l')
%         Tau_pass_opt_GC(:,Tau_pass_opt_inv) = Tau_pass_opt_GC(:,:);
%     end
    
    % Create .mot file for OpenSim GUI
    q_opt_GUI_GC = zeros(2*N,1+nq.all+2);
    q_opt_GUI_GC(1:N-IC1i+1,1) = tgrid(:,IC1i:end-1)';
    q_opt_GUI_GC(N-IC1i+2:N-IC1i+1+N,1)  = tgrid(:,1:end-1)' + tgrid(end);
    q_opt_GUI_GC(N-IC1i+2+N:2*N,1) = tgrid(:,1:IC1i-1)' + 2*tgrid(end);    
    q_opt_GUI_GC(:,2:end-2) = Qs_opt;
    q_opt_GUI_GC(:,end-1:end) = 1.51*180/pi*ones(2*N,2); % pro_sup (locked)
    q_opt_GUI_GC(:,1) = q_opt_GUI_GC(:,1)-q_opt_GUI_GC(1,1);   
    if writeIKmotion
        pathOpenSim = [pathRepo,'/OpenSim'];
        addpath(genpath(pathOpenSim));
        JointAngle.labels = {'time','pelvis_tilt','pelvis_list',...
        'pelvis_rotation','pelvis_tx','pelvis_ty','pelvis_tz',...
        'hip_flexion_l','hip_adduction_l','hip_rotation_l',...
        'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
        'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
        'subtalar_angle_l','subtalar_angle_r',...
        'lumbar_extension','lumbar_bending','lumbar_rotation',...
        'arm_flex_l','arm_add_l','arm_rot_l',...
        'arm_flex_r','arm_add_r','arm_rot_r',...
        'elbow_flex_l','elbow_flex_r',...
        'pro_sup_l','pro_sup_r'};        
        % Two gait cycles
        % Joint angles
        q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
        q_opt_GUI_GC_2(2*N+1:4*N,1) = q_opt_GUI_GC_2(2*N+1:4*N,1) + ...
            q_opt_GUI_GC_2(end,1) + ...
            q_opt_GUI_GC_2(end,1)-q_opt_GUI_GC_2(end-1,1);
        q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) = ...
            q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) + ...
            2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);
        % Muscle activations (to have muscles turning red when activated).
        Acts_opt_GUI = [Acts_opt;Acts_opt];
        % Combine data joint angles and muscle activations
        JointAngleMuscleAct.data = [q_opt_GUI_GC_2,Acts_opt_GUI];
        % Get muscle labels
        muscleNamesAll = cell(1,NMuscle);
        for i = 1:NMuscle/2
            muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
            muscleNamesAll{i+NMuscle/2} = [muscleNames{i}(1:end-2),'_r'];
        end  
        % Combine labels joint angles and muscle activations
        JointAngleMuscleAct.labels = JointAngle.labels;
        for i = 1:NMuscle
            JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC_2,2)} = ...
                [muscleNamesAll{i},'/activation'];
        end
        filenameJointAngles = [pathRepo,'/Results/',namescript,...
                '/IK',savename,'.mot'];
        write_motionFile(JointAngleMuscleAct, filenameJointAngles)
    end
        
    %% Metabolic cost of transport
%     Qs_opt_rad = Qs_opt;
%     Qs_opt_rad(:,roti) = Qs_opt_rad(:,roti).*pi/180;
%     qdot_opt_GC_rad = Qdots_opt;
%     qdot_opt_GC_rad(:,roti)= qdot_opt_GC_rad(:,roti).*pi/180;   
%     % Pre-allocations
%     e_mo_opt = zeros(2*N,1);     
%     e_mo_optb = zeros(2*N,1);    
%     vMtilde_opt_all = zeros(2*N, NMuscle);
%     lMtilde_opt_all = zeros(2*N, NMuscle);
%     for nn = 1:2*N
%         % Get muscle-tendon lengths, velocities, moment arms
%         % Left leg
%         qin_l_opt = [Qs_opt_rad(nn,jointi.hip_flex.l),...
%             Qs_opt_rad(nn,jointi.hip_add.l), ...
%             Qs_opt_rad(nn,jointi.hip_rot.l), ...
%             Qs_opt_rad(nn,jointi.knee.l), ...
%             Qs_opt_rad(nn,jointi.ankle.l),...
%             Qs_opt_rad(nn,jointi.subt.l),...
%             Qs_opt_rad(nn,jointi.trunk.ext),...
%             Qs_opt_rad(nn,jointi.trunk.ben),...
%             Qs_opt_rad(nn,jointi.trunk.rot)];  
%         qdotin_l_opt = [qdot_opt_GC_rad(nn,jointi.hip_flex.l),...
%             qdot_opt_GC_rad(nn,jointi.hip_add.l),...
%             qdot_opt_GC_rad(nn,jointi.hip_rot.l),...
%             qdot_opt_GC_rad(nn,jointi.knee.l),...
%             qdot_opt_GC_rad(nn,jointi.ankle.l),...
%             qdot_opt_GC_rad(nn,jointi.subt.l),...
%             qdot_opt_GC_rad(nn,jointi.trunk.ext),...
%             qdot_opt_GC_rad(nn,jointi.trunk.ben),...
%             qdot_opt_GC_rad(nn,jointi.trunk.rot)];  
%         [lMTk_l_opt,vMTk_l_opt,~] = ...
%             f_lMT_vMT_dM(qin_l_opt,qdotin_l_opt);    
%         % Right leg
%         qin_r_opt = [Qs_opt_rad(nn,jointi.hip_flex.r),...
%             Qs_opt_rad(nn,jointi.hip_add.r),...
%             Qs_opt_rad(nn,jointi.hip_rot.r),...
%             Qs_opt_rad(nn,jointi.knee.r),...
%             Qs_opt_rad(nn,jointi.ankle.r),...
%             Qs_opt_rad(nn,jointi.subt.r),...
%             Qs_opt_rad(nn,jointi.trunk.ext),...
%             Qs_opt_rad(nn,jointi.trunk.ben),...
%             Qs_opt_rad(nn,jointi.trunk.rot)];  
%         qdotin_r_opt = [qdot_opt_GC_rad(nn,jointi.hip_flex.r),...
%             qdot_opt_GC_rad(nn,jointi.hip_add.r),...
%             qdot_opt_GC_rad(nn,jointi.hip_rot.r),...
%             qdot_opt_GC_rad(nn,jointi.knee.r),...
%             qdot_opt_GC_rad(nn,jointi.ankle.r),...
%             qdot_opt_GC_rad(nn,jointi.subt.r),...
%             qdot_opt_GC_rad(nn,jointi.trunk.ext),...
%             qdot_opt_GC_rad(nn,jointi.trunk.ben),...
%             qdot_opt_GC_rad(nn,jointi.trunk.rot)];      
%         [lMTk_r_opt,vMTk_r_opt,~] = ...
%             f_lMT_vMT_dM(qin_r_opt,qdotin_r_opt);
%         % Both legs
%         lMTk_lr_opt     = [lMTk_l_opt([1:43,47:49],1);lMTk_r_opt(1:46,1)];
%         vMTk_lr_opt     = [vMTk_l_opt([1:43,47:49],1);vMTk_r_opt(1:46,1)];    
%         [Hilldiff_optt,FT_optt,Fce_optt,Fpass_optt,Fiso_optt,...
%             vMmax_optt,massM_optt] = ...
%                 f_forceEquilibrium_FtildeState_all_tendon(...
%                 Acts_opt(nn,:)',FTtilde_opt_GC(nn,:)',...
%                 dFTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),...
%                 full(vMTk_lr_opt),tensions,aTendon,shift);             
%         [~,lMtilde_opt] = f_FiberLength_TendonForce_tendon(...
%             FTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),aTendon,shift);         
%         lMtilde_opt_all(nn,:) = full(lMtilde_opt)';   
%         [vM_opt,vMtilde_opt] = ...
%             f_FiberVelocity_TendonForce_tendon(FTtilde_opt_GC(nn,:)',...
%             dFTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),full(vMTk_lr_opt),...
%             aTendon,shift);        
%         vMtilde_opt_all(nn,:) = full(vMtilde_opt)';    
%         if mE == 0 % Bhargava et al. (2004)
%         [~,~,~,~,~,e_mot] = ...
%             fgetMetabolicEnergySmooth2004all(Acts_opt(nn,:)',...
%                 Acts_opt(nn,:)',full(lMtilde_opt),full(vM_opt),...
%                 full(Fce_optt),full(Fpass_optt),full(massM_optt),pctsts,...
%                 full(Fiso_optt),MTparameters_m(1,:)',body_mass,10);           
%         elseif mE == 1 % Umberger et al. (2003)
%             % vMtilde defined for this model as vM/lMopt
%             vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
%             [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2003all(...
%                 Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
%                 vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
%                 full(massM_optt)',pctsts,full(vMmax_optt)',...
%                 full(Fiso_optt)',body_mass,10); 
%         elseif mE == 2 % Umberger (2010)
%             % vMtilde defined for this model as vM/lMopt
%             vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
%             [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2010all(...
%                 Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
%                 vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
%                 full(massM_optt)',pctsts,full(vMmax_optt)',...
%                 full(Fiso_optt)',body_mass,10); 
%         elseif mE == 3 % Uchida et al. (2016)
%             % vMtilde defined for this model as vM/lMopt
%             vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
%             [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2016all(...
%                 Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
%                 vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
%                 full(massM_optt)',pctsts,full(vMmax_optt)',...
%                 full(Fiso_optt)',body_mass,10); 
%         elseif mE == 4 % Umberger (2010) treating muscle lengthening 
%             % heat rate as Umberger et al. (2003)
%             % vMtilde defined for this model as vM/lMopt
%             vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
%             [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2010all_hl(...
%                 Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
%                 vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
%                 full(massM_optt)',pctsts,full(vMmax_optt)',...
%                 full(Fiso_optt)',body_mass,10); 
%         elseif mE == 5 % Umberger (2010) treating negative mechanical 
%             % work as Umberger et al. (2003)
%             % vMtilde defined for this model as vM/lMopt
%             vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
%             [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2010all_neg(...
%                 Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
%                 vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
%                 full(massM_optt)',pctsts,full(vMmax_optt)',...
%                 full(Fiso_optt)',body_mass,10);             
%         end   
%         e_mo_opt(nn,:) = full(e_mot)';   
%         % Bhargava et al. (2004)
%         [~,~,~,~,~,e_motb] = ...
%             fgetMetabolicEnergySmooth2004all(Acts_opt(nn,:)',...
%             Acts_opt(nn,:)',full(lMtilde_opt),full(vM_opt),...
%             full(Fce_optt),full(Fpass_optt),full(massM_optt),pctsts,...
%             full(Fiso_optt),MTparameters_m(1,:)',body_mass,10);
%         e_mo_optb(nn,:) = full(e_motb)'; 
%     end       
%     % Get COT
%     dist_trav_opt_GC = Qs_opt_rad(end,jointi.pelvis.tx) - ...
%         Qs_opt_rad(1,jointi.pelvis.tx); % distance traveled    
%     time_GC = q_opt_GUI_GC(:,1);
%     e_mo_opt_tr = trapz(time_GC,e_mo_opt);
%     e_mo_opt_trb = trapz(time_GC,e_mo_optb);
%     % Cost of transport: J/kg/m
%     % Energy model used for optimization
%     COT_opt_InOpt = e_mo_opt_tr/body_mass/dist_trav_opt_GC; 
%     % Energy model from Bhargava et al. (2004)
%     COT_opt = e_mo_opt_trb/body_mass/dist_trav_opt_GC; 
    
    %% Optimal cost and CPU time
%     pathDiary = [pathresults,'/',namescript,'/D',savename];
%     [CPU_IPOPT,CPU_NLP,~,Cost,~,~,~,~,OptSol] = readDiary(pathDiary);

    %% Save results       
    if saveResults
        if (exist([pathresults,'/',namescript,...
                '/Results_all.mat'],'file')==2) 
            load([pathresults,'/',namescript,'/Results_all.mat']);
        else
            Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
                (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
                (['W_MuscleActivity_',num2str(W.A)]). ...            
                (['W_JointAcceleration_',num2str(W.Ak)]). ...
                (['W_PassiveTorque_',num2str(W.passMom)]). ...
                (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
                (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
                (['InitialGuessType_',num2str(IGsel)]). ...
                (['InitialGuessMode_',num2str(IGm)]). ...
                (['InitialGuessCase_',num2str(IGcase)]). ... 
                (['WeaknessHipActuators_',num2str(h_weak)]). ...
                (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
                (['MetabolicEnergyModel_',num2str(mE)]). ...            
                (['ContactModel_',num2str(cm)]). ...  
                (['Number_MeshIntervals_',num2str(N)]). ...
                (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
                (['CoContraction_',num2str(coCont)]) = ...
                struct('Qs_opt',[]);
        end    
        % Structure Results_all
        Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['CoContraction_',num2str(coCont)]).Qs_opt = ...
            Qs_opt;
        Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['CoContraction_',num2str(coCont)]). ...
            Qdots_opt = Qdots_opt;
%         Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
%             (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
%             (['W_MuscleActivity_',num2str(W.A)]). ...            
%             (['W_JointAcceleration_',num2str(W.Ak)]). ...
%             (['W_PassiveTorque_',num2str(W.passMom)]). ...
%             (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
%             (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
%             (['InitialGuessType_',num2str(IGsel)]). ...
%             (['InitialGuessMode_',num2str(IGm)]). ...
%             (['InitialGuessCase_',num2str(IGcase)]). ... 
%             (['WeaknessHipActuators_',num2str(h_weak)]). ...
%             (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
%             (['MetabolicEnergyModel_',num2str(mE)]). ...            
%             (['ContactModel_',num2str(cm)]). ...  
%             (['Number_MeshIntervals_',num2str(N)]). ...
%             (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
%             (['CoContraction_',num2str(coCont)]).GRFs_opt =... 
%             GRFs_opt;
%         Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
%             (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
%             (['W_MuscleActivity_',num2str(W.A)]). ...            
%             (['W_JointAcceleration_',num2str(W.Ak)]). ...
%             (['W_PassiveTorque_',num2str(W.passMom)]). ...
%             (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
%             (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
%             (['InitialGuessType_',num2str(IGsel)]). ...
%             (['InitialGuessMode_',num2str(IGm)]). ...
%             (['InitialGuessCase_',num2str(IGcase)]). ... 
%             (['WeaknessHipActuators_',num2str(h_weak)]). ...
%             (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
%             (['MetabolicEnergyModel_',num2str(mE)]). ...            
%             (['ContactModel_',num2str(cm)]). ...  
%             (['Number_MeshIntervals_',num2str(N)]). ...
%             (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
%             (['CoContraction_',num2str(coCont)]).Ts_opt = ...
%             Ts_opt;
        Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['CoContraction_',num2str(coCont)]).Acts_opt =...
            Acts_opt;
%         Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
%             (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
%             (['W_MuscleActivity_',num2str(W.A)]). ...            
%             (['W_JointAcceleration_',num2str(W.Ak)]). ...
%             (['W_PassiveTorque_',num2str(W.passMom)]). ...
%             (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
%             (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
%             (['InitialGuessType_',num2str(IGsel)]). ...
%             (['InitialGuessMode_',num2str(IGm)]). ...
%             (['InitialGuessCase_',num2str(IGcase)]). ... 
%             (['WeaknessHipActuators_',num2str(h_weak)]). ...
%             (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
%             (['MetabolicEnergyModel_',num2str(mE)]). ...            
%             (['ContactModel_',num2str(cm)]). ...  
%             (['Number_MeshIntervals_',num2str(N)]). ...
%             (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
%             (['CoContraction_',num2str(coCont)]).COT_opt = ...
%             COT_opt;
        Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['CoContraction_',num2str(coCont)]). ...
            StrideLength_opt = StrideLength_opt;
        Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['CoContraction_',num2str(coCont)]). ...
            StepWidth_opt = stride_width_mean;
%         Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
%             (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
%             (['W_MuscleActivity_',num2str(W.A)]). ...            
%             (['W_JointAcceleration_',num2str(W.Ak)]). ...
%             (['W_PassiveTorque_',num2str(W.passMom)]). ...
%             (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
%             (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
%             (['InitialGuessType_',num2str(IGsel)]). ...
%             (['InitialGuessMode_',num2str(IGm)]). ...
%             (['InitialGuessCase_',num2str(IGcase)]). ... 
%             (['WeaknessHipActuators_',num2str(h_weak)]). ...
%             (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
%             (['MetabolicEnergyModel_',num2str(mE)]). ...            
%             (['ContactModel_',num2str(cm)]). ...  
%             (['Number_MeshIntervals_',num2str(N)]). ...
%             (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
%             (['CoContraction_',num2str(coCont)]). ...
%             CPU_IPOPT = CPU_IPOPT;
%         Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
%             (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
%             (['W_MuscleActivity_',num2str(W.A)]). ...            
%             (['W_JointAcceleration_',num2str(W.Ak)]). ...
%             (['W_PassiveTorque_',num2str(W.passMom)]). ...
%             (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
%             (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
%             (['InitialGuessType_',num2str(IGsel)]). ...
%             (['InitialGuessMode_',num2str(IGm)]). ...
%             (['InitialGuessCase_',num2str(IGcase)]). ... 
%             (['WeaknessHipActuators_',num2str(h_weak)]). ...
%             (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
%             (['MetabolicEnergyModel_',num2str(mE)]). ...            
%             (['ContactModel_',num2str(cm)]). ...  
%             (['Number_MeshIntervals_',num2str(N)]). ...
%             (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
%             (['CoContraction_',num2str(coCont)]).CPU_NLP = ...
%             CPU_NLP;
%         Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
%             (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
%             (['W_MuscleActivity_',num2str(W.A)]). ...            
%             (['W_JointAcceleration_',num2str(W.Ak)]). ...
%             (['W_PassiveTorque_',num2str(W.passMom)]). ...
%             (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
%             (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
%             (['InitialGuessType_',num2str(IGsel)]). ...
%             (['InitialGuessMode_',num2str(IGm)]). ...
%             (['InitialGuessCase_',num2str(IGcase)]). ... 
%             (['WeaknessHipActuators_',num2str(h_weak)]). ...
%             (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
%             (['MetabolicEnergyModel_',num2str(mE)]). ...            
%             (['ContactModel_',num2str(cm)]). ...  
%             (['Number_MeshIntervals_',num2str(N)]). ...
%             (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
%             (['CoContraction_',num2str(coCont)]).Cost = ...
%             Cost;    
%         Results_all.(['Speed_',num2str(v_tgt_id*100)]). ...  
%             (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
%             (['W_MuscleActivity_',num2str(W.A)]). ...            
%             (['W_JointAcceleration_',num2str(W.Ak)]). ...
%             (['W_PassiveTorque_',num2str(W.passMom)]). ...
%             (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
%             (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
%             (['InitialGuessType_',num2str(IGsel)]). ...
%             (['InitialGuessMode_',num2str(IGm)]). ...
%             (['InitialGuessCase_',num2str(IGcase)]). ... 
%             (['WeaknessHipActuators_',num2str(h_weak)]). ...
%             (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
%             (['MetabolicEnergyModel_',num2str(mE)]). ...            
%             (['ContactModel_',num2str(cm)]). ...  
%             (['Number_MeshIntervals_',num2str(N)]). ...
%             (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
%             (['CoContraction_',num2str(coCont)]).OptSol = ...
%             OptSol;   
        Results_all.colheaders.joints = joints;
        Results_all.colheaders.GRF = {'fore_aft_r','vertical_r',...
        'lateral_r','fore_aft_l','vertical_l','lateral_l'};
        for i = 1:NMuscle/2
            Results_all.colheaders.muscles{i} = ...
                [muscleNames{i}(1:end-2),'_l'];
            Results_all.colheaders.muscles{i+NMuscle/2} = ...
                [muscleNames{i}(1:end-2),'_r'];
        end       
        % Save data
        save([pathresults,'/',namescript,...
            '/Results_all_opti_int_NLPSol_order.mat'], 'Results_all');
    end    
end
end
