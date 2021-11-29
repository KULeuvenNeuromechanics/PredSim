% This file contains the settings used for the predictive simulations and 
% is loaded in the different files that are processing and analyzing the
% results

settings = [    
    % A. Varying prescribed gait speed
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0;    % 1  
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0;    % 2
    % C. Subject-specific contact model
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0;    % 3
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0;    % 4   
    % C. 2nd Subject-specific contact model
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0;    % 5
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 3, 1, 0, 0, 0, 0, 0, 0;    % 6  
    % D. 3nd Subject-specific contact model
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 4, 0, 0, 0, 0, 0, 0, 0;    % 7
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 4, 1, 0, 0, 0, 0, 0, 0;    % 8 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Added constraint on distance between toes origins
    % D. 3nd Subject-specific contact model
    % Quasi-random initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 1, 4, 0, 0, 0, 0, 0, 0, 0;    % 9
    % Data-informed (walking) initial guess
    1.33, 4, 50, 500, 50000, 1000000, 1000, 2000, 2, 2, 4, 1, 0, 0, 0, 0, 0, 0;    % 10 
];