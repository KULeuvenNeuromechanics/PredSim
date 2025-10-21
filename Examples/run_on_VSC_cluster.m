% --------------------------------------------------------------------------
% run_on_VSC_cluster
%   Run PredSim on the VSC cluster. KU Leuven provides compute resources to 
%   researchers in the High Performance Computing service. The HPC clusters 
%   of KU Leuven are part of the Vlaams Supercomputer Centrum  (VSC).
% 
%   https://icts.kuleuven.be/sc/onderzoeksgegevens/hpc
%   https://www.vscentrum.be/
%
% Original author: Lars D'Hondt
% Original date: 01/October/2025
% --------------------------------------------------------------------------

clear
close all
clc


% Check BLAS/LAPACK version; add functions from LinearAlgebra subdirectory
% to path in case Intel is *not* used
blas_version = version('-blas')
lapack_version = version('-lapack')
if ~startsWith(lapack_version, 'Intel')
    addpath(fullfile(getenv('PWD'), 'LinearAlgebra'))
end

[pathExDir,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathExDir);
[pathRepoFolder,~,~] = fileparts(pathRepo);

addpath(fullfile(pathRepo,'DefaultSettings'))
addpath(pathRepo)

% if the OpenSim module is loaded, make its Java library available
if isenv('EBROOTOPENSIM')
    javaclasspath(fullfile(getenv('EBROOTOPENSIM'), 'sdk', 'Java', 'org-opensim-modeling.jar'));
end

% if the CasADi-MATLAB module is loaded, expose its matlab bindings
if isenv('EBROOTCASADIMINMATLAB')
    addpath(fullfile(getenv('EBROOTCASADIMINMATLAB'), 'matlab'))
end

%% Initialize S

[S] = initializeSettings('Falisse_et_al_2022');

%% Settings

% name of the subject
S.subject.name = 'Falisse_et_al_2022';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;
% S.solver.IG_selection = 'quasi-random';

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);


%% Run predictive simulations

[savename] = runPredSim(S, osim_path);

