function [varargout] = PredSim_wrapper_for_app_v2(U,sf)
% --------------------------------------------------------------------------
% PredSim_wrapper_for_app
%   This functions replaces main.m when using the Vitruvian_Man apps
% 
% INPUT:
%   - U -
%   * struct with fields
%       U.savefolder    path to high-level folder to save results
%       U.GroupName     name of subfolder to save results
%       U.ModelName     name of the model, 
%       U.Height        heigt of the model (in m)
%       U.Mass          mass of the model (in kg)
%       U.Speed         target speed, or -1 to maximise speed
%
%   - sf -
%   * struct with scale factors
% 
%
% OUTPUT:
%   - (no outputs) -
% 
% Original author: Lars D'Hondt
% Original date: October/2022
%
% Last edit by: Ellis Van Can
% Last edit date: 30/09/2025
% --------------------------------------------------------------------------


%% set paths
[pathApp,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathApp);
cd(pathRepo);

modelpath = [pathRepo, '\Subjects\gait1018\gait1018.osim'];
output_dir = fullfile(pathRepo, 'Subjects',U.GroupName, U.ModelName);
osim_output_name = [U.ModelName '.osim'];

% Scale + adapt gravity
% scaleOsim(pathRepo, U, sf);


%% changes to the model
[S] = initializeSettings();

%%% scaling strength
S.subject.muscle_strength = {{'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
                                    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r',...
                                    'hamstrings_l','bifemsh_l','glut_max_l','iliopsoas_l',...
                                    'rect_fem_l','vasti_l','gastroc_l','soleus_l','tib_ant_l'},U.Force_sf};

%%% scaling length + contactspheres
scaleOsimModel_2025(pathRepo, U)

import org.opensim.modeling.*
model = Model(fullfile(output_dir,osim_output_name));

%%% scaling max isometric force
if U.Mass ~= 62
muscleSet = model.getMuscles;
old_model_name = char(model.getName.toString);
    if ~contains(old_model_name,'_sf')
    model.setName([old_model_name '_sf'])
    % Loop through all muscles
        for m=1:muscleSet.getSize
            clear muscleName
            muscleName = char(muscleSet.get(m-1));
            model = scaleFMO(muscleName,model,62,U.Mass);
        end
    end
end

%%% scaling bounds+limittorques+cost function weights
% S = scale4PredSim_2025(S,U,62,160);

%%% adapt gravity
gravity_model =  U.sf_Gravity*-9.8066499999999994;
model.setGravity( Vec3(0, gravity_model, 0) );

if ~isfolder(output_dir)
    mkdir(output_dir)
end

%%% save all changes to the model
model.print(fullfile(output_dir,osim_output_name));

%% Inputs PredSim
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings)

S.misc.main_path = pathRepo;

addpath([S.misc.main_path '\VariousFunctions'])
addpath([S.misc.main_path '\AdaptOpenSimModel'])

% name of the subject
S.subject.name = U.ModelName;

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(U.savefolder,U.GroupName,U.ModelName); 

% give the path to the osim model of your subject
osim_path = fullfile(output_dir,osim_output_name);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

% casadi folder
S.solver.CasADi_path    = U.PathCasadi;

S.solver.tol_ipopt = 3;
S.solver.max_iter = 1000;
% S.OpenSimADOptions.compiler = 'Visual Studio 14 2015 Win64';

%%% speed
S.subject.v_pelvis_x_trgt = U.Speed;
%% IG
if strcmp(U.ModelName,'Jullie_lengte')
    S.subject.IG_selection = "C:\GBW_MyPrograms\KinderuniversiteitApp\OCP\Jullie_lengte_v1.mot";
    S.subject.IG_selection_gaitCyclePercent = 200;
elseif strcmp(U.ModelName,'Kleuter')
    S.subject.IG_selection = "C:\GBW_MyPrograms\KinderuniversiteitApp\OCP\Kleuter_v1.mot";
    S.subject.IG_selection_gaitCyclePercent = 200;
elseif strcmp(U.ModelName,'Volwassenen')
    S.subject.IG_selection = "C:\GBW_MyPrograms\KinderuniversiteitApp\OCP\Volwassenen_v1.mot";
    S.subject.IG_selection_gaitCyclePercent = 200;
else
    S.subject.IG_selection = "C:\GBW_MyPrograms\KinderuniversiteitApp\OCP\Ik_Guess_2D.mot";
    S.subject.IG_selection_gaitCyclePercent = 200;
end

% S.subject.adapt_IG_pelvis_y = 1;
% S.subject.IG_pelvis_y = U.IG_pelvis_y;

%%
% S.misc.visualize_bounds    = 1;
S.OpenSimADOptions.verbose_mode = false;

%% Run predictive simulations
savename = run_pred_sim(S,osim_path);

savepath = fullfile(S.subject.save_folder,savename);

if nargout>=1
    varargout{1} = savepath;
end

end