function [MuscleData, varargout] = muscleAnalysis(S,osim_path,model_info,varargin)
% --------------------------------------------------------------------------
% muscleAnalysis
%   Analyse the musculoskeletal geometry of the given osim model.
%       1) Create a dummy motion
%       2) Run the OpenSim muscle analysis tool to calculate the muscle-tendon
%       lengths and moment arms in the model for this dummy motion.
%       3) Import the muscle-tendon lengths and moment arms, and store them
%       along with the dummy motion.
%
%   note 1: The boundaries of the dummy motion are hard-coded.
%
%   note 2: This function is replaced by muscleAnalysisAPI, because using
%   the OpenSim API to evaluate muscle-tendon lengths and momentarms is 
%   approximately 3 times faster than running a muscle analysis. The
%   analysis brings the system to dynamic stage, while the API only
%   requires position stage. In "\Tests\test_muscleAnalysisAPI.m", this
%   function is used to perform a unit test of muscleAnalysisAPI.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - n_data_points (optional input) -
%   * number of data points for the dummy motion. Default is 5000.
%
% OUTPUT:
%   - MuscleData -
%   * Structure with joint angles and according muscle-tendon lengths and
%   momentarms. This will be used to fit expressions to the muscle-tendon 
%   lengths and momentarms.
%
%   - Qs (optional output, for testing purpose) -
%   * Coordinate values of the dummy motion
% 
% Original author: Lars D'Hondt
% Original date: 17/January/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

MA_path = fullfile(S.misc.subject_path, '/MuscleAnalysis');
if ~isfolder(MA_path)
    mkdir(MA_path);
end

%% Create dummy motion

% bounds of training dataset (in degrees)
% hard-coded for now, make input later
% Bounds(1,:) = {'hip_flex',[-50 50]};
% Bounds(2,:) = {'hip_add',[-30 30]};
% Bounds(3,:) = {'hip_rot',[-30 30]};
% Bounds(4,:) = {'knee',[-90 0]};
% Bounds(5,:) = {'ankle',[-30 30]};
% Bounds(6,:) = {'subt',[-30 30]};
% Bounds(7,:) = {'mtj',[-30 30]};
% Bounds(8,:) = {'mtp',[-20 50]};
% Bounds(9,:) = {'lumbar_ext',[-30 30]};
% Bounds(10,:) = {'lumbar_bend',[-30 30]};
% Bounds(11,:) = {'lumbar_rot',[-30 30]};


% coordinate names
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);

% muscle names
muscle_names = model_info.muscle_info.muscle_names;

% sizes and indices
n_coord = length(coordinate_names);
% if length(varargin)>=1
%     n_data_points = varargin{1};
% else
%     n_data_points = 5000;
% end

n_data_points = S.misc.msk_geom_n_samples;

% % default bounds
% Q_bounds = [-30;30]*ones(1,n_coord);
% 
% % adjust bounds
% for i=1:size(Bounds,1)
%     idx = find(contains(coordinate_names,Bounds{i,1}));
%     if ~isempty(idx)
%         Q_bounds(:,idx) = reshape(Bounds{i,2},2,1)*ones(1,length(idx));
%     end
% end
% 
% % construct scale from bounds
% Q_scale = diff(Q_bounds);
% 
% % generate random joint angles (range 0-1)
% Qs = lhsdesign(n_data_points,n_coord);
% 
% % scale and offset random joint angles to fit bounds
% Qs = Qs.*(ones(n_data_points,1)*Q_scale) + ones(n_data_points,1)*Q_bounds(1,:);
% 
% % base does not have to move for analysis
% Qs(:,model_info.ExtFunIO.jointi.floating_base) = 0;
% 
% % angles in degrees
% Qs(:,model_info.ExtFunIO.jointi.rotations) = Qs(:,model_info.ExtFunIO.jointi.rotations)*pi/180;

Qs = generate_dummy_motion(S,model_info,n_data_points);

% time vector
time = (1:n_data_points)./100;

% column labels for dummy motion
q_in.labels = [{'time'},coordinate_names'];

% data for dummy motion
q_in.data = [time',Qs];

% angles are in radians
q_in.inDeg = 'no';

% path to save dummy motion
pathDummyMotion = fullfile(MA_path,'dummy_motion.mot');

% generate motion file
% write_motionFile_v40(q_in,pathDummyMotion)

q_out = read_motionFile_v40(pathDummyMotion);
Qs = q_out.data(:,2:end);

%% Run analysis
% Narrow down coordinates to analyse, to save time. The floating base dofs
% should never have muscles crossing.
idx_coord_analyse = setdiff(1:n_coord,model_info.ExtFunIO.jointi.floating_base);

% import the opensim API
import org.opensim.modeling.*

% run muscle analysis
OpenSim_Muscle_Analysis(pathDummyMotion, osim_path, MA_path, [time(1) time(end)],...
    coordinate_names(idx_coord_analyse));


% import the muscle analysis data
lMT = importdata([MA_path,'/dummy_motion_MuscleAnalysis_Length.sto']);
for i=1:n_coord
    coord_name_i = coordinate_names{i};
    if sum(idx_coord_analyse==i)==1
%         disp(coord_name_i)
        MA.(coord_name_i) = importdata([MA_path,'/dummy_motion_MuscleAnalysis_MomentArm_' coord_name_i '.sto']);
    end
end

%% Store analysis results
% structure and fieldnames of MuscleData are to be consistent with
% PolynomialFit.m

% coordinate names
MuscleData.dof_names = coordinate_names;

% muscle names
MuscleData.muscle_names = muscle_names;

% joint angles training data
MuscleData.q = Qs;

for m = 1:length(muscle_names)
    idx = strcmp(lMT.colheaders,muscle_names{m});

    % muscle-tendon lengths
    MuscleData.lMT(:,m) = lMT.data(:,idx);

    % moment arms
    for i=1:n_coord
        if sum(idx_coord_analyse==i)==1
            MuscleData.dM(:,m,i) = MA.(coordinate_names{i}).data(:,idx);
        else
            MuscleData.dM(:,m,i) = zeros(n_data_points,1);
        end

    end
end

muscle_spanning_joint_info = squeeze( max(abs(MuscleData.dM),[],1) );
muscle_spanning_joint_info(muscle_spanning_joint_info~=0) = 1;

if nargout>=2
    varargout{1} = Qs;
end

if nargout>=2
    varargout{2} = muscle_spanning_joint_info;
end

end % end of function
