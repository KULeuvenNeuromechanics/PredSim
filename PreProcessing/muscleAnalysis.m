%
% TO DO: finish this


function [] = muscleAnalysis(S,osimPath,model_info)

MA_path=fullfile(S.subject.name,'MuscleAnalysis');
if ~isfolder(MA_path)
    mkdir(MA_path);
end

%% Create dummy motion

% Training dataset for fitting has same bounds of simulation
[bounds,scaling] = getBounds(S,model_info);
% coordinate names
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);
q_in.labels = [{'time'},coordinate_names'];
% angles are in radians
q_in.inDeg = 'no';
% path with dummy motion
pathDummyMotion = fullfile(MA_path,'dummy_motion.mot');
% motion data
idx_not_base = setdiff((1:length(coordinate_names)),model_info.ExtFunIO.jointi.floating_base);
n = 5000;
p = length(idx_not_base);
X = lhsdesign(n,p);
X_scale = scaling.Qs(idx_not_base);
X_min = bounds.Qs.lower(idx_not_base).*scaling.Qs(idx_not_base);
Angles=X.*(ones(n,1)*X_scale)+(ones(n,1)*X_min);

time=(1:n)./100;
data=zeros(length(time),length(q_in.labels));
data(:,1)=time;
data(:,idx_not_base) = Angles;


q_in.data = data;

% generate motion file
write_motionFile_v40(q_in,pathDummyMotion)

%% Run analysis
% Import the opensim API
import org.opensim.modeling.*

