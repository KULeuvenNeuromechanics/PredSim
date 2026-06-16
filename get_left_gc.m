%Code to get the left side gait cycle 
% Load the data 

results_folder =  "results_path";
SubjectName = 'Subject_name';
resultname = [SubjectName, '.mat'];
mat = fullfile(results_folder, resultname);
result = load(mat);

% Read results
R_right = result.R;
model_info = result.model_info;
dist_trav_opt = R_right.spatiotemp.dist_trav;
Qs_GC_right = R_right.kinematics.Qs;

%Extract GRF
GRFs = R_right.ground_reaction;   
GRF_right = GRFs.GRF_r;              
GRF_left  = GRFs.GRF_l;             

% swap to get left-side gait cycle
GRFk_opt = [GRF_left, GRF_right];
[idx_GC,idx_GC_base_forward_offset,~,threshold] = getStancePhaseSimulation(GRFk_opt,model_info.mass/3);

%% copy R structure then overwrite with left-side GC variable

%kinematics
R_left.kinematics.Qs = R_right.kinematics.Qs(idx_GC, :);
R_left.kinematics.Qdots    = R_right.kinematics.Qdots(idx_GC, :);
R_left.kinematics.Qddots = R_right.kinematics.Qddots(idx_GC, :);
R_left.kinematics.Qs_rad = R_right.kinematics.Qs_rad(idx_GC, :);
R_left.kinematics.Qdots_rad = R_right.kinematics.Qdots_rad(idx_GC, :);
R_left.kinematics.Qddots_rad = R_right.kinematics.Qddots_rad(idx_GC, :);

% adjust forward position to be continuous and start at 0 - only for Qs!
R_left.kinematics.Qs(idx_GC_base_forward_offset,model_info.ExtFunIO.jointi.base_forward) = R_left.kinematics.Qs(idx_GC_base_forward_offset,model_info.ExtFunIO.jointi.base_forward) + dist_trav_opt;
R_left.kinematics.Qs(:,model_info.ExtFunIO.jointi.base_forward) = R_left.kinematics.Qs(:,model_info.ExtFunIO.jointi.base_forward) - R_left.kinematics.Qs(1,model_info.ExtFunIO.jointi.base_forward);

% GRFs
R_left.ground_reaction.GRF_r = GRFs.GRF_r(idx_GC, :);
R_left.ground_reaction.GRF_l = GRFs.GRF_l(idx_GC, :);

% Spatiotemp
R_left.spatiotemp = R_right.spatiotemp;

% Save to a new .mat file
save(mat, 'R_left', '-append');