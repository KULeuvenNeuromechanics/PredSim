function [] = muscleAnalysis(S,osimPath,IO)

MA_path=fullfile(SubjFolder,'MuscleAnalysis');
if ~isfolder(MA_path)
    mkdir(MA_path);
end

muscles = IO.muscle.params.names;
jointi

[bounds,~] = getBounds_all_mtp(Qs_walk,NMuscle,nq,jointi,S.v_tgt);


% Import the opensim API
import org.opensim.modeling.*

