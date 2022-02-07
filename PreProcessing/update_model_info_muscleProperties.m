function [model_info] = update_model_info_muscleProperties(model_info)

NMuscle = model_info.muscle_info.NMuscle;

model_info.muscle_info.aTendon = 35*ones(NMuscle,1);
model_info.muscle_info.aTendon(model_info.muscle_info.IndexCalf) = 20; % Why different k for calf muslces?
model_info.muscle_info.shift = getShift(model_info.muscle_info.aTendon);

muscleNames = model_info.muscle_info.muscle_names;

model_info.muscle_info.tensions = getSpecificTensions(muscleNames);
model_info.muscle_info.pctsts = getSlowTwitchRatios(muscleNames);

model_info.muscle_info.muscle_mass = GetMuscleMass(muscleNames,model_info.muscle_info.params);