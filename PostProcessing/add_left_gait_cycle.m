function R = add_left_gait_cycle(R, model_info)
% add_left_gait_cycle
%   Adds one left gait cycle to the PredSim results structure R, using
%   GRF-based gait event detection. The GRF order is swapped to detect the
%   left stance phase, the corresponding time-varying results are extracted,
%   and the forward position is made continuous and starts at zero.
%
%   The left cycle is stored in R.left with the same layout as the
%   (right) fields in R, so plot functions can use
%   R.left.<category>.<field> exactly like R.<category>.<field>.
%
% INPUT:
%   R          - PredSim results structure (right gait cycle)
%   model_info - PredSim model information structure
%
% OUTPUT:
%   R          - PredSim results structure with added field R.left:
%                  R.left.idx_GC           frame indices of the left cycle
%                  R.left.kinematics       Qs, Qdots, Qddots (+ _rad)
%                  R.left.kinetics         T_ID, ...
%                  R.left.ground_reaction  GRF_r, GRF_l, ...
%                  R.left.muscles          (if R.muscles exists)
%   Fields that do not change between cycles (colheaders, spatiotemp)
%   are not duplicated; use R.colheaders and R.spatiotemp.
%
% Original author: Alexandra Duarte Fernandes
% Original date: June 16, 2026
%
% Last edit by: Ellis Van Can
% Last edit date: October 20, 2026

%% Number of frames in the (right) gait cycle
N = size(R.kinematics.Qs,1);

%% Detect left gait cycle
% Swap GRF order so the algorithm detects the left stance phase
GRFk_opt = [R.ground_reaction.GRF_l, R.ground_reaction.GRF_r];
[idx_GC, idx_GC_base_forward_offset, ~, ~] = ...
    getStancePhaseSimulation(GRFk_opt, model_info.mass/3);

R.left = struct();
R.left.idx_GC = idx_GC;

%% Copy all time-varying fields (rows = frames)
categories = {'kinematics','kinetics','ground_reaction','muscles','metabolics'};
for i = 1:numel(categories)
    cat_i = categories{i};
    if ~isfield(R,cat_i) || ~isstruct(R.(cat_i))
        continue
    end
    fn = fieldnames(R.(cat_i));
    for j = 1:numel(fn)
        val = R.(cat_i).(fn{j});
        if isnumeric(val) && size(val,1) == N
            % e.g. R.kinematics.Qs
            R.left.(cat_i).(fn{j}) = val(idx_GC,:,:);
    
        elseif isstruct(val) && isscalar(val)
            % e.g. R.metabolics.Bhargava2004
            fn_sub = fieldnames(val);
            for k = 1:numel(fn_sub)
                val_sub = val.(fn_sub{k});
                if isnumeric(val_sub) && size(val_sub,1) == N
                    R.left.(cat_i).(fn{j}).(fn_sub{k}) = val_sub(idx_GC,:,:);
                end
            end
        end
    end
end
%% Stance phase indices (R.ground_reaction.idx_stance_r / _l)
% These are frame numbers of the right cycle, not one row per frame, so
% they are not copied above. Convert them to frame numbers of the left
% cycle: new frame k corresponds to old frame idx_GC(k).
sides = {'r','l'};
for s = 1:numel(sides)
    fld = ['idx_stance_' sides{s}];
    if isfield(R.ground_reaction,fld)
        idx_old = R.ground_reaction.(fld);
        if islogical(idx_old)
            % logical mask with one element per frame
            R.left.ground_reaction.(fld) = idx_old(idx_GC);
        else
            % list of frame numbers
            R.left.ground_reaction.(fld) = find(ismember(idx_GC(:), idx_old(:)));
        end
    end
end

%% Correct forward position
% Make forward position continuous and start at zero
i_fw = model_info.ExtFunIO.jointi.base_forward;
dist_trav = R.spatiotemp.dist_trav;
q_fields = {'Qs','Qs_rad'};
for k = 1:numel(q_fields)
    if isfield(R.left.kinematics, q_fields{k})
        q = R.left.kinematics.(q_fields{k});
        q(idx_GC_base_forward_offset,i_fw) = ...
            q(idx_GC_base_forward_offset,i_fw) + dist_trav;
        q(:,i_fw) = q(:,i_fw) - q(1,i_fw);
        R.left.kinematics.(q_fields{k}) = q;
    end
end
%% Correct forward position of the centre of pressure (x = column 1)
% Same correction as base_forward, only for frames in stance of that foot
shift = R.kinematics.Qs(idx_GC(1),i_fw);   % pelvis position at start of left cycle
rows_offset = false(N,1);
rows_offset(idx_GC_base_forward_offset) = true;
for s = 1:numel(sides)
    fld_cop = ['COP_' sides{s}];
    fld_stance = ['idx_stance_' sides{s}];
    if ~isfield(R.left.ground_reaction,fld_cop)
        continue
    end
    cop = R.left.ground_reaction.(fld_cop);
    % frames in stance of this foot
    in_contact = false(N,1);
    if isfield(R.left.ground_reaction,fld_stance)
        in_contact(R.left.ground_reaction.(fld_stance)) = true;
    else
        in_contact = any(cop ~= 0, 2);
    end
    cop(rows_offset & in_contact,1) = cop(rows_offset & in_contact,1) + dist_trav;
    cop(in_contact,1) = cop(in_contact,1) - shift;
    R.left.ground_reaction.(fld_cop) = cop;
end

end