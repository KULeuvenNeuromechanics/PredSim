function [f_TrackSynW] = createCasadi_TrackSynW(S,model_info)
% --------------------------------------------------------------------------
% meaningful_function_name
%   Explanation of what this function does. Length depends on function 
%   complexity. Include relevant citations. If applicable, provide an
%   example and/or refer to a unit test. 
% 
%
%   References
%   [1] cite a relevant paper if applicable
%
%   See also related_function_name
%
% INPUT:
%   - input_1 -
%   * brief description of input_1, including data type.
% 
%   - input_2 -
%   * brief description of input_2. For more complex inputs, such as structs,
%   include an example input.
%
%   - input_3 - (optional) Default: default_value
%   * brief description of input_3, including data type.
%
% OUTPUT:
%   - output_1 -
%   * brief description of output_1
% 
% Original author: (First name Last name)
% Original date: (Using "30/May/2022" format avoids confusion)
% --------------------------------------------------------------------------

import casadi.*

idx_m_r = model_info.muscle_info.idx_right;
idx_m_l = model_info.muscle_info.idx_left;

muscleNames_r = model_info.muscle_info.muscle_names(idx_m_r);
muscleNames_l = model_info.muscle_info.muscle_names(idx_m_l);

SynW_rk         = SX.sym('SynW_rk',length(idx_m_r),S.subject.NSyn_r);
SynW_lk         = SX.sym('SynW_lk',length(idx_m_l),S.subject.NSyn_l);

J = 0;


if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle') % same weights right and left
    knownSynW_r_all = unpack_name_value_combinations(S.subject.knownSynW_r,muscleNames_r,[S.subject.TrackSynW_NSyn_r]);
    k = 1;
    for i = 1:length(idx_m_r)
        if ~isnan(knownSynW_r_all(1,i))
            knownSynW_idx(k) = i;
            knownSynW(k,:) = knownSynW_r_all(:,i);
            k = k + 1;
        end
    end
    SynW_rk_sel = SynW_rk(knownSynW_idx,:);
    Jtemp = 0;
    for i = 1:S.subject.TrackSynW_NSyn_r
        for k = 1:length(knownSynW_idx)
            Jtemp = Jtemp + (SynW_rk_sel(k,i)-knownSynW(k,i)).^2;
        end
    end
    Jtemp = Jtemp/(length(knownSynW_idx));
    J = J + Jtemp;
elseif strcmp(S.misc.gaitmotion_type,'FullGaitCycle') % Options: Track only right, only left, or both
    switch S.subject.TrackSynW_side
        case 'onlyRight'
            knownSynW_r_all = unpack_name_value_combinations(S.subject.knownSynW_r,muscleNames_r,[S.subject.TrackSynW_NSyn_r]);
            k = 1;
            for i = 1:length(idx_m_r)
                if ~isnan(knownSynW_r_all(1,i))
                    knownSynW_idx(k) = i;
                    knownSynW(k,:) = knownSynW_r_all(:,i);
                    k = k + 1;
                end
            end
            SynW_rk_sel = SynW_rk(knownSynW_idx,:);
            Jtemp = 0;
            for i = 1:S.subject.TrackSynW_NSyn_r
                for k = 1:length(knownSynW_idx)
                    Jtemp = Jtemp + (SynW_rk_sel(k,i)-knownSynW(k,i)).^2;
                end
            end
            Jtemp = Jtemp/(length(knownSynW_idx));
            J = J + Jtemp;
        case 'onlyLeft'
            knownSynW_l_all = unpack_name_value_combinations(S.subject.knownSynW_l,muscleNames_l,[S.subject.TrackSynW_NSyn_l]);
            k = 1;
            for i = 1:length(idx_m_l)
                if ~isnan(knownSynW_l_all(1,i))
                    knownSynW_idx(k) = i;
                    knownSynW(k,:) = knownSynW_l_all(:,i);
                    k = k + 1;
                end
            end
            SynW_lk_sel = SynW_lk(knownSynW_idx,:);
            Jtemp = 0;
            for i = 1:S.subject.TrackSynW_NSyn_l
                for k = 1:length(knownSynW_idx)
                    Jtemp = Jtemp + (SynW_lk_sel(k,i)-knownSynW(k,i)).^2;
                end
            end
            Jtemp = Jtemp/(length(knownSynW_idx));
            J = J + Jtemp;
        case 'RightLeft'
            % right
            knownSynW_idx = [];
            knownSynW = [];
            knownSynW_r_all = unpack_name_value_combinations(S.subject.knownSynW_r,muscleNames_r,[S.subject.TrackSynW_NSyn_r]);
            k = 1;
            for i = 1:length(idx_m_r)
                if ~isnan(knownSynW_r_all(1,i))
                    knownSynW_idx(k) = i;
                    knownSynW(k,:) = knownSynW_r_all(:,i);
                    k = k + 1;
                end
            end
            SynW_rk_sel = SynW_rk(knownSynW_idx,:);
            Jtemp = 0;
            for i = 1:S.subject.TrackSynW_NSyn_r
                for k = 1:length(knownSynW_idx)
                    Jtemp = Jtemp + (SynW_rk_sel(k,i)-knownSynW(k,i)).^2;
                end
            end
            Jtemp = Jtemp/(length(knownSynW_idx));
            J = J + Jtemp;
            knownSynW_idx = [];
            knownSynW = [];
            % left
            knownSynW_l_all = unpack_name_value_combinations(S.subject.knownSynW_l,muscleNames_l,[S.subject.TrackSynW_NSyn_l]);
            k = 1;
            for i = 1:length(idx_m_l)
                if ~isnan(knownSynW_l_all(1,i))
                    knownSynW_idx(k) = i;
                    knownSynW(k,:) = knownSynW_l_all(:,i);
                    k = k + 1;
                end
            end
            SynW_lk_sel = SynW_lk(knownSynW_idx,:);
            Jtemp = 0;
            for i = 1:S.subject.TrackSynW_NSyn_l
                for k = 1:length(knownSynW_idx)
                    Jtemp = Jtemp + (SynW_lk_sel(k,i)-knownSynW(k,i)).^2;
                end
            end
            Jtemp = Jtemp/(length(knownSynW_idx));
            J = J + Jtemp;
    end
end


f_TrackSynW = Function('f_TrackSynW', {SynW_rk, SynW_lk}, {J},...
    {'SynW_rk', 'SynW_lk'}, {'SynW_track_cost'});

end