if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle') % same weights right and left
    knownSynW_r_all = unpack_name_value_combinations(S.subject.knownSynW_r,muscleNames_r,[S.subject.TrackSynW_NSyn_r]);
    k = 1;
    for i = 1:NMuscle/2
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
    J = J + W.TrackSynW * Jtemp;
elseif strcmp(S.misc.gaitmotion_type,'FullGaitCycle') % Options: Track only right, only left, or both
    switch S.subject.TrackSynW_side
        case 'onlyRight'
            knownSynW_r_all = unpack_name_value_combinations(S.subject.knownSynW_r,muscleNames_r,[S.subject.TrackSynW_NSyn_r]);
            k = 1;
            for i = 1:NMuscle/2
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
            J = J + W.TrackSynW * Jtemp;
        case 'onlyLeft'
            knownSynW_l_all = unpack_name_value_combinations(S.subject.knownSynW_l,muscleNames_l,[S.subject.TrackSynW_NSyn_l]);
            k = 1;
            for i = 1:NMuscle/2
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
            J = J + W.TrackSynW * Jtemp;
        case 'RightLeft'
            % right
            knownSynW_idx = [];
            knownSynW = [];
            knownSynW_r_all = unpack_name_value_combinations(S.subject.knownSynW_r,muscleNames_r,[S.subject.TrackSynW_NSyn_r]);
            k = 1;
            for i = 1:NMuscle/2
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
            J = J + W.TrackSynW * Jtemp;
            knownSynW_idx = [];
            knownSynW = [];
            % left
            knownSynW_l_all = unpack_name_value_combinations(S.subject.knownSynW_l,muscleNames_l,[S.subject.TrackSynW_NSyn_l]);
            k = 1;
            for i = 1:NMuscle/2
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
            J = J + W.TrackSynW * Jtemp;
    end
end