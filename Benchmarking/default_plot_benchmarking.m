function [] = default_plot_benchmarking(Dat, dofs_plot, study_name, ...
    msim, Lsim, bool_rot_grf)
%Default plot function
% input argument:
%   Dat =  matlab structure with all sim results. This should contain
%       - field R: (simulation results)
%       - field benchmarking:
%   dofs_plot = data structure with coord names to plot

nsim = length(Dat);
nsubplot_dofs = length(dofs_plot);
if nsim == 0
    return
end
Colours = lines(nsim);
g = 9.81;
for isim = 1:nsim
    for field = {'ik','id','grf_r','stride_frequency','Pmetab_mean'}
        if ~isfield(Dat(isim).benchmark,field{1})
            Dat(isim).benchmark.(field{1}) = [];
        end
    end
end

for isim = 1:nsim
    [~, folder_file, ~] = fileparts(Dat(isim).R.S.misc.save_folder);
    headers{isim} = folder_file;
end

if any(arrayfun(@(d) ~isempty(d.benchmark.ik),Dat))
    figure('Name',[study_name ': kinematics'],'Color',[1 1 1]);
    t = tiledlayout(2,nsubplot_dofs,'TileSpacing','compact','Padding','compact');

    for idof = 1:length(dofs_plot)
        for isim = 1:length(Dat)

            Cs = Colours(isim,:);

            % plot experimental data
            tile_number = idof;
            nexttile(tile_number);
            angles = benchmark_coordinate(Dat(isim).benchmark.ik,dofs_plot{idof});
            if isfield(Dat(isim).benchmark,'study') && ...
                    strcmp(Dat(isim).benchmark.study,'vanderzee2022') && ...
                    ~ismember(dofs_plot{idof},{'pelvis_tx','pelvis_ty','pelvis_tz'})
                angles = rad2deg(angles);
            end
            plot(angles,...
                'Color',Cs); hold on;
            % plot simulation data
            tile_number = nsubplot_dofs + idof;
            nexttile(tile_number);
            dsel_int = simulation_coordinate(Dat(isim).R.kinematics.Qs,...
                Dat(isim).R.colheaders.coordinates,dofs_plot{idof});
            legs(isim) = plot(dsel_int,'Color',Cs); hold on;
        end
    end


    hL = legend(legs,headers,'NumColumns',3,'Box','off', ...
        'FontSize',12,'Interpreter','none');
    clear legs
    % Move the legend to the right side of the figure
    hL.Layout.Tile = 'North';

    for isubpl =1:length(dofs_plot)*2
        nexttile(isubpl);
        set(gca,'box','off');
        set(gca,'FontSize',10);
        if isubpl<=length(dofs_plot)
            title(dofs_plot{isubpl},'interpreter','none');
        else
            xlabel('% gait cycle');
        end
        if isubpl == 1
            ylabel({'experiment','joint angle [deg]'});
        elseif isubpl == (length(dofs_plot)+1)
            ylabel({'simulation','joint angle [deg]'});
        end
    end
end

% Plot joint moments
if any(arrayfun(@(d) ~isempty(d.benchmark.id),Dat))
    figure('Name',[study_name ': kinetics'],'Color',[1 1 1]);
    t = tiledlayout(2,nsubplot_dofs,'TileSpacing','compact','Padding','compact');
    for idof = 1:length(dofs_plot)
        for isim = 1:length(Dat)
            % select color
            Cs = Colours(isim,:);

            % plot experimental data
            nexttile(idof);
            id_exp = benchmark_coordinate(Dat(isim).benchmark.id,dofs_plot{idof});
            id_exp = id_exp.*(msim*g*Lsim); % scale to subject
            plot(id_exp,'Color',Cs); hold on;


            % plot simulation data
            nexttile(idof+nsubplot_dofs);
            dsel_int = simulation_coordinate(Dat(isim).R.kinetics.T_ID,...
                Dat(isim).R.colheaders.coordinates,dofs_plot{idof});
            legs(isim) = plot(dsel_int,'Color',Cs); hold on;
        end
    end
    hL = legend(legs,headers,'NumColumns',3,'Box','off', ...
        'FontSize',12,'Interpreter','none');
    clear legs
    hL.Layout.Tile = 'North';

    for isubpl =1:length(dofs_plot)*2
        nexttile(isubpl);
        set(gca,'box','off');
        set(gca,'FontSize',10);
        if isubpl<=length(dofs_plot)
            title(dofs_plot{isubpl},'interpreter','none');
        else
            xlabel('% gait cycle');
        end
        if isubpl == 1
            ylabel({'experiment','joint moment [Nm]'});
        elseif isubpl == (length(dofs_plot)+1)
            ylabel({'simulation','joint moment [Nm]'});
        end
    end
end

% Plot ground reaction forces
if any(arrayfun(@(d) ~isempty(d.benchmark.grf_r),Dat))
    figure('Name',[study_name ': grf'],'Color',[1 1 1]);
    t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    grf_headers = {'Fx','Fy','Fz'};
    for coord = 1:3
        for isim = 1:nsim
            % select color
            Cs = Colours(isim,:);
            nexttile(coord);
            % experimental grf
            if isempty(Dat(isim).benchmark.grf_r)
                Fsel = NaN;
            elseif isnumeric(Dat(isim).benchmark.grf_r)
                Fsel = Dat(isim).benchmark.grf_r(:,coord);
            else
                Fsel = Dat(isim).benchmark.grf_r.(grf_headers{coord});
            end
            Fsel = Fsel*msim*g;
            plot(Fsel,'Color',Cs);hold on;

            % simulated grf
            nexttile(coord+3);
            forces = Dat(isim).R.ground_reaction.GRF_r;
            if bool_rot_grf
                if isfield(Dat(isim).model_info,'slope')
                    fi = atan(Dat(isim).model_info.slope);
                    Rotm = benchmark_rotation_z(fi);
                    forces = forces*Rotm(1:3,1:3)';
                else
                    disp(['warning could not rotate forces for slope walking']);
                end
            end
            dsel = forces(:,coord);
            dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
            legs(isim) =plot(dsel_int,'Color',Cs);hold on;
        end
    end
    hL = legend(legs,headers,'NumColumns',3,'Box','off', ...
        'FontSize',12,'Interpreter','none');
    clear legs
    hL.Layout.Tile = 'North';
    title_grf = {'GRFx','GRFy','GRFz'};
    for isubpl =1:6
        nexttile(isubpl)
        set(gca,'box','off');
        set(gca,'FontSize',10);
        if isubpl<=3
            title(title_grf{isubpl},'interpreter','none');
        else
            xlabel('% gait cycle');
        end
        if isubpl == 1
            ylabel({'experiment','force [N]'});
        elseif isubpl == 4
            ylabel({'simulation','force [N]'});
        end
    end
end

% plot stride frequency
if any(arrayfun(@(d) ~isempty(d.benchmark.stride_frequency),Dat))
    figure('Name',[study_name ': stride frequency'],'Color',[1 1 1]);

    % experimental stride frequency
    exp_freq = nan(nsim,1);
    sim_freq = nan(nsim,1);
    for isim = 1:nsim
        if ~isempty(Dat(isim).benchmark.stride_frequency)
            exp_freq(isim) = Dat(isim).benchmark.stride_frequency .* (sqrt(g/Lsim));
        end
        sim_freq(isim) = Dat(isim).R.spatiotemp.stride_freq;
    end
    Cs = [0 0 0];
    mk = 4;
    plot([min(exp_freq) max(exp_freq)], [min(exp_freq) max(exp_freq)],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    plot(exp_freq,sim_freq,'ok','Color',Cs,'MarkerFaceColor',Cs,...
        'MarkerSize',mk);
    set(gca,'box','off');
    set(gca,'FontSize',10);
    xlabel('measured stride frequency');
    ylabel('simulated stride frequency');
end



% plot metabolic power
% if ~isempty(Dat(1).benchmark.Pmetab_mean)
figure('Name',[study_name ': metabolic power'],'Color',[1 1 1]);

% experimental stride frequency
exp_metab = nan(nsim,1);
sim_metab = nan(nsim,1);
for isim = 1:nsim
    % measured metabolic power
    if ~isempty(Dat(isim).benchmark.Pmetab_mean)
        exp_metab(isim) = Dat(isim).benchmark.Pmetab_mean ...
            .* (msim*sqrt(Lsim)*g^1.5);
    end

    % simulated metabolic power
    sim_metab(isim) = benchmark_mean_metabolic_power(Dat(isim).R);
end
Cs = [0 0 0];
mk = 4;
plot([min(exp_metab) max(exp_metab)], [min(exp_metab) max(exp_metab)],...
    '--','Color',[0 0 0],'LineWidth',1.3); hold on;
plot(exp_metab,sim_metab,'ok','Color',Cs,'MarkerFaceColor',Cs,...
    'MarkerSize',mk);
set(gca,'box','off');
set(gca,'FontSize',10);
xlabel('measured metab. power');
ylabel('simulated metab. power');

% end




end


function values = simulation_coordinate(data,names,name)
% Keep a missing model coordinate blank without aborting the other figures.
column = strcmp(name,names);
if ~any(column)
    warning('PredSim:MissingPlotCoordinate','No simulated coordinate %s; leaving its curve blank.',name);
    values = nan(1,100);
    return
end
values = interp1(1:size(data,1),data(:,column),linspace(1,size(data,1),100));
end

function values = benchmark_coordinate(data,name)
% Some studies store bilateral names, others omit the side suffix.
if isempty(data)
    values = NaN;
    return
end
if istable(data)
    names = data.Properties.VariableNames;
else
    names = fieldnames(data);
end
if ~ismember(name,names)
    name = regexprep(name,'_[rl]$','');
end
if ismember(name,names)
    values = data.(name);
else
    values = NaN;
end
end
