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
Colours = linspecer(nsim);
g = 9.81;

for isim = 1:nsim
    [~, folder_file, ~] = fileparts(Dat(isim).R.S.misc.save_folder);
    headers{isim} = folder_file;
end

if ~isempty(Dat(1).benchmark.ik)
    figure('Name',[study_name ': kinematics'],'Color',[1 1 1]);
    t = tiledlayout(2,nsubplot_dofs,'TileSpacing','compact','Padding','compact');

    for idof = 1:length(dofs_plot)
        for isim = 1:length(Dat)

            Cs = Colours(isim,:);

            % plot experimental data
            tile_number = tilenum(t, 1, idof);
            nexttile(tile_number);
            plot(Dat(isim).benchmark.ik.(dofs_plot{idof}(1:end-2)),...
                'Color',Cs); hold on;
            % plot simulation data
            tile_number = tilenum(t, 2, idof);
            nexttile(tile_number);
            icol = strcmp(dofs_plot{idof},Dat(isim).R.colheaders.coordinates);
            dsel = Dat(isim).R.kinematics.Qs(:,icol);
            dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
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
if ~isempty(Dat(1).benchmark.id)
    figure('Name',[study_name ': kinetics'],'Color',[1 1 1]);
    t = tiledlayout(2,nsubplot_dofs,'TileSpacing','compact','Padding','compact');
    for idof = 1:length(dofs_plot)
        for isim = 1:length(Dat)
            % select color
            Cs = Colours(isim,:);

            % plot experimental data
            nexttile(idof);
            id_exp = Dat(isim).benchmark.id.(dofs_plot{idof}(1:end-2));
            id_exp = id_exp.*(msim*g*Lsim); % scale to subject
            plot(id_exp,'Color',Cs); hold on;


            % plot simulation data
            nexttile(idof+nsubplot_dofs);
            icol = strcmp(dofs_plot{idof},Dat(isim).R.colheaders.coordinates);
            dsel = Dat(isim).R.kinetics.T_ID(:,icol);
            dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
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
            ylabel({'simulation','joint angle [Nm]'});
        end
    end
end

% Plot ground reaction forces
if ~isempty(Dat(1).benchmark.grf)
    figure('Name',[study_name ': grf'],'Color',[1 1 1]);
    t = tiledlayout(2,nsubplot_dofs,'TileSpacing','compact','Padding','compact');
    grf_headers = {'Fx','Fy','Fz'};
    for coord = 1:3
        for isim = 1:nsim
            % select color
            Cs = Colours(isim,:);
            nexttile(coord);
            % experimental grf
            Fsel = Dat(isim).benchmark.grf.(grf_headers{coord});
            Fsel = Fsel*msim*g;
            plot(Fsel,'Color',Cs);hold on;

            % simulated grf
            nexttile(coord+3);
            if bool_rot_grf
                dsel = Dat(isim).R.ground_reaction.GRF_r_rot(:,coord);
            else
                dsel = Dat(isim).R.ground_reaction.GRF_r(:,coord);
            end
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
if ~isempty(Dat(1).benchmark.stride_frequency)
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
    t = Dat(isim).R.time.mesh_GC;
    dt = t(end)-t(1);
    Pmetab = Dat(isim).R.metabolics.Bhargava2004.Edot_gait;
    metab_work  = trapz(t(1:end-1)',Pmetab);
    P_mean = sum(metab_work)./dt;
    sim_metab(isim) = P_mean;
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