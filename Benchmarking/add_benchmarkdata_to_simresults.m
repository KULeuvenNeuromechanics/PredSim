function [] = add_benchmarkdata_to_simresults(benchmarking_folder,varargin)
%add_benchmarkdata_to_simresults Compares simulations to experiments
%   input arguments:
%       (1) benchmarking_folder = folder with benchmarking results. this is
%       the S_benchmark.out_folder folder when running benchmark_predsim
%       (2) optional input arguments:
%           - 'BoolPlot', true, makes some default plots with benchmarking results
%           - 'dofs', {'ankle_angle_r','knee_angle_r','hip_flexion_r'}:
%           plots these dofs
%           - studies: studies you want to include in plotting


%% load benchmarking settings

load(fullfile(benchmarking_folder,'benchmark_settings.mat'),...
    'S','osim_path','S_benchmark');

%% input parser
p = inputParser;
addParameter(p, 'BoolPlot', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'dofs', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'studies', {}, @(x) iscellstr(x) || isstring(x));
parse(p, varargin{:});
BoolPlot = logical(p.Results.BoolPlot);
dofs_plot = cellstr(p.Results.dofs);
studies_plot = cellstr(p.Results.studies);

% defaults settings
if isempty(dofs_plot)
    dofs_plot = {'ankle_angle_r','knee_angle_r','hip_flexion_r'};
end
nsubplot_dofs = length(dofs_plot);

if isempty(studies_plot)
    studies_plot = S_benchmark.studies;
end


%% Download benchmarking data if desired

bool_overwrite = true; % overwrite existing data (if true it downloads it every time)
[datafolder] = download_benchmarkdata(bool_overwrite);

% add data processing functions to matlab path
addpath(fullfile(datafolder,'functions'));

% read all experimental data
[data,studyList, identifierList, slopeList, speedList] = get_all_benchmarkdata();


%% Read all simulation results and add experimental data

% find all .mat files in benchmarking_folder
mat_files = dir(fullfile(benchmarking_folder, '**', '*.mat'));

for i =1:length(mat_files)
    % get current filename
    filename = fullfile(mat_files(i).folder, mat_files(i).name);
    % check if this is a simultion results file
    vars = who('-file', filename);
    if ismember('R', vars)
        % load the .mat file
        load(filename,'R');
        % find data with same ID
        if isfield(R.S.misc,'benchmark_id') && ~isempty(R.S.misc.benchmark_id)
            % find id in exp datalist
            try
                id_exp = strcmp(R.S.misc.benchmark_id,identifierList);
            catch
                disp('error'); % stupid bug we have to solve in inputs
            end
            if ~isempty(id_exp) && sum(id_exp) == 1
                benchmark = data{id_exp};
                save(filename, 'benchmark', "-append");
            elseif sum(id_exp)>1
                disp(['import warning ! I found the id ' R.S.misc.benchmark_id,...
                    ' ' num2str(sum(id_exp)) ' times in the experimental dataset' ]);
            end
        end
    end
end


%% Plotting

%% Default plots
if BoolPlot
    
    % hard coded model properties for now
    msim = 62; % should read this from the model in the future
    Height_Sim = 1.70; % I should read this from the model in the future
    prop_leg_length = 0.5;
    Lsim = Height_Sim.*prop_leg_length;
    g = 9.81;

    % allowed to adapt figure defaults ?
    set(0,'defaultLineLineWidth',1.6);
    set(0,'defaultAxesLineWidth',1.2);
    set(0,'defaultAxesFontSize',12);
    set(0,'defaultFigureColor',[1 1 1]);
    set(0,'defaultLineMarkerSize',4);
    set(0,'defaultAxesBox','off');

    % the idea is here to make a function that works for all studies
    % to do so I need to provide some input arguments such as
    %  1. how to loop over different simulations / gait conditions ?
    %  2. study name as input
    %  3. angles in rad to deg ?
    %  4. add options to handle ik, id and grf input as imported mot files
    %     or matlab tables
    %   => ToDo: adapt custom code for vanderzee2022 to this default
    %   function

    %% Plot figure with stride frequency and metabolic power for all datapoints
    % the approach is quite simple here, we load all .mat files with
    % simulations results and plot the datapoints.

    % also store table with all results
    headers_table = {'sim_stride_frequency','sim_metabolic_power','speed',...
        'slope','id_study','exp_stride_frequency','exp_metabolic_power'};
    data_table = nan(1000, length(headers_table));
    ct_sim = 1;
    study_id_header = {};
    for istudy = 1:length(studies_plot)
        % find all mat files in this folder
        study_name = studies_plot{istudy};
        study_id_header{istudy} = study_name;
        mat_files = dir(fullfile(benchmarking_folder,study_name, '**', '*.mat'));
        for ifile = 1:length(mat_files)
            % get current filename
            filename = fullfile(mat_files(ifile).folder, mat_files(ifile).name);
            % check if this is a simultion results file
            vars = who('-file', filename);
            if ismember('R', vars)
                sim_res_folder = mat_files(ifile).folder;
                id_study = istudy;
                [data_table, ct_sim] = add_benchmark_to_table(ct_sim,...
                    sim_res_folder, headers_table, data_table, id_study,...
                    Lsim, msim);
            end
        end
    end

    data_table(ct_sim:end,:) = [];
    table_all = array2table(data_table,...
        'VariableNames',headers_table);

    % plot all data on one graph
    h_figallp = figure('Name','All data','Color',[1 1 1]);
    t_layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % get study ids and assign colors
    study_ids = unique(table_all.id_study);
    n_studies = length(study_ids);
    cols_sel = linspecer(n_studies);
    mk = 4;

    % plot stride frequency
    nexttile(1);
    plot([min(table_all.exp_stride_frequency) max(table_all.exp_stride_frequency)],...
        [min(table_all.exp_stride_frequency) max(table_all.exp_stride_frequency)],...
        '--','Color',[0 0 0],'LineWidth',1.3); hold on;
    for istudy = 1:n_studies
        Cs = cols_sel(istudy,:);
        rows_sel = table_all.id_study == study_ids(istudy);
        plot(table_all.exp_stride_frequency(rows_sel),...
            table_all.sim_stride_frequency(rows_sel),...
            'ok','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk)
    end
    set(gca,'box','off')
    set(gca,'FontSize',10);
    xlabel('measured stride frequency');
    ylabel('simulated stride frequency');

    % plot metabolic power
    nexttile(2);
    plot([min(table_all.exp_metabolic_power) max(table_all.exp_metabolic_power)],...
        [min(table_all.exp_metabolic_power) max(table_all.exp_metabolic_power)],...
        '--','Color',[0 0 0],'LineWidth',1.3); hold on;
    legs = [];
    for istudy = 1:n_studies
        Cs = cols_sel(istudy,:);
        rows_sel = table_all.id_study == study_ids(istudy);
        legs(istudy) = plot(table_all.exp_metabolic_power(rows_sel),...
            table_all.sim_metabolic_power(rows_sel),...
            'ok','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    end
    set(gca,'box','off')
    set(gca,'FontSize',10);
    xlabel('measured metab. power');
    ylabel('simulated metab. power');


    hL = legend(legs,study_id_header(study_ids),'NumColumns',5,'Box','off', ...
        'FontSize',10,'Interpreter','none');
    hL.Layout.Tile = 'North';

    %% Van Der Zee 2022: custom function here
    % mainly useful as an example on how to compare experiments and
    % simulations using custom code
    ct_sim = 1;
    if any(strcmp(studies_plot,'vanderzee2022'))
        % default function to plot van der zee results
        % to do: add this to a function
        %
        speeds = S_benchmark.vanderzee.gait_speeds;
        [speeds_sort,isort] = sort(speeds);
        for i= 1:length(speeds_sort)
            WalkSpeed_legend{i} = [num2str(speeds_sort(i)) 'ms^{-1}'];
        end
        % Colors_Speeds = [146,133,101;... % 0.7
        %     126,108,62;... % 0.9
        %     238,202,102;...% 1.1
        %     125,162,197;...% 1.6
        %     238,202,102;...% 1.8
        %     125,162,197;...% 2.0
        %     238,202,102]./255 % 1.4; % these colors are pretty bad, also in paper ?
        % Colors_Speeds = flipud(Colors_Speeds);
        Colors_Speeds = copper(length(speeds)+2);
        % load data of all simulation conditions first
        clear Dat
        for ispeed_sort = 1:length(speeds_sort)
            ispeed = isort(ispeed_sort);
            sim_res_folder = fullfile(benchmarking_folder,'vanderzee2022',...
                S_benchmark.vanderzee.names{ispeed});
            mat_files = dir(fullfile(sim_res_folder,'*.mat'));
            if length(mat_files) ~= 1
                disp('warning mutiple mat files in folder')
                disp(sim_res_folder);
                disp([' assumes that file ' mat_files(1).name, ...
                    'contains the simulation results'])
            end
            sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
            load(sim_res_file,'R','benchmark','model_info');
            try
                Dat(ispeed_sort).benchmark = benchmark;
            catch
                disp('error');
            end
            Dat(ispeed_sort).R = R;
            Dat(ispeed_sort).model_info = model_info;
        end

        % plot kinematics
        if ~isempty(Dat(1).benchmark.ik)
            figure('Name','vanderzee: kinematics','Color',[1 1 1]);

            for ispeed = 1: length(speeds_sort)
                for idof = 1:length(dofs_plot)
                    % plot experimental data
                    subplot(2, nsubplot_dofs, idof)
                    plot(Dat(ispeed).benchmark.ik.(dofs_plot{idof})*180/pi,...
                        'Color',Colors_Speeds(ispeed,:)); hold on;
                    % plot simulation data
                    subplot(2, nsubplot_dofs, idof+nsubplot_dofs)
                    icol = strcmp(dofs_plot{idof},Dat(ispeed).R.colheaders.coordinates);
                    dsel = Dat(ispeed).R.kinematics.Qs(:,icol);
                    dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
                    l(ispeed) = plot(dsel_int,'Color',Colors_Speeds(ispeed,:)); hold on;
                end
            end
            legend(l,WalkSpeed_legend,'Interpreter','tex');
            legend boxoff;
            for isubpl =1:length(dofs_plot)*2
                subplot(2, nsubplot_dofs,isubpl)
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

        % plot joint moments
        if ~isempty(Dat(1).benchmark.id)
            figure('Name','vanderzee: kinetics','Color',[1 1 1]);

            for ispeed = 1: length(speeds_sort)
                for idof = 1:length(dofs_plot)
                    % plot experimental data
                    subplot(2, nsubplot_dofs, idof)
                    id_exp = Dat(ispeed).benchmark.id.(dofs_plot{idof});
                    id_exp = id_exp.*(msim*g*Lsim); % scale to subject
                    plot(id_exp,'Color',Colors_Speeds(ispeed,:)); hold on;


                    % plot simulation data
                    subplot(2, nsubplot_dofs, idof+nsubplot_dofs)
                    icol = strcmp(dofs_plot{idof},Dat(ispeed).R.colheaders.coordinates);
                    dsel = Dat(ispeed).R.kinetics.T_ID(:,icol);
                    dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
                    l(ispeed) = plot(dsel_int,'Color',Colors_Speeds(ispeed,:)); hold on;
                end
            end
            legend(l,WalkSpeed_legend,'Interpreter','tex');
            legend boxoff;
            for isubpl =1:length(dofs_plot)*2
                subplot(2, nsubplot_dofs,isubpl)
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


        % plot ground reaction forces
        if ~isempty(Dat(1).benchmark.grf_r)
            figure('Name','vanderzee: grf','Color',[1 1 1]);
            for ispeed = 1: length(speeds_sort)
                for coord = 1:3
                    subplot(2,3,coord)
                    % experimental grf
                    Fsel = Dat(ispeed).benchmark.grf_r(:,coord);
                    Fsel = Fsel.*msim*g;
                    plot(Fsel,'Color',Colors_Speeds(ispeed,:));hold on;

                    % simulated grf
                    subplot(2,3,coord+3)
                    dsel = Dat(ispeed).R.ground_reaction.GRF_r(:,coord);
                    dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
                    l(ispeed) =plot(dsel_int,'Color',Colors_Speeds(ispeed,:));hold on;


                end
            end
            legend(l,WalkSpeed_legend,'Interpreter','tex');
            legend boxoff;
            title_grf = {'grfx','grfy','grfz'};
            for isubpl =1:6
                subplot(2, 3,isubpl)
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
            figure('Name','vanderzee: stride frequency','Color',[1 1 1]);

            % experimental stride frequency
            exp_freq = nan(length(speeds),1);
            sim_freq = nan(length(speeds),1);
            for ispeed = 1:length(speeds)
                exp_freq(ispeed) = Dat(ispeed).benchmark.stride_frequency .* (sqrt(g/Lsim));
                sim_freq(ispeed) = Dat(ispeed).R.spatiotemp.stride_freq;
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
        if ~isempty(Dat(1).benchmark.Pmetab_mean)
            figure('Name','vanderzee: metabolic power','Color',[1 1 1]);

            % experimental stride frequency
            exp_metab = nan(length(speeds),1);
            sim_metab = nan(length(speeds),1);
            for ispeed = 1:length(speeds)
                % measured metabolic power
                exp_metab(ispeed) = Dat(ispeed).benchmark.Pmetab_mean ...
                    .* (m*sqrt(Lsim)*g^1.5);

                % simulated metabolic power
                t = Dat(ispeed).R.time.mesh_GC;
                dt = t(end)-t(1);
                Pmetab = Dat(ispeed).R.metabolics.Bhargava2004.Edot_gait;
                metab_work  = trapz(t(1:end-1)',Pmetab);
                P_mean = sum(metab_work)./dt;
                sim_metab(ispeed) = P_mean;
            end
            Cs = [0 0 0];
            mk = 4;
            plot([min(exp_metab) max(exp_metab)], [min(exp_metab) max(exp_metab)],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
            plot(exp_metab,sim_metab,'ok','Color',Cs,'MarkerFaceColor',Cs,...
                'MarkerSize',mk);
            set(gca,'box','off');
            set(gca,'FontSize',10);
            xlabel('measured metab. power');
            ylabel('simulated metab. power');
        end
        clear Dat
    end

    %% Koelewijn2019: default functions
    % example how to compare experiments and simulations using a
    % generic/default function (this saves you quite a bit of coding, but ugly)
    if any(strcmp(studies_plot,'koelewijn2019'))
        % read all the data
        Dat = [];
        for isim =1:6%length(S_benchmark.koelewijn.names)
            for idof = 1:length(dofs_plot)
                % load sim file
                sim_res_folder = fullfile(benchmarking_folder,'koelewijn2019',...
                    S_benchmark.koelewijn.names{isim});
                mat_files = dir(fullfile(sim_res_folder,'*.mat'));
                if ~isempty(mat_files)
                    if length(mat_files) > 1
                        disp('warning mutiple mat files in folder')
                        disp(sim_res_folder);
                        disp([' assumes that file ' mat_files(1).name, ...
                            'contains the simulation results'])
                    end
                    sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
                    load(sim_res_file,'R','benchmark','model_info');
                    Dat(isim).benchmark = benchmark;
                    Dat(isim).R = R;
                    Dat(isim).model_info = model_info;
                end
            end
        end

        % default plot function for koelewijn
        if ~isempty(Dat)
            bool_rot_grf = true;
            default_plot_benchmarking(Dat, dofs_plot, 'Koelewijn',msim,Lsim,bool_rot_grf);            
        end
        clear Dat
    end

    %% Browning2008: default functions
    % example how to compare experiments and simulations using a
    % generic/default function (this saves you quite a bit of coding, but ugly)
    if any(strcmp(studies_plot,'browning2008'))
        % read all the data
        nsim = length(S_benchmark.browning.names);
        for isim =1:nsim
            for idof = 1:length(dofs_plot)
                % load sim file
                sim_res_folder = fullfile(benchmarking_folder,'browning2008',...
                    S_benchmark.browning.names{isim});
                mat_files = dir(fullfile(sim_res_folder,'*.mat'));
                if length(mat_files) ~= 1
                    disp('warning mutiple mat files in folder')
                    disp(sim_res_folder);
                    disp([' assumes that file ' mat_files(1).name, ...
                        'contains the simulation results'])
                end
                sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
                load(sim_res_file,'R','benchmark','model_info');
                Dat(isim).benchmark = benchmark;
                Dat(isim).R = R;
                Dat(isim).model_info = model_info;
            end
        end
        % default plot for Koelewijn
        bool_rot_grf = false;
        default_plot_benchmarking(Dat, dofs_plot, 'Browning',msim,Lsim,bool_rot_grf);
        clear Dat
    end
end

end












