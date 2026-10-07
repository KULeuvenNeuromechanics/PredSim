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
%           - SubjectMass: unloaded model mass in kg (default: model/settings)
%           - LegLength: normalization length in m (default: 0.85 for Falisse)


%% load benchmarking settings

load(fullfile(benchmarking_folder,'benchmark_settings.mat'),...
    'S','osim_path','S_benchmark');

%% input parser
p = inputParser;
addParameter(p, 'BoolPlot', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'dofs', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'studies', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'OverwriteData', false, @(x) islogical(x) && isscalar(x));
positive_scalar = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0;
addParameter(p, 'SubjectMass', [], @(x) isempty(x) || positive_scalar(x));
addParameter(p, 'LegLength', 0.85, positive_scalar);
parse(p, varargin{:});
BoolPlot = logical(p.Results.BoolPlot);
dofs_plot = cellstr(p.Results.dofs);
studies_plot = cellstr(p.Results.studies);

% defaults settings
if isempty(dofs_plot)
    dofs_plot = {'ankle_angle_r','knee_angle_r','hip_flexion_r'};
end

if isempty(studies_plot)
    studies_plot = S_benchmark.studies;
end

if BoolPlot
    msim = p.Results.SubjectMass;
    if isempty(msim)
        if isfield(S,'subject') && isfield(S.subject,'mass') && ~isempty(S.subject.mass)
            msim = S.subject.mass;
        elseif isfile(osim_path)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'VariousFunctions'));
            msim = getModelMass(osim_path);
        else
            error('PredSim:MissingBenchmarkMass',...
                'The original model is unavailable. Supply SubjectMass (unloaded mass in kg) for plotting.');
        end
    end
    assert(positive_scalar(msim),'PredSim:InvalidBenchmarkMass',...
        'SubjectMass must be a positive finite scalar in kg.');
    Lsim = p.Results.LegLength;
end


%% Download benchmarking data if desired

bool_overwrite = p.Results.OverwriteData;
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
        loaded = load(filename,'R');
        R = loaded.R;
        benchmark = [];
        % find data with same ID
        if isfield(R.S.misc,'benchmark_id') && ~isempty(R.S.misc.benchmark_id)
            % find id in exp datalist
            id_exp = match_benchmark_id(R.S.misc.benchmark_id,identifierList);
            if numel(id_exp) == 1
                benchmark = data{id_exp};
            elseif numel(id_exp)>1
                disp(['import warning ! I found the id ' R.S.misc.benchmark_id,...
                    ' ' num2str(numel(id_exp)) ' times in the experimental dataset' ]);
            elseif ~startsWith(R.S.misc.benchmark_id,'gait_speeds_')
                warning('PredSim:MissingBenchmark','No experimental match for %s in %s.',...
                    R.S.misc.benchmark_id,filename);
            end
        end
        % Clear stale attachments if the current ID has no unique match.
        save(filename, 'benchmark', '-append');
    end
end


%% Plotting

%% Default plots
if BoolPlot
    


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
                sim_res_folder = filename;
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
    if isempty(table_all)
        warning('PredSim:NoBenchmarkResults','No completed results to plot in %s.',benchmarking_folder);
        return
    end

    % plot all data on one graph
    h_figallp = figure('Name','All data','Color',[1 1 1]);
    t_layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % get study ids and assign colors
    study_ids = unique(table_all.id_study);
    n_studies = length(study_ids);
    cols_sel = lines(n_studies);
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

    % Use one result-loading and plotting path for all selected studies.
    for istudy = 1:numel(studies_plot)
        study = studies_plot{istudy};
        files = dir(fullfile(benchmarking_folder,study,'**','*.mat'));
        Dat = struct('R',{},'benchmark',{},'model_info',{});
        for k = 1:numel(files)
            filename = fullfile(files(k).folder,files(k).name);
            [R,benchmark,model_info] = load_sim_file(filename);
            if isempty(R) || isempty(benchmark)
                continue
            end
            if ~all(isfield(R,{'kinematics','kinetics','ground_reaction','spatiotemp','metabolics','time'}))
                warning('PredSim:IncompleteBenchmarkResult','Skipping incomplete result %s.',filename);
                continue
            end
            Dat(end+1) = struct('R',R,'benchmark',benchmark,'model_info',model_info);
        end
        if ~isempty(Dat)
            speeds = arrayfun(@(d) d.R.S.misc.forward_velocity,Dat);
            [~,order] = sort(speeds);
            default_plot_benchmarking(Dat(order),dofs_plot,study,msim,Lsim,true);
        end
    end
end
end
