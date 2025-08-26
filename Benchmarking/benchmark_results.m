function [] = benchmark_results(benchmarking_folder)
%benchmark_results Compares simulations to experiments
%   input arguments:
%       (1) benchmarking_folder = folder with benchmarking results. this is
%       the S_benchmark.out_folder folder when running benchmark_predsim

% Note:
% this is work in progress. the plan is to work out a minimal example with
% benchmarking the study of Van Der Zee with variations in gait speed. I
% will discuss this example at a PredSim meeting

%% load benchmarking settings

load(fullfile(benchmarking_folder,'benchmark_settings.mat'),'S','osim_path','S_benchmark');


%% load experimental data

% DOI 10.5281/zenodo.14524619
data_exp_path = fullfile(benchmarking_folder,'exp_data');
if ~isfolder(data_exp_path)
    mkdir(data_exp_path);
end

% general data
% check if we have to download the data
if ~isfolder(fullfile(benchmarking_folder,'exp_data','benchmarking_folder'))
    % download data as in Afschrift et al. 2025
    LinkZipFile = ['https://zenodo.org/records/14524620/files/',...
        'PredSim_gait_conditions-master_19_12_2024.zip?download=1'];
    websave(fullfile(data_exp_path,'githubrepo_data.zip'),LinkZipFile);
    unzip(fullfile(data_exp_path,'githubrepo_data.zip'),...
        fullfile(data_exp_path,'github_repo_temp'));
    delete(fullfile(data_exp_path,'githubrepo_data.zip'));
    %       change data format a bit here, we don't need everything
    copyfile(fullfile(data_exp_path,'github_repo_temp',...
        'PredSim_gait_conditions-master','ExperimentalData'),...
        fullfile(data_exp_path,'general_exp_data'))
    rmdir(fullfile(data_exp_path,'github_repo_temp'),'s');
end

% Van Der Zee kinematic data
if isfield(S_benchmark,'studies') && ~isempty(S_benchmark.studies)
    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        if ~isfolder(fullfile(benchmarking_folder,'exp_data','vanderzee2022'))
            % url to datafile
            LinkZipFile = ['https://zenodo.org/records/14524620/files/',...
                'VanDerZee_Averages_Addbiomech.zip?download=1'];
            % download data
            websave(fullfile(data_exp_path,'vanderzee2022.zip'),LinkZipFile);
            % unzip and delete zip file
            unzip(fullfile(data_exp_path,'vanderzee2022.zip'),...
                fullfile(data_exp_path,'vanderzee2022'));
            delete(fullfile(data_exp_path,'vanderzee2022.zip'));
        end
    end
end

%% Attach experimental data to simulation results

if isfield(S_benchmark,'studies') && ~isempty(S_benchmark.studies)
    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        % loop over simulations
        n_sim = length(S_benchmark.vanderzee.names);
        for isim = 1:n_sim
            % path to experimental data
            speed = S_benchmark.vanderzee.gait_speeds(isim);
            grf_datafile = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_GRF.csv']);
            grf_datafile_std = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_GRF_std.csv']);
            id_datafile = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_ID.csv']);
            id_datafile_std = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_ID_std.csv']);
            ik_datafile = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_IK.csv']);
            ik_datafile_std = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_IK_std.csv']);

            % read GRF, ID and IK data
            grf_dat = readtable(grf_datafile);
            grf_dat_std = readtable(grf_datafile_std);
            id_dat = readtable(id_datafile);
            id_dat_std = readtable(id_datafile_std);
            ik_dat= readtable(ik_datafile);
            ik_dat_std = readtable(ik_datafile_std);

            % attach to simulation results
            sim_res_file = fullfile()




        end
    end
end








outputArg1 = inputArg1;
outputArg2 = inputArg2;
end