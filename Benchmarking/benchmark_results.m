function [] = benchmark_results(benchmarking_folder,varargin)
%benchmark_results Compares simulations to experiments
%   input arguments:
%       (1) benchmarking_folder = folder with benchmarking results. this is
%       the S_benchmark.out_folder folder when running benchmark_predsim
%       (2) optional input arguments:
%           - 'BoolPlot', true, makes some default plots with benchmarking results


% Note:
% this is work in progress. the plan is to work out a minimal example with
% benchmarking the study of Van Der Zee with variations in gait speed. I
% will discuss this example at a PredSim meeting

% ToDO:
%
%       Change datastorage to osf instead of zenodo:
% I decided to move to osf framework to store data. The main reason is that
% the store of data of vanderzee is now a bit weird on zenode (I did the
% normalisation wrong, forgot to divided by g and so. I don't want this
% here because it is ugly and unclear)


%% input parser
p = inputParser;
addParameter(p, 'BoolPlot', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'dofs', {}, @(x) iscellstr(x) || isstring(x));
parse(p, varargin{:});
BoolPlot = logical(p.Results.BoolPlot);
dofs_plot = cellstr(p.Results.dofs);

% defaults settings
if isempty(dofs_plot)
    dofs_plot = {'ankle_angle_r','knee_angle_r','hip_flexion_r'};
end
nsubplot_dofs = length(dofs_plot);


%% load benchmarking settings

load(fullfile(benchmarking_folder,'benchmark_settings.mat'),...
    'S','osim_path','S_benchmark');


%% load experimental data

% DOI 10.5281/zenodo.14524619
data_exp_path = fullfile(benchmarking_folder,'exp_data');
if ~isfolder(data_exp_path)
    mkdir(data_exp_path);
end

% general data
% check if we have to download the data
% a bit unclear how we want to handle this.
% osf doesnt work anymore so I think that I will create a zenodo project
% once I'm happy with all the data from my benchmarking paper.
if ~isfolder(fullfile(benchmarking_folder,'exp_data','general_exp_data'))
    % download data as in Afschrift et al. 2025
    LinkZipFile = ['https://zenodo.org/records/16995320/files/' ...
        'general_exp_data.zip?download=1'];
    websave(fullfile(data_exp_path,'general_exp_data.zip'),LinkZipFile);
    unzip(fullfile(data_exp_path,'general_exp_data.zip'),...
        fullfile(data_exp_path,'general_exp_data'));
    delete(fullfile(data_exp_path,'general_exp_data.zip'));
end

% Van Der Zee kinematic data
if isfield(S_benchmark,'studies') && ~isempty(S_benchmark.studies)
    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        if ~isfolder(fullfile(benchmarking_folder,'exp_data','vanderzee2022'))
            % url to datafile
            LinkZipFile = ['https://zenodo.org/records/16995320/files/' ...
                'vanderzee2022_nondim.zip?download=1'];
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

    % subject properties for all experimental studies
    Studyname = {'Gomenuka','Browning','Schertzer','Huang','Koelewijn',...
        'vanderzee2022','McDonald','Abe2015','Strutzenberger'};
    prop_leg_length = 0.5; % assumes leg length is half of model
    Mass = [71.6, 74.16, 74.88, 71.1, 70, 73.5, 69.6, 58.9, 73.1];
    Height = [1.78, 1.82, 1.78, 0.99/prop_leg_length, 1.73, 1.76, 1.70, 1.70, 1.77];
    LegLength = Height.*prop_leg_length;
    Exp_SubjProp = table(Studyname',Mass',Height',LegLength','VariableNames',{'Study','mass','height','LegLength'});
    msim = 62; % should read this from the model in the future
    Height_Sim = 1.70; % I should read this from the model in the future
    Lsim = Height_Sim.*prop_leg_length;
    g = 9.81;


    % attach to simulation results
    % store everything in a matlab datastruture with the fields
    % everything should be nondim
    %   - grf           OR: GRF_r and GRF_l
    %   - grf_std       OR: GRF_r_std and GRF_l_std
    %   - ik
    %   - ik_std
    %   - id
    %   - id_std
    %   - stride_frequency
    %   - Pmetab_mean
    %   - subject_height: (average) height of participant
    %   - subject_mass: (average) mass of participant
    % note that for ik and id it is important that the same model is used
    % for experiment and predictive simulation. The generalized coordinates
    % should at least have the same name between both models if you want that
    % the code runs.
    %
    % to nondim:
    %   frequency:  sqrt(g/l)
    %   moments:    m*g*l
    %   forces:     m*g
    %   Pmetab:     m*g^1.5*sqrt(l)


    % ----- Van Der Zee
    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        % load the sf datafile, this contains unfortunatly all conditions
        % todo: see in script of benchmarking paper to select correct ones
        sf_datafile = fullfile(data_exp_path,'vanderzee2022','stridefreq.csv');
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

            sim_res_folder = fullfile(benchmarking_folder,'vanderzee2022',...
                S_benchmark.vanderzee.names{isim});
            mat_files = dir(fullfile(sim_res_folder,'*.mat'));
            if length(mat_files) ~= 1
                disp('warning mutiple mat files in folder')
                disp(sim_res_folder);
                disp(['assumes that file ' mat_files(1).name, ...
                    'contains the simulation results'])
            end
            sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);

            % get average subject properties
            % we can compute average subject properties here from a file
            % Lexp = Exp_SubjProp.height(strcmp(Exp_SubjProp.Study,'vanderzee2022'));
            % mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'vanderzee2022'));
            file_subj_prop = fullfile(data_exp_path,'vanderzee2022',...
                'subj_prop.csv');
            subj_prop_all = readtable(file_subj_prop);
            mexp = nanmean(subj_prop_all.mass);
            Lexp = nanmean(subj_prop_all.height);

            % compute average stride frequency
            sf_all  = readtable(sf_datafile);
            sf_mean = nanmean(sf_all.stridefreq(sf_all.walkspeed == speed));

            % add all data to benchmark structure
            % note that we expact that all data in benchmark is nondim
            benchmark.subject_height    = Lexp;
            benchmark.subject_mass      = mexp; % to do check in publication
            benchmark.prop_leg_length   = prop_leg_length;
            benchmark.leglength         = benchmark.subject_height *  benchmark.prop_leg_length;
            benchmark.grf               = readtable(grf_datafile);
            benchmark.grf_std           = readtable(grf_datafile_std);
            benchmark.id                = readtable(id_datafile); % check if already norm to mass
            benchmark.id_std            = readtable(id_datafile_std);
            benchmark.ik                = readtable(ik_datafile);
            benchmark.ik_std            = readtable(ik_datafile_std);
            benchmark.Pmetab_mean       = []; % no data
            benchmark.stride_frequency  = sf_mean;

            % note I decide here on grf headers like this:
            %  Flx         Fly        Flz         Frx       Fry        Frz
            % I might want to do this with the same headers as in predsim
            % grf_headers = GRF_r: [100×3 double] and GRF_l: [100×3 double]
            % ToDo: discuss this in predsim meeting. I think I'm in favour of adapting
            % experimental data such that it has the same structure as the predsim results
            benchmark.GRF_r = [-benchmark.grf.Flx, benchmark.grf.Fly, -benchmark.grf.Flz]; % have to switch l and r and rotate here for some reason
            benchmark.GRF_l = [-benchmark.grf.Frx, benchmark.grf.Fry, -benchmark.grf.Frz]; % have to switch l and r and rotate here for some reason
            benchmark.GRF_r_std = [benchmark.grf_std.Frx, benchmark.grf_std.Fry, benchmark.grf_std.Frz];
            benchmark.GRF_l_std = [benchmark.grf_std.Flx, benchmark.grf_std.Fly, benchmark.grf_std.Flz];

            % add benchmark to .mat file
            save(sim_res_file, 'benchmark', "-append");
            clear benchmark;

        end
    end
    
    % attach data Koelewijn
    if any(strcmp(S_benchmark.studies,'koelewijn2019'))

        % load matlab structure with generic experimental data
        if ~exist('ExperimentalData','var')
            load(fullfile(benchmarking_folder,'exp_data',...
                'ExperimentalData.mat'),'ExperimentalData');
        end
        sExp = 11; % selected subject Koelewijn
        % names of experimental data associated with simulations
        % be carefull that the order matters here
        Names_ExpData = {['Koel2019_S' num2str(sExp) '_downhill_08'],...
            ['Koel2019_S' num2str(sExp) '_downhill_13'],...
            ['Koel2019_S' num2str(sExp) '_uphill_08'],...
            ['Koel2019_S' num2str(sExp) '_uphill_13'],...
            ['Koel2019_S' num2str(sExp) '_level_08'],...
            ['Koel2019_S' num2str(sExp) '_Level_13']};

        % read table with experimental data metabolic energy
        % experimental data
        Koel.Exp = importdata(fullfile(benchmarking_folder,'exp_data',...
            'general_exp_data','Koelewijn2019','DataTableKoelewijn.mat'));
        slopes_sim = [-8 -8, 8, 8, 0, 0];


        for isim =1:6%length(S_benchmark.koelewijn.names)
            % load sim file
            sim_res_folder = fullfile(benchmarking_folder,'koelewijn2019',...
                S_benchmark.koelewijn.names{isim});
            mat_files = dir(fullfile(sim_res_folder,'*.mat'));
            if length(mat_files) ~= 1
                disp('warning mutiple mat files in folder')
                disp(sim_res_folder);
                disp(['assumes that file ' mat_files(1).name, ...
                    'contains the simulation results'])
            end
            sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
            SimRes = load(sim_res_file);
            % attach experimental data
            name_exp_sel = Names_ExpData{isim};
            benchmark.subject_height    = Exp_SubjProp.height(strcmp(Exp_SubjProp.Study,'Koelewijn'));
            benchmark.subject_mass      = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn')); % to do check in publication
            benchmark.prop_leg_length   = prop_leg_length;
            benchmark.leglength         = benchmark.subject_height *  benchmark.prop_leg_length;
            benchmark.grf = array2table(ExperimentalData.GRFs.(name_exp_sel).mean,...
                'VariableNames',{'Fx','Fy','Fz'});
            benchmark.grf_std = array2table(ExperimentalData.GRFs.(name_exp_sel).std,...
                'VariableNames',{'Fx','Fy','Fz'});
            benchmark.id = array2table(ExperimentalData.Torques.(name_exp_sel).mean,...
                'VariableNames',ExperimentalData.Torques.(name_exp_sel).colheaders);
            benchmark.id_std = array2table(ExperimentalData.Torques.(name_exp_sel).std,...
                'VariableNames',ExperimentalData.Torques.(name_exp_sel).colheaders);;
            benchmark.ik = array2table(ExperimentalData.Q.(name_exp_sel).Qs.mean,...
                'VariableNames',ExperimentalData.Q.(name_exp_sel).Qs.colheaders);
            benchmark.ik_std = array2table(ExperimentalData.Q.(name_exp_sel).Qs.std,...
                'VariableNames',ExperimentalData.Q.(name_exp_sel).Qs.colheaders);            
            benchmark.stride_frequency  = ExperimentalData.StrideFreq.(name_exp_sel).mean;

            % metabolic energy
            %   find correct data
            metab = squeeze(Koel.Exp.DatKoel(:,strcmp(Koel.Exp.HeaderKoel,'Pmetab'),:));
            speed = squeeze(Koel.Exp.DatKoel(:,strcmp(Koel.Exp.HeaderKoel,'speed'),1));
            slope = squeeze(Koel.Exp.DatKoel(:,strcmp(Koel.Exp.HeaderKoel,'slope'),1));
            irow = speed == SimRes.R.S.misc.forward_velocity & slope == slopes_sim(isim);
            cot_mean = nanmean(metab,2); % this is J/kg/m
            metab= cot_mean(irow) * SimRes.R.S.misc.forward_velocity * benchmark.subject_mass ;
            benchmark.Pmetab_mean       = metab./(benchmark.subject_mass* ...
                9.81^1.5*sqrt(benchmark.leglength));

            % nondim id
            benchmark.id.Variables = benchmark.id.Variables / ...
                (benchmark.leglength * 9.81 * benchmark.subject_mass);
            benchmark.id_std.Variables = benchmark.id_std.Variables / ...
                (benchmark.leglength  * 9.81 * benchmark.subject_mass);

            % nondim sf
            benchmark.stride_frequency = benchmark.stride_frequency./ ...
                sqrt(9.81/benchmark.leglength);

            % non grf: weird unit in Antoines experimental data structure
            %   it is expressed there as a % of body weight
            %   todo: adapt this in the paper (adapt unit of ylabel)
            benchmark.grf.Variables = benchmark.grf.Variables/100;
            benchmark.grf_std.Variables = benchmark.grf_std.Variables/100;

            % adapt simulated GRF such that y-axis points in vertical
            % direction
            % problem is that like this you can only do this once. not
            % desired
            slope =slopes_sim(isim)/100;
            fi = atan(slope);
            Rotm = rotz(fi);
            GRF_rot = SimRes.R.ground_reaction.GRF_r*Rotm(1:3,1:3)';
            SimRes.R.ground_reaction.GRF_r_rot =  GRF_rot;
            GRF_rot = SimRes.R.ground_reaction.GRF_l*Rotm(1:3,1:3)';
            SimRes.R.ground_reaction.GRF_l_rot =  GRF_rot;
            

            % add benchmark to .mat file
            R = SimRes.R;
            save(sim_res_file, 'benchmark','R', "-append");
            clear benchmark;

        end
    end

    % -- Browning 2008
    if any(strcmp(S_benchmark.studies,'browning2008'))
        % no experimental kinematics and kinetics
        % load experimental data from table
        % create table here manually, in same order as in the predictive
        % simulations
        tab_browning = {
            'femur4', NaN, NaN;...
            'femur8', NaN, 3.04;...
            'femur16', 0.89, 3.56;...
            'foot4', 0.84, 3.15;...
            'foot8', 0.80, 3.94;...
            'pelvis4', NaN, 2.53;...
            'pelvis8', NaN, 2.67; ...
            'pelvis12', NaN, 2.67; ...
            'pelvis16', 0.88, 3.10;...
            'tibia4', NaN, 2.73;...
            'tibia8', 0.87, 2.89;...     
            'level', 0.88, 2.34};
        browning_walkspeed = 1.25;
        % temp fix for a bug
        if isfield(S_benchmark,'browning2008')
            S_benchmark.browning.names = S_benchmark.browning2008.names;
        end
        nsim = length(S_benchmark.browning.names);
        for isim = 1:nsim
            sim_res_folder = fullfile(benchmarking_folder,'browning2008',...
                S_benchmark.browning.names{isim});
            mat_files = dir(fullfile(sim_res_folder,'*.mat'));
            if length(mat_files) ~= 1
                disp('warning mutiple mat files in folder')
                disp(sim_res_folder);
                disp(['assumes that file ' mat_files(1).name, ...
                    'contains the simulation results'])
            end
            sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
            SimRes = load(sim_res_file);
            % attach data to benchmarking

            benchmark.subject_height    = Exp_SubjProp.height(strcmp(Exp_SubjProp.Study,'Browning'));
            benchmark.subject_mass      = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning')); % to do check in publication
            benchmark.prop_leg_length   = prop_leg_length;
            benchmark.leglength         = benchmark.subject_height *  benchmark.prop_leg_length;
            benchmark.grf               = [];
            benchmark.grf_std           = [];
            benchmark.id                = []; % check if already norm to mass
            benchmark.id_std            = [];
            benchmark.ik                = [];
            benchmark.ik_std            = [];
            % metabolic power
            benchmark.Pmetab_mean       = tab_browning{isim,3}*benchmark.subject_mass./...
                (benchmark.subject_mass*9.81^1.5*sqrt(benchmark.leglength));
            % stride frequency
            benchmark.stride_frequency  =  tab_browning{isim,2}./(sqrt(g/benchmark.leglength));
            save(sim_res_file, 'benchmark', "-append");
            clear benchmark;

        end
    end
end


%% Default plots
if BoolPlot
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

    ct_sim = 1;

    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
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
        for ispeed_sort = 1:length(speeds_sort)
            ispeed = isort(ispeed_sort);
            sim_res_folder = fullfile(benchmarking_folder,'vanderzee2022',...
                S_benchmark.vanderzee.names{ispeed});
            mat_files = dir(fullfile(sim_res_folder,'*.mat'));
            if length(mat_files) ~= 1
                disp('warning mutiple mat files in folder')
                disp(sim_res_folder);
                disp(['assumes that file ' mat_files(1).name, ...
                    'contains the simulation results'])
            end
            sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
            load(sim_res_file,'R','benchmark');
            Dat(ispeed_sort).benchmark = benchmark;
            Dat(ispeed_sort).R = R;
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
        if ~isempty(Dat(1).benchmark.GRF_r)
            figure('Name','vanderzee: grf','Color',[1 1 1]);
            for ispeed = 1: length(speeds_sort)
                for coord = 1:3
                    subplot(2,3,coord)
                    % experimental grf
                    Fsel = Dat(ispeed).benchmark.GRF_r(:,coord);
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
            title_grf = {'GRFx','GRFy','GRFz'};
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
        
        % store results in a table
        



        clear Dat
    end


    % Plot results Koelewijn
    if any(strcmp(S_benchmark.studies,'koelewijn2019'))
        % read all the data
        for isim =1:6%length(S_benchmark.koelewijn.names)
            for idof = 1:length(dofs_plot)
                % load sim file
                sim_res_folder = fullfile(benchmarking_folder,'koelewijn2019',...
                    S_benchmark.koelewijn.names{isim});
                mat_files = dir(fullfile(sim_res_folder,'*.mat'));
                if length(mat_files) ~= 1
                    disp('warning mutiple mat files in folder')
                    disp(sim_res_folder);
                    disp(['assumes that file ' mat_files(1).name, ...
                        'contains the simulation results'])
                end
                sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
                load(sim_res_file,'R','benchmark');
                Dat(isim).benchmark = benchmark;
                Dat(isim).R = R;
            end
        end

        % default plot function for koelewijn
        bool_rot_grf = true;
        default_plot_benchmarking(Dat, dofs_plot, 'Koelewijn',msim,Lsim,bool_rot_grf);        
        clear Dat
    end

    % Plot results Koelewijn
    if any(strcmp(S_benchmark.studies,'browning2008'))


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
                    disp(['assumes that file ' mat_files(1).name, ...
                        'contains the simulation results'])
                end
                sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
                load(sim_res_file,'R','benchmark');
                Dat(isim).benchmark = benchmark;
                Dat(isim).R = R;
            end
        end

        % default plot for Browning
        bool_rot_grf = false;
        default_plot_benchmarking(Dat, dofs_plot, 'Browning',msim,Lsim,bool_rot_grf);        
        clear Dat    
    end


    %% Plot figure with stride frequency and metabolic power for all datapoints
    % the approach is quite simple here, we load all .mat files with
    % simulations results and plot the datapoints.
    
    disp('work in progress');

        % also store table with all results
    headers_table = {'sim_stride_frequency','sim_metabolic_power','speed',...
        'slope','id_study','exp_stride_frequency','exp_metabolic_power'};
    data_table = nan(1000, length(headers_table));
    ct_sim = 1;
    % id_study: 1 = vanderzee
    %           2 = koelewijn
    %           3 = browning
    if any(strcmp(S_benchmark.studies,'browning2008'))
        % read all the data
        nsim = length(S_benchmark.browning.names);
        for isim =1:nsim
            % load sim file
            sim_res_folder = fullfile(benchmarking_folder,'browning2008',...
                S_benchmark.browning.names{isim});
            slope = 0;
            id_study = 3;
            [data_table, ct_sim] = add_benchmark_to_table(ct_sim,...
                sim_res_folder, headers_table, data_table, slope, id_study,...
                Lsim, msim);
        end
    end
    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        % read all the data
        nsim = length(S_benchmark.vanderzee.names);
        for isim =1:nsim
            % load sim file
            sim_res_folder = fullfile(benchmarking_folder,'vanderzee2022',...
                S_benchmark.vanderzee.names{isim});
            slope = 0;
            id_study = 1;
            [data_table, ct_sim] = add_benchmark_to_table(ct_sim,...
                sim_res_folder, headers_table, data_table, slope, id_study,...
                Lsim, msim);
        end
    end

    % Plot results Koelewijn
    if any(strcmp(S_benchmark.studies,'koelewijn2019'))
        % read all the data
        slopes = [-8, -8, 8, 8, 0, 0];
        for isim =1:6%length(S_benchmark.koelewijn.names)
            % load sim file
            sim_res_folder = fullfile(benchmarking_folder,'koelewijn2019',...
                S_benchmark.koelewijn.names{isim});
            slope = slopes(isim);
            id_study = 2;
            [data_table, ct_sim] = add_benchmark_to_table(ct_sim,...
                sim_res_folder, headers_table, data_table, slope, id_study,...
                Lsim, msim);
        end
    end

    data_table(ct_sim,:) = [];
    table_all = array2table(data_table,...
    'VariableNames',headers_table);

    
    







end




end