classdef Vitruvian_Man_NL_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        UITable                       matlab.ui.control.Table
        PeesstijfheidEditField        matlab.ui.control.NumericEditField
        PeesstijfheidEditFieldLabel   matlab.ui.control.Label
        Label_3                       matlab.ui.control.Label
        msLabel                       matlab.ui.control.Label
        SnelheidSpinner               matlab.ui.control.Spinner
        SnelheidSpinnerLabel          matlab.ui.control.Label
        GroepEditFieldLabel_2         matlab.ui.control.Label
        KlaarvoorsimulatieLamp        matlab.ui.control.Lamp
        KlaarvoorsimulatieLampLabel   matlab.ui.control.Label
        AanhetsimulerenLamp           matlab.ui.control.Lamp
        AanhetsimulerenLampLabel      matlab.ui.control.Label
        PlantairflexielimitEditField  matlab.ui.control.NumericEditField
        PlantairflexielimitEditFieldLabel  matlab.ui.control.Label
        Label_2                       matlab.ui.control.Label
        SluitvideoButton              matlab.ui.control.Button
        SpeelvideoButton              matlab.ui.control.Button
        SpierkrachtEditField          matlab.ui.control.NumericEditField
        SpierkrachtEditFieldLabel     matlab.ui.control.Label
        Label                         matlab.ui.control.Label
        MassaEditField                matlab.ui.control.NumericEditField
        MassaEditFieldLabel           matlab.ui.control.Label
        kgLabel                       matlab.ui.control.Label
        NaamEditField                 matlab.ui.control.EditField
        StartsimulatieButton          matlab.ui.control.Button
        Onderzoeksstage2023PHABLabel  matlab.ui.control.Label
        cmLabel                       matlab.ui.control.Label
        LichaamslengteEditField       matlab.ui.control.NumericEditField
        LichaamslengteEditFieldLabel  matlab.ui.control.Label
        GroepEditField                matlab.ui.control.EditField
        GroepEditFieldLabel           matlab.ui.control.Label
    end

    
    properties (Access = public)
        % Layout
        ink_colour = [139,69,19]/256; % RGB code for "ink" i.e. text and figures
        paper_colour = [235,222,173]/256; % RGB code for "paper" i.e. background

        % storage
        path_repo
        path_results
        path_savefolder
        path_geom
        path_casadi
        ModelName
        GroupName
        tbldata = cell(1,4);
        idx_sel_data = 1;
        sel_mot_file
        sel_osim_file
        t_name
        t_dist
        t_vel
        Table

        t_mass
        t_musStr
        t_pflim
        t_pftenstiff
        t_subj_name


    end
    
    methods (Access = private)
        
        % load results and put them in table
        function loadResultsTable(app)
            % search folder
            MatFiles = dir(fullfile(app.path_savefolder,app.GroupName,'*.mat'));
            % Empty display table if folder is empty
            if isempty(MatFiles)
                app.NaamAfstandSnelheidListBox.Items = {};
                app.NaamAfstandSnelheidListBox.ItemsData = {};
                return
            end
%             t = cell(length(MatFiles), 3);
            % loop over files
            for i=1:length(MatFiles)
                % calc distance walked for given energy budget and get
                % other variables
                [distance_i, avg_v_i, body_mass_i, muscle_strength_i,...
                    plant_flex_lim_i, pf_tend_stiff_i,subject_name_i] = getInfoFromMat(app,char(fullfile(MatFiles(i).folder,MatFiles(i).name)));
                name_tmp = MatFiles(i).name(1:end-4);
                app.t_name{i,1} = name_tmp;
                app.t_dist(i,1) = distance_i;
                app.t_vel(i,1) = avg_v_i;
                app.t_mass(i,1) = body_mass_i;
                app.t_musStr(i,1) = muscle_strength_i*100;
                app.t_pflim(i,1) = plant_flex_lim_i;
                app.t_pftenstiff(i,1) = pf_tend_stiff_i*100;
                app.t_subj_name{i,1} = subject_name_i;

            end

            % populate table
            app.Table = table(app.t_name, app.t_dist, app.t_vel,app.t_mass,...
                app.t_musStr,app.t_pflim,app.t_pftenstiff,app.t_subj_name);
            app.UITable.Data = app.Table;
            app.UITable.ColumnName = {'Filename', 'Afgelegde afstand (m)', 'Snelheid (m/s)',...
                'Massa (kg)', 'Spierkracht (%)', 'Plantair flexie limit (°)',...
                'Peesstijfheid (%)','Modelnaam'};

        end
        

        function [total_distance, avg_vel, body_mass, muscle_strength,...
                    plant_flex_lim, pf_tend_stiff,subject_name] = getInfoFromMat(app,savepath)
            % get values
            E_metab = 778e3; % 778 kJ = 1 pancake
            load(savepath,'R');
            COT = R.metabolics.Bhargava2004.COT;
            total_distance = round(E_metab/(COT*R.misc.body_mass));
            avg_vel = R.kinematics.avg_velocity;

            body_mass = R.misc.body_mass;
            muscle_strength = R.S.subject.MT_params{1, 3};
            plant_flex_lim = R.S.subject.plant_flex_lim_deg;
            pf_tend_stiff = R.S.subject.tendon_stiff_scale{1, 2};
            subject_name = R.S.subject.name;
            clear('R');

        end
        
        function setPaths(app)
            [init_path_savefolder,init_path_geom,init_path_casadi] = getPaths();
            % test given paths
            if ~exist(init_path_savefolder,'dir')
                error('Please set "init_path_savefolder"')
            end
            if ~exist(init_path_geom,'dir')
                error('Please set "init_path_geom"')
            end
            if ~exist(init_path_casadi,'dir')
                error('Please set "init_path_casadi"')
            end

            % set properties
            app.path_savefolder = init_path_savefolder;
            app.path_geom = init_path_geom;
            app.path_casadi = init_path_casadi;

        end
        

    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            % path repo
            [pathApp,~,~] = fileparts(mfilename('fullpath'));
            [pathRepo,~,~] = fileparts(pathApp);
            app.path_repo = pathRepo;

            addpath([pathRepo '\VariousFunctions'])
            addpath([pathRepo '\App'])

            app.sel_mot_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.mot');

            app.sel_osim_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.osim');

            setPaths(app)

        end

        % Button pushed function: StartsimulatieButton
        function StartsimulatieButtonPushed(app, event)

            if isempty(app.ModelName)
                warning('Geef een naam.')
                return
            end
            if isempty(app.GroupName)
                warning('Kies een groepsnaam.')
                return
            end

            app.StartsimulatieButton.Enable = 'off';
            app.AanhetsimulerenLamp.Color = [0.93,0.69,0.13]; % orange
            app.KlaarvoorsimulatieLamp.Color = [1 1 1]; % white
            pause(0.5); % pause added to allow for color update

            U.savefolder = app.path_savefolder;
            U.ModelName = app.ModelName;
            U.GroupName = app.GroupName;
            U.Height = app.LichaamslengteEditField.Value/100;
            U.Mass = app.MassaEditField.Value;
            U.Force_sf = app.SpierkrachtEditField.Value/100;
            U.plant_flex_lim = app.PlantairflexielimitEditField.Value/180*pi;
            U.Speed = app.SnelheidSpinner.Value;
            U.pf_stiff_scale = app.PeesstijfheidEditField.Value/100;
            U.PathCasadi = app.path_casadi;

            sf.foot = 1;
            sf.low_leg = 1;
            sf.upp_leg = 1;
            sf.upp_arm = 1;
            sf.low_arm = 1;
            sf.shoulder = 1;
            sf.torso = 1;

            % start simulation
            PredSim_wrapper_for_app(U,sf);

            app.AanhetsimulerenLamp.Color = [1 1 1]; % white
            app.KlaarvoorsimulatieLamp.Color = [0 1 0]; % green

            % add result to table
            loadResultsTable(app)

            app.StartsimulatieButton.Enable = 'on';

        end

        % Value changed function: GroepEditField
        function GroepEditFieldValueChanged(app, event)
            if isempty(app.GroepEditField.Value)
                return
            end
            % avoid issues with paths
            GroupName_tmp = app.GroepEditField.Value;
            GroupName_tmp(strfind(GroupName_tmp,' ')) = '_';
            GroupName_tmp(strfind(GroupName_tmp,'/')) = '_';
            GroupName_tmp(strfind(GroupName_tmp,'\')) = '_';
            
            while strcmp(GroupName_tmp(1),'_') && length(GroupName_tmp)>1
                GroupName_tmp = GroupName_tmp(2:end);
            end
            
            while strcmp(GroupName_tmp(end),'_') && length(GroupName_tmp)>1
                GroupName_tmp = GroupName_tmp(1:end-1);
            end
            app.GroupName = GroupName_tmp;

            % load table with results for this group
            loadResultsTable(app);
            
        end

        % Value changed function: NaamEditField
        function NaamEditFieldValueChanged(app, event)
            if isempty(app.NaamEditField.Value)
                return
            end
            ModelName_tmp = app.NaamEditField.Value;
            ModelName_tmp(strfind(ModelName_tmp,' ')) = '_';
            ModelName_tmp(strfind(ModelName_tmp,'/')) = '_';
            ModelName_tmp(strfind(ModelName_tmp,'\')) = '_';
            
            while strcmp(ModelName_tmp(1),'_') && length(ModelName_tmp)>1
                ModelName_tmp = ModelName_tmp(2:end);
            end
            
            while strcmp(ModelName_tmp(end),'_') && length(ModelName_tmp)>1
                ModelName_tmp = ModelName_tmp(1:end-1);
            end
            app.ModelName = ModelName_tmp;

        end

        % Button pushed function: SpeelvideoButton
        function SpeelvideoButtonPushed(app, event)

            flag = 0;
            if ~exist(app.sel_osim_file,'file')
                flag = 1;
                disp(['Could not find "' app.sel_osim_file '"'])
            end
            if ~exist(app.sel_mot_file,'file')
                flag = 1;
                disp(['Could not find "' app.sel_mot_file '"'])
            end
            if flag
                app.sel_osim_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.osim');
                app.sel_mot_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.mot');
            end

            % function
            try
                playVideo(app.sel_osim_file,app.sel_mot_file,app.path_geom);
            catch
            end

            % make close button active
            app.SpeelvideoButton.Enable = 'off';
            app.SpeelvideoButton.Visible = 'off';

            app.SluitvideoButton.Enable = 'on';
            app.SluitvideoButton.Visible = 'on';

        end

        % Button pushed function: SluitvideoButton
        function SluitvideoButtonPushed(app, event)
            try
                [~,~] = system('taskkill /IM simbody-visualizer.exe');
            catch
            end

            app.SluitvideoButton.Enable = 'off';
            app.SluitvideoButton.Visible = 'off';

            app.SpeelvideoButton.Enable = 'on';
            app.SpeelvideoButton.Visible = 'on';

        end

        % Value changed function: LichaamslengteEditField
        function LichaamslengteEditFieldValueChanged(app, event)
            app.NaamEditField.Value = '';    
            app.NaamEditField.Placeholder = '(geef nieuwe naam)';
        end

        % Value changed function: MassaEditField
        function MassaEditFieldValueChanged(app, event)
            app.NaamEditField.Value = '';    
            app.NaamEditField.Placeholder = '(geef nieuwe naam)';
        end

        % Cell selection callback: UITable
        function UITableCellSelection(app, event)
            % Get the selected row index
            selectedRow = event.Indices(1);
            
            % Find the index of the 'Filename' column
            filenameColumnIndex = find(strcmp(app.UITable.ColumnName, 'Filename'));
            modelnameColumnIndex = find(strcmp(app.UITable.ColumnName, 'Modelnaam'));
            
            % Retrieve the value from the 'Filename' column of the selected row
            selectedFilename = app.UITable.Data{selectedRow, filenameColumnIndex};
            selectedModelname = app.UITable.Data{selectedRow, modelnameColumnIndex};

            % generate correct file paths
            app.sel_mot_file = fullfile(app.path_savefolder,app.GroupName,[char(selectedFilename) '.mot']);
            osim_tmp = fullfile(app.path_repo,'Subjects',char(selectedModelname),[char(selectedModelname) '.osim']);

            if exist(osim_tmp,'file')
                app.sel_osim_file = osim_tmp;
            else
                res_tmp = app.sel_mot_file;
                res_tmp(end-1) = 'a'; % .mot -> .mat
                load (res_tmp,'model_info');
                app.sel_osim_file = model_info.osim_path;
            end
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Color = [0.9216 0.8706 0.6706];
            app.UIFigure.Position = [1 41 1602 968];
            app.UIFigure.Name = 'MATLAB App';

            % Create GroepEditFieldLabel
            app.GroepEditFieldLabel = uilabel(app.UIFigure);
            app.GroepEditFieldLabel.FontName = 'Bahnschrift';
            app.GroepEditFieldLabel.FontSize = 30;
            app.GroepEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.GroepEditFieldLabel.Position = [88 767 86 42];
            app.GroepEditFieldLabel.Text = 'Groep';

            % Create GroepEditField
            app.GroepEditField = uieditfield(app.UIFigure, 'text');
            app.GroepEditField.ValueChangedFcn = createCallbackFcn(app, @GroepEditFieldValueChanged, true);
            app.GroepEditField.HorizontalAlignment = 'center';
            app.GroepEditField.FontName = 'Bahnschrift';
            app.GroepEditField.FontSize = 22;
            app.GroepEditField.FontColor = [0.5412 0.2706 0.0706];
            app.GroepEditField.BackgroundColor = [1 0.9725 0.8627];
            app.GroepEditField.Placeholder = '(Naam groep)';
            app.GroepEditField.Position = [206 771 306 31];

            % Create LichaamslengteEditFieldLabel
            app.LichaamslengteEditFieldLabel = uilabel(app.UIFigure);
            app.LichaamslengteEditFieldLabel.FontName = 'Bahnschrift';
            app.LichaamslengteEditFieldLabel.FontSize = 30;
            app.LichaamslengteEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.LichaamslengteEditFieldLabel.Position = [88 627 220 37];
            app.LichaamslengteEditFieldLabel.Text = 'Lichaamslengte';

            % Create LichaamslengteEditField
            app.LichaamslengteEditField = uieditfield(app.UIFigure, 'numeric');
            app.LichaamslengteEditField.Limits = [80 300];
            app.LichaamslengteEditField.ValueDisplayFormat = '%111g';
            app.LichaamslengteEditField.ValueChangedFcn = createCallbackFcn(app, @LichaamslengteEditFieldValueChanged, true);
            app.LichaamslengteEditField.HorizontalAlignment = 'center';
            app.LichaamslengteEditField.FontName = 'Bahnschrift';
            app.LichaamslengteEditField.FontSize = 30;
            app.LichaamslengteEditField.FontColor = [0.5412 0.2706 0.0706];
            app.LichaamslengteEditField.BackgroundColor = [1 0.9725 0.8627];
            app.LichaamslengteEditField.Position = [388 626 74 38];
            app.LichaamslengteEditField.Value = 180;

            % Create cmLabel
            app.cmLabel = uilabel(app.UIFigure);
            app.cmLabel.FontName = 'Bahnschrift';
            app.cmLabel.FontSize = 30;
            app.cmLabel.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel.Position = [487 627 46 37];
            app.cmLabel.Text = 'cm';

            % Create Onderzoeksstage2023PHABLabel
            app.Onderzoeksstage2023PHABLabel = uilabel(app.UIFigure);
            app.Onderzoeksstage2023PHABLabel.HorizontalAlignment = 'center';
            app.Onderzoeksstage2023PHABLabel.FontName = 'Bahnschrift';
            app.Onderzoeksstage2023PHABLabel.FontSize = 80;
            app.Onderzoeksstage2023PHABLabel.FontColor = [0.5412 0.2706 0.0706];
            app.Onderzoeksstage2023PHABLabel.Position = [154 844 1295 111];
            app.Onderzoeksstage2023PHABLabel.Text = 'Onderzoeksstage 2023 PH AB';

            % Create StartsimulatieButton
            app.StartsimulatieButton = uibutton(app.UIFigure, 'push');
            app.StartsimulatieButton.ButtonPushedFcn = createCallbackFcn(app, @StartsimulatieButtonPushed, true);
            app.StartsimulatieButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.StartsimulatieButton.FontName = 'Bahnschrift';
            app.StartsimulatieButton.FontSize = 30;
            app.StartsimulatieButton.FontColor = [0.5412 0.2706 0.0706];
            app.StartsimulatieButton.Position = [118 75 314 61];
            app.StartsimulatieButton.Text = 'Start simulatie';

            % Create NaamEditField
            app.NaamEditField = uieditfield(app.UIFigure, 'text');
            app.NaamEditField.ValueChangedFcn = createCallbackFcn(app, @NaamEditFieldValueChanged, true);
            app.NaamEditField.HorizontalAlignment = 'center';
            app.NaamEditField.FontName = 'Bahnschrift';
            app.NaamEditField.FontSize = 22;
            app.NaamEditField.FontColor = [0.5412 0.2706 0.0706];
            app.NaamEditField.BackgroundColor = [1 0.9725 0.8627];
            app.NaamEditField.Placeholder = '(Modelnaam)';
            app.NaamEditField.Position = [269 727 244 31];

            % Create kgLabel
            app.kgLabel = uilabel(app.UIFigure);
            app.kgLabel.FontName = 'Bahnschrift';
            app.kgLabel.FontSize = 30;
            app.kgLabel.FontColor = [0.5412 0.2706 0.0706];
            app.kgLabel.Position = [487 562 38 37];
            app.kgLabel.Text = 'kg';

            % Create MassaEditFieldLabel
            app.MassaEditFieldLabel = uilabel(app.UIFigure);
            app.MassaEditFieldLabel.FontName = 'Bahnschrift';
            app.MassaEditFieldLabel.FontSize = 30;
            app.MassaEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.MassaEditFieldLabel.Position = [88 557 124 42];
            app.MassaEditFieldLabel.Text = 'Massa';

            % Create MassaEditField
            app.MassaEditField = uieditfield(app.UIFigure, 'numeric');
            app.MassaEditField.Limits = [20 300];
            app.MassaEditField.ValueDisplayFormat = '%111g';
            app.MassaEditField.ValueChangedFcn = createCallbackFcn(app, @MassaEditFieldValueChanged, true);
            app.MassaEditField.HorizontalAlignment = 'center';
            app.MassaEditField.FontName = 'Bahnschrift';
            app.MassaEditField.FontSize = 30;
            app.MassaEditField.FontColor = [0.5412 0.2706 0.0706];
            app.MassaEditField.BackgroundColor = [1 0.9725 0.8627];
            app.MassaEditField.Position = [388 561 74 38];
            app.MassaEditField.Value = 75;

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.FontName = 'Bahnschrift';
            app.Label.FontSize = 30;
            app.Label.FontColor = [0.5412 0.2706 0.0706];
            app.Label.Position = [487 487 26 37];
            app.Label.Text = '%';

            % Create SpierkrachtEditFieldLabel
            app.SpierkrachtEditFieldLabel = uilabel(app.UIFigure);
            app.SpierkrachtEditFieldLabel.FontName = 'Bahnschrift';
            app.SpierkrachtEditFieldLabel.FontSize = 30;
            app.SpierkrachtEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.SpierkrachtEditFieldLabel.Position = [88 487 163 37];
            app.SpierkrachtEditFieldLabel.Text = 'Spierkracht';

            % Create SpierkrachtEditField
            app.SpierkrachtEditField = uieditfield(app.UIFigure, 'numeric');
            app.SpierkrachtEditField.Limits = [20 500];
            app.SpierkrachtEditField.ValueDisplayFormat = '%111g';
            app.SpierkrachtEditField.HorizontalAlignment = 'center';
            app.SpierkrachtEditField.FontName = 'Bahnschrift';
            app.SpierkrachtEditField.FontSize = 30;
            app.SpierkrachtEditField.FontColor = [0.5412 0.2706 0.0706];
            app.SpierkrachtEditField.BackgroundColor = [1 0.9725 0.8627];
            app.SpierkrachtEditField.Position = [388 486 74 38];
            app.SpierkrachtEditField.Value = 100;

            % Create SpeelvideoButton
            app.SpeelvideoButton = uibutton(app.UIFigure, 'push');
            app.SpeelvideoButton.ButtonPushedFcn = createCallbackFcn(app, @SpeelvideoButtonPushed, true);
            app.SpeelvideoButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.SpeelvideoButton.FontName = 'Bahnschrift';
            app.SpeelvideoButton.FontSize = 30;
            app.SpeelvideoButton.FontColor = [0.5412 0.2706 0.0706];
            app.SpeelvideoButton.Position = [1220 75 314 61];
            app.SpeelvideoButton.Text = 'Speel video';

            % Create SluitvideoButton
            app.SluitvideoButton = uibutton(app.UIFigure, 'push');
            app.SluitvideoButton.ButtonPushedFcn = createCallbackFcn(app, @SluitvideoButtonPushed, true);
            app.SluitvideoButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.SluitvideoButton.FontName = 'Bahnschrift';
            app.SluitvideoButton.FontSize = 30;
            app.SluitvideoButton.FontColor = [0.5412 0.2706 0.0706];
            app.SluitvideoButton.Enable = 'off';
            app.SluitvideoButton.Visible = 'off';
            app.SluitvideoButton.Position = [1220 75 314 61];
            app.SluitvideoButton.Text = 'Sluit video';

            % Create Label_2
            app.Label_2 = uilabel(app.UIFigure);
            app.Label_2.FontName = 'Bahnschrift';
            app.Label_2.FontSize = 30;
            app.Label_2.FontColor = [0.5412 0.2706 0.0706];
            app.Label_2.Position = [487 419 25 37];
            app.Label_2.Text = '°';

            % Create PlantairflexielimitEditFieldLabel
            app.PlantairflexielimitEditFieldLabel = uilabel(app.UIFigure);
            app.PlantairflexielimitEditFieldLabel.FontName = 'Bahnschrift';
            app.PlantairflexielimitEditFieldLabel.FontSize = 30;
            app.PlantairflexielimitEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.PlantairflexielimitEditFieldLabel.Position = [88 419 263 37];
            app.PlantairflexielimitEditFieldLabel.Text = 'Plantair flexie limit';

            % Create PlantairflexielimitEditField
            app.PlantairflexielimitEditField = uieditfield(app.UIFigure, 'numeric');
            app.PlantairflexielimitEditField.Limits = [-35 35];
            app.PlantairflexielimitEditField.ValueDisplayFormat = '%111g';
            app.PlantairflexielimitEditField.HorizontalAlignment = 'center';
            app.PlantairflexielimitEditField.FontName = 'Bahnschrift';
            app.PlantairflexielimitEditField.FontSize = 30;
            app.PlantairflexielimitEditField.FontColor = [0.5412 0.2706 0.0706];
            app.PlantairflexielimitEditField.BackgroundColor = [1 0.9725 0.8627];
            app.PlantairflexielimitEditField.Position = [388 418 74 38];
            app.PlantairflexielimitEditField.Value = -30;

            % Create AanhetsimulerenLampLabel
            app.AanhetsimulerenLampLabel = uilabel(app.UIFigure);
            app.AanhetsimulerenLampLabel.FontName = 'Bahnschrift';
            app.AanhetsimulerenLampLabel.FontSize = 30;
            app.AanhetsimulerenLampLabel.FontColor = [0.5412 0.2706 0.0706];
            app.AanhetsimulerenLampLabel.Position = [104 216 254 37];
            app.AanhetsimulerenLampLabel.Text = 'Aan het simuleren';

            % Create AanhetsimulerenLamp
            app.AanhetsimulerenLamp = uilamp(app.UIFigure);
            app.AanhetsimulerenLamp.Position = [420 216 42 42];
            app.AanhetsimulerenLamp.Color = [1 1 1];

            % Create KlaarvoorsimulatieLampLabel
            app.KlaarvoorsimulatieLampLabel = uilabel(app.UIFigure);
            app.KlaarvoorsimulatieLampLabel.FontName = 'Bahnschrift';
            app.KlaarvoorsimulatieLampLabel.FontSize = 30;
            app.KlaarvoorsimulatieLampLabel.FontColor = [0.5412 0.2706 0.0706];
            app.KlaarvoorsimulatieLampLabel.Position = [104 169 281 37];
            app.KlaarvoorsimulatieLampLabel.Text = 'Klaar voor simulatie';

            % Create KlaarvoorsimulatieLamp
            app.KlaarvoorsimulatieLamp = uilamp(app.UIFigure);
            app.KlaarvoorsimulatieLamp.Position = [420 169 42 42];

            % Create GroepEditFieldLabel_2
            app.GroepEditFieldLabel_2 = uilabel(app.UIFigure);
            app.GroepEditFieldLabel_2.FontName = 'Bahnschrift';
            app.GroepEditFieldLabel_2.FontSize = 30;
            app.GroepEditFieldLabel_2.FontColor = [0.5412 0.2706 0.0706];
            app.GroepEditFieldLabel_2.Position = [92 717 160 42];
            app.GroepEditFieldLabel_2.Text = 'Modelnaam';

            % Create SnelheidSpinnerLabel
            app.SnelheidSpinnerLabel = uilabel(app.UIFigure);
            app.SnelheidSpinnerLabel.FontName = 'Bahnschrift';
            app.SnelheidSpinnerLabel.FontSize = 30;
            app.SnelheidSpinnerLabel.FontColor = [0.5412 0.2706 0.0706];
            app.SnelheidSpinnerLabel.Position = [88 354 122 37];
            app.SnelheidSpinnerLabel.Text = 'Snelheid';

            % Create SnelheidSpinner
            app.SnelheidSpinner = uispinner(app.UIFigure);
            app.SnelheidSpinner.Step = 0.01;
            app.SnelheidSpinner.FontName = 'Bahnschrift';
            app.SnelheidSpinner.FontSize = 30;
            app.SnelheidSpinner.FontColor = [0.5412 0.2706 0.0706];
            app.SnelheidSpinner.BackgroundColor = [1 0.9686 0.8588];
            app.SnelheidSpinner.Position = [388 353 100 38];
            app.SnelheidSpinner.Value = 1.33;

            % Create msLabel
            app.msLabel = uilabel(app.UIFigure);
            app.msLabel.FontName = 'Bahnschrift';
            app.msLabel.FontSize = 30;
            app.msLabel.FontColor = [0.5412 0.2706 0.0706];
            app.msLabel.Position = [503 357 58 37];
            app.msLabel.Text = 'm/s';

            % Create Label_3
            app.Label_3 = uilabel(app.UIFigure);
            app.Label_3.FontName = 'Bahnschrift';
            app.Label_3.FontSize = 30;
            app.Label_3.FontColor = [0.5412 0.2706 0.0706];
            app.Label_3.Position = [487 286 26 37];
            app.Label_3.Text = '%';

            % Create PeesstijfheidEditFieldLabel
            app.PeesstijfheidEditFieldLabel = uilabel(app.UIFigure);
            app.PeesstijfheidEditFieldLabel.FontName = 'Bahnschrift';
            app.PeesstijfheidEditFieldLabel.FontSize = 30;
            app.PeesstijfheidEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.PeesstijfheidEditFieldLabel.Position = [88 286 178 37];
            app.PeesstijfheidEditFieldLabel.Text = 'Peesstijfheid';

            % Create PeesstijfheidEditField
            app.PeesstijfheidEditField = uieditfield(app.UIFigure, 'numeric');
            app.PeesstijfheidEditField.Limits = [20 500];
            app.PeesstijfheidEditField.ValueDisplayFormat = '%111g';
            app.PeesstijfheidEditField.HorizontalAlignment = 'center';
            app.PeesstijfheidEditField.FontName = 'Bahnschrift';
            app.PeesstijfheidEditField.FontSize = 30;
            app.PeesstijfheidEditField.FontColor = [0.5412 0.2706 0.0706];
            app.PeesstijfheidEditField.BackgroundColor = [1 0.9725 0.8627];
            app.PeesstijfheidEditField.Position = [388 285 74 38];
            app.PeesstijfheidEditField.Value = 100;

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.ColumnName = {'Filename'; 'Afgelegde afstand (m)'; 'Snelheid (m/s)'; 'Massa (kg)'; 'Spierkracht (%)'; 'Plantair flexie limit (°)'; 'Peesstijfheid (%)'; 'Modelnaam'};
            app.UITable.RowName = {};
            app.UITable.CellSelectionCallback = createCallbackFcn(app, @UITableCellSelection, true);
            app.UITable.Position = [607 169 963 632];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Vitruvian_Man_NL_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end