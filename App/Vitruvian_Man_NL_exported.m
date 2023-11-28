classdef Vitruvian_Man_NL_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        KlaarvoorsimulatieLamp        matlab.ui.control.Lamp
        KlaarvoorsimulatieLampLabel   matlab.ui.control.Label
        AanhetsimulerenLamp           matlab.ui.control.Lamp
        AanhetsimulerenLampLabel      matlab.ui.control.Label
        ErrorLamp                     matlab.ui.control.Lamp
        ErrorLampLabel                matlab.ui.control.Label
        PlantairflexielimitEditField  matlab.ui.control.NumericEditField
        PlantairflexielimitEditFieldLabel  matlab.ui.control.Label
        Label_2                       matlab.ui.control.Label
        kmuLabel                      matlab.ui.control.Label
        SnelheidSlider                matlab.ui.control.Slider
        SnelheidSliderLabel           matlab.ui.control.Label
        NaamAfstandSnelheidListBox    matlab.ui.control.ListBox
        NaamAfstandLabel              matlab.ui.control.Label
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
            % loop over files
            for i=1:length(MatFiles)
                % calc distance walked for given energy budget
                [distance_i,avg_v_i] = calcDistance(app,char(fullfile(MatFiles(i).folder,MatFiles(i).name)));
                name_tmp = MatFiles(i).name(1:end-4);
                name_tmp2 = getName(app,name_tmp);
            
                % add name, distance, mot file to array
                app.tbldata{i,1} = [name_tmp2 ': ' num2str(distance_i*1e-3,'%.1f') ' km, aan ' num2str(avg_v_i) ' km/h'];
                app.tbldata{i,2} = fullfile(MatFiles(i).folder,[name_tmp '.mot']);
                app.tbldata{i,3} = distance_i;
                app.tbldata{i,3} = avg_v_i;
                
            end

            if ~isempty(app.tbldata{1,1})
                app.NaamAfstandSnelheidListBox.Items = app.tbldata(:,1);
                app.NaamAfstandSnelheidListBox.ItemsData = app.tbldata(:,2);
            end

        end
        
        function [total_distance, avg_vel] = calcDistance(app,savepath)
            %% get distance
            E_metab = 778e3; % 778 kJ = 1 pancake
            load(savepath,'R');
            COT = R.metabolics.Bhargava2004.COT;
            total_distance = round(E_metab/(COT*R.misc.body_mass));
            avg_vel = round(R.kinematics.avg_velocity*3.6);
            clear('R');

        end

        function name_tmp2 = getName(app,name_tmp)
            if strcmp(name_tmp(end-2:end-1),'_v')
                name_tmp2 = name_tmp(1:end-3);
            elseif strcmp(name_tmp(end-3:end-2),'_v')
                name_tmp2 = name_tmp(1:end-4);
            else
                name_tmp2 = name_tmp;
            end
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

        % Callback function
        function MaaktekeningButtonPushed(app, event)

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

            app.ErrorLamp.Color = [1 1 1]; % white
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
            U.Speed = app.SnelheidSlider.Value/3.6;
            U.PathCasadi = app.path_casadi;

            sf.foot = 1;
            sf.low_leg = 1;
            sf.upp_leg = 1;
            sf.upp_arm = 1;
            sf.low_arm = 1;
            sf.shoulder = 1;
            sf.torso = 1;

            % start simulation
            resultpath = PredSim_wrapper_for_app(U,sf);

            app.ErrorLamp.Color = [1 1 1]; % white
            app.AanhetsimulerenLamp.Color = [1 1 1]; % white
            app.KlaarvoorsimulatieLamp.Color = [0 1 0]; % green

            % add result to table
            loadResultsTable(app)

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

        % Value changed function: NaamAfstandSnelheidListBox
        function NaamAfstandSnelheidListBoxValueChanged(app, event)
            app.sel_mot_file = app.NaamAfstandSnelheidListBox.Value;

            [~,name_tmp,~] = fileparts(app.sel_mot_file);
            name_tmp2 = getName(app,name_tmp);
            osim_tmp = fullfile(app.path_repo,'Subjects',name_tmp2,[name_tmp2 '.osim']);
            if exist(osim_tmp,'file')
                app.sel_osim_file = osim_tmp;
            else
                res_tmp = app.sel_mot_file;
                res_tmp(end-1) = 'a'; % .mot -> .mat
                load (res_tmp,'model_info');
                app.sel_osim_file = model_info.osim_path;
            end
            

        end

        % Value changed function: SnelheidSlider
        function SnelheidSliderValueChanged(app, event)

        end

        % Value changed function: LichaamslengteEditField
        function LichaamslengteEditFieldValueChanged(app, event)

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
            app.UIFigure.Position = [1 41 1920 963];
            app.UIFigure.Name = 'MATLAB App';

            % Create GroepEditFieldLabel
            app.GroepEditFieldLabel = uilabel(app.UIFigure);
            app.GroepEditFieldLabel.FontName = 'Bahnschrift';
            app.GroepEditFieldLabel.FontSize = 30;
            app.GroepEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.GroepEditFieldLabel.Position = [88 798 86 42];
            app.GroepEditFieldLabel.Text = 'Groep';

            % Create GroepEditField
            app.GroepEditField = uieditfield(app.UIFigure, 'text');
            app.GroepEditField.ValueChangedFcn = createCallbackFcn(app, @GroepEditFieldValueChanged, true);
            app.GroepEditField.FontName = 'Bahnschrift';
            app.GroepEditField.FontSize = 22;
            app.GroepEditField.FontColor = [0.5412 0.2706 0.0706];
            app.GroepEditField.BackgroundColor = [1 0.9725 0.8627];
            app.GroepEditField.Placeholder = '(Geef de naam van jouw groepje)';
            app.GroepEditField.Position = [173 802 344 31];

            % Create LichaamslengteEditFieldLabel
            app.LichaamslengteEditFieldLabel = uilabel(app.UIFigure);
            app.LichaamslengteEditFieldLabel.FontName = 'Bahnschrift';
            app.LichaamslengteEditFieldLabel.FontSize = 30;
            app.LichaamslengteEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.LichaamslengteEditFieldLabel.Position = [88 688 220 37];
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
            app.LichaamslengteEditField.Position = [388 687 74 38];
            app.LichaamslengteEditField.Value = 180;

            % Create cmLabel
            app.cmLabel = uilabel(app.UIFigure);
            app.cmLabel.FontName = 'Bahnschrift';
            app.cmLabel.FontSize = 30;
            app.cmLabel.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel.Position = [487 688 46 37];
            app.cmLabel.Text = 'cm';

            % Create Onderzoeksstage2023PHABLabel
            app.Onderzoeksstage2023PHABLabel = uilabel(app.UIFigure);
            app.Onderzoeksstage2023PHABLabel.HorizontalAlignment = 'center';
            app.Onderzoeksstage2023PHABLabel.FontName = 'Bahnschrift';
            app.Onderzoeksstage2023PHABLabel.FontSize = 80;
            app.Onderzoeksstage2023PHABLabel.FontColor = [0.5412 0.2706 0.0706];
            app.Onderzoeksstage2023PHABLabel.Position = [335 832 1295 111];
            app.Onderzoeksstage2023PHABLabel.Text = 'Onderzoeksstage 2023 PH AB';

            % Create StartsimulatieButton
            app.StartsimulatieButton = uibutton(app.UIFigure, 'push');
            app.StartsimulatieButton.ButtonPushedFcn = createCallbackFcn(app, @StartsimulatieButtonPushed, true);
            app.StartsimulatieButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.StartsimulatieButton.FontName = 'Bahnschrift';
            app.StartsimulatieButton.FontSize = 30;
            app.StartsimulatieButton.FontColor = [0.5412 0.2706 0.0706];
            app.StartsimulatieButton.Position = [1022 59 314 61];
            app.StartsimulatieButton.Text = 'Start simulatie';

            % Create NaamEditField
            app.NaamEditField = uieditfield(app.UIFigure, 'text');
            app.NaamEditField.ValueChangedFcn = createCallbackFcn(app, @NaamEditFieldValueChanged, true);
            app.NaamEditField.HorizontalAlignment = 'center';
            app.NaamEditField.FontName = 'Bahnschrift';
            app.NaamEditField.FontSize = 30;
            app.NaamEditField.FontColor = [0.5412 0.2706 0.0706];
            app.NaamEditField.BackgroundColor = [1 0.9725 0.8627];
            app.NaamEditField.Placeholder = '(Jouw naam)';
            app.NaamEditField.Position = [145 743 244 41];

            % Create kgLabel
            app.kgLabel = uilabel(app.UIFigure);
            app.kgLabel.FontName = 'Bahnschrift';
            app.kgLabel.FontSize = 30;
            app.kgLabel.FontColor = [0.5412 0.2706 0.0706];
            app.kgLabel.Position = [487 623 38 37];
            app.kgLabel.Text = 'kg';

            % Create MassaEditFieldLabel
            app.MassaEditFieldLabel = uilabel(app.UIFigure);
            app.MassaEditFieldLabel.FontName = 'Bahnschrift';
            app.MassaEditFieldLabel.FontSize = 30;
            app.MassaEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.MassaEditFieldLabel.Position = [88 618 124 42];
            app.MassaEditFieldLabel.Text = 'Massa';

            % Create MassaEditField
            app.MassaEditField = uieditfield(app.UIFigure, 'numeric');
            app.MassaEditField.Limits = [20 300];
            app.MassaEditField.ValueDisplayFormat = '%111g';
            app.MassaEditField.HorizontalAlignment = 'center';
            app.MassaEditField.FontName = 'Bahnschrift';
            app.MassaEditField.FontSize = 30;
            app.MassaEditField.FontColor = [0.5412 0.2706 0.0706];
            app.MassaEditField.BackgroundColor = [1 0.9725 0.8627];
            app.MassaEditField.Position = [388 622 74 38];
            app.MassaEditField.Value = 75;

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.FontName = 'Bahnschrift';
            app.Label.FontSize = 30;
            app.Label.FontColor = [0.5412 0.2706 0.0706];
            app.Label.Position = [487 548 26 37];
            app.Label.Text = '%';

            % Create SpierkrachtEditFieldLabel
            app.SpierkrachtEditFieldLabel = uilabel(app.UIFigure);
            app.SpierkrachtEditFieldLabel.FontName = 'Bahnschrift';
            app.SpierkrachtEditFieldLabel.FontSize = 30;
            app.SpierkrachtEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.SpierkrachtEditFieldLabel.Position = [88 548 163 37];
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
            app.SpierkrachtEditField.Position = [388 547 74 38];
            app.SpierkrachtEditField.Value = 100;

            % Create SpeelvideoButton
            app.SpeelvideoButton = uibutton(app.UIFigure, 'push');
            app.SpeelvideoButton.ButtonPushedFcn = createCallbackFcn(app, @SpeelvideoButtonPushed, true);
            app.SpeelvideoButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.SpeelvideoButton.FontName = 'Bahnschrift';
            app.SpeelvideoButton.FontSize = 30;
            app.SpeelvideoButton.FontColor = [0.5412 0.2706 0.0706];
            app.SpeelvideoButton.Position = [1481 59 314 61];
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
            app.SluitvideoButton.Position = [1481 59 314 61];
            app.SluitvideoButton.Text = 'Sluit video';

            % Create NaamAfstandLabel
            app.NaamAfstandLabel = uilabel(app.UIFigure);
            app.NaamAfstandLabel.BackgroundColor = [0.9216 0.8706 0.6706];
            app.NaamAfstandLabel.FontName = 'Bahnschrift';
            app.NaamAfstandLabel.FontSize = 22;
            app.NaamAfstandLabel.FontColor = [0.5412 0.2706 0.0706];
            app.NaamAfstandLabel.Position = [1430 798 499 42];
            app.NaamAfstandLabel.Text = 'Naam              Afstand              Snelheid';

            % Create NaamAfstandSnelheidListBox
            app.NaamAfstandSnelheidListBox = uilistbox(app.UIFigure);
            app.NaamAfstandSnelheidListBox.Items = {};
            app.NaamAfstandSnelheidListBox.ValueChangedFcn = createCallbackFcn(app, @NaamAfstandSnelheidListBoxValueChanged, true);
            app.NaamAfstandSnelheidListBox.FontName = 'Blackadder ITC';
            app.NaamAfstandSnelheidListBox.FontSize = 30;
            app.NaamAfstandSnelheidListBox.FontColor = [0.5412 0.2706 0.0706];
            app.NaamAfstandSnelheidListBox.BackgroundColor = [0.9216 0.8706 0.6706];
            app.NaamAfstandSnelheidListBox.Position = [1430 157 384 644];
            app.NaamAfstandSnelheidListBox.Value = {};

            % Create SnelheidSliderLabel
            app.SnelheidSliderLabel = uilabel(app.UIFigure);
            app.SnelheidSliderLabel.FontName = 'Bahnschrift';
            app.SnelheidSliderLabel.FontSize = 25;
            app.SnelheidSliderLabel.FontColor = [0.5412 0.2706 0.0706];
            app.SnelheidSliderLabel.Position = [88 88 123 42];
            app.SnelheidSliderLabel.Text = 'Snelheid';

            % Create SnelheidSlider
            app.SnelheidSlider = uislider(app.UIFigure);
            app.SnelheidSlider.Limits = [1 10];
            app.SnelheidSlider.MajorTicks = [1 5 10 15 20];
            app.SnelheidSlider.MajorTickLabels = {'1', '5', '10', '15', 'max'};
            app.SnelheidSlider.ValueChangedFcn = createCallbackFcn(app, @SnelheidSliderValueChanged, true);
            app.SnelheidSlider.MinorTicks = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
            app.SnelheidSlider.FontName = 'Bahnschrift';
            app.SnelheidSlider.FontSize = 25;
            app.SnelheidSlider.FontColor = [0.5412 0.2706 0.0706];
            app.SnelheidSlider.Position = [226 119 224 3];
            app.SnelheidSlider.Value = 5;

            % Create kmuLabel
            app.kmuLabel = uilabel(app.UIFigure);
            app.kmuLabel.FontName = 'Bahnschrift';
            app.kmuLabel.FontSize = 30;
            app.kmuLabel.FontColor = [0.5412 0.2706 0.0706];
            app.kmuLabel.Position = [474 91 75 37];
            app.kmuLabel.Text = 'km/u';

            % Create Label_2
            app.Label_2 = uilabel(app.UIFigure);
            app.Label_2.FontName = 'Bahnschrift';
            app.Label_2.FontSize = 30;
            app.Label_2.FontColor = [0.5412 0.2706 0.0706];
            app.Label_2.Position = [487 480 25 37];
            app.Label_2.Text = '°';

            % Create PlantairflexielimitEditFieldLabel
            app.PlantairflexielimitEditFieldLabel = uilabel(app.UIFigure);
            app.PlantairflexielimitEditFieldLabel.FontName = 'Bahnschrift';
            app.PlantairflexielimitEditFieldLabel.FontSize = 30;
            app.PlantairflexielimitEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.PlantairflexielimitEditFieldLabel.Position = [88 480 263 37];
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
            app.PlantairflexielimitEditField.Position = [388 479 74 38];
            app.PlantairflexielimitEditField.Value = -30;

            % Create ErrorLampLabel
            app.ErrorLampLabel = uilabel(app.UIFigure);
            app.ErrorLampLabel.FontName = 'Bahnschrift';
            app.ErrorLampLabel.FontSize = 30;
            app.ErrorLampLabel.FontColor = [0.5412 0.2706 0.0706];
            app.ErrorLampLabel.Position = [92 403 78 37];
            app.ErrorLampLabel.Text = 'Error';

            % Create ErrorLamp
            app.ErrorLamp = uilamp(app.UIFigure);
            app.ErrorLamp.Position = [408 403 42 42];
            app.ErrorLamp.Color = [1 1 1];

            % Create AanhetsimulerenLampLabel
            app.AanhetsimulerenLampLabel = uilabel(app.UIFigure);
            app.AanhetsimulerenLampLabel.FontName = 'Bahnschrift';
            app.AanhetsimulerenLampLabel.FontSize = 30;
            app.AanhetsimulerenLampLabel.FontColor = [0.5412 0.2706 0.0706];
            app.AanhetsimulerenLampLabel.Position = [92 353 254 37];
            app.AanhetsimulerenLampLabel.Text = 'Aan het simuleren';

            % Create AanhetsimulerenLamp
            app.AanhetsimulerenLamp = uilamp(app.UIFigure);
            app.AanhetsimulerenLamp.Position = [408 353 42 42];
            app.AanhetsimulerenLamp.Color = [1 1 1];

            % Create KlaarvoorsimulatieLampLabel
            app.KlaarvoorsimulatieLampLabel = uilabel(app.UIFigure);
            app.KlaarvoorsimulatieLampLabel.FontName = 'Bahnschrift';
            app.KlaarvoorsimulatieLampLabel.FontSize = 30;
            app.KlaarvoorsimulatieLampLabel.FontColor = [0.5412 0.2706 0.0706];
            app.KlaarvoorsimulatieLampLabel.Position = [92 306 281 37];
            app.KlaarvoorsimulatieLampLabel.Text = 'Klaar voor simulatie';

            % Create KlaarvoorsimulatieLamp
            app.KlaarvoorsimulatieLamp = uilamp(app.UIFigure);
            app.KlaarvoorsimulatieLamp.Position = [408 306 42 42];

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