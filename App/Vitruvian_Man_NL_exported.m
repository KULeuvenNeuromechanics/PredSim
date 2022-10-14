classdef Vitruvian_Man_NL_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        Image2                          matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        kmuLabel                        matlab.ui.control.Label
        SnelheidSlider                  matlab.ui.control.Slider
        SnelheidSliderLabel             matlab.ui.control.Label
        NaamAfstandSnelheidListBox      matlab.ui.control.ListBox
        NaamAfstandLabel                matlab.ui.control.Label
        SluitvideoButton                matlab.ui.control.Button
        SpeelvideoButton                matlab.ui.control.Button
        SpierkrachtEditField            matlab.ui.control.NumericEditField
        SpierkrachtEditFieldLabel       matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        MassaEditField                  matlab.ui.control.NumericEditField
        MassaEditFieldLabel             matlab.ui.control.Label
        kgLabel                         matlab.ui.control.Label
        NaamEditField                   matlab.ui.control.EditField
        StartsimulatieButton            matlab.ui.control.Button
        DeManvanVitruviusLabel          matlab.ui.control.Label
        MaaktekeningButton              matlab.ui.control.Button
        LengtevanvoetEditField          matlab.ui.control.NumericEditField
        LengtevanvoetEditFieldLabel     matlab.ui.control.Label
        cmLabel_9                       matlab.ui.control.Label
        cmLabel_6                       matlab.ui.control.Label
        cmLabel_5                       matlab.ui.control.Label
        cmLabel_4                       matlab.ui.control.Label
        cmLabel_3                       matlab.ui.control.Label
        cmLabel_2                       matlab.ui.control.Label
        cmLabel                         matlab.ui.control.Label
        AfstandtussenschoudersEditField  matlab.ui.control.NumericEditField
        AfstandtussenschoudersEditFieldLabel  matlab.ui.control.Label
        AfstandvangrondtotknieEditField  matlab.ui.control.NumericEditField
        AfstandvangrondtotknieEditFieldLabel  matlab.ui.control.Label
        AfstandvanknietotheupEditField  matlab.ui.control.NumericEditField
        AfstandvanknietotheupEditFieldLabel  matlab.ui.control.Label
        AfstandvanelleboogtotschouderEditField  matlab.ui.control.NumericEditField
        AfstandvanelleboogtotschouderEditFieldLabel  matlab.ui.control.Label
        AfstandvanvingertoptotelleboogEditField  matlab.ui.control.NumericEditField
        AfstandvanvingertoptotelleboogLabel  matlab.ui.control.Label
        HoogteEditField                 matlab.ui.control.NumericEditField
        HoogteEditFieldLabel            matlab.ui.control.Label
        GroepEditField                  matlab.ui.control.EditField
        GroepEditFieldLabel             matlab.ui.control.Label
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        % Layout
        ink_colour = [139,69,19]/256; % RGB code for "ink" i.e. text and figures
        paper_colour = [235,222,173]/256; % RGB code for "paper" i.e. background

        % Default values
        default_height = 1.8; % [m]
        default_ratio_fingertip_elbow = 1/4;
        default_ratio_elbow_shoulder = 1/8;
        default_ratio_shoulder_width = 1/4;
        default_ratio_hip_knee = 1/4;
        default_ratio_knee_ground = 1/4;
        default_ratio_foot_length = 1/6;

        % Allocate user input values
        usr_height = 1.8;
        usr_fingertip_elbow = 1/4*1.8;
        usr_elbow_shoulder = 1/8*1.8;
        usr_shoulder_width = 1/4*1.8;
        usr_hip_knee = 1/4*1.8;
        usr_knee_ground = 1/4*1.8;
        usr_foot_length = 1/6*1.8;

        usr_speed = 1.2;

        % scale factors
        scale_factors

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

        % user inputs with default values
        function [] = updateUserInput(app)
%             app.usr_height = app.default_height;
            app.usr_fingertip_elbow = app.default_ratio_fingertip_elbow*app.usr_height;
            app.usr_elbow_shoulder = app.default_ratio_elbow_shoulder*app.usr_height;
            app.usr_shoulder_width = app.default_ratio_shoulder_width*app.usr_height;
            app.usr_hip_knee = app.default_ratio_hip_knee*app.usr_height;
            app.usr_knee_ground = app.default_ratio_knee_ground*app.usr_height;
            app.usr_foot_length = app.default_ratio_foot_length*app.usr_height;
        end

        % read user inputs
        function [] = readUserInput(app)
            app.usr_height = app.HoogteEditField.Value;
            app.usr_fingertip_elbow = app.AfstandvanvingertoptotelleboogEditField.Value;
            app.usr_elbow_shoulder = app.AfstandvanelleboogtotschouderEditField.Value;
            app.usr_shoulder_width = app.AfstandtussenschoudersEditField.Value;
            app.usr_hip_knee = app.AfstandvanknietotheupEditField.Value;
            app.usr_knee_ground = app.AfstandvangrondtotknieEditField.Value;
            app.usr_foot_length = app.LengtevanvoetEditField.Value;
        end

        % set default user inputs
        function [] = writeDefaultUserInput(app)
            app.AfstandvanvingertoptotelleboogEditField.Value = app.usr_fingertip_elbow;
            app.AfstandvanelleboogtotschouderEditField.Value = app.usr_elbow_shoulder;
            app.AfstandtussenschoudersEditField.Value = app.usr_shoulder_width;
            app.AfstandvanknietotheupEditField.Value = app.usr_hip_knee;
            app.AfstandvangrondtotknieEditField.Value = app.usr_knee_ground;
            app.LengtevanvoetEditField.Value = app.usr_foot_length;
        end
        
        
        % calls updateDrawing
        function updateDrawingWrapper(app)
            % read user inputs
            readUserInput(app)
            % call function to update drawing
            app.scale_factors = updateDrawing(app.usr_height,app.usr_fingertip_elbow,app.usr_elbow_shoulder,...
            app.usr_shoulder_width,app.usr_hip_knee,app.usr_knee_ground,app.usr_foot_length,...
            app.UIAxes,app.ink_colour,app.paper_colour);
        end
        
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
            % fill these in based on your computer
            init_path_savefolder = 'C:\Users\u0150099\OneDrive - KU Leuven\Resultaten_KinderUniversiteit';
            init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.3\Geometry';
            init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';

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
            % make default drawing
            updateDrawingWrapper(app);

            % path repo
            [pathApp,~,~] = fileparts(mfilename('fullpath'));
            [pathRepo,~,~] = fileparts(pathApp);
            app.path_repo = pathRepo;

            addpath([pathRepo '\VariousFunctions'])

            app.sel_mot_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.mot');

            app.sel_osim_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.osim');

            setPaths(app)
            
        end

        % Button pushed function: MaaktekeningButton
        function MaaktekeningButtonPushed(app, event)
           updateDrawingWrapper(app);

        end

        % Button pushed function: StartsimulatieButton
        function StartsimulatieButtonPushed(app, event)
            % update drawing and get scale factors
            updateDrawingWrapper(app);

            if isempty(app.ModelName)
                warning('Geef een naam.')
                return
            end
            if isempty(app.GroupName)
                warning('Kies een groepsnaam.')
                return
            end

            app.Image.Enable = 'on';
            app.Image.Visible = 'on';
            app.Image2.Enable = 'on';
            app.Image2.Visible = 'on';

            U.savefolder = app.path_savefolder;
            U.ModelName = app.ModelName;
            U.GroupName = app.GroupName;
            U.Height = app.usr_height;
            U.Mass = app.MassaEditField.Value;
            U.Force_sf = app.SpierkrachtEditField.Value/100;
            U.Speed = app.usr_speed;
            U.PathCasadi = app.path_casadi;


            % start simulation
            PredSim_wrapper_for_app(U,app.scale_factors);

            app.Image.Enable = 'off';
            app.Image.Visible = 'off';
            app.Image2.Enable = 'off';
            app.Image2.Visible = 'off';

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
            v = app.SnelheidSlider.Value;
            if v < app.SnelheidSlider.Limits(2)
                app.usr_speed = round(app.SnelheidSlider.Value)/3.6;
            else
                app.usr_speed = -1;
            end
            
        end

        % Value changed function: HoogteEditField
        function HoogteEditFieldValueChanged(app, event)
%             value = app.HoogteEditField.Value;
            readUserInput(app)
            updateUserInput(app)
            writeDefaultUserInput(app)
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

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            app.UIAxes.DataAspectRatio = [1 1 1];
            app.UIAxes.XColor = 'none';
            app.UIAxes.YColor = 'none';
            app.UIAxes.Color = 'none';
            app.UIAxes.NextPlot = 'add';
            app.UIAxes.Position = [589 157 846 627];

            % Create GroepEditFieldLabel
            app.GroepEditFieldLabel = uilabel(app.UIFigure);
            app.GroepEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.GroepEditFieldLabel.FontSize = 30;
            app.GroepEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.GroepEditFieldLabel.Position = [88 800 83 40];
            app.GroepEditFieldLabel.Text = 'Groep';

            % Create GroepEditField
            app.GroepEditField = uieditfield(app.UIFigure, 'text');
            app.GroepEditField.ValueChangedFcn = createCallbackFcn(app, @GroepEditFieldValueChanged, true);
            app.GroepEditField.FontName = 'Edwardian Script ITC';
            app.GroepEditField.FontSize = 30;
            app.GroepEditField.FontColor = [0.5412 0.2706 0.0706];
            app.GroepEditField.BackgroundColor = [1 0.9725 0.8627];
            app.GroepEditField.Position = [172 799 344 41];

            % Create HoogteEditFieldLabel
            app.HoogteEditFieldLabel = uilabel(app.UIFigure);
            app.HoogteEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.HoogteEditFieldLabel.FontSize = 30;
            app.HoogteEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.HoogteEditFieldLabel.Position = [88 685 124 40];
            app.HoogteEditFieldLabel.Text = 'Hoogte';

            % Create HoogteEditField
            app.HoogteEditField = uieditfield(app.UIFigure, 'numeric');
            app.HoogteEditField.Limits = [80 300];
            app.HoogteEditField.ValueDisplayFormat = '%111g';
            app.HoogteEditField.ValueChangedFcn = createCallbackFcn(app, @HoogteEditFieldValueChanged, true);
            app.HoogteEditField.HorizontalAlignment = 'center';
            app.HoogteEditField.FontName = 'Edwardian Script ITC';
            app.HoogteEditField.FontSize = 30;
            app.HoogteEditField.FontColor = [0.5412 0.2706 0.0706];
            app.HoogteEditField.BackgroundColor = [1 0.9725 0.8627];
            app.HoogteEditField.Position = [388 687 74 38];
            app.HoogteEditField.Value = 180;

            % Create AfstandvanvingertoptotelleboogLabel
            app.AfstandvanvingertoptotelleboogLabel = uilabel(app.UIFigure);
            app.AfstandvanvingertoptotelleboogLabel.FontName = 'Edwardian Script ITC';
            app.AfstandvanvingertoptotelleboogLabel.FontSize = 30;
            app.AfstandvanvingertoptotelleboogLabel.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvanvingertoptotelleboogLabel.Position = [88 555 285 40];
            app.AfstandvanvingertoptotelleboogLabel.Text = 'Afstand van vingertop tot elleboog';

            % Create AfstandvanvingertoptotelleboogEditField
            app.AfstandvanvingertoptotelleboogEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvanvingertoptotelleboogEditField.Limits = [10 100];
            app.AfstandvanvingertoptotelleboogEditField.ValueDisplayFormat = '%111g';
            app.AfstandvanvingertoptotelleboogEditField.HorizontalAlignment = 'center';
            app.AfstandvanvingertoptotelleboogEditField.FontName = 'Edwardian Script ITC';
            app.AfstandvanvingertoptotelleboogEditField.FontSize = 30;
            app.AfstandvanvingertoptotelleboogEditField.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvanvingertoptotelleboogEditField.BackgroundColor = [1 0.9725 0.8627];
            app.AfstandvanvingertoptotelleboogEditField.Position = [388 557 74 38];
            app.AfstandvanvingertoptotelleboogEditField.Value = 45;

            % Create AfstandvanelleboogtotschouderEditFieldLabel
            app.AfstandvanelleboogtotschouderEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandvanelleboogtotschouderEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.AfstandvanelleboogtotschouderEditFieldLabel.FontSize = 30;
            app.AfstandvanelleboogtotschouderEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvanelleboogtotschouderEditFieldLabel.Position = [88 490 279 40];
            app.AfstandvanelleboogtotschouderEditFieldLabel.Text = 'Afstand van elleboog tot schouder';

            % Create AfstandvanelleboogtotschouderEditField
            app.AfstandvanelleboogtotschouderEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvanelleboogtotschouderEditField.Limits = [5 100];
            app.AfstandvanelleboogtotschouderEditField.ValueDisplayFormat = '%111g';
            app.AfstandvanelleboogtotschouderEditField.HorizontalAlignment = 'center';
            app.AfstandvanelleboogtotschouderEditField.FontName = 'Edwardian Script ITC';
            app.AfstandvanelleboogtotschouderEditField.FontSize = 30;
            app.AfstandvanelleboogtotschouderEditField.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvanelleboogtotschouderEditField.BackgroundColor = [1 0.9725 0.8627];
            app.AfstandvanelleboogtotschouderEditField.Position = [388 492 74 38];
            app.AfstandvanelleboogtotschouderEditField.Value = 22.5;

            % Create AfstandvanknietotheupEditFieldLabel
            app.AfstandvanknietotheupEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandvanknietotheupEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.AfstandvanknietotheupEditFieldLabel.FontSize = 30;
            app.AfstandvanknietotheupEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvanknietotheupEditFieldLabel.Position = [88 295 218 40];
            app.AfstandvanknietotheupEditFieldLabel.Text = 'Afstand van knie tot heup';

            % Create AfstandvanknietotheupEditField
            app.AfstandvanknietotheupEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvanknietotheupEditField.Limits = [10 100];
            app.AfstandvanknietotheupEditField.ValueDisplayFormat = '%111g';
            app.AfstandvanknietotheupEditField.HorizontalAlignment = 'center';
            app.AfstandvanknietotheupEditField.FontName = 'Edwardian Script ITC';
            app.AfstandvanknietotheupEditField.FontSize = 30;
            app.AfstandvanknietotheupEditField.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvanknietotheupEditField.BackgroundColor = [1 0.9725 0.8627];
            app.AfstandvanknietotheupEditField.Position = [388 297 74 38];
            app.AfstandvanknietotheupEditField.Value = 45;

            % Create AfstandvangrondtotknieEditFieldLabel
            app.AfstandvangrondtotknieEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandvangrondtotknieEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.AfstandvangrondtotknieEditFieldLabel.FontSize = 30;
            app.AfstandvangrondtotknieEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvangrondtotknieEditFieldLabel.Position = [88 360 234 40];
            app.AfstandvangrondtotknieEditFieldLabel.Text = 'Afstand van grond tot knie';

            % Create AfstandvangrondtotknieEditField
            app.AfstandvangrondtotknieEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvangrondtotknieEditField.Limits = [10 100];
            app.AfstandvangrondtotknieEditField.ValueDisplayFormat = '%111g';
            app.AfstandvangrondtotknieEditField.HorizontalAlignment = 'center';
            app.AfstandvangrondtotknieEditField.FontName = 'Edwardian Script ITC';
            app.AfstandvangrondtotknieEditField.FontSize = 30;
            app.AfstandvangrondtotknieEditField.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandvangrondtotknieEditField.BackgroundColor = [1 0.9725 0.8627];
            app.AfstandvangrondtotknieEditField.Position = [388 362 74 38];
            app.AfstandvangrondtotknieEditField.Value = 45;

            % Create AfstandtussenschoudersEditFieldLabel
            app.AfstandtussenschoudersEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandtussenschoudersEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.AfstandtussenschoudersEditFieldLabel.FontSize = 30;
            app.AfstandtussenschoudersEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandtussenschoudersEditFieldLabel.Position = [88 425 209 40];
            app.AfstandtussenschoudersEditFieldLabel.Text = 'Afstand tussen schouders';

            % Create AfstandtussenschoudersEditField
            app.AfstandtussenschoudersEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandtussenschoudersEditField.Limits = [10 100];
            app.AfstandtussenschoudersEditField.ValueDisplayFormat = '%111g';
            app.AfstandtussenschoudersEditField.HorizontalAlignment = 'center';
            app.AfstandtussenschoudersEditField.FontName = 'Edwardian Script ITC';
            app.AfstandtussenschoudersEditField.FontSize = 30;
            app.AfstandtussenschoudersEditField.FontColor = [0.5412 0.2706 0.0706];
            app.AfstandtussenschoudersEditField.BackgroundColor = [1 0.9725 0.8627];
            app.AfstandtussenschoudersEditField.Position = [388 427 74 38];
            app.AfstandtussenschoudersEditField.Value = 45;

            % Create cmLabel
            app.cmLabel = uilabel(app.UIFigure);
            app.cmLabel.FontName = 'Edwardian Script ITC';
            app.cmLabel.FontSize = 30;
            app.cmLabel.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel.Position = [487 685 29 40];
            app.cmLabel.Text = 'cm';

            % Create cmLabel_2
            app.cmLabel_2 = uilabel(app.UIFigure);
            app.cmLabel_2.FontName = 'Edwardian Script ITC';
            app.cmLabel_2.FontSize = 30;
            app.cmLabel_2.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel_2.Position = [487 556 29 40];
            app.cmLabel_2.Text = 'cm';

            % Create cmLabel_3
            app.cmLabel_3 = uilabel(app.UIFigure);
            app.cmLabel_3.FontName = 'Edwardian Script ITC';
            app.cmLabel_3.FontSize = 30;
            app.cmLabel_3.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel_3.Position = [487 490 29 40];
            app.cmLabel_3.Text = 'cm';

            % Create cmLabel_4
            app.cmLabel_4 = uilabel(app.UIFigure);
            app.cmLabel_4.FontName = 'Edwardian Script ITC';
            app.cmLabel_4.FontSize = 30;
            app.cmLabel_4.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel_4.Position = [487 425 29 40];
            app.cmLabel_4.Text = 'cm';

            % Create cmLabel_5
            app.cmLabel_5 = uilabel(app.UIFigure);
            app.cmLabel_5.FontName = 'Edwardian Script ITC';
            app.cmLabel_5.FontSize = 30;
            app.cmLabel_5.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel_5.Position = [487 360 29 40];
            app.cmLabel_5.Text = 'cm';

            % Create cmLabel_6
            app.cmLabel_6 = uilabel(app.UIFigure);
            app.cmLabel_6.FontName = 'Edwardian Script ITC';
            app.cmLabel_6.FontSize = 30;
            app.cmLabel_6.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel_6.Position = [487 296 29 40];
            app.cmLabel_6.Text = 'cm';

            % Create cmLabel_9
            app.cmLabel_9 = uilabel(app.UIFigure);
            app.cmLabel_9.FontName = 'Edwardian Script ITC';
            app.cmLabel_9.FontSize = 30;
            app.cmLabel_9.FontColor = [0.5412 0.2706 0.0706];
            app.cmLabel_9.Position = [487 232 29 40];
            app.cmLabel_9.Text = 'cm';

            % Create LengtevanvoetEditFieldLabel
            app.LengtevanvoetEditFieldLabel = uilabel(app.UIFigure);
            app.LengtevanvoetEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.LengtevanvoetEditFieldLabel.FontSize = 30;
            app.LengtevanvoetEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.LengtevanvoetEditFieldLabel.Position = [88 231 140 40];
            app.LengtevanvoetEditFieldLabel.Text = 'Lengte van voet';

            % Create LengtevanvoetEditField
            app.LengtevanvoetEditField = uieditfield(app.UIFigure, 'numeric');
            app.LengtevanvoetEditField.Limits = [10 100];
            app.LengtevanvoetEditField.ValueDisplayFormat = '%111g';
            app.LengtevanvoetEditField.HorizontalAlignment = 'center';
            app.LengtevanvoetEditField.FontName = 'Edwardian Script ITC';
            app.LengtevanvoetEditField.FontSize = 30;
            app.LengtevanvoetEditField.FontColor = [0.5412 0.2706 0.0706];
            app.LengtevanvoetEditField.BackgroundColor = [1 0.9725 0.8627];
            app.LengtevanvoetEditField.Position = [388 233 74 38];
            app.LengtevanvoetEditField.Value = 30;

            % Create MaaktekeningButton
            app.MaaktekeningButton = uibutton(app.UIFigure, 'push');
            app.MaaktekeningButton.ButtonPushedFcn = createCallbackFcn(app, @MaaktekeningButtonPushed, true);
            app.MaaktekeningButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.MaaktekeningButton.FontName = 'Edwardian Script ITC';
            app.MaaktekeningButton.FontSize = 30;
            app.MaaktekeningButton.FontColor = [0.5412 0.2706 0.0706];
            app.MaaktekeningButton.Position = [590 59 314 61];
            app.MaaktekeningButton.Text = 'Maak tekening';

            % Create DeManvanVitruviusLabel
            app.DeManvanVitruviusLabel = uilabel(app.UIFigure);
            app.DeManvanVitruviusLabel.HorizontalAlignment = 'center';
            app.DeManvanVitruviusLabel.FontName = 'Edwardian Script ITC';
            app.DeManvanVitruviusLabel.FontSize = 80;
            app.DeManvanVitruviusLabel.FontColor = [0.5412 0.2706 0.0706];
            app.DeManvanVitruviusLabel.Position = [313 839 1295 104];
            app.DeManvanVitruviusLabel.Text = 'De Man van Vitruvius';

            % Create StartsimulatieButton
            app.StartsimulatieButton = uibutton(app.UIFigure, 'push');
            app.StartsimulatieButton.ButtonPushedFcn = createCallbackFcn(app, @StartsimulatieButtonPushed, true);
            app.StartsimulatieButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.StartsimulatieButton.FontName = 'Edwardian Script ITC';
            app.StartsimulatieButton.FontSize = 30;
            app.StartsimulatieButton.FontColor = [0.5412 0.2706 0.0706];
            app.StartsimulatieButton.Position = [1064 59 314 61];
            app.StartsimulatieButton.Text = 'Start simulatie';

            % Create NaamEditField
            app.NaamEditField = uieditfield(app.UIFigure, 'text');
            app.NaamEditField.ValueChangedFcn = createCallbackFcn(app, @NaamEditFieldValueChanged, true);
            app.NaamEditField.HorizontalAlignment = 'center';
            app.NaamEditField.FontName = 'Edwardian Script ITC';
            app.NaamEditField.FontSize = 30;
            app.NaamEditField.FontColor = [0.5412 0.2706 0.0706];
            app.NaamEditField.BackgroundColor = [1 0.9725 0.8627];
            app.NaamEditField.Position = [839 799 244 41];
            app.NaamEditField.Value = '(Kies een naam)';

            % Create kgLabel
            app.kgLabel = uilabel(app.UIFigure);
            app.kgLabel.FontName = 'Edwardian Script ITC';
            app.kgLabel.FontSize = 30;
            app.kgLabel.FontColor = [0.5412 0.2706 0.0706];
            app.kgLabel.Position = [487 620 25 40];
            app.kgLabel.Text = 'kg';

            % Create MassaEditFieldLabel
            app.MassaEditFieldLabel = uilabel(app.UIFigure);
            app.MassaEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.MassaEditFieldLabel.FontSize = 30;
            app.MassaEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.MassaEditFieldLabel.Position = [88 620 124 40];
            app.MassaEditFieldLabel.Text = 'Massa';

            % Create MassaEditField
            app.MassaEditField = uieditfield(app.UIFigure, 'numeric');
            app.MassaEditField.Limits = [20 300];
            app.MassaEditField.ValueDisplayFormat = '%111g';
            app.MassaEditField.HorizontalAlignment = 'center';
            app.MassaEditField.FontName = 'Edwardian Script ITC';
            app.MassaEditField.FontSize = 30;
            app.MassaEditField.FontColor = [0.5412 0.2706 0.0706];
            app.MassaEditField.BackgroundColor = [1 0.9725 0.8627];
            app.MassaEditField.Position = [388 622 74 38];
            app.MassaEditField.Value = 75;

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.FontName = 'Edwardian Script ITC';
            app.Label.FontSize = 30;
            app.Label.FontColor = [0.5412 0.2706 0.0706];
            app.Label.Position = [487 165 25 40];
            app.Label.Text = '%';

            % Create SpierkrachtEditFieldLabel
            app.SpierkrachtEditFieldLabel = uilabel(app.UIFigure);
            app.SpierkrachtEditFieldLabel.FontName = 'Edwardian Script ITC';
            app.SpierkrachtEditFieldLabel.FontSize = 30;
            app.SpierkrachtEditFieldLabel.FontColor = [0.5412 0.2706 0.0706];
            app.SpierkrachtEditFieldLabel.Position = [88 165 124 40];
            app.SpierkrachtEditFieldLabel.Text = 'Spierkracht';

            % Create SpierkrachtEditField
            app.SpierkrachtEditField = uieditfield(app.UIFigure, 'numeric');
            app.SpierkrachtEditField.Limits = [20 500];
            app.SpierkrachtEditField.ValueDisplayFormat = '%111g';
            app.SpierkrachtEditField.HorizontalAlignment = 'center';
            app.SpierkrachtEditField.FontName = 'Edwardian Script ITC';
            app.SpierkrachtEditField.FontSize = 30;
            app.SpierkrachtEditField.FontColor = [0.5412 0.2706 0.0706];
            app.SpierkrachtEditField.BackgroundColor = [1 0.9725 0.8627];
            app.SpierkrachtEditField.Position = [388 167 74 38];
            app.SpierkrachtEditField.Value = 100;

            % Create SpeelvideoButton
            app.SpeelvideoButton = uibutton(app.UIFigure, 'push');
            app.SpeelvideoButton.ButtonPushedFcn = createCallbackFcn(app, @SpeelvideoButtonPushed, true);
            app.SpeelvideoButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.SpeelvideoButton.FontName = 'Edwardian Script ITC';
            app.SpeelvideoButton.FontSize = 30;
            app.SpeelvideoButton.FontColor = [0.5412 0.2706 0.0706];
            app.SpeelvideoButton.Position = [1538 59 314 61];
            app.SpeelvideoButton.Text = 'Speel video';

            % Create SluitvideoButton
            app.SluitvideoButton = uibutton(app.UIFigure, 'push');
            app.SluitvideoButton.ButtonPushedFcn = createCallbackFcn(app, @SluitvideoButtonPushed, true);
            app.SluitvideoButton.BackgroundColor = [0.9216 0.8706 0.6706];
            app.SluitvideoButton.FontName = 'Edwardian Script ITC';
            app.SluitvideoButton.FontSize = 30;
            app.SluitvideoButton.FontColor = [0.5412 0.2706 0.0706];
            app.SluitvideoButton.Enable = 'off';
            app.SluitvideoButton.Visible = 'off';
            app.SluitvideoButton.Position = [1538 59 314 61];
            app.SluitvideoButton.Text = 'Sluit video';

            % Create NaamAfstandLabel
            app.NaamAfstandLabel = uilabel(app.UIFigure);
            app.NaamAfstandLabel.BackgroundColor = [0.9216 0.8706 0.6706];
            app.NaamAfstandLabel.FontName = 'Edwardian Script ITC';
            app.NaamAfstandLabel.FontSize = 30;
            app.NaamAfstandLabel.FontColor = [0.5412 0.2706 0.0706];
            app.NaamAfstandLabel.Position = [1467 800 384 40];
            app.NaamAfstandLabel.Text = 'Naam              Afstand         Snelheid';

            % Create NaamAfstandSnelheidListBox
            app.NaamAfstandSnelheidListBox = uilistbox(app.UIFigure);
            app.NaamAfstandSnelheidListBox.Items = {};
            app.NaamAfstandSnelheidListBox.ValueChangedFcn = createCallbackFcn(app, @NaamAfstandSnelheidListBoxValueChanged, true);
            app.NaamAfstandSnelheidListBox.FontName = 'Edwardian Script ITC';
            app.NaamAfstandSnelheidListBox.FontSize = 30;
            app.NaamAfstandSnelheidListBox.FontColor = [0.5412 0.2706 0.0706];
            app.NaamAfstandSnelheidListBox.BackgroundColor = [0.9216 0.8706 0.6706];
            app.NaamAfstandSnelheidListBox.Position = [1467 157 384 644];
            app.NaamAfstandSnelheidListBox.Value = {};

            % Create SnelheidSliderLabel
            app.SnelheidSliderLabel = uilabel(app.UIFigure);
            app.SnelheidSliderLabel.FontName = 'Edwardian Script ITC';
            app.SnelheidSliderLabel.FontSize = 30;
            app.SnelheidSliderLabel.FontColor = [0.5412 0.2706 0.0706];
            app.SnelheidSliderLabel.Position = [88 90 123 40];
            app.SnelheidSliderLabel.Text = 'Snelheid';

            % Create SnelheidSlider
            app.SnelheidSlider = uislider(app.UIFigure);
            app.SnelheidSlider.Limits = [1 20];
            app.SnelheidSlider.MajorTicks = [1 5 10 15 20];
            app.SnelheidSlider.MajorTickLabels = {'1', '5', '10', '15', 'max'};
            app.SnelheidSlider.ValueChangedFcn = createCallbackFcn(app, @SnelheidSliderValueChanged, true);
            app.SnelheidSlider.MinorTicks = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
            app.SnelheidSlider.FontName = 'Edwardian Script ITC';
            app.SnelheidSlider.FontSize = 20;
            app.SnelheidSlider.FontColor = [0.5412 0.2706 0.0706];
            app.SnelheidSlider.Position = [226 119 224 3];
            app.SnelheidSlider.Value = 4;

            % Create kmuLabel
            app.kmuLabel = uilabel(app.UIFigure);
            app.kmuLabel.FontName = 'Edwardian Script ITC';
            app.kmuLabel.FontSize = 30;
            app.kmuLabel.FontColor = [0.5412 0.2706 0.0706];
            app.kmuLabel.Position = [474 88 56 40];
            app.kmuLabel.Text = 'km/u';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Enable = 'off';
            app.Image.Visible = 'off';
            app.Image.Position = [1082 687 254 210];
            app.Image.ImageSource = 'balloon.png';

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.Enable = 'off';
            app.Image2.Visible = 'off';
            app.Image2.Position = [1140 744 138 100];
            app.Image2.ImageSource = 'allgifs.gif';

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