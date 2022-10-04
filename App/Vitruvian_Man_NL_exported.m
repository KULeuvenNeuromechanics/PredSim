classdef Vitruvian_Man_NL_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
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
        NaamvangroepEditField           matlab.ui.control.EditField
        NaamvangroepEditFieldLabel      matlab.ui.control.Label
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
        default_ratio_foot_length = 1/7;

        % Allocate user input values
        usr_height = 1.8;
        usr_fingertip_elbow = 1/4*1.8;
        usr_elbow_shoulder = 1/8*1.8;
        usr_shoulder_width = 1/4*1.8;
        usr_hip_knee = 1/4*1.8;
        usr_knee_ground = 1/4*1.8;
        usr_foot_length = 1/7*1.8;

    end
    
    methods (Access = private)

        % initialise
        function [] = initDrawing(app)

            set(app.UIFigure,'defaultAxesFontName','Edwardian Script ITC');

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
        
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: MaaktekeningButton
        function MaaktekeningButtonPushed(app, event)
            % read inputs
            readUserInput(app);
            % call function to update drawing
            updateDrawing(app.usr_height,app.usr_fingertip_elbow,app.usr_elbow_shoulder,...
                app.usr_shoulder_width,app.usr_hip_knee,app.usr_knee_ground,app.usr_foot_length,...
                app.UIAxes,app.ink_colour,app.paper_colour);

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1083 701];
            app.UIFigure.Name = 'De Man van Vitruvius';
            app.UIFigure.Color = app.paper_colour;

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            app.UIAxes.Position = [441 39 592 539];

            % Create NaamvangroepEditFieldLabel
            app.NaamvangroepEditFieldLabel = uilabel(app.UIFigure);
            app.NaamvangroepEditFieldLabel.HorizontalAlignment = 'right';
            app.NaamvangroepEditFieldLabel.Position = [35 577 94 22];
            app.NaamvangroepEditFieldLabel.Text = 'Naam van groep';

            % Create NaamvangroepEditField
            app.NaamvangroepEditField = uieditfield(app.UIFigure, 'text');
            app.NaamvangroepEditField.Position = [140 577 169 22];

            % Create HoogteEditFieldLabel
            app.HoogteEditFieldLabel = uilabel(app.UIFigure);
            app.HoogteEditFieldLabel.Position = [34 446 124 22];
            app.HoogteEditFieldLabel.Text = 'Hoogte';

            % Create HoogteEditField
            app.HoogteEditField = uieditfield(app.UIFigure, 'numeric');
            app.HoogteEditField.Position = [243 446 66 22];
            app.HoogteEditField.Value = 180;

            % Create AfstandvanvingertoptotelleboogLabel
            app.AfstandvanvingertoptotelleboogLabel = uilabel(app.UIFigure);
            app.AfstandvanvingertoptotelleboogLabel.Position = [34 392 187 22];
            app.AfstandvanvingertoptotelleboogLabel.Text = 'Afstand van vingertop tot elleboog';

            % Create AfstandvanvingertoptotelleboogEditField
            app.AfstandvanvingertoptotelleboogEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvanvingertoptotelleboogEditField.Position = [243 392 66 22];
            app.AfstandvanvingertoptotelleboogEditField.Value = 45;

            % Create AfstandvanelleboogtotschouderEditFieldLabel
            app.AfstandvanelleboogtotschouderEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandvanelleboogtotschouderEditFieldLabel.Position = [34 338 187 22];
            app.AfstandvanelleboogtotschouderEditFieldLabel.Text = 'Afstand van elleboog tot schouder';

            % Create AfstandvanelleboogtotschouderEditField
            app.AfstandvanelleboogtotschouderEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvanelleboogtotschouderEditField.Position = [243 338 66 22];
            app.AfstandvanelleboogtotschouderEditField.Value = 22.5;

            % Create AfstandvanknietotheupEditFieldLabel
            app.AfstandvanknietotheupEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandvanknietotheupEditFieldLabel.Position = [34 176 141 22];
            app.AfstandvanknietotheupEditFieldLabel.Text = 'Afstand van knie tot heup';

            % Create AfstandvanknietotheupEditField
            app.AfstandvanknietotheupEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvanknietotheupEditField.Position = [243 176 66 22];
            app.AfstandvanknietotheupEditField.Value = 45;

            % Create AfstandvangrondtotknieEditFieldLabel
            app.AfstandvangrondtotknieEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandvangrondtotknieEditFieldLabel.Position = [34 230 145 22];
            app.AfstandvangrondtotknieEditFieldLabel.Text = 'Afstand van grond tot knie';

            % Create AfstandvangrondtotknieEditField
            app.AfstandvangrondtotknieEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandvangrondtotknieEditField.Position = [243 230 66 22];
            app.AfstandvangrondtotknieEditField.Value = 45;

            % Create AfstandtussenschoudersEditFieldLabel
            app.AfstandtussenschoudersEditFieldLabel = uilabel(app.UIFigure);
            app.AfstandtussenschoudersEditFieldLabel.Position = [34 284 144 22];
            app.AfstandtussenschoudersEditFieldLabel.Text = 'Afstand tussen schouders';

            % Create AfstandtussenschoudersEditField
            app.AfstandtussenschoudersEditField = uieditfield(app.UIFigure, 'numeric');
            app.AfstandtussenschoudersEditField.Position = [243 284 66 22];
            app.AfstandtussenschoudersEditField.Value = 45;

            % Create cmLabel
            app.cmLabel = uilabel(app.UIFigure);
            app.cmLabel.Position = [321 446 25 22];
            app.cmLabel.Text = 'cm';

            % Create cmLabel_2
            app.cmLabel_2 = uilabel(app.UIFigure);
            app.cmLabel_2.Position = [321 392 25 22];
            app.cmLabel_2.Text = 'cm';

            % Create cmLabel_3
            app.cmLabel_3 = uilabel(app.UIFigure);
            app.cmLabel_3.Position = [321 338 25 22];
            app.cmLabel_3.Text = 'cm';

            % Create cmLabel_4
            app.cmLabel_4 = uilabel(app.UIFigure);
            app.cmLabel_4.Position = [321 284 25 22];
            app.cmLabel_4.Text = 'cm';

            % Create cmLabel_5
            app.cmLabel_5 = uilabel(app.UIFigure);
            app.cmLabel_5.Position = [321 230 25 22];
            app.cmLabel_5.Text = 'cm';

            % Create cmLabel_6
            app.cmLabel_6 = uilabel(app.UIFigure);
            app.cmLabel_6.Position = [321 176 25 22];
            app.cmLabel_6.Text = 'cm';

            % Create cmLabel_9
            app.cmLabel_9 = uilabel(app.UIFigure);
            app.cmLabel_9.Position = [322 123 25 22];
            app.cmLabel_9.Text = 'cm';

            % Create LengtevanvoetEditFieldLabel
            app.LengtevanvoetEditFieldLabel = uilabel(app.UIFigure);
            app.LengtevanvoetEditFieldLabel.Position = [35 122 91 22];
            app.LengtevanvoetEditFieldLabel.Text = 'Lengte van voet';

            % Create LengtevanvoetEditField
            app.LengtevanvoetEditField = uieditfield(app.UIFigure, 'numeric');
            app.LengtevanvoetEditField.Position = [244 122 66 22];
            app.LengtevanvoetEditField.Value = 25.7142857142857;

            % Create MaaktekeningButton
            app.MaaktekeningButton = uibutton(app.UIFigure, 'push');
            app.MaaktekeningButton.ButtonPushedFcn = createCallbackFcn(app, @MaaktekeningButtonPushed, true);
            app.MaaktekeningButton.Position = [32 37 314 61];
            app.MaaktekeningButton.Text = 'Maak tekening';

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

            initDrawing(app);
            
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