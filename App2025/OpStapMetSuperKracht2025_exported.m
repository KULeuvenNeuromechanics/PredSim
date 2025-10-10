classdef OpStapMetSuperKracht2025_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        LengteMenu                      matlab.ui.container.Menu
        KrachtMenu                      matlab.ui.container.Menu
        ZwaartekrachtMenu               matlab.ui.container.Menu
        CreerjezelfalsmodelMenu         matlab.ui.container.Menu
        GIF                             matlab.ui.control.Image
        Ballon                          matlab.ui.control.Image
        CheckBox_Maui                   matlab.ui.control.CheckBox
        CheckBox_MrImpossible           matlab.ui.control.CheckBox
        CheckBox_Jullie_strength        matlab.ui.control.CheckBox
        CheckBox_Spongebob              matlab.ui.control.CheckBox
        CheckBox_Olaf                   matlab.ui.control.CheckBox
        Button                          matlab.ui.control.Button
        Maui                            matlab.ui.control.Image
        MrIncredible                    matlab.ui.control.Image
        Jullie_strength                 matlab.ui.control.Image
        Spongebob                       matlab.ui.control.Image
        Olaf                            matlab.ui.control.Image
        GroepEditField                  matlab.ui.control.EditField
        SpeelvideoButton                matlab.ui.control.Button
        SluitvideoButton                matlab.ui.control.Button
        StartsimulatieButton            matlab.ui.control.Button
        OpstapmetsuperkrachtWhite       matlab.ui.control.Label
        OpstapmetsuperkrachtBlue        matlab.ui.control.Label
        BackgroundImageSmurfReusKracht  matlab.ui.control.Image
        earth                           matlab.ui.control.Image
        Mars                            matlab.ui.control.Image
        Jupiter                         matlab.ui.control.Image
        Spaceman                        matlab.ui.control.Image
        spaceBackground                 matlab.ui.control.Image
        Image2                          matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        ZwaartekrachtLabelEigenModel    matlab.ui.control.Label
        SpierkrachtLabelEigenModel      matlab.ui.control.Label
        cmLabelEigenModel               matlab.ui.control.Label
        LengteLabelEigenModel           matlab.ui.control.Label
        LengteSliderEigenModel          matlab.ui.control.RangeSlider
        gLabelEigenModel                matlab.ui.control.Label
        ZwaartekrachtEigenModel         matlab.ui.control.NumericEditField
        PlanetXblack                    matlab.ui.control.Image
        NaamModelEditField              matlab.ui.control.EditField
        MassaEditFieldEigenModel        matlab.ui.control.NumericEditField
        MassaEditFieldLabel             matlab.ui.control.Label
        kgLabelEigenModel               matlab.ui.control.Label
        SpierkrachtEigenModel           matlab.ui.control.NumericEditField
        procentlabelEigenModel          matlab.ui.control.Label
        maan                            matlab.ui.control.Image
        GIFsmurf                        matlab.ui.control.Image
        DENKsmurf                       matlab.ui.control.Image
        GIFreus                         matlab.ui.control.Image
        DENKreus                        matlab.ui.control.Image
        smurf                           matlab.ui.control.Image
        giant                           matlab.ui.control.Image
        ZwaartekrachtEditFieldMaan      matlab.ui.control.NumericEditField
        ZwaartekrachtEditFieldPlanetX   matlab.ui.control.NumericEditField
        ZwaartekrachtEditFieldLabelPlanetX  matlab.ui.control.Label
        gLabel_3                        matlab.ui.control.Label
        ZwaartekrachtEditFieldLabelMaan  matlab.ui.control.Label
        MaanLabel                       matlab.ui.control.Label
        gLabel                          matlab.ui.control.Label
        PlaneetXLabel                   matlab.ui.control.Label
        LengteSliderReusSmurf           matlab.ui.control.RangeSlider
        PlaneetX                        matlab.ui.control.Image
        GroepEditFieldLabel             matlab.ui.control.Label
        Background_creerjeeigenmodel    matlab.ui.control.Image
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

        usr_speed = 1.38;

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
    properties (Access = private)
    CurrentMenu string = "Lengte";   % startwaarde
    end
    
    methods (Access = private)

        % Menu selected function: Lengte
        function LengteMenuSelected(app,event)
            app.CurrentMenu = "Lengte";  
            app.BackgroundImageSmurfReusKracht.Visible = 'on';
            app.OpstapmetsuperkrachtBlue.Visible = 'on';
            app.OpstapmetsuperkrachtWhite.Visible = 'off';

            % related to smurf reus
            app.giant.Visible = 'on';
            app.smurf.Visible = 'on';
            app.LengteSliderReusSmurf.Visible = 'on';

            % related to kracht 
            app.Spongebob.Visible = 'off';
            app.CheckBox_Spongebob.Visible = 'off';

            app.Maui.Visible = 'off';
            app.CheckBox_Maui.Visible = 'off';

            app.MrIncredible.Visible = 'off';
            app.CheckBox_MrImpossible.Visible = 'off';


            app.Jullie_strength.Visible  = 'off';
            app.CheckBox_Jullie_strength.Visible = 'off';
            
            app.Olaf.Visible  = 'off';
            app.CheckBox_Olaf.Visible = 'off';



            % space
            app.Spaceman.Visible = 'off';
            app.spaceBackground.Visible = 'off';
            app.maan.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.PlaneetXLabel.Visible = 'off';
            app.gLabel.Visible = 'off';
            app.ZwaartekrachtEditFieldMaan.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelMaan.Visible = 'off';
            app.gLabel_3.Visible = 'off';
            app.MaanLabel.Visible = 'off';
           app.ZwaartekrachtEditFieldPlanetX.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelPlanetX.Visible = 'off';
            app.earth.Visible ='off';
            app.Mars.Visible ='off';
            app.Jupiter.Visible ='off';
            


             app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
            app.NaamModelEditField.Visible = 'off';
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.MassaEditFieldLabel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';

             app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off';            
        end


        % Menu selected function: Kracht
        function KrachtMenuSelected(app,event)
            app.CurrentMenu = "Kracht";  
            app.BackgroundImageSmurfReusKracht.Visible = 'on';
            app.OpstapmetsuperkrachtBlue.Visible = 'on';
            app.OpstapmetsuperkrachtWhite.Visible = 'off';

            % related to kracht 
            app.Spongebob.Visible = 'on';
            app.CheckBox_Spongebob.Visible  = 'on';

            app.Maui.Visible = 'on';
            app.CheckBox_Maui.Visible = 'on';

            app.MrIncredible.Visible = 'on';
            app.CheckBox_MrImpossible.Visible = 'on';


            app.Jullie_strength.Visible  = 'on';
            app.CheckBox_Jullie_strength.Visible = 'on';
            
            app.Olaf.Visible  = 'on';
            app.CheckBox_Olaf.Visible = 'on';

            % related to smurf reus
            app.giant.Visible = 'off';
            app.smurf.Visible = 'off';
            app.LengteSliderReusSmurf.Visible = 'off';


            % space
            app.Spaceman.Visible = 'off';
            app.spaceBackground.Visible = 'off';
            app.maan.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.PlaneetXLabel.Visible = 'off';
            app.gLabel.Visible = 'off';
            app.ZwaartekrachtEditFieldMaan.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelMaan.Visible = 'off';
            app.gLabel_3.Visible = 'off';
            app.MaanLabel.Visible = 'off';
           app.ZwaartekrachtEditFieldPlanetX.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelPlanetX.Visible = 'off';
            app.earth.Visible ='off';
            app.Mars.Visible ='off';
            app.Jupiter.Visible ='off';


            app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
            app.NaamModelEditField.Visible = 'off';
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.MassaEditFieldLabel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';

             app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off'; 



        end


        % Menu selected function: ZwaarteKracht
        function ZwaartekrachtMenuSelectedFix(app,event)
                app.CurrentMenu = "Zwaartekracht";  


            % related to smurf reus kracht
            app.BackgroundImageSmurfReusKracht.Visible = 'off';
            app.OpstapmetsuperkrachtBlue.Visible = 'off';

            % related to kracht 
            app.Spongebob.Visible = 'off';
            app.CheckBox_Spongebob.Visible = 'off';

            app.Maui.Visible = 'off';
            app.CheckBox_Maui.Visible = 'off';

            app.MrIncredible.Visible = 'off';
            app.CheckBox_MrImpossible.Visible = 'off';


            app.Jullie_strength.Visible  = 'off';
            app.CheckBox_Jullie_strength.Visible = 'off';
            
            app.Olaf.Visible  = 'off';
            app.CheckBox_Olaf.Visible = 'off';

            % related to smurf reus
            app.giant.Visible = 'off';
            app.smurf.Visible = 'off';
            app.LengteSliderReusSmurf.Visible = 'off';

            % space
            app.OpstapmetsuperkrachtWhite.Visible = 'on';
            app.Spaceman.Visible = 'on';
            app.spaceBackground.Visible = 'on';
            app.maan.Visible = 'on';
            app.PlaneetX.Visible = 'on';
            app.PlaneetXLabel.Visible = 'on';
            app.gLabel.Visible = 'on';
            app.ZwaartekrachtEditFieldMaan.Visible = 'on';
            app.ZwaartekrachtEditFieldLabelMaan.Visible = 'on';
            app.gLabel_3.Visible = 'on';
            app.MaanLabel.Visible = 'on';
           app.ZwaartekrachtEditFieldPlanetX.Visible = 'on';
            app.ZwaartekrachtEditFieldLabelPlanetX.Visible = 'on';
            app.earth.Visible ='on';
            app.Mars.Visible ='on';
            app.Jupiter.Visible ='on';

            app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
            app.NaamModelEditField.Visible = 'off';
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.MassaEditFieldLabel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';


             app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off';

        end


        function CreerjezelfalsmodelMenuSelected(app,event)
            app.CurrentMenu = "EigenModel";  
            % related to smurf reus kracht
            app.BackgroundImageSmurfReusKracht.Visible = 'off';
            app.OpstapmetsuperkrachtBlue.Visible = 'off';

            % related to kracht 
              app.Spongebob.Visible = 'off';
            app.CheckBox_Spongebob.Visible = 'off';

            app.Maui.Visible = 'off';
            app.CheckBox_Maui.Visible = 'off';

            app.MrIncredible.Visible = 'off';
            app.CheckBox_MrImpossible.Visible = 'off';


            app.Jullie_strength.Visible  = 'off';
            app.CheckBox_Jullie_strength.Visible = 'off';
            
            app.Olaf.Visible  = 'off';
            app.CheckBox_Olaf.Visible = 'off';
           

            % related to smurf reus
            app.giant.Visible = 'off';
            app.smurf.Visible = 'off';
            app.LengteSliderReusSmurf.Visible = 'off';


            % space eigen model
            app.OpstapmetsuperkrachtWhite.Visible = 'on';
            
            % space
            app.Spaceman.Visible = 'off';
            app.spaceBackground.Visible = 'off';
            app.maan.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.PlaneetXLabel.Visible = 'off';
            app.gLabel.Visible = 'off';
            app.ZwaartekrachtEditFieldMaan.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelMaan.Visible = 'off';
            app.gLabel_3.Visible = 'off';
            app.MaanLabel.Visible = 'off';
           app.ZwaartekrachtEditFieldPlanetX.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelPlanetX.Visible = 'off';
           app.earth.Visible ='off';
            app.Mars.Visible ='off';
            app.Jupiter.Visible ='off';

            app.cmLabelEigenModel.Visible = 'on';
            app.LengteLabelEigenModel.Visible = 'on';
            app.LengteSliderEigenModel.Visible = 'on';
            app.gLabelEigenModel.Visible = 'on';
            app.PlanetXblack.Visible = 'on';
            app.NaamModelEditField.Visible = 'on';
            app.MassaEditFieldEigenModel.Visible = 'on';
            app.MassaEditFieldLabel.Visible = 'on';
            app.kgLabelEigenModel.Visible = 'on';
            app.procentlabelEigenModel.Visible = 'on';
            app.SpierkrachtEigenModel.Visible = 'on';
            app.Background_creerjeeigenmodel.Visible = 'on';
            app.ZwaartekrachtEigenModel.Visible = 'on';
            app.ZwaartekrachtLabelEigenModel.Visible = 'on';
            app.SpierkrachtLabelEigenModel.Visible = 'on';


            % eigen model

        end


  



% 
%         % user inputs with default values
%         function [] = updateUserInput(app)
% %             app.usr_height = app.default_height;
%             app.usr_fingertip_elbow = app.default_ratio_fingertip_elbow*app.usr_height;
%             app.usr_elbow_shoulder = app.default_ratio_elbow_shoulder*app.usr_height;
%             app.usr_shoulder_width = app.default_ratio_shoulder_width*app.usr_height;
%             app.usr_hip_knee = app.default_ratio_hip_knee*app.usr_height;
%             app.usr_knee_ground = app.default_ratio_knee_ground*app.usr_height;
%             app.usr_foot_length = app.default_ratio_foot_length*app.usr_height;
%         end
% 
%         % read user inputs
%         function [] = readUserInput(app)
%             % app.usr_height = app.LichaamslengteEditField.Value;
%             app.usr_fingertip_elbow = app.AfstandelleboogtotvingertopEditField.Value;
%             app.usr_elbow_shoulder = app.AfstandschoudertotelleboogEditField.Value;
%             app.usr_shoulder_width = app.AfstandtussenschoudersEditField.Value;
%             app.usr_hip_knee = app.AfstandheuptotknieEditField.Value;
%             app.usr_knee_ground = app.AfstandknietotgrondEditField.Value;
%             app.usr_foot_length = app.LengtevanvoetEditField.Value;
%         end
% 




        
        
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
            name = getenv('COMPUTERNAME');
% 
%             % fill these in based on your computer
            if strcmp(name,'GBW-L-W2195')
                init_path_savefolder = 'C:\GBW_MyPrograms\KinderuniversiteitApp\Resultaten';
                init_path_geom = 'C:\OpenSim 4.4\Geometry';
                init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
%             elseif strcmp(name, 'GBW-L-W2075')
%                 init_path_savefolder = 'C:\Users\u0138016\OneDrive - KU Leuven\Outreach\Kinderuniversiteit\2022\Resultaten';
%                 init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.3\Geometry';
%                 init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
%             else
%                 init_path_savefolder = 'C:\Users\u0138016\OneDrive - KU Leuven\Outreach\Kinderuniversiteit\2022\Resultaten';
%                 init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.3\Geometry';
%                 init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
            end


            % [init_path_savefolder,init_path_geom,init_path_casadi] = getPaths();
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
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % make default drawing

            % path repo
            [pathApp,~,~] = fileparts(mfilename('fullpath'));
            [pathRepo,~,~] = fileparts(pathApp);
            app.path_repo = pathRepo;

            addpath([pathRepo '\VariousFunctions'])
            addpath([pathRepo '\App'])

            app.sel_mot_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.mot');

            app.sel_osim_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.osim');

            setPaths(app)

            % menu
           app.LengteMenu.MenuSelectedFcn = @(src, event) LengteMenuSelected(app,event);
           app.KrachtMenu.MenuSelectedFcn = @(src, event) KrachtMenuSelected(app,event);
           app.ZwaartekrachtMenu.MenuSelectedFcn = @(src,event) ZwaartekrachtMenuSelectedFix(app,event); 
           app.CreerjezelfalsmodelMenu.MenuSelectedFcn =  @(src,event) CreerjezelfalsmodelMenuSelected(app,event); 
            
           app.CurrentMenu = "Lengte";  
          
            app.Spongebob.Visible = 'off';
            app.CheckBox_Spongebob.Visible = 'off';

            app.Maui.Visible = 'off';
            app.CheckBox_Maui.Visible = 'off';

            app.MrIncredible.Visible = 'off';
            app.CheckBox_MrImpossible.Visible = 'off';


            app.Jullie_strength.Visible  = 'off';
            app.CheckBox_Jullie_strength.Visible = 'off';
            
            app.Olaf.Visible  = 'off';
            app.CheckBox_Olaf.Visible = 'off';
              
         
           app.OpstapmetsuperkrachtWhite.Visible = 'off';
            app.Spaceman.Visible = 'off';
            app.spaceBackground.Visible = 'off';
            app.maan.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.PlaneetXLabel.Visible = 'off';
            app.gLabel.Visible = 'off';
            app.ZwaartekrachtEditFieldMaan.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelMaan.Visible  = 'off';
            app.gLabel_3.Visible = 'off';
            app.MaanLabel.Visible = 'off';
           app.ZwaartekrachtEditFieldPlanetX.Visible = 'off';
            app.ZwaartekrachtEditFieldLabelPlanetX.Visible = 'off';
            app.earth.Visible ='off';
            app.Mars.Visible ='off';
            app.Jupiter.Visible ='off';


             app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
            app.NaamModelEditField.Visible = 'off';
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.MassaEditFieldLabel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';
            app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off';
            app.OpstapmetsuperkrachtBlue.Visible = 'on';
            

        end

        % Callback function
        function MaaktekeningButtonPushed(app, event)

        end

        % Button pushed function: StartsimulatieButton
        function StartsimulatieButtonPushed(app, event)

            if isempty(app.GroepEditField)
                warning('Kies een groepsnaam.')
                return
            end

            if strcmp(app.CurrentMenu,'Kracht')
                allChecks = [app.CheckBox_Olaf.Value, app.CheckBox_Spongebob.Value, ...
                app.CheckBox_Jullie_strength.Value, app.CheckBox_MrImpossible.Value, ...
                 app.CheckBox_Maui.Value];
                modelNames = {'Olaf','SpongeBob','Jullie','MrImpossible','Maui'};
                sf_strength = [0.1 0.5 1 2 5];
                U.Force_sf = allChecks.*sf_strength;
                i_checks = logical(allChecks);
                U.ModelName = modelNames{i_checks};
                xpos_ballon_all = [265 508 693 1009 1336];
                ypos_ballon_all = [413 413 441 492 432];
                xpos_ballon_current = xpos_ballon_all(i_checks);
                ypos_ballon_current = ypos_ballon_all(i_checks);

                xpos_gif_all = [313 556 741 1057 1384];
                ypos_gif_all = [470 470 498 549 489];
                xpos_gif_current = xpos_gif_all(i_checks);
                ypos_gif_current = ypos_gif_all(i_checks);

                app.Ballon.Position = [xpos_ballon_current ypos_ballon_current 214 177];
                app.GIF.Position = [xpos_gif_current ypos_gif_current 98 71];
            
            else
                U.Force_sf = 1;
            end

              app.Ballon.Enable = 'on';
              app.Ballon.Visible = 'on';
              app.GIF.Enable = 'on';
              app.GIF.Visible = 'on';

            U.savefolder = app.path_savefolder;
            U.GroupName = app.GroupName;
            U.PathCasadi = app.path_casadi;
            U.Height = 180;
            U.Speed = 1.33;
            U.Mass = 62;
            U.sf_Gravity = 1;
            sf.foot = 1;
            sf.upp_leg = 1;
            sf.low_leg = 1;
            sf.torso = 1;
            sf.shoulder = 1;
            sf.low_arm = 1;
            sf.upp_arm = 1;

            % start simulation
            resultpath = PredSim_wrapper_for_app(U,sf);


              app.Ballon.Enable = 'off';
              app.Ballon.Visible = 'off';
              app.GIF.Enable = 'off';
              app.GIF.Visible = 'off';


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

        % Callback function
        function NaamEditFieldValueChanged(app, event)
            if isempty(app.NaamModelEditField.Value)
                return
            end
            ModelName_tmp = app.NaamModelEditField.Value;
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

        % Callback function
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

        % Callback function
        function KrachtSliderValueChanged(app, event)
            v = app.KrachtSlider.Value;
            if v < app.KrachtSlider.Limits(2)
                app.usr_speed = round(app.KrachtSlider.Value)/3.6;
            else
                app.usr_speed = -1;
            end
            
        end

        % Callback function
        function LichaamslengteEditFieldValueChanged(app, event)
            readUserInput(app)
            updateUserInput(app)
            writeDefaultUserInput(app)
            updateInputLimits(app)
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Color = [0.302 0.7451 0.9333];
            app.UIFigure.Position = [1 41 1545 763];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.WindowState = 'fullscreen';

            % Create LengteMenu
            app.LengteMenu = uimenu(app.UIFigure);
            app.LengteMenu.ForegroundColor = [0 0.4471 0.7412];
            app.LengteMenu.Text = '         Wat gebeurt er als je heel klein bent of juist heel groot?         ';

            % Create KrachtMenu
            app.KrachtMenu = uimenu(app.UIFigure);
            app.KrachtMenu.ForegroundColor = [0 0.4471 0.7412];
            app.KrachtMenu.Text = '         Wat als je spieren extra sterk of juist heel zwak zijn?          ';

            % Create ZwaartekrachtMenu
            app.ZwaartekrachtMenu = uimenu(app.UIFigure);
            app.ZwaartekrachtMenu.ForegroundColor = [0 0.4471 0.7412];
            app.ZwaartekrachtMenu.Text = '         Hoe beweeg je als de zwaartekracht anders is, zoals op de maan of op een andere planeet?         ';

            % Create CreerjezelfalsmodelMenu
            app.CreerjezelfalsmodelMenu = uimenu(app.UIFigure);
            app.CreerjezelfalsmodelMenu.ForegroundColor = [0 0.4471 0.7412];
            app.CreerjezelfalsmodelMenu.Text = 'Creëer jezelf  als model!';

            % Create Background_creerjeeigenmodel
            app.Background_creerjeeigenmodel = uiimage(app.UIFigure);
            app.Background_creerjeeigenmodel.Position = [-177 -67 1953 940];
            app.Background_creerjeeigenmodel.ImageSource = fullfile(pathToMLAPP, 'createyourownmodel_v2.png');

            % Create GroepEditFieldLabel
            app.GroepEditFieldLabel = uilabel(app.UIFigure);
            app.GroepEditFieldLabel.FontName = 'Footlight MT Light';
            app.GroepEditFieldLabel.FontSize = 30;
            app.GroepEditFieldLabel.FontColor = [1 0.9686 0.8588];
            app.GroepEditFieldLabel.Position = [123 588 83 42];
            app.GroepEditFieldLabel.Text = 'Groep';

            % Create PlaneetX
            app.PlaneetX = uiimage(app.UIFigure);
            app.PlaneetX.Position = [797 102 846 581];
            app.PlaneetX.ImageSource = fullfile(pathToMLAPP, 'PlanetX_purple.png');

            % Create LengteSliderReusSmurf
            app.LengteSliderReusSmurf = uislider(app.UIFigure, 'range');
            app.LengteSliderReusSmurf.MajorTicks = [0 50 100];
            app.LengteSliderReusSmurf.MajorTickLabels = {'Smurf', 'Mens', 'Reus'};
            app.LengteSliderReusSmurf.MinorTicks = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100];
            app.LengteSliderReusSmurf.Step = 50;
            app.LengteSliderReusSmurf.FontName = 'Footlight MT Light';
            app.LengteSliderReusSmurf.FontSize = 36;
            app.LengteSliderReusSmurf.FontColor = [0 0.4471 0.7412];
            app.LengteSliderReusSmurf.Position = [436 148 704 3];

            % Create PlaneetXLabel
            app.PlaneetXLabel = uilabel(app.UIFigure);
            app.PlaneetXLabel.FontName = 'Footlight MT Light';
            app.PlaneetXLabel.FontSize = 48;
            app.PlaneetXLabel.FontColor = [1 0.9686 0.8588];
            app.PlaneetXLabel.Position = [1137 413 187 63];
            app.PlaneetXLabel.Text = 'Planeet X';

            % Create gLabel
            app.gLabel = uilabel(app.UIFigure);
            app.gLabel.FontName = 'Footlight MT Light';
            app.gLabel.FontSize = 30;
            app.gLabel.FontColor = [1 0.9686 0.8588];
            app.gLabel.Position = [1338 370 33 42];
            app.gLabel.Text = 'g';

            % Create MaanLabel
            app.MaanLabel = uilabel(app.UIFigure);
            app.MaanLabel.FontName = 'Footlight MT Light';
            app.MaanLabel.FontSize = 48;
            app.MaanLabel.FontColor = [0.1098 0.1804 0.2];
            app.MaanLabel.Position = [312 471 119 63];
            app.MaanLabel.Text = 'Maan';

            % Create ZwaartekrachtEditFieldLabelMaan
            app.ZwaartekrachtEditFieldLabelMaan = uilabel(app.UIFigure);
            app.ZwaartekrachtEditFieldLabelMaan.FontName = 'Footlight MT Light';
            app.ZwaartekrachtEditFieldLabelMaan.FontSize = 30;
            app.ZwaartekrachtEditFieldLabelMaan.FontColor = [0.1098 0.1804 0.2];
            app.ZwaartekrachtEditFieldLabelMaan.Position = [210 435 186 42];
            app.ZwaartekrachtEditFieldLabelMaan.Text = 'Zwaartekracht';

            % Create gLabel_3
            app.gLabel_3 = uilabel(app.UIFigure);
            app.gLabel_3.FontName = 'Footlight MT Light';
            app.gLabel_3.FontSize = 30;
            app.gLabel_3.FontColor = [0.1098 0.1804 0.2];
            app.gLabel_3.Position = [499 433 33 42];
            app.gLabel_3.Text = 'g';

            % Create ZwaartekrachtEditFieldLabelPlanetX
            app.ZwaartekrachtEditFieldLabelPlanetX = uilabel(app.UIFigure);
            app.ZwaartekrachtEditFieldLabelPlanetX.FontName = 'Footlight MT Light';
            app.ZwaartekrachtEditFieldLabelPlanetX.FontSize = 30;
            app.ZwaartekrachtEditFieldLabelPlanetX.FontColor = [1 0.9686 0.8588];
            app.ZwaartekrachtEditFieldLabelPlanetX.Position = [1049 364 186 42];
            app.ZwaartekrachtEditFieldLabelPlanetX.Text = 'Zwaartekracht';

            % Create ZwaartekrachtEditFieldPlanetX
            app.ZwaartekrachtEditFieldPlanetX = uieditfield(app.UIFigure, 'numeric');
            app.ZwaartekrachtEditFieldPlanetX.Limits = [0.165 9];
            app.ZwaartekrachtEditFieldPlanetX.ValueDisplayFormat = '%111g';
            app.ZwaartekrachtEditFieldPlanetX.HorizontalAlignment = 'center';
            app.ZwaartekrachtEditFieldPlanetX.FontName = 'Footlight MT Light';
            app.ZwaartekrachtEditFieldPlanetX.FontSize = 30;
            app.ZwaartekrachtEditFieldPlanetX.FontColor = [0.1098 0.1804 0.2];
            app.ZwaartekrachtEditFieldPlanetX.BackgroundColor = [1 0.9686 0.8588];
            app.ZwaartekrachtEditFieldPlanetX.Position = [1247 364 63 41];
            app.ZwaartekrachtEditFieldPlanetX.Value = 1;

            % Create ZwaartekrachtEditFieldMaan
            app.ZwaartekrachtEditFieldMaan = uieditfield(app.UIFigure, 'numeric');
            app.ZwaartekrachtEditFieldMaan.Limits = [0.165 0.165];
            app.ZwaartekrachtEditFieldMaan.ValueDisplayFormat = '%111g';
            app.ZwaartekrachtEditFieldMaan.HorizontalAlignment = 'center';
            app.ZwaartekrachtEditFieldMaan.FontName = 'Footlight MT Light';
            app.ZwaartekrachtEditFieldMaan.FontSize = 30;
            app.ZwaartekrachtEditFieldMaan.FontColor = [0.1098 0.1804 0.2];
            app.ZwaartekrachtEditFieldMaan.BackgroundColor = [1 0.9686 0.8588];
            app.ZwaartekrachtEditFieldMaan.Position = [397 436 99 41];
            app.ZwaartekrachtEditFieldMaan.Value = 0.165;

            % Create giant
            app.giant = uiimage(app.UIFigure);
            app.giant.Position = [986 1 548 552];
            app.giant.ImageSource = fullfile(pathToMLAPP, 'Giant_nobackground.png');

            % Create smurf
            app.smurf = uiimage(app.UIFigure);
            app.smurf.Position = [219 61 237 204];
            app.smurf.ImageSource = fullfile(pathToMLAPP, 'Smurf_nobackground.png');

            % Create DENKreus
            app.DENKreus = uiimage(app.UIFigure);
            app.DENKreus.Enable = 'off';
            app.DENKreus.Visible = 'off';
            app.DENKreus.Position = [1316 480 214 177];
            app.DENKreus.ImageSource = 'balloon.png';

            % Create GIFreus
            app.GIFreus = uiimage(app.UIFigure);
            app.GIFreus.Enable = 'off';
            app.GIFreus.Visible = 'off';
            app.GIFreus.Position = [1369 533 98 71];
            app.GIFreus.ImageSource = 'allgifs.gif';

            % Create DENKsmurf
            app.DENKsmurf = uiimage(app.UIFigure);
            app.DENKsmurf.Enable = 'off';
            app.DENKsmurf.Visible = 'off';
            app.DENKsmurf.Position = [405 248 214 177];
            app.DENKsmurf.ImageSource = 'balloon.png';

            % Create GIFsmurf
            app.GIFsmurf = uiimage(app.UIFigure);
            app.GIFsmurf.Enable = 'off';
            app.GIFsmurf.Visible = 'off';
            app.GIFsmurf.Position = [455 301 98 71];
            app.GIFsmurf.ImageSource = 'allgifs.gif';

            % Create maan
            app.maan = uiimage(app.UIFigure);
            app.maan.Position = [143 306 457 377];
            app.maan.ImageSource = fullfile(pathToMLAPP, 'moon_drawing.png');

            % Create procentlabelEigenModel
            app.procentlabelEigenModel = uilabel(app.UIFigure);
            app.procentlabelEigenModel.FontName = 'Footlight MT Light';
            app.procentlabelEigenModel.FontSize = 30;
            app.procentlabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.procentlabelEigenModel.Position = [1111 304 27 42];
            app.procentlabelEigenModel.Text = '%';

            % Create SpierkrachtEigenModel
            app.SpierkrachtEigenModel = uieditfield(app.UIFigure, 'numeric');
            app.SpierkrachtEigenModel.Limits = [20 500];
            app.SpierkrachtEigenModel.ValueDisplayFormat = '%111g';
            app.SpierkrachtEigenModel.HorizontalAlignment = 'center';
            app.SpierkrachtEigenModel.FontName = 'Footlight MT Light';
            app.SpierkrachtEigenModel.FontSize = 30;
            app.SpierkrachtEigenModel.FontColor = [0 0.4471 0.7412];
            app.SpierkrachtEigenModel.BackgroundColor = [1 0.9725 0.8627];
            app.SpierkrachtEigenModel.Position = [1003 306 74 41];
            app.SpierkrachtEigenModel.Value = 100;

            % Create kgLabelEigenModel
            app.kgLabelEigenModel = uilabel(app.UIFigure);
            app.kgLabelEigenModel.FontName = 'Footlight MT Light';
            app.kgLabelEigenModel.FontSize = 30;
            app.kgLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.kgLabelEigenModel.Position = [1466 302 33 42];
            app.kgLabelEigenModel.Text = 'kg';

            % Create MassaEditFieldLabel
            app.MassaEditFieldLabel = uilabel(app.UIFigure);
            app.MassaEditFieldLabel.FontName = 'Footlight MT Light';
            app.MassaEditFieldLabel.FontSize = 30;
            app.MassaEditFieldLabel.FontColor = [0 0.4471 0.7412];
            app.MassaEditFieldLabel.Position = [1235 302 124 42];
            app.MassaEditFieldLabel.Text = 'Massa';

            % Create MassaEditFieldEigenModel
            app.MassaEditFieldEigenModel = uieditfield(app.UIFigure, 'numeric');
            app.MassaEditFieldEigenModel.Limits = [20 300];
            app.MassaEditFieldEigenModel.ValueDisplayFormat = '%111g';
            app.MassaEditFieldEigenModel.HorizontalAlignment = 'center';
            app.MassaEditFieldEigenModel.FontName = 'Footlight MT Light';
            app.MassaEditFieldEigenModel.FontSize = 30;
            app.MassaEditFieldEigenModel.FontColor = [0 0.4471 0.7412];
            app.MassaEditFieldEigenModel.BackgroundColor = [1 0.9725 0.8627];
            app.MassaEditFieldEigenModel.Position = [1370 303 74 41];
            app.MassaEditFieldEigenModel.Value = 75;

            % Create NaamModelEditField
            app.NaamModelEditField = uieditfield(app.UIFigure, 'text');
            app.NaamModelEditField.FontName = 'Footlight MT Light';
            app.NaamModelEditField.FontSize = 30;
            app.NaamModelEditField.BackgroundColor = [1 0.9686 0.8588];
            app.NaamModelEditField.Placeholder = '(Geef de naam van jouw groepje)';
            app.NaamModelEditField.Position = [814 469 677 44];
            app.NaamModelEditField.Value = '(Geef de naam van het model)';

            % Create PlanetXblack
            app.PlanetXblack = uiimage(app.UIFigure);
            app.PlanetXblack.Position = [953 79 244 244];
            app.PlanetXblack.ImageSource = fullfile(pathToMLAPP, 'PlanetX_dark.png');

            % Create ZwaartekrachtEigenModel
            app.ZwaartekrachtEigenModel = uieditfield(app.UIFigure, 'numeric');
            app.ZwaartekrachtEigenModel.Limits = [0.1 5];
            app.ZwaartekrachtEigenModel.ValueDisplayFormat = '%111g';
            app.ZwaartekrachtEigenModel.HorizontalAlignment = 'center';
            app.ZwaartekrachtEigenModel.FontName = 'Footlight MT Light';
            app.ZwaartekrachtEigenModel.FontSize = 30;
            app.ZwaartekrachtEigenModel.BackgroundColor = [1 0.9686 0.8588];
            app.ZwaartekrachtEigenModel.Position = [1025 224 74 41];
            app.ZwaartekrachtEigenModel.Value = 1;

            % Create gLabelEigenModel
            app.gLabelEigenModel = uilabel(app.UIFigure);
            app.gLabelEigenModel.HorizontalAlignment = 'center';
            app.gLabelEigenModel.FontName = 'Footlight MT Light';
            app.gLabelEigenModel.FontSize = 30;
            app.gLabelEigenModel.FontColor = [1 0.9686 0.8588];
            app.gLabelEigenModel.Position = [1101 227 43 42];
            app.gLabelEigenModel.Text = 'g';

            % Create LengteSliderEigenModel
            app.LengteSliderEigenModel = uislider(app.UIFigure, 'range');
            app.LengteSliderEigenModel.Limits = [80 180];
            app.LengteSliderEigenModel.MajorTicks = [80 100 120 140 160 180];
            app.LengteSliderEigenModel.MajorTickLabels = {'80', '100', '120', '140', '160', '180'};
            app.LengteSliderEigenModel.MinorTicks = [80 85 90.95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175180];
            app.LengteSliderEigenModel.Step = 1;
            app.LengteSliderEigenModel.FontName = 'Footlight MT Light';
            app.LengteSliderEigenModel.FontSize = 30;
            app.LengteSliderEigenModel.FontColor = [0 0.4471 0.7412];
            app.LengteSliderEigenModel.Position = [929 426 495 3];
            app.LengteSliderEigenModel.Value = [80 180];

            % Create LengteLabelEigenModel
            app.LengteLabelEigenModel = uilabel(app.UIFigure);
            app.LengteLabelEigenModel.FontName = 'Footlight MT Light';
            app.LengteLabelEigenModel.FontSize = 30;
            app.LengteLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.LengteLabelEigenModel.Position = [814 403 84 39];
            app.LengteLabelEigenModel.Text = 'Lengte';

            % Create cmLabelEigenModel
            app.cmLabelEigenModel = uilabel(app.UIFigure);
            app.cmLabelEigenModel.FontName = 'Footlight MT Light';
            app.cmLabelEigenModel.FontSize = 30;
            app.cmLabelEigenModel.Position = [1466 405 43 42];
            app.cmLabelEigenModel.Text = 'cm';

            % Create SpierkrachtLabelEigenModel
            app.SpierkrachtLabelEigenModel = uilabel(app.UIFigure);
            app.SpierkrachtLabelEigenModel.FontName = 'Footlight MT Light';
            app.SpierkrachtLabelEigenModel.FontSize = 30;
            app.SpierkrachtLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.SpierkrachtLabelEigenModel.Position = [813 308 146 39];
            app.SpierkrachtLabelEigenModel.Text = 'Spierkracht';

            % Create ZwaartekrachtLabelEigenModel
            app.ZwaartekrachtLabelEigenModel = uilabel(app.UIFigure);
            app.ZwaartekrachtLabelEigenModel.FontName = 'Footlight MT Light';
            app.ZwaartekrachtLabelEigenModel.FontSize = 30;
            app.ZwaartekrachtLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.ZwaartekrachtLabelEigenModel.Position = [810 224 186 39];
            app.ZwaartekrachtLabelEigenModel.Text = 'Zwaartekracht';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Enable = 'off';
            app.Image.Visible = 'off';
            app.Image.Position = [1313 -511 214 177];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'balloon_reflected.png');

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.Enable = 'off';
            app.Image2.Visible = 'off';
            app.Image2.Position = [1369 -469 98 71];
            app.Image2.ImageSource = 'allgifs.gif';

            % Create spaceBackground
            app.spaceBackground = uiimage(app.UIFigure);
            app.spaceBackground.Position = [-229 -188 2004 1032];
            app.spaceBackground.ImageSource = fullfile(pathToMLAPP, 'Space_drawing.png');

            % Create Spaceman
            app.Spaceman = uiimage(app.UIFigure);
            app.Spaceman.Position = [482 39 591 505];
            app.Spaceman.ImageSource = fullfile(pathToMLAPP, 'Spaceman_nobackground.png');

            % Create Jupiter
            app.Jupiter = uiimage(app.UIFigure);
            app.Jupiter.Position = [19 -46 377 393];
            app.Jupiter.ImageSource = fullfile(pathToMLAPP, 'Jupiter_noBackground.png');

            % Create Mars
            app.Mars = uiimage(app.UIFigure);
            app.Mars.Position = [958 379 317 317];
            app.Mars.ImageSource = fullfile(pathToMLAPP, 'Mars.png');

            % Create earth
            app.earth = uiimage(app.UIFigure);
            app.earth.Position = [1196 10 306 306];
            app.earth.ImageSource = fullfile(pathToMLAPP, 'Earth_nobackground.png');

            % Create BackgroundImageSmurfReusKracht
            app.BackgroundImageSmurfReusKracht = uiimage(app.UIFigure);
            app.BackgroundImageSmurfReusKracht.ScaleMethod = 'fill';
            app.BackgroundImageSmurfReusKracht.HandleVisibility = 'callback';
            app.BackgroundImageSmurfReusKracht.Position = [-15 -19 1586 844];
            app.BackgroundImageSmurfReusKracht.ImageSource = fullfile(pathToMLAPP, 'Background_landscape_path.png');

            % Create OpstapmetsuperkrachtBlue
            app.OpstapmetsuperkrachtBlue = uilabel(app.UIFigure);
            app.OpstapmetsuperkrachtBlue.HorizontalAlignment = 'center';
            app.OpstapmetsuperkrachtBlue.FontName = 'Footlight MT Light';
            app.OpstapmetsuperkrachtBlue.FontSize = 80;
            app.OpstapmetsuperkrachtBlue.FontColor = [0 0.4471 0.7412];
            app.OpstapmetsuperkrachtBlue.Position = [378 595 858 197];
            app.OpstapmetsuperkrachtBlue.Text = 'Op stap met superkracht!';

            % Create OpstapmetsuperkrachtWhite
            app.OpstapmetsuperkrachtWhite = uilabel(app.UIFigure);
            app.OpstapmetsuperkrachtWhite.HorizontalAlignment = 'center';
            app.OpstapmetsuperkrachtWhite.FontName = 'Footlight MT Light';
            app.OpstapmetsuperkrachtWhite.FontSize = 80;
            app.OpstapmetsuperkrachtWhite.FontColor = [1 0.9686 0.8588];
            app.OpstapmetsuperkrachtWhite.Position = [377 595 858 197];
            app.OpstapmetsuperkrachtWhite.Text = 'Op stap met superkracht!';

            % Create StartsimulatieButton
            app.StartsimulatieButton = uibutton(app.UIFigure, 'push');
            app.StartsimulatieButton.ButtonPushedFcn = createCallbackFcn(app, @StartsimulatieButtonPushed, true);
            app.StartsimulatieButton.IconAlignment = 'center';
            app.StartsimulatieButton.BackgroundColor = [1 0.9686 0.8588];
            app.StartsimulatieButton.FontName = 'Footlight MT Light';
            app.StartsimulatieButton.FontSize = 30;
            app.StartsimulatieButton.FontColor = [0 0.451 0.7412];
            app.StartsimulatieButton.Position = [673 569 314 61];
            app.StartsimulatieButton.Text = 'Start simulatie';

            % Create SluitvideoButton
            app.SluitvideoButton = uibutton(app.UIFigure, 'push');
            app.SluitvideoButton.ButtonPushedFcn = createCallbackFcn(app, @SluitvideoButtonPushed, true);
            app.SluitvideoButton.BackgroundColor = [1 0.9686 0.8588];
            app.SluitvideoButton.FontName = 'Footlight MT Light';
            app.SluitvideoButton.FontSize = 30;
            app.SluitvideoButton.FontColor = [0 0.4471 0.7412];
            app.SluitvideoButton.Enable = 'off';
            app.SluitvideoButton.Visible = 'off';
            app.SluitvideoButton.Position = [1056 569 314 61];
            app.SluitvideoButton.Text = 'Sluit video';

            % Create SpeelvideoButton
            app.SpeelvideoButton = uibutton(app.UIFigure, 'push');
            app.SpeelvideoButton.ButtonPushedFcn = createCallbackFcn(app, @SpeelvideoButtonPushed, true);
            app.SpeelvideoButton.BackgroundColor = [1 0.9686 0.8588];
            app.SpeelvideoButton.FontName = 'Footlight MT Light';
            app.SpeelvideoButton.FontSize = 30;
            app.SpeelvideoButton.FontColor = [0 0.4471 0.7412];
            app.SpeelvideoButton.Position = [1056 569 314 61];
            app.SpeelvideoButton.Text = 'Speel video';

            % Create GroepEditField
            app.GroepEditField = uieditfield(app.UIFigure, 'text');
            app.GroepEditField.ValueChangedFcn = createCallbackFcn(app, @GroepEditFieldValueChanged, true);
            app.GroepEditField.FontName = 'Footlight MT Light';
            app.GroepEditField.FontSize = 30;
            app.GroepEditField.FontColor = [0 0.4471 0.7412];
            app.GroepEditField.BackgroundColor = [1 0.9725 0.8627];
            app.GroepEditField.Placeholder = '(Geef de naam van jouw groepje)';
            app.GroepEditField.Position = [207 569 412 61];

            % Create Olaf
            app.Olaf = uiimage(app.UIFigure);
            app.Olaf.Position = [53 172 307 307];
            app.Olaf.ImageSource = fullfile(pathToMLAPP, 'Olaf_nobackground.png');

            % Create Spongebob
            app.Spongebob = uiimage(app.UIFigure);
            app.Spongebob.Position = [330 176 247 247];
            app.Spongebob.ImageSource = fullfile(pathToMLAPP, 'Spongebob.png');

            % Create Jullie_strength
            app.Jullie_strength = uiimage(app.UIFigure);
            app.Jullie_strength.Position = [546 147 323 323];
            app.Jullie_strength.ImageSource = fullfile(pathToMLAPP, 'Jullie_nobackground.png');

            % Create MrIncredible
            app.MrIncredible = uiimage(app.UIFigure);
            app.MrIncredible.Position = [789 139 377 377];
            app.MrIncredible.ImageSource = fullfile(pathToMLAPP, 'incredibles2-mrincredible2.png');

            % Create Maui
            app.Maui = uiimage(app.UIFigure);
            app.Maui.Position = [1079 84 432 472];
            app.Maui.ImageSource = fullfile(pathToMLAPP, 'maui_nobackground_v2.png');

            % Create Button
            app.Button = uibutton(app.UIFigure, 'push');
            app.Button.Icon = fullfile(pathToMLAPP, 'Olaf_nobackground.png');
            app.Button.Position = [-465 -1 330 594];
            app.Button.Text = '';

            % Create CheckBox_Olaf
            app.CheckBox_Olaf = uicheckbox(app.UIFigure);
            app.CheckBox_Olaf.Text = '';
            app.CheckBox_Olaf.Position = [197 100 77 62];

            % Create CheckBox_Spongebob
            app.CheckBox_Spongebob = uicheckbox(app.UIFigure);
            app.CheckBox_Spongebob.Text = '';
            app.CheckBox_Spongebob.Position = [405 100 77 62];

            % Create CheckBox_Jullie_strength
            app.CheckBox_Jullie_strength = uicheckbox(app.UIFigure);
            app.CheckBox_Jullie_strength.Text = '';
            app.CheckBox_Jullie_strength.Position = [676 100 77 62];

            % Create CheckBox_MrImpossible
            app.CheckBox_MrImpossible = uicheckbox(app.UIFigure);
            app.CheckBox_MrImpossible.Text = '';
            app.CheckBox_MrImpossible.Position = [958 100 77 62];

            % Create CheckBox_Maui
            app.CheckBox_Maui = uicheckbox(app.UIFigure);
            app.CheckBox_Maui.Text = '';
            app.CheckBox_Maui.Position = [1273 100 77 62];

            % Create Ballon
            app.Ballon = uiimage(app.UIFigure);
            app.Ballon.Enable = 'off';
            app.Ballon.Visible = 'off';
            app.Ballon.Position = [265 413 214 177];
            app.Ballon.ImageSource = 'balloon.png';

            % Create GIF
            app.GIF = uiimage(app.UIFigure);
            app.GIF.Enable = 'off';
            app.GIF.Visible = 'off';
            app.GIF.Position = [313 470 98 71];
            app.GIF.ImageSource = 'allgifs.gif';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = OpStapMetSuperKracht2025_exported

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