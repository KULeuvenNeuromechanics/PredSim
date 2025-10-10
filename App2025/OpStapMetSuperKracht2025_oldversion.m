classdef OpStapMetSuperKracht2025_FINAL_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        LengteMenu                      matlab.ui.container.Menu
        KrachtMenu                      matlab.ui.container.Menu
        ZwaartekrachtMenu               matlab.ui.container.Menu
        CreerjezelfalsmodelMenu         matlab.ui.container.Menu
        ProbeeropnieuwButton            matlab.ui.control.Button
        SimulatieNietGeluktWarning      matlab.ui.control.Label
        LengteSliderEigenModel          matlab.ui.control.Slider
        gLabelEigenModel                matlab.ui.control.Label
        CheckBox_Reus                   matlab.ui.control.CheckBox
        smurf                           matlab.ui.control.Image
        Kleuters                        matlab.ui.control.Image
        CheckBox_Volwassenen            matlab.ui.control.CheckBox
        CheckBox_Jullie_lengte          matlab.ui.control.CheckBox
        CheckBox_Kleuter                matlab.ui.control.CheckBox
        CheckBox_Smurf                  matlab.ui.control.CheckBox
        CheckBoxJupiter                 matlab.ui.control.CheckBox
        GIF                             matlab.ui.control.Image
        Ballon                          matlab.ui.control.Image
        GroepEditField                  matlab.ui.control.EditField
        GroepEditFieldLabel             matlab.ui.control.Label
        CheckBox_Maui                   matlab.ui.control.CheckBox
        CheckBox_MrImpossible           matlab.ui.control.CheckBox
        CheckBox_Jullie_strength        matlab.ui.control.CheckBox
        CheckBox_Spongebob              matlab.ui.control.CheckBox
        CheckBox_Olaf                   matlab.ui.control.CheckBox
        Spongebob                       matlab.ui.control.Image
        Olaf                            matlab.ui.control.Image
        SpeelvideoButton                matlab.ui.control.Button
        SluitvideoButton                matlab.ui.control.Button
        StartsimulatieButton            matlab.ui.control.Button
        OpstapmetsuperkrachtWhite       matlab.ui.control.Label
        OpstapmetsuperkrachtBlue        matlab.ui.control.Label
        CheckBoxEarth                   matlab.ui.control.CheckBox
        CheckBoxMars                    matlab.ui.control.CheckBox
        CheckBoxMaan                    matlab.ui.control.CheckBox
        earth                           matlab.ui.control.Image
        Jupiter                         matlab.ui.control.Image
        maan                            matlab.ui.control.Image
        CheckBoxPluto                   matlab.ui.control.CheckBox
        Image2                          matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        giant                           matlab.ui.control.Image
        ZwaartekrachtLabelEigenModel    matlab.ui.control.Label
        ZwaartekrachtEditFieldEigenModel  matlab.ui.control.NumericEditField
        Jullie_lengte                   matlab.ui.control.Image
        Mars                            matlab.ui.control.Image
        Volwassenen                     matlab.ui.control.Image
        MrIncredible                    matlab.ui.control.Image
        Jullie_strength                 matlab.ui.control.Image
        Spaceman                        matlab.ui.control.Image
        PlanetXblack                    matlab.ui.control.Image
        SpierkrachtLabelEigenModel      matlab.ui.control.Label
        GewichtLabelEigenModel          matlab.ui.control.Label
        cmLabelEigenModel               matlab.ui.control.Label
        LengteLabelEigenModel           matlab.ui.control.Label
        ZwaartekrachtEigenModel         matlab.ui.control.NumericEditField
        kgLabelEigenModel               matlab.ui.control.Label
        SpierkrachtEigenModel           matlab.ui.control.NumericEditField
        procentlabelEigenModel          matlab.ui.control.Label
        GIFsmurf                        matlab.ui.control.Image
        DENKsmurf                       matlab.ui.control.Image
        GIFreus                         matlab.ui.control.Image
        DENKreus                        matlab.ui.control.Image
        PlaneetX                        matlab.ui.control.Image
        Maui                            matlab.ui.control.Image
        NaamModelEditField              matlab.ui.control.EditField
        ModelEditFieldLabel             matlab.ui.control.Label
        MassaEditFieldEigenModel        matlab.ui.control.NumericEditField
        Background_creerjeeigenmodel    matlab.ui.control.Image
        spaceBackground                 matlab.ui.control.Image
        BackgroundImageSmurfReusKracht  matlab.ui.control.Image
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
            app.CheckBox_Reus.Visible = 'on';
            app.smurf.Visible = 'on';
            app.CheckBox_Smurf.Visible = 'on';
            app.Kleuters.Visible = 'on';
            app.CheckBox_Kleuter.Visible = 'on';
            app.Jullie_lengte.Visible = 'on';
            app.CheckBox_Jullie_lengte.Visible = 'on';
            app.Volwassenen.Visible = 'on';
            app.CheckBox_Volwassenen.Visible = 'on';
            %app.LengteSliderReusSmurf.Visible = 'on';
            app.GroepEditField.Visible = 'on';

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
            app.CheckBoxMaan.Visible = 'off';
            app.earth.Visible ='off';
            app.CheckBoxEarth.Visible = 'off';
            app.Mars.Visible ='off';
            app.CheckBoxMars.Visible = 'off';
            app.Jupiter.Visible ='off';
            app.CheckBoxJupiter.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.CheckBoxPluto.Visible = 'off'; 
            
            % Eigen model
            app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.GewichtLabelEigenModel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';

            app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off';   
            app.GewichtLabelEigenModel.Visible = 'off';
            app.ZwaartekrachtEditFieldEigenModel.Visible = 'off';
            app.NaamModelEditField.Visible = 'off';
            app.ModelEditFieldLabel.Visible = 'off';
            
        
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
            app.GroepEditField.Visible = 'on';

            % related to smurf reus
            app.giant.Visible = 'off';
            app.CheckBox_Reus.Visible = 'off';
            app.smurf.Visible = 'off';
            app.CheckBox_Smurf.Visible = 'off';
            app.Kleuters.Visible = 'off';
            app.CheckBox_Kleuter.Visible = 'off';
            app.Jullie_lengte.Visible = 'off';
            app.CheckBox_Jullie_lengte.Visible = 'off';
            app.Volwassenen.Visible = 'off';
            app.CheckBox_Volwassenen.Visible = 'off';
            %app.LengteSliderReusSmurf.Visible = 'on';


            % space
            app.Spaceman.Visible = 'off';
            app.spaceBackground.Visible = 'off';
            app.maan.Visible= 'off'; 
            app.CheckBoxMaan.Visible = 'off';
            app.earth.Visible ='off';
            app.CheckBoxEarth.Visible = 'off';
            app.Mars.Visible ='off';
            app.CheckBoxMars.Visible = 'off';
            app.Jupiter.Visible ='off';
            app.CheckBoxJupiter.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.CheckBoxPluto.Visible = 'off'; 
            
            % Eigen model
            app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
           
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.GewichtLabelEigenModel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';

            app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off'; 
            app.GewichtLabelEigenModel.Visible = 'off';
            app.ZwaartekrachtEditFieldEigenModel.Visible = 'off';
             app.NaamModelEditField.Visible = 'off';
             app.ModelEditFieldLabel.Visible = 'off';

    

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
            app.CheckBox_Reus.Visible = 'off';
            app.smurf.Visible = 'off';
            app.CheckBox_Smurf.Visible = 'off';
            app.Kleuters.Visible = 'off';
            app.CheckBox_Kleuter.Visible = 'off';
            app.Jullie_lengte.Visible = 'off';
            app.CheckBox_Jullie_lengte.Visible = 'off';
            app.Volwassenen.Visible = 'off';
            app.CheckBox_Volwassenen.Visible = 'off';
            %app.LengteSliderReusSmurf.Visible = 'on';

            % Space
            app.Spaceman.Visible = 'on';
            app.spaceBackground.Visible = 'on';
            app.maan.Visible= 'on'; 
            app.CheckBoxMaan.Visible = 'on';
            app.earth.Visible ='on';
            app.CheckBoxEarth.Visible = 'on';
            app.Mars.Visible ='on';
            app.CheckBoxMars.Visible = 'on';
            app.Jupiter.Visible ='on';
            app.CheckBoxJupiter.Visible = 'on';
            app.PlaneetX.Visible = 'on';
            app.CheckBoxPluto.Visible = 'on'; 
            app.GroepEditField.Visible = 'on';

            % Eigen model
            app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
        
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.GewichtLabelEigenModel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';


            app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off';
            app.GewichtLabelEigenModel.Visible = 'off';
            app.ZwaartekrachtEditFieldEigenModel.Visible = 'off';
             app.NaamModelEditField.Visible = 'off';
             app.ModelEditFieldLabel.Visible = 'off';

    

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
            app.CheckBox_Reus.Visible = 'off';
            app.smurf.Visible = 'off';
            app.CheckBox_Smurf.Visible = 'off';
            app.Kleuters.Visible = 'off';
            app.CheckBox_Kleuter.Visible = 'off';
            app.Jullie_lengte.Visible = 'off';
            app.CheckBox_Jullie_lengte.Visible = 'off';
            app.Volwassenen.Visible = 'off';
            app.CheckBox_Volwassenen.Visible = 'off';
            %app.LengteSliderReusSmurf.Visible = 'on';

            % space eigen model
            app.OpstapmetsuperkrachtWhite.Visible = 'on';
            
            % Space
            app.Spaceman.Visible = 'off';
            app.spaceBackground.Visible = 'off';
            app.maan.Visible= 'off'; 
            app.CheckBoxMaan.Visible = 'off';
            app.earth.Visible ='off';
            app.CheckBoxEarth.Visible = 'off';
            app.Mars.Visible ='off';
            app.CheckBoxMars.Visible = 'off';
            app.Jupiter.Visible ='off';
            app.CheckBoxJupiter.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.CheckBoxPluto.Visible = 'off'; 

            % Eigen model    
            app.cmLabelEigenModel.Visible = 'on';
            app.LengteLabelEigenModel.Visible = 'on';
            app.LengteSliderEigenModel.Visible = 'on';
            app.gLabelEigenModel.Visible = 'on';
            app.PlanetXblack.Visible = 'on';
            app.GroepEditField.Visible = 'on';
            app.MassaEditFieldEigenModel.Visible = 'on';
            app.GewichtLabelEigenModel.Visible = 'on';
            app.kgLabelEigenModel.Visible = 'on';
            app.procentlabelEigenModel.Visible = 'on';
            app.SpierkrachtEigenModel.Visible = 'on';
            app.Background_creerjeeigenmodel.Visible = 'on';
            app.ZwaartekrachtEigenModel.Visible = 'on';
            app.ZwaartekrachtLabelEigenModel.Visible = 'on';
            app.SpierkrachtLabelEigenModel.Visible = 'on';
            app.GewichtLabelEigenModel.Visible = 'on';
            app.ZwaartekrachtEditFieldEigenModel.Visible = 'on';
             app.NaamModelEditField.Visible = 'on';
             app.ModelEditFieldLabel.Visible = 'on';



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

 

         function initPaths(app)
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
            % make default drawing

            % path repo
            [pathApp,~,~] = fileparts(mfilename('fullpath'));
            [pathRepo,~,~] = fileparts(pathApp);
            app.path_repo = pathRepo;

            addpath([pathRepo '\VariousFunctions'])
            addpath([pathRepo '\App2025'])

            % app.sel_mot_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.mot');
            % 
            % app.sel_osim_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.osim');
            % 

            % setPaths(app)
            initPaths(app)

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
            app.GroepEditField.Visible = 'on';
              
            
            app.OpstapmetsuperkrachtWhite.Visible = 'off';
            app.Spaceman.Visible = 'off';
            app.spaceBackground.Visible = 'off';
            app.maan.Visible= 'off'; 
            app.CheckBoxMaan.Visible = 'off';
            app.earth.Visible ='off';
            app.CheckBoxEarth.Visible = 'off';
            app.Mars.Visible ='off';
            app.CheckBoxMars.Visible = 'off';
            app.Jupiter.Visible ='off';
            app.CheckBoxJupiter.Visible = 'off';
            app.PlaneetX.Visible = 'off';
            app.CheckBoxPluto.Visible = 'off'; 


            app.cmLabelEigenModel.Visible = 'off';
            app.LengteLabelEigenModel.Visible = 'off';
            app.LengteSliderEigenModel.Visible = 'off';
            app.gLabelEigenModel.Visible = 'off';
            app.PlanetXblack.Visible = 'off';
            app.GroepEditField.Visible = 'off';
            app.MassaEditFieldEigenModel.Visible = 'off';
            app.GewichtLabelEigenModel.Visible = 'off';
            app.kgLabelEigenModel.Visible = 'off';
            app.procentlabelEigenModel.Visible = 'off';
            app.SpierkrachtEigenModel.Visible = 'off';
            app.Background_creerjeeigenmodel.Visible = 'off';
            app.ZwaartekrachtEigenModel.Visible = 'off';
            app.ZwaartekrachtLabelEigenModel.Visible = 'off';
            app.SpierkrachtLabelEigenModel.Visible = 'off';
            app.OpstapmetsuperkrachtBlue.Visible = 'on';
            app.GewichtLabelEigenModel.Visible = 'off';
            app.ZwaartekrachtEditFieldEigenModel.Visible = 'off';
            app.NaamModelEditField.Visible = 'off';
            app.ModelEditFieldLabel.Visible = 'off';



            

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

            try
            if strcmp(app.CurrentMenu,'Kracht')
                allChecks = [app.CheckBox_Olaf.Value, app.CheckBox_Spongebob.Value, ...
                app.CheckBox_Jullie_strength.Value, app.CheckBox_MrImpossible.Value, ...
                 app.CheckBox_Maui.Value];
                modelNames = {'Olaf','SpongeBob','Jullie','MrImpossible','Maui'};
                sf_strength = [0.1 0.5 1 2 5];
                % i_checks = logical(allChecks);
                U.Force_sf = sum(allChecks.*sf_strength);
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
            
            if strcmp(app.CurrentMenu,'Zwaartekracht')
                allChecks = [app.CheckBoxMaan.Value, app.CheckBoxMars.Value, ...
                app.CheckBoxEarth.Value, app.CheckBoxJupiter.Value, app.CheckBoxPluto.Value];
                modelNames = {'Maan','Mars','Aarde','Jupiter','Pluto'};
                sf_gravity = [0.17 0.38 1 2.53 0.06]; % units = g
                % i_checks = logical(allChecks);
                U.sf_Gravity = sum(allChecks.*sf_gravity);
                i_checks = logical(allChecks);
                U.ModelName = modelNames{i_checks};
                % xpos_ballon_all = [265 508 693 1009 1336];
                % ypos_ballon_all = [413 413 441 492 432];
                % xpos_ballon_current = xpos_ballon_all(i_checks);
                % ypos_ballon_current = ypos_ballon_all(i_checks);
                % 
                % xpos_gif_all = [313 556 741 1057 1384];
                % ypos_gif_all = [470 470 498 549 489];
                % xpos_gif_current = xpos_gif_all(i_checks);
                % ypos_gif_current = ypos_gif_all(i_checks);

                app.Ballon.Position = [693 441 214 177];
                app.GIF.Position = [741 498 98 71];
            
            else
                U.sf_Gravity = 1;
            end

            if strcmp(app.CurrentMenu,'Lengte')
                allChecks = [app.CheckBox_Smurf.Value, app.CheckBox_Kleuter.Value, ...
                app.CheckBox_Jullie_lengte.Value, app.CheckBox_Volwassenen.Value, ...
                 app.CheckBox_Reus.Value];
                modelNames = {'Smurf','Kleuter','Jullie_lengte','Volwassenen','Reus'};
                v_Length = [20 80 120 160 400];
                v_Mass = [2 12 35 62 200];
                v_Speed = [0.20 0.53 0.85 1.20 2.7];
                %v_IG_pelvis_y = [0.115 0.460 0.690 1.02 2.3];
                %i_checks = logical(allChecks);
                U.Height = sum(allChecks.*v_Length);
                U.Mass = sum(allChecks.*v_Mass);
                U.Speed = sum(allChecks.*v_Speed);
                
                %U.IG_pelvis_y = sum(allChecks.*v_IG_pelvis_y);
                i_checks = logical(allChecks);
                U.ModelName = modelNames{i_checks};
                xpos_ballon_all = [265 508 693 1009 1387];
                ypos_ballon_all = [413 413 441 492 567];
                xpos_ballon_current = xpos_ballon_all(i_checks);
                ypos_ballon_current = ypos_ballon_all(i_checks);

                xpos_gif_all = [313 556 741 1057 1445];
                ypos_gif_all = [470 470 498 549 620];
                xpos_gif_current = xpos_gif_all(i_checks);
                ypos_gif_current = ypos_gif_all(i_checks);

                app.Ballon.Position = [xpos_ballon_current ypos_ballon_current 214 177];
                app.GIF.Position = [xpos_gif_current ypos_gif_current 98 71];
            
            else
                U.Height = 160;
                U.Speed = 1.33;
                U.Mass = 62;
                %U.IG_pelvis_y = [];
                
            end

            if strcmp(app.CurrentMenu,'EigenModel')
                Length_val = app.LengteSliderEigenModel.Value;
                Strength_val = app.SpierkrachtEigenModel.Value;
                Mass_val = app.MassaEditFieldEigenModel.Value;
                g_val = app.ZwaartekrachtEditFieldEigenModel.Value;


                U.Height = Length_val;
                U.Mass = Mass_val;
                U.Force_sf = Strength_val/100;
                U.sf_Gravity = g_val;
                U.Speed = 1.33;


                U.ModelName = app.NaamModelEditField.Value;

                app.Ballon.Position = [576 511 214 177];
                app.GIF.Position = [627 569 98 71];

            else

            end

              app.Ballon.Enable = 'on';
              app.Ballon.Visible = 'on';
              app.GIF.Enable = 'on';
              app.GIF.Visible = 'on';

            U.savefolder = app.path_savefolder;
            U.GroupName = app.GroupName;
            U.PathCasadi = app.path_casadi;
    
            sf_gen = U.Height/160;
            sf.foot = sf_gen;
            sf.upp_leg = sf_gen;
            sf.low_leg = sf_gen;
            sf.torso = sf_gen;
            sf.shoulder = sf_gen;
            sf.low_arm = sf_gen;
            sf.upp_arm = sf_gen;

            % start simulation
            if strcmp(app.CurrentMenu,'Lengte') && strcmp(U.ModelName ,'Reus') ||...
                    strcmp(app.CurrentMenu,'EigenModel')
                resultpath = PredSim_wrapper_for_app_Old(U,sf);
            else
                resultpath = PredSim_wrapper_for_app(U,sf);
            end

            catch ME

            if strcmp(app.CurrentMenu,'EigenModel')
            app.ProbeeropnieuwButton.Visible = 'on';  
            app.SimulatieNietGeluktWarning.Visible = 'on';  
            app.SimulatieNietGeluktWarning.Enable = 'on';  
            app.ProbeeropnieuwButton.Enable = 'on';  
            end


            app.Ballon.Enable = 'off';
            app.Ballon.Visible = 'off';
            app.GIF.Enable = 'off';
            app.GIF.Visible = 'off';     
            end


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
            % loadResultsTable(app);
            
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
            app.ModelName = app.ModelName_tmp;

        end

        % Button pushed function: SpeelvideoButton
        function SpeelvideoButtonPushed(app, event)
            if isempty(app.ModelName)
            app.ModelName = app.NaamModelEditField.Value;
            end
            if isempty(app.GroupName)
            app.GroupName = app.GroepEditField.Value;
            end
             flag = 0;
             app.sel_osim_file = fullfile(app.path_repo,'Subjects',app.GroupName,app.ModelName,[app.ModelName '.osim']);
             mot_folder = fullfile(app.path_savefolder ,app.GroupName);
                files = dir(fullfile(mot_folder, [app.ModelName '_v*.mot']));
                if isempty(files), disp('No mot files found.'); end
                

                nums = cellfun(@(s) sscanf(s, [app.ModelName '_v%d.mot']), {files.name});
                [~, idx] = max(nums);

                app.sel_mot_file = fullfile(mot_folder, files(idx).name);
                
            if ~exist(app.sel_osim_file,'file')
                flag = 1;
                disp(['Could not find "' app.sel_osim_file '"'])
                app.sel_osim_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.osim');
                app.sel_mot_file = fullfile(app.path_repo,'Subjects','Vitruvian_Man','Vitruvian_Man.mot');

            end
            if ~exist(app.sel_mot_file,'file')
                flag = 1;
                disp(['Could not find "' app.sel_mot_file '"'])
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

        % Button pushed function: ProbeeropnieuwButton
        function ProbeeropnieuwButtonPushed(app, event)
           app.ProbeeropnieuwButton.Visible = 'off';  
            app.SimulatieNietGeluktWarning.Visible = 'off';  
            app.SimulatieNietGeluktWarning.Enable = 'off';  
            app.ProbeeropnieuwButton.Enable = 'off';  
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

            % Create BackgroundImageSmurfReusKracht
            app.BackgroundImageSmurfReusKracht = uiimage(app.UIFigure);
            app.BackgroundImageSmurfReusKracht.ScaleMethod = 'fill';
            app.BackgroundImageSmurfReusKracht.HandleVisibility = 'callback';
            app.BackgroundImageSmurfReusKracht.Position = [-15 -19 1586 844];
            app.BackgroundImageSmurfReusKracht.ImageSource = fullfile(pathToMLAPP, 'Background_landscape_path.png');

            % Create spaceBackground
            app.spaceBackground = uiimage(app.UIFigure);
            app.spaceBackground.Position = [-229 -188 2004 1032];
            app.spaceBackground.ImageSource = fullfile(pathToMLAPP, 'Space_drawing.png');

            % Create Background_creerjeeigenmodel
            app.Background_creerjeeigenmodel = uiimage(app.UIFigure);
            app.Background_creerjeeigenmodel.Position = [-177 -67 1953 940];
            app.Background_creerjeeigenmodel.ImageSource = fullfile(pathToMLAPP, 'createyourownmodel_v2.png');

            % Create MassaEditFieldEigenModel
            app.MassaEditFieldEigenModel = uieditfield(app.UIFigure, 'numeric');
            app.MassaEditFieldEigenModel.Limits = [20 300];
            app.MassaEditFieldEigenModel.ValueDisplayFormat = '%111g';
            app.MassaEditFieldEigenModel.HorizontalAlignment = 'center';
            app.MassaEditFieldEigenModel.FontName = 'Footlight MT Light';
            app.MassaEditFieldEigenModel.FontSize = 30;
            app.MassaEditFieldEigenModel.FontColor = [0 0.451 0.7412];
            app.MassaEditFieldEigenModel.BackgroundColor = [1 0.9686 0.8588];
            app.MassaEditFieldEigenModel.Position = [1370 303 74 41];
            app.MassaEditFieldEigenModel.Value = 75;

            % Create ModelEditFieldLabel
            app.ModelEditFieldLabel = uilabel(app.UIFigure);
            app.ModelEditFieldLabel.FontName = 'Footlight MT Light';
            app.ModelEditFieldLabel.FontSize = 30;
            app.ModelEditFieldLabel.FontColor = [0 0.451 0.7412];
            app.ModelEditFieldLabel.Position = [813 466 83 42];
            app.ModelEditFieldLabel.Text = 'Model';

            % Create NaamModelEditField
            app.NaamModelEditField = uieditfield(app.UIFigure, 'text');
            app.NaamModelEditField.HorizontalAlignment = 'center';
            app.NaamModelEditField.FontName = 'Footlight MT Light';
            app.NaamModelEditField.FontSize = 30;
            app.NaamModelEditField.FontColor = [0 0.451 0.7412];
            app.NaamModelEditField.BackgroundColor = [1 0.9686 0.8588];
            app.NaamModelEditField.Placeholder = '(Geef de naam van het model)';
            app.NaamModelEditField.Position = [936 463 555 61];

            % Create Maui
            app.Maui = uiimage(app.UIFigure);
            app.Maui.Position = [1079 84 432 472];
            app.Maui.ImageSource = fullfile(pathToMLAPP, 'maui_nobackground_v2.png');

            % Create PlaneetX
            app.PlaneetX = uiimage(app.UIFigure);
            app.PlaneetX.Position = [1096 261 218 203];
            app.PlaneetX.ImageSource = fullfile(pathToMLAPP, 'PlanetX_purple.png');

            % Create DENKreus
            app.DENKreus = uiimage(app.UIFigure);
            app.DENKreus.Enable = 'off';
            app.DENKreus.Visible = 'off';
            app.DENKreus.Position = [1387 567 214 177];
            app.DENKreus.ImageSource = 'balloon.png';

            % Create GIFreus
            app.GIFreus = uiimage(app.UIFigure);
            app.GIFreus.Enable = 'off';
            app.GIFreus.Visible = 'off';
            app.GIFreus.Position = [1445 620 98 71];
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

            % Create LengteLabelEigenModel
            app.LengteLabelEigenModel = uilabel(app.UIFigure);
            app.LengteLabelEigenModel.FontName = 'Footlight MT Light';
            app.LengteLabelEigenModel.FontSize = 30;
            app.LengteLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.LengteLabelEigenModel.Position = [813 399 84 39];
            app.LengteLabelEigenModel.Text = 'Lengte';

            % Create cmLabelEigenModel
            app.cmLabelEigenModel = uilabel(app.UIFigure);
            app.cmLabelEigenModel.FontName = 'Footlight MT Light';
            app.cmLabelEigenModel.FontSize = 30;
            app.cmLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.cmLabelEigenModel.Position = [1466 405 43 42];
            app.cmLabelEigenModel.Text = 'cm';

            % Create GewichtLabelEigenModel
            app.GewichtLabelEigenModel = uilabel(app.UIFigure);
            app.GewichtLabelEigenModel.FontName = 'Footlight MT Light';
            app.GewichtLabelEigenModel.FontSize = 30;
            app.GewichtLabelEigenModel.FontColor = [0 0.451 0.7412];
            app.GewichtLabelEigenModel.Position = [1242 306 106 42];
            app.GewichtLabelEigenModel.Text = 'Gewicht';

            % Create SpierkrachtLabelEigenModel
            app.SpierkrachtLabelEigenModel = uilabel(app.UIFigure);
            app.SpierkrachtLabelEigenModel.FontName = 'Footlight MT Light';
            app.SpierkrachtLabelEigenModel.FontSize = 30;
            app.SpierkrachtLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.SpierkrachtLabelEigenModel.Position = [813 308 146 39];
            app.SpierkrachtLabelEigenModel.Text = 'Spierkracht';

            % Create PlanetXblack
            app.PlanetXblack = uiimage(app.UIFigure);
            app.PlanetXblack.Position = [953 79 244 244];
            app.PlanetXblack.ImageSource = fullfile(pathToMLAPP, 'PlanetX_dark.png');

            % Create Spaceman
            app.Spaceman = uiimage(app.UIFigure);
            app.Spaceman.Position = [358 -25 591 505];
            app.Spaceman.ImageSource = fullfile(pathToMLAPP, 'Spaceman_nobackground.png');

            % Create Jullie_strength
            app.Jullie_strength = uiimage(app.UIFigure);
            app.Jullie_strength.Position = [546 147 323 323];
            app.Jullie_strength.ImageSource = fullfile(pathToMLAPP, 'Jullie_nobackground.png');

            % Create MrIncredible
            app.MrIncredible = uiimage(app.UIFigure);
            app.MrIncredible.Position = [789 139 377 377];
            app.MrIncredible.ImageSource = fullfile(pathToMLAPP, 'incredibles2-mrincredible2.png');

            % Create Volwassenen
            app.Volwassenen = uiimage(app.UIFigure);
            app.Volwassenen.Position = [848 154 386 378];
            app.Volwassenen.ImageSource = fullfile(pathToMLAPP, 'Volwassenen.png');

            % Create Mars
            app.Mars = uiimage(app.UIFigure);
            app.Mars.Position = [708 383 317 317];
            app.Mars.ImageSource = fullfile(pathToMLAPP, 'Mars.png');

            % Create Jullie_lengte
            app.Jullie_lengte = uiimage(app.UIFigure);
            app.Jullie_lengte.Position = [626 162 324 328];
            app.Jullie_lengte.ImageSource = fullfile(pathToMLAPP, 'Jullie_nobackground.png');

            % Create ZwaartekrachtEditFieldEigenModel
            app.ZwaartekrachtEditFieldEigenModel = uieditfield(app.UIFigure, 'numeric');
            app.ZwaartekrachtEditFieldEigenModel.Limits = [1 9];
            app.ZwaartekrachtEditFieldEigenModel.ValueDisplayFormat = '%111g';
            app.ZwaartekrachtEditFieldEigenModel.HorizontalAlignment = 'center';
            app.ZwaartekrachtEditFieldEigenModel.FontName = 'Footlight MT Light';
            app.ZwaartekrachtEditFieldEigenModel.FontSize = 30;
            app.ZwaartekrachtEditFieldEigenModel.FontColor = [0 0.4471 0.7412];
            app.ZwaartekrachtEditFieldEigenModel.BackgroundColor = [1 0.9686 0.8588];
            app.ZwaartekrachtEditFieldEigenModel.Position = [1025 221 74 41];
            app.ZwaartekrachtEditFieldEigenModel.Value = 1;

            % Create ZwaartekrachtLabelEigenModel
            app.ZwaartekrachtLabelEigenModel = uilabel(app.UIFigure);
            app.ZwaartekrachtLabelEigenModel.FontName = 'Footlight MT Light';
            app.ZwaartekrachtLabelEigenModel.FontSize = 30;
            app.ZwaartekrachtLabelEigenModel.FontColor = [0 0.4471 0.7412];
            app.ZwaartekrachtLabelEigenModel.Position = [810 224 186 39];
            app.ZwaartekrachtLabelEigenModel.Text = 'Zwaartekracht';

            % Create giant
            app.giant = uiimage(app.UIFigure);
            app.giant.Position = [1096 46 548 552];
            app.giant.ImageSource = fullfile(pathToMLAPP, 'Giant_nobackground.png');

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

            % Create CheckBoxPluto
            app.CheckBoxPluto = uicheckbox(app.UIFigure);
            app.CheckBoxPluto.Text = '';
            app.CheckBoxPluto.Position = [1196 344 25 22];

            % Create maan
            app.maan = uiimage(app.UIFigure);
            app.maan.Position = [143 306 457 377];
            app.maan.ImageSource = fullfile(pathToMLAPP, 'moon_drawing.png');

            % Create Jupiter
            app.Jupiter = uiimage(app.UIFigure);
            app.Jupiter.Position = [19 -46 377 393];
            app.Jupiter.ImageSource = fullfile(pathToMLAPP, 'Jupiter_noBackground.png');

            % Create earth
            app.earth = uiimage(app.UIFigure);
            app.earth.Position = [1196 10 306 306];
            app.earth.ImageSource = fullfile(pathToMLAPP, 'Earth_nobackground.png');

            % Create CheckBoxMaan
            app.CheckBoxMaan = uicheckbox(app.UIFigure);
            app.CheckBoxMaan.Text = '';
            app.CheckBoxMaan.Position = [358 476 25 22];

            % Create CheckBoxMars
            app.CheckBoxMars = uicheckbox(app.UIFigure);
            app.CheckBoxMars.Text = '';
            app.CheckBoxMars.Position = [854 530 25 22];

            % Create CheckBoxEarth
            app.CheckBoxEarth = uicheckbox(app.UIFigure);
            app.CheckBoxEarth.Text = '';
            app.CheckBoxEarth.Position = [1356 140 25 22];

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

            % Create Olaf
            app.Olaf = uiimage(app.UIFigure);
            app.Olaf.Position = [53 172 307 307];
            app.Olaf.ImageSource = fullfile(pathToMLAPP, 'Olaf_nobackground.png');

            % Create Spongebob
            app.Spongebob = uiimage(app.UIFigure);
            app.Spongebob.Position = [330 176 247 247];
            app.Spongebob.ImageSource = fullfile(pathToMLAPP, 'Spongebob.png');

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

            % Create GroepEditFieldLabel
            app.GroepEditFieldLabel = uilabel(app.UIFigure);
            app.GroepEditFieldLabel.FontName = 'Footlight MT Light';
            app.GroepEditFieldLabel.FontSize = 30;
            app.GroepEditFieldLabel.FontColor = [0 0.451 0.7412];
            app.GroepEditFieldLabel.Position = [123 588 83 42];
            app.GroepEditFieldLabel.Text = 'Groep';

            % Create GroepEditField
            app.GroepEditField = uieditfield(app.UIFigure, 'text');
            app.GroepEditField.ValueChangedFcn = createCallbackFcn(app, @GroepEditFieldValueChanged, true);
            app.GroepEditField.FontName = 'Footlight MT Light';
            app.GroepEditField.FontSize = 30;
            app.GroepEditField.FontColor = [0 0.451 0.7412];
            app.GroepEditField.BackgroundColor = [1 0.9686 0.8588];
            app.GroepEditField.Placeholder = '(Geef de naam van jouw groepje)';
            app.GroepEditField.Position = [207 569 421 61];

            % Create Ballon
            app.Ballon = uiimage(app.UIFigure);
            app.Ballon.Enable = 'off';
            app.Ballon.Visible = 'off';
            app.Ballon.Position = [576 511 214 177];
            app.Ballon.ImageSource = 'balloon.png';

            % Create GIF
            app.GIF = uiimage(app.UIFigure);
            app.GIF.Enable = 'off';
            app.GIF.Visible = 'off';
            app.GIF.Position = [627 569 98 71];
            app.GIF.ImageSource = 'allgifs.gif';

            % Create CheckBoxJupiter
            app.CheckBoxJupiter = uicheckbox(app.UIFigure);
            app.CheckBoxJupiter.Text = '';
            app.CheckBoxJupiter.Position = [194 118 25 22];

            % Create CheckBox_Smurf
            app.CheckBox_Smurf = uicheckbox(app.UIFigure);
            app.CheckBox_Smurf.Text = '';
            app.CheckBox_Smurf.Position = [277 139 25 22];

            % Create CheckBox_Kleuter
            app.CheckBox_Kleuter = uicheckbox(app.UIFigure);
            app.CheckBox_Kleuter.Text = '';
            app.CheckBox_Kleuter.Position = [523 137 25 22];

            % Create CheckBox_Jullie_lengte
            app.CheckBox_Jullie_lengte = uicheckbox(app.UIFigure);
            app.CheckBox_Jullie_lengte.Text = '';
            app.CheckBox_Jullie_lengte.Position = [765 137 25 22];

            % Create CheckBox_Volwassenen
            app.CheckBox_Volwassenen = uicheckbox(app.UIFigure);
            app.CheckBox_Volwassenen.Text = '';
            app.CheckBox_Volwassenen.Position = [1033 139 25 22];

            % Create Kleuters
            app.Kleuters = uiimage(app.UIFigure);
            app.Kleuters.Position = [377 158 310 280];
            app.Kleuters.ImageSource = fullfile(pathToMLAPP, 'Kleuters.png');

            % Create smurf
            app.smurf = uiimage(app.UIFigure);
            app.smurf.Position = [172 160 237 204];
            app.smurf.ImageSource = fullfile(pathToMLAPP, 'Smurf_nobackground.png');

            % Create CheckBox_Reus
            app.CheckBox_Reus = uicheckbox(app.UIFigure);
            app.CheckBox_Reus.Text = '';
            app.CheckBox_Reus.Position = [1317 139 25 22];

            % Create gLabelEigenModel
            app.gLabelEigenModel = uilabel(app.UIFigure);
            app.gLabelEigenModel.FontName = 'Footlight MT Light';
            app.gLabelEigenModel.FontSize = 30;
            app.gLabelEigenModel.FontColor = [1 0.9686 0.8588];
            app.gLabelEigenModel.Position = [1109 224 25 39];
            app.gLabelEigenModel.Text = 'g';

            % Create LengteSliderEigenModel
            app.LengteSliderEigenModel = uislider(app.UIFigure);
            app.LengteSliderEigenModel.Limits = [80 180];
            app.LengteSliderEigenModel.MinorTicks = [80 85 90.95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180];
            app.LengteSliderEigenModel.FontName = 'Footlight MT Light';
            app.LengteSliderEigenModel.FontSize = 30;
            app.LengteSliderEigenModel.FontColor = [0 0.451 0.7412];
            app.LengteSliderEigenModel.Position = [936 425 467 3];
            app.LengteSliderEigenModel.Value = 120;

            % Create SimulatieNietGeluktWarning
            app.SimulatieNietGeluktWarning = uilabel(app.UIFigure);
            app.SimulatieNietGeluktWarning.BackgroundColor = [1 0.9686 0.8588];
            app.SimulatieNietGeluktWarning.FontName = 'Footlight MT Light';
            app.SimulatieNietGeluktWarning.FontSize = 30;
            app.SimulatieNietGeluktWarning.FontColor = [0 0.451 0.7412];
            app.SimulatieNietGeluktWarning.Visible = 'off';
            app.SimulatieNietGeluktWarning.Position = [676 564 694 71];
            app.SimulatieNietGeluktWarning.Text = '             Simulatie niet gelukt 😔';

            % Create ProbeeropnieuwButton
            app.ProbeeropnieuwButton = uibutton(app.UIFigure, 'push');
            app.ProbeeropnieuwButton.ButtonPushedFcn = createCallbackFcn(app, @ProbeeropnieuwButtonPushed, true);
            app.ProbeeropnieuwButton.BackgroundColor = [0 0.4471 0.7412];
            app.ProbeeropnieuwButton.FontName = 'Footlight MT Light';
            app.ProbeeropnieuwButton.FontSize = 20;
            app.ProbeeropnieuwButton.FontColor = [1 0.9686 0.8588];
            app.ProbeeropnieuwButton.Enable = 'off';
            app.ProbeeropnieuwButton.Visible = 'off';
            app.ProbeeropnieuwButton.Position = [1137 579 172 40];
            app.ProbeeropnieuwButton.Text = 'Probeer opnieuw';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = OpStapMetSuperKracht2025_FINAL_exported

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