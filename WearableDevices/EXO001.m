classdef EXO001 < Orthosis
% --------------------------------------------------------------------------
% EXO001
%   Interface to define custom controllers for an ankle exoskeleton
%   emulator. The emulator is modelled after the EXO001 by Humotech.
%
%   See also Orthosis
% 
% Original author: Lars D'Hondt
% Original date: January 2024
% --------------------------------------------------------------------------

    properties(Access = private)
        % indicate left or right side
        side = 'right'

        % dimentions of exo frame
        L_shank = 0.3; % length ankle-shank strap
        H_ff; % heigth ankle-mtp
        L_ff; % length ankle-mtp
        r_pos_ff_force = 0.3; % distance to ankle relative to ankle-mtp distance

        % reference points in OpenSim model
        osimKneeCentreBody = 'tibia_r';
        osimKneeCentrePosInBody = [0;0;0];
        osimAnkleCentreBody = 'talus_r';
        osimAnkleCentrePosInBody = [0;0;0];
        osimMTPCentreBody = 'toes_r';
        osimMTPCentrePosInBody = [0;0;0];

        % OpenSim bodies on which the exo transmits force
        osimShankForceBody = 'tibia_r';
        osimHeelForceBody = 'calcn_r';
        osimFootForceBody = 'midfoot_r';

    end

    methods

    %% constructor
        function self = EXO0001(name,shoeSizeEU,init,isTimeVarying)
            arguments
                name char
                shoeSizeEU (1,1) double
                init struct = [];
                isTimeVarying (1,1) logical = false;
            end
            % Create a new EXO001
            %   Constructor should be used inside a function.
            %   
            %   See also parametricAFO, ankleExoZhang2017 
            %
            % INPUT:
            %   - name - [char]
            %   * User-defined name. 
            % 
            %   - shoeSizeEU - [double]
            %   * Size of the shoe.
            %
            %   - init - [struct]
            %   * Information used to initialise the Orthosis. 
            %   This is be the 1st input to the function where the EXO001 
            %   is defined, and should be passed to the constructor.
            %
            %   - isTimeVarying - [boolean] (optional) Default: false
            %   * Boolean indicating whether the expression that determines
            %   the generalised forces is time-varying. E.g. a predefined
            %   torque profile.
            %
            % OUTPUT:
            %   - self - [EXO001]
            %   * new EXO001
            
            self.name = name;

            dself.setForefootDimentions(shoeSizeEU)

            if ~isempty(init)
                if isfield(init,'Nmesh')
                    self.Nmesh = init.Nmesh;
                end
                if isfield(init,'osimPath')
                    self.osimPath = init.osimPath;
                    updatePropertiesFromOsimModel(self);
                end
            end

            if ~isTimeVarying
                self.Nmesh = 1;
            end
        end

    %% getters and setters
        function setForefootDimentions(self,shoeSizeEU)
            % Set the dimentions of the forefoot frame based on shoe size

            switch true
                case 34 < shoeSizeEU && shoeSizeEU <= 37 % 36
%                     self.H_ff = 
%                     self.L_ff = 
                case 37 < shoeSizeEU && shoeSizeEU <= 41 % 38
%                     self.H_ff = 
%                     self.L_ff = 
                case 41 < shoeSizeEU && shoeSizeEU <= 45 % 43
                    self.H_ff = 0.14;
                    self.L_ff = 0.18;
                case 45 < shoeSizeEU && shoeSizeEU <= 47 % 46
%                     self.H_ff = 
%                     self.L_ff = 
            end
        end
    
        function setSide(self, side)
            arguments
                self
                side char {mustBeMember(side,{'left','right','l','r'})}
            end
            % Exo is on left or right side?

            sideProperties = ["osimKneeCentreBody", "osimAnkleCentreBody",...
                "osimMTPCentreBody", "osimShankForceBody", "osimHeelForceBody",...
                "osimFootForceBody"];

            if side(1)=='l'
                self.side = 'left';
                for prop=sideProperties
                    [self.(prop),~] = mirrorName(self.(prop));
                end

            elseif side(1)=='r'
                self.side = 'right';
                for prop=sideProperties
                    [~,self.(prop)] = mirrorName(self.(prop));
                end

            end
        end % end of setSide
    
        function setKneeJointCentre(self, osim_body_name, location_in_body)
            arguments
                self
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
            end
            % Set the knee joint centre. Used for alignment of the exo frame.

            self.osimKneeCentreBody = osim_body_name;
            self.osimKneeCentrePosInBody = location_in_body;    
        end

        function setAnkleJointCentre(self, osim_body_name, location_in_body)
            arguments
                self
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
            end
            % Set the ankle joint centre. Used for alignment of the exo frame.

            self.osimAnkleCentreBody = osim_body_name;
            self.osimAnkleCentrePosInBody = location_in_body;    
        end

        function setMTPJointCentre(self, osim_body_name, location_in_body)
            arguments
                self
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
            end
            % Set the MTP joint centre. Used for alignment of the exo frame.

            self.osimMTPCentreBody = osim_body_name;
            self.osimMTPCentrePosInBody = location_in_body;    
        end



    end % end of methods

    methods(Hidden=true)

        function [] = addCoordForce(self,varargin)
            % 

        end

    end
end