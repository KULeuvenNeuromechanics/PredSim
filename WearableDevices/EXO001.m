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
        % general properties
        side = 'right' % indicate left or right side
        PFMoment; % plantarflexion moment

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
        osimShankForcePosInBody
        osimHeelForceBody = 'calcn_r';
        osimHeelForcePosInBody
        osimFootForceBody = 'midfoot_r';
        osimFootForcePosInBody

    end

    methods

    %% constructor
        function self = EXO001(name,shoeSizeEU,init,isTimeVarying)
            arguments
                name char
                shoeSizeEU (1,1) double
                init struct = [];
                isTimeVarying (1,1) logical = false;
            end
            % Create a new EXO001
            %   Constructor should be used inside a function.
            %   
            %   See also 
            %
            % INPUT:
            %   - name - [char]
            %   * User-defined name. 
            % 
            %   - shoeSizeEU - [double]
            %   * Size of the shoe that is mounted.
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

            self.setForefootDimentions(shoeSizeEU)

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
            arguments
                self EXO001
                shoeSizeEU (1,1) double
            end
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
                self EXO001
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
                self EXO001
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
            end
            % Set the knee joint centre. Used for alignment of the exo frame.

            self.osimKneeCentreBody = osim_body_name;
            self.osimKneeCentrePosInBody = location_in_body;    
        end

        function setAnkleJointCentre(self, osim_body_name, location_in_body)
            arguments
                self EXO001
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
            end
            % Set the ankle joint centre. Used for alignment of the exo frame.

            self.osimAnkleCentreBody = osim_body_name;
            self.osimAnkleCentrePosInBody = location_in_body;    
        end

        function setMTPJointCentre(self, osim_body_name, location_in_body)
            arguments
                self EXO001
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
            end
            % Set the MTP joint centre. Used for alignment of the exo frame.

            self.osimMTPCentreBody = osim_body_name;
            self.osimMTPCentrePosInBody = location_in_body;    
        end

        function setPlantarFlexionMoment(self,M)
            arguments
                self EXO001
                M (1,:)
            end
            % Set ankle plantarflexion moment provided by exoskeleton

            self.PFMoment = M;
        end
        
    %% CasADi
        function [] = createCasadiFunction(self)
            % create CasADi Function

            % create function for geometry
            f_geometry = self.geometryFunctionFactory();

            % define variables for point positions needed by geometry function
            knee_pos = self.var_point(['knee_pos_',self.side],...
                self.osimKneeCentreBody,self.osimKneeCentrePosInBody,'pos');

            ankle_pos = self.var_point(['ankle_pos_',self.side],...
                self.osimAnkleCentreBody,self.osimAnkleCentrePosInBody,'pos');

            mtp_pos = self.var_point(['mtp_pos_',self.side],...
                self.osimMTPCentreBody,self.osimMTPCentrePosInBody,'pos');

            % evaluate geometry function to get force vectors
            [F_shank,F_heel,F_foot] = f_geometry(knee_pos, ankle_pos, mtp_pos,...
                self.PFMoment);

            % create bodyForces from force vectors
            addBodyForce(self,F_shank,['F_shank_',self.side],...
                self.osimShankForceBody, self.osimShankForcePosInBody,...
                self.osimShankForceBody);

            addBodyForce(self,F_heel,['F_heel_',self.side],...
                self.osimHeelForceBody, self.osimHeelForcePosInBody,...
                self.osimShankForceBody);

            addBodyForce(self,F_foot,['F_foot_',self.side],...
                self.osimFootForceBody, self.osimFootForcePosInBody,...
                self.osimShankForceBody);

            % create the casadi function
            createCasadiFunction@Orthosis(self)
        end

    end % end of methods

    methods(Hidden=true)
        % Override methods of Orthosis that are not supported here and set
        % them to 'Hidden' so they do not show up in documentation.
        function [] = addCoordForce(self,varargin)
            addCoordForce@Orthosis(self,varargin)
        end
        function [] = addBodyForce(self,varargin)
            addBodyForce@Orthosis(self,varargin)
        end
        function [] = addBodyMoment(self,varargin)
            addBodyMoment@Orthosis(self,varargin)
        end
    end

    methods(Access=protected) % protected — Access from methods in class or subclasses

        function [f_geometry] = geometryFunctionFactory(self)
            % Create a casadi Function with the geometric relationship
            % between exo moment and forces acting on the wearer.

            % positions of joint centres wrt ground frame
            pos_knee = SX.sym('pos_knee',3,1);
            pos_ankle = SX.sym('pos_ankle',3,1);
            pos_mtp = SX.sym('pos_mtp',3,1);
            
            % moment setpoint from controller
            M = SX.sym('M',1,1);
            
            % orientations of exo frames
            vec_shank = pos_knee - pos_ankle;
            vec_foot = pos_mtp - pos_ankle;
            
            % angle between exo frames
            q_exo_hw_0 = (90 + atan(self.H_ff/self.L_ff)) *pi/180;
            q_exo_hw = acos( dot(vec_shank/norm(vec_shank), vec_foot/norm(vec_foot)));
            q_exo = q_exo_hw - q_exo_hw_0;
            
            % angle of heel rope wrt shank
            theta_0 = 30 *pi/180;
            theta = theta_0 + q_exo;
            
            % angle of forefoot force wrt shank
            phi_0 = atan(self.H_ff/self.L_ff);
            phi = phi_0 + q_exo;
            
            % momentarm forefoot force wrt ankle centre
            D_ff_force = sqrt(self.H_ff^2 + self.L_ff^2)*self.r_pos_ff_force;
            
            % Moment equilibrium of shank frame (quasi-static)
            % M - L_shank * F_shank_force = 0
            F_shank_force = M/self.L_shank;
            
            % Moment equilibrium of foot frame (quasi-static)
            % M - D_ff_force * F_ff_force = 0
            F_ff_force = M/D_ff_force;
            
            % Vertical force equilibrium of exo frame (quasi-static)
            % F_ff_force*cos(phi) - F_heel_force*cos(theta) = 0
            %   assume that forces on toe plate, from shoe and from ground, cancel out
            F_heel_force = F_ff_force*cos(phi)/cos(theta);
            
            
            % assuming origin of talus is ankle joint, and exo hinge is aligned with
            % ankle joint centre
            fr_ankle = model.get_BodySet().get(self.osimAnkleCentreBody);
            pos_ankle_0 = fr_ankle.findStationLocationInGround(state,...
                Vec3.createFromMat(self.osimAnkleCentrePosInBody)).getAsMat;
            
            fr_ground = model.getGround();
            
            % positions of exo force acting on tibia
            fr_shank = model.get_BodySet().get(self.osimShankForceBody);
            self.osimShankForcePosInBody = fr_ground.findStationLocationInAnotherFrame(state,...
                Vec3.createFromMat(pos_ankle_0 + [0;self.L_shank;0]), fr_shank).getAsMat;
            
            % positions of exo force acting on midfoot
            fr_foot = model.get_BodySet().get(self.osimFootForceBody);
            self.osimFootForcePosInBody = fr_ground.findStationLocationInAnotherFrame(state,...
                Vec3.createFromMat(pos_ankle_0 + [self.L_ff;-self.H_ff;0]*self.r_pos_ff_force),...
                fr_foot).getAsMat;
            
            % positions of exo force acting on heel
            fr_heel = model.get_BodySet().get(self.osimHeelForceBody);
            self.osimHeelForcePosInBody = fr_ground.findStationLocationInAnotherFrame(state,...
                Vec3.createFromMat(pos_ankle_0 + [0;0;0]), fr_heel).getAsMat;
            
            
            % exo forces acting on body
            %   all vectors are expressed in shank frame
            F_shank = [-F_shank_force; 0; 0];
            
            F_heel = [sin(theta); cos(theta); 0]*F_heel_force;
            
            F_foot = [-sin(phi); -cos(phi); 0]*F_ff_force;
            
            
            f_geometry = Function('f_geometry',{pos_knee,pos_ankle,pos_mtp,M},...
                {F_shank, F_heel, F_foot});

        end % end of geometryFunctionFactory

    end % end of protected methods

end