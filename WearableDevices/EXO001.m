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

    properties(Access = protected)
        % general properties
        side = 'r' % indicate left or right side
        PFMoment; % plantarflexion moment

        % dimentions of exo frame
        L_shank = 0.3; % length ankle-shank strap
        H_ff; % heigth ankle-mtp
        L_ff; % length ankle-mtp
        r_pos_ff_force = 1; % distance to ankle relative to ankle-mtp distance
        D_toe_force;

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
        osimFootForceBody = 'forefoot_r';
        osimFootForcePosInBody
        osimToeForcePosInBody

        % Casadi Functions
        f_encoderPos = [];
            
        % Internal states
        encoderPos = [];
        encoderVel = [];
        encoder_offset = [];

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
            
            % call constructor of superclass (Orthosis)
            self@Orthosis(name,init,isTimeVarying);

            % set size property values
            self.setForefootDimentions(shoeSizeEU)

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
                self.side = 'l';
                for prop=sideProperties
                    [self.(prop),~] = mirrorName(self.(prop));
                end

            elseif side(1)=='r'
                self.side = 'r';
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

        function setShankForceBody(self, osim_body_name)
            arguments
                self EXO001
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
            end
            % Set the body on which the shank force acts.

            self.osimShankForceBody = osim_body_name;  
        end

        function setHeelForceBody(self, osim_body_name)
            arguments
                self EXO001
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
            end
            % Set the body on which the heel force acts.

            self.osimHeelForceBody = osim_body_name;  
        end

        function setFootForceBody(self, osim_body_name)
            arguments
                self EXO001
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
            end
            % Set the body on which the foot force acts.

            self.osimFootForceBody = osim_body_name;  
        end

        function setPlantarFlexionMoment(self,M)
            arguments
                self EXO001
                M (1,:)
            end
            % Set ankle plantarflexion moment provided by exoskeleton

            self.PFMoment = M;
        end
        
    %% Create a variable
    function enc = var_encoder(self,pos_vel)
            arguments
                self EXO001
                pos_vel char {mustBeMember(pos_vel,{'pos','vel'})} = 'pos';
            end
            % Create a variable for the ankle encoder position or velocity

            switch pos_vel
                case 'pos'
                    if isempty(self.encoderPos)
                        calcEncoder(self);
                    end
                    enc = self.encoderPos;
%                 case 'vel'
%                     if isempty(self.encoderVel)
%                         calcEncoder(self);
%                     end
%                     enc = self.encoderVel;
            end
    end % end of var_encoder

    %% CasADi
        function [] = createCasadiFunction(self)
            % create CasADi Function

            import casadi.*

            % get encoder position
            if isempty(self.encoderPos)
                calcEncoder(self);
            end

            % create function for geometry
            f_geometry = self.geometryFunctionFactory();
            
            q_encoder = self.encoderPos;

            % evaluate geometry function to get force vectors
            [F_shank,F_heel,F_foot,F_toe] = f_geometry(q_encoder, self.PFMoment);

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

            addBodyForce(self,F_toe,['F_toe_',self.side],...
                self.osimFootForceBody, self.osimToeForcePosInBody,...
                self.osimShankForceBody);

            self.addVarToPostProcessing(F_shank,'F_shank');
            self.addVarToPostProcessing(F_heel,'F_heel');
            self.addVarToPostProcessing(F_foot,'F_foot');
            self.addVarToPostProcessing(F_toe,'F_toe');

            % create the casadi function
            createCasadiFunction@Orthosis(self)
        end

    end % end of methods

    methods(Hidden=true)
        % Override methods of Orthosis that are not supported here and set
        % them to 'Hidden' so they do not show up in documentation.
        function [] = addCoordForce(self,varargin)
            addCoordForce@Orthosis(self,varargin{:})
        end
        function [] = addBodyForce(self,varargin)
            addBodyForce@Orthosis(self,varargin{:})
        end
        function [] = addBodyMoment(self,varargin)
            addBodyMoment@Orthosis(self,varargin{:})
        end
    end

    methods(Access=protected) % protected — Access from methods in class or subclasses

        function [] = encoderPosFunctionFactory(self)
            % Create a casadi Function to calculate the encoder position.

            import casadi.*

            if isempty(self.encoder_offset)
                self.zeroEncoder();
            end

            % positions of joint centres wrt ground frame
            pos_knee = SX.sym('pos_knee',3,1);
            pos_ankle = SX.sym('pos_ankle',3,1);
            pos_mtp = SX.sym('pos_mtp',3,1);

            % orientations of exo frames
            vec_shank = pos_knee - pos_ankle;
            vec_foot = pos_mtp - pos_ankle;
            
            % angle between exo frames
            q_exo_hw = acos( dot(vec_shank, vec_foot)./...
                ( norm(vec_foot)*norm(vec_shank) ) );

            % account for encoder offset
            q_exo = q_exo_hw - self.encoder_offset;

            self.f_encoderPos = Function('f_encoderPos',...
                {pos_knee,pos_ankle,pos_mtp},{q_exo});
        end % end of encoderPosFunctionFactory

        function [] = zeroEncoder(self)

            import org.opensim.modeling.*

            model = Model(self.osimPath);
            state = model.initSystem();

            %
            fr_knee = model.get_BodySet().get(self.osimKneeCentreBody);
            pos_knee= fr_knee.findStationLocationInGround(state,...
                Vec3.createFromMat(self.osimKneeCentrePosInBody)).getAsMat;
            
            fr_ankle = model.get_BodySet().get(self.osimAnkleCentreBody);
            pos_ankle= fr_ankle.findStationLocationInGround(state,...
                Vec3.createFromMat(self.osimAnkleCentrePosInBody)).getAsMat;

            fr_mtp = model.get_BodySet().get(self.osimMTPCentreBody);
            pos_mtp= fr_mtp.findStationLocationInGround(state,...
                Vec3.createFromMat(self.osimMTPCentrePosInBody)).getAsMat;

            % orientations of exo frames
            vec_shank = pos_knee - pos_ankle;
            vec_foot = pos_mtp - pos_ankle;

            % encoder angle at neutral position
            enc_offset = acos( dot(vec_shank, vec_foot)./...
                ( norm(vec_foot)*norm(vec_shank) ) );

             self.encoder_offset = enc_offset;

             self.D_toe_force = norm(vec_foot);
             self.H_ff = abs(vec_foot(2));
             self.L_ff = abs(vec_foot(1));

        end % end of zeroEncoder

        function [] = calcEncoder(self)
            % Calculate the encoder position and velocity

            import casadi.*

            % define variables for point positions needed by geometry function
            knee_pos = self.var_point(['knee_pos_',self.side],...
                self.osimKneeCentreBody,self.osimKneeCentrePosInBody,'pos');

            ankle_pos = self.var_point(['ankle_pos_',self.side],...
                self.osimAnkleCentreBody,self.osimAnkleCentrePosInBody,'pos');

            mtp_pos = self.var_point(['mtp_pos_',self.side],...
                self.osimMTPCentreBody,self.osimMTPCentrePosInBody,'pos');

            % add point positions to post-processing
            self.addVarToPostProcessing(knee_pos,'knee_pos');
            self.addVarToPostProcessing(ankle_pos,'ankle_pos');
            self.addVarToPostProcessing(mtp_pos,'mtp_pos');

%             knee_vel = self.var_point(['knee_pos_',self.side],...
%                 self.osimKneeCentreBody,self.osimKneeCentrePosInBody,'pvel');
% 
%             ankle_vel = self.var_point(['ankle_pos_',self.side],...
%                 self.osimAnkleCentreBody,self.osimAnkleCentrePosInBody,'vel');
% 
%             mtp_vel = self.var_point(['mtp_pos_',self.side],...
%                 self.osimMTPCentreBody,self.osimMTPCentrePosInBody,'vel');

            if isempty(self.f_encoderPos)
                encoderPosFunctionFactory(self);
            end

            self.encoderPos = self.f_encoderPos(knee_pos,ankle_pos,mtp_pos);

            self.addVarToPostProcessing(self.encoderPos,'encoderPos');

        end % end of calcEncoder

        function [f_geometry] = geometryFunctionFactory(self)
            % Create a casadi Function with the geometric relationship
            % between exo moment and forces acting on the wearer.

            import casadi.*
            import org.opensim.modeling.*

            model = Model(self.osimPath);
            state = model.initSystem();

            % encoder position
            q_encoder = SX.sym('q_encoder',1,1);
            % moment setpoint from controller
            M = SX.sym('M',1,1);
            
            % angle of heel rope wrt shank
            theta_0 = 30 *pi/180;
            theta = theta_0 + q_encoder;
            
            % angle of forefoot force wrt shank
            phi_0 = atan(self.H_ff/self.L_ff);
            phi = phi_0 + q_encoder;
            
            % momentarm forefoot force wrt ankle centre
            D_ff_force = sqrt(self.H_ff^2 + self.L_ff^2)*self.r_pos_ff_force;
            
            % Moment equilibrium of shank frame (quasi-static)
            % M - L_shank * F_shank_force = 0
            F_S = M/self.L_shank;
            
%             % Moment equilibrium of foot frame (quasi-static)
%             % M - D_ff_force * F_ff_force = 0
%             F_ff_force = M/D_ff_force;
%             
%             % Vertical force equilibrium of exo frame (quasi-static)
%             % F_ff_force*cos(phi) - F_heel_force*cos(theta) = 0
%             %   assume that forces on toe plate, from shoe and from ground, cancel out
%             F_heel_force = F_ff_force*cos(phi)/cos(theta);
            

%             F_F = ( F_S*(cos(q_encoder) - sin(q_encoder)*tan(theta_0)) +...
%                 M/self.D_toe_force ) /...
%                 (cos(phi_0)*tan(theta_0) + D_ff_force/self.D_toe_force - sin(phi_0));
% 
%             F_T = (F_F*D_ff_force - M)/self.D_toe_force;
% 
%             F_H = (F_S*sin(q_encoder) + F_F*cos(phi_0))/cos(theta_0);


            F_F = M/D_ff_force;

            F_H = (F_F + F_S*sin(phi_0+q_encoder))/cos(theta_0-phi_0);

            F_T = F_H*sin(theta_0-phi_0) - F_S*cos(phi_0+q_encoder);

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
            
            self.osimToeForcePosInBody = fr_ground.findStationLocationInAnotherFrame(state,...
                Vec3.createFromMat(pos_ankle_0 + [self.L_ff;-self.D_toe_force;0]),...
                fr_foot).getAsMat;

            % positions of exo force acting on heel
            fr_heel = model.get_BodySet().get(self.osimHeelForceBody);
            self.osimHeelForcePosInBody = fr_ground.findStationLocationInAnotherFrame(state,...
                Vec3.createFromMat(pos_ankle_0 + [0;0;0]), fr_heel).getAsMat;
            
            % exo forces acting on body
            %   all vectors are expressed in shank frame
            F_shank = [-F_S; 0; 0];
            
            F_heel = [sin(theta); cos(theta); 0]*F_H;
            
            F_foot = [-sin(phi); -cos(phi); 0]*F_F;

            F_toe = [-cos(phi); sin(phi); 0]*F_T;
            
            
            f_geometry = Function('f_geometry',{q_encoder,M},...
                {F_shank, F_heel, F_foot, F_toe});

        end % end of geometryFunctionFactory

    end % end of protected methods

end