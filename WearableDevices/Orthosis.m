classdef Orthosis < handle
% --------------------------------------------------------------------------
% Orthosis
%   Interface to define custom orthosis or exoskeleton devices to be used
%   in predictive simulations. For examples on how to use this, see
%   parametricAFO, ankleExoZhang2017.
%
%   General workflow:
%   1) Create a MATLAB function describing an orthosis
%       a) function header
%           function [Orthosis_object] = exampleOrthosis(init, settings_orthosis)
%       
%       b) create Orthosis object
%           Orthosis_object = Orthosis('exampleName', init);
%
%       c) read parameter values
%           param_1 = settings_orthosis.param_1;
%           param_2 = settings_orthosis.param_2;
%
%       d) use methods to create variables
%           see Orthosis/var, Orthosis/var_coord, Orthosis/var_point,
%           Orthosis/var_GRF, Orthosis/var_muscle
%
%       e) calculate forces exerted by the orthosis
%           be creative ;)
%
%       f) use methods to add the forces to the musculoskeletal model
%           see Orthosis/addCoordForce, Orthosis/addBodyForce,
%           Orthosis/addBodyMoment
%
%   2) Add a custom orthosis to the predictive simulations via main.m
%       a) use the function we just created
%           S.orthosis.settings{1}.function_name = 'exampleOrthosis';
%
%       b) set the parameter values
%           S.orthosis.settings{1}.param_1 = 42;
%           S.orthosis.settings{1}.param_2 = "whatever";
%
%
%   See also parametricAFO, ankleExoZhang2017 
% 
% Original author: Lars D'Hondt
% Original date: January 2024
% --------------------------------------------------------------------------

    properties(Access = protected)

        name = []; % name of the orthosis
        Nmesh = 1; % number of meshpoints used to describe the time-varying behaviour of the orthosis. Set to 1 if time-independent.


        arg = {}; % input arguments of CasADi Function describing orthosis mechanics
        res = {}; % output arguments of CasADi Function describing orthosis mechanics
        res_pp = {}; % output arguments for post-processing function
        names_arg = {}; % names of input arguments of CasADi Function describing orthosis mechanics
        names_res = {}; % names of output arguments of CasADi Function describing orthosis mechanics
        names_res_pp = {}; % names of post-processing output arguments
        meta_arg = []; % metadata of input arguments of CasADi Function describing orthosis mechanics
        meta_res = []; % metadata of output arguments of CasADi Function describing orthosis mechanics
        fun = []; % handle of CasADi Function


        BodyForces = {}; % input3DBodyForces for OpenSimAD
        BodyMoments = {}; % input3DBodyMoments for OpenSimAD
        PointPositions = {}; % export3DPositions for OpenSimAD
        PointVelocities = {}; % export3DVelocities for OpenSimAD
        

        osimPath = []; % path to model file
        warningOsimPathNotSet = false; % only warn once that osimPath was not set.

        osimCoordsAll = {}; % coordinate names
        osimBodiesAll = {}; % body names
        osimContactsAll = {}; % contact sphere names
        osimMusclesAll = {}; % muscle names

        osimCoordsUsed = {}; % coordinate names used by orthosis
        osimBodiesUsed = {}; % body names used by orthosis
        osimContactsUsed = {}; % contact sphere names used by orthosis
        osimMusclesUsed = {}; % muscle names used by orthosis
    end

    methods

    %% constructor
        function self = Orthosis(name,init,isTimeVarying)
            arguments
                name char
                init struct = [];
                isTimeVarying (1,1) logical = false;
            end
            % Create a new Orthosis.
            %   Constructor should be used inside a function.
            %   
            %   See also parametricAFO, ankleExoZhang2017 
            %
            % INPUT:
            %   - name - [char]
            %   * User-defined name. 
            % 
            %   - init - [struct]
            %   * Information used to initialise the Orthosis. 
            %   This is be the 1st input to the function where the orthosis 
            %   is defined, and should be passed to the constructor.
            %
            %   - isTimeVarying - [boolean] (optional) Default: false
            %   * Boolean indicating whether the expression that determines
            %   the generalised forces is time-varying. E.g. a predefined
            %   torque profile.
            %
            % OUTPUT:
            %   - orthosis_obj - [Orthosis]
            %   * new Orthosis (object)
            %
            
            self.name = name;

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
        function name = getName(self)
            % Get the name of the Orthosis
            name = self.name;
        end

        function N_mesh = getNmesh(self)
            % Get the number of mesh intervals of the simulation
            N_mesh = self.Nmesh;
        end

        function fun = getFunction(self)
            % Get the CasADi Function describing the Orthosis
            if isempty(self.fun)
                createCasadiFunction(self);
            end
            fun = self.fun;
        end

        function BodyForces = getBodyForces(self)
            % Get cell array of structs describing body forces in OpenSimAD format
            BodyForces = self.BodyForces;
        end

        function BodyMoments = getBodyMoments(self)
            % Get cell array of structs describing body moments in OpenSimAD format
            BodyMoments = self.BodyMoments;
        end

        function PointPositions = getPointPositions(self)
            % Get cell array of structs describing point positions in OpenSimAD format
            PointPositions = self.PointPositions;
        end

        function PointVelocities = getPointVelocities(self)
            % Get cell array of structs describing point velocities in OpenSimAD format
            PointVelocities = self.PointPositions;
        end

        function [] = setOsimPath(self,osimPath)
            arguments
                self Orthosis
                osimPath char {mustBeFile}
            end
            % Set the path to the OpenSim model file associated with the Orthosis
            self.osimPath = osimPath;
            updatePropertiesFromOsimModel(self);
        end

    %% create a variable 
        function coord = var_coord(self,osim_coord_name,pos_vel_acc)
            arguments
                self Orthosis
                osim_coord_name char {inputExistsInOsimModel(self,osim_coord_name,'coord')}
                pos_vel_acc char {mustBeMember(pos_vel_acc,{'pos','vel','acc'})} = 'pos';
            end
            % Create a variable for a coordinate position, velocity or acceleration
            %
            % EXAMPLE
            %   Position of the right ankle:
            %   q_ankle_r = orthosis_obj.var_coord('ankle_angle_r');
            %   
            %   Angular acceleration of the left knee:
            %   ddq_knee_l = orthosis_obj.var_coord('knee_angle_l','acc');
            %
            % INPUT:
            %   - osim_coord_name - [char]
            %   * Name of a coordinate in the OpenSim model
            %
            %   - pos_vel_acc - [char] (optional) Default: 'pos'
            %   * Use position, velocity or acceleration of the coordinate
            %
            % OUTPUT:
            %   - coord - [1x1 variable]
            %   * Variable for the position, velocity or acceleration of a
            %   coordinate. This can be used to describe the Orthosis.
            %

            var_name = ['coord_',osim_coord_name,'_',pos_vel_acc];
            coord = casadi.SX.sym(var_name,1,self.Nmesh);
            self.arg{end+1} = coord;
            self.names_arg{end+1} = var_name;
            self.osimCoordsUsed{end+1} = osim_coord_name;
            self.meta_arg(end+1).name = osim_coord_name;
            self.meta_arg(end).type = 'coordi';
            switch pos_vel_acc
                case 'pos'
                    self.meta_arg(end).subtype = 'q';
                case 'vel'
                    self.meta_arg(end).subtype = 'qdot';
                case 'acc'
                    self.meta_arg(end).subtype = 'qddot';
            end
        end % end of var_coord

        function point = var_point(self,point_name,osim_body_name,location_in_body,pos_vel)
            arguments
                self Orthosis
                point_name char
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
                pos_vel char {mustBeMember(pos_vel,{'pos','vel'})} = 'pos';
            end
            % Create a variable for the position or velocity (w.r.t. ground) of a point
            %
            % EXAMPLE
            %   Position of the right hip joint centre:
            %   pos_hip_r = orthosis_obj.var_point('hip_r','femur_r');
            %   
            %   Velocity of the head:
            %   vel_head = orthosis_obj.var_point('head','torso',[0.06, 0.54, 0],'vel');
            %
            % INPUT:
            %   - point_name - [char]
            %   * Name of the point. Point names, defined here and in
            %   S.bounds.points, shoule be unique.
            %
            %   - osim_body_name - [char]
            %   * Name of a body in the OpenSim model, on which the point
            %   is located
            %
            %   - location in body - [1x3 double] (optional) Default: [0, 0, 0]
            %   * Location of the point in the body.
            %
            %   - pos_vel - [char] (optional) Default: 'pos'
            %   * Use position or velocity of the point
            %
            % OUTPUT:
            %   - point - [3x1 variable]
            %   * Variable for the position or velocity of a point. This
            %   can be used to describe the Orthosis. 
            %

            var_name = ['point_',point_name,'_',pos_vel];
            point = casadi.SX.sym(var_name,3,self.Nmesh);
            self.arg{end+1} = point;
            self.names_arg{end+1} = var_name;
            self.osimBodiesUsed{end+1} = osim_body_name;
            self.meta_arg(end+1).name = point_name;
            self.meta_arg(end).subtype = 'fromExtFun';
            switch pos_vel
                case 'pos'
                    self.PointPositions(end+1).body = osim_body_name;
                    self.PointPositions(end).point_in_body = location_in_body;
                    self.PointPositions(end).name = point_name;
                    self.meta_arg(end).type = 'position';
                case 'vel'
                    self.PointVelocities(end+1).body = osim_body_name;
                    self.PointVelocities(end).point_in_body = location_in_body;
                    self.PointVelocities(end).name = point_name;
                    self.meta_arg(end).type = 'velocity';
            end
        end % end of var_point

        function GRF = var_GRF(self,osim_contact_name,F_d)
            arguments
                self Orthosis
                osim_contact_name char {inputExistsInOsimModel(self,osim_contact_name,'contact')}
                F_d char {mustBeMember(F_d,{'Force','indentation'})} = 'Force';
            end
            % Create a variable for a ground reaction force or contact indentation
            %
            % EXAMPLE
            %   GRF under the right heel:
            %   GRF_heel_r = orthosis_obj.var_GRF('heel_r');
            %
            % INPUT:
            %   - osim_contact_name - [char]
            %   * Name of a contact sphere in the OpenSim model.
            %
            %   - F_d - ['Force' or 'indentation'] (optional) Default: 'Force'
            %   * Use force or indentation. 
            %
            % OUTPUT:
            %   - GRF - [3x1 or 1x1 variable]
            %   * F_d = 'Force': vector of xyz ground reaction forces.
            %     F_d = 'indentation': indentation of the contact sphere
            %

            var_name = ['GRF_',osim_contact_name,'_',F_d];
            switch F_d
                case 'Force'
                    GRF = casadi.SX.sym(var_name,3,self.Nmesh);
                    self.arg{end+1} = GRF;
                    self.names_arg{end+1} = var_name;
                    self.osimContactsUsed{end+1} = osim_contact_name;
                    self.meta_arg(end+1).name = contact_name;
                    self.meta_arg(end).subtype = 'fromExtFun';
                    self.meta_arg(end).type = 'GRFs';
                case 'indentation'
                    GRF = calcContactIndentation(self,osim_contact_name);
                    % Metadata for contact indentation is added when
                    % creating the point position variable.
            end
        end % end of var_GRF

        function act = var_muscle(self,osim_muscle_name)
            arguments
                self Orthosis
                osim_muscle_name char {inputExistsInOsimModel(self,osim_muscle_name,'muscle')}
            end
            % Create a variable for the activation of a muscle
            %
            % EXAMPLE
            %   Activation of the right soleus:
            %   act_sol_r = orthosis_obj.var_muscle('soleus_r');
            %
            % INPUT:
            %   - muscle_name - [char]
            %   * Name of a muscle in the OpenSim model.
            %
            % OUTPUT:
            %   - act - [1x1 variable]
            %   * Variable for the activation of a muscle.
            %

            % note: can be extended to include e.g. fibre length

            var_name = ['muscle_',osim_muscle_name,'_act'];
            act = casadi.SX.sym(var_name,1,self.Nmesh);
            self.arg{end+1} = act;
            self.names_arg{end+1} = var_name;
            self.osimMusclesUsed{end+1} = osim_muscle_name;
            self.meta_arg(end+1).name = osim_muscle_name;
            self.meta_arg(end).type = 'muscle';
            self.meta_arg(end).subtype = 'act';
        end % end of var_muscle

        function var = var(self,varargin)
            % Create a variable. 
            % See also var_coord, var_point, var_GRF, var_muscle

            try
                if ismember(varargin{1},self.osimCoordsAll)
                    var = var_coord(varargin{:});
                elseif ismember(varargin{2},self.osimBodiesAll)
                    var = var_point(varargin{:});
                elseif ismember(varargin{1},self.osimContactsAll)
                    var = var_GRF(varargin{:});
                elseif ismember(varargin{1},self.osimMusclesAll)
                    var = var_muscle(varargin{:});
                else
                    error('')
                end

            catch
                error(['Unable to use %s.var with the provided arguments.',...
                    ' Try using %s.var_coord, %s.var_point, %s.var_GRF, or %s.var_mus'],...
                    inputname(1),inputname(1),inputname(1),inputname(1),inputname(1));
            end
        end % end of var


    %% add a force/moment
        function [] = addCoordForce(self,value,osim_coord_name)
            arguments
                self Orthosis
                value (1,:)
                osim_coord_name char {inputExistsInOsimModel(self,osim_coord_name,'coord')}
            end
            % Add a force or moment acting on a coordinate.
            %
            % EXAMPLE
            %   Add moment to the right ankle:
            %   orthosis_obj.addCoordForce(5,'ankle_angle_r');
            %
            % INPUT:
            %   - value - [1x1 variable or double]
            %   * Value of the force or moment (in N or Nm)
            %
            %   - osim_coord_name - [char]
            %   * Name of a coordinate in the OpenSim model
            %

            % value should be a row vector with Nmesh elements
            if size(value,1)~=1 || size(value,2)~=self.Nmesh
                error('Expected "%s" to have size %ix%i, but the size is %ix%i.',...
                    inputname(2),1,self.Nmesh,size(value,1),size(value,2));
            end

            F_name = ['CoordForce_',osim_coord_name];

            idx = find(cellfun(@(x)strcmp(x,F_name), self.names_res));
            if isempty(idx)
                self.res{end+1} = value;
                self.names_res{end+1} = F_name;
                self.osimCoordsUsed{end+1} = osim_coord_name;
                self.meta_res(end+1).name = osim_coord_name;
                self.meta_res(end).type = 'coordi';
                self.meta_res(end).subtype = 'Mcoord';
            else
                self.res{idx} = self.res{idx} + value;
            end
        end % end of addCoordForce

        function [] = addBodyForce(self,value,force_name,osim_body_name,...
                location_in_body,reference_frame)
            arguments
                self Orthosis
                value (3,:)
                force_name char
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
                reference_frame char {inputExistsInOsimModel(self,reference_frame,'frame')} = osim_body_name;
            end
            % Add a force vector acting on a body.
            %
            % EXAMPLE
            %   Add a force on the pelvis, with the vector given in pelvis reference frame:
            %   orthosis_obj.addBodyForce(F_push, 'F_push', 'pelvis');
            %
            %   Add a force on the pelvis, with the vector given in ground reference frame:
            %   orthosis_obj.addBodyForce(F_push, 'F_push', 'pelvis', [0,0,0], 'ground');
            %
            % INPUT:
            %   - value - [3x1 variable or double]
            %   * Value of the force (in N)
            %
            %   - force_name - [char]
            %   * Name of the force.
            %
            %   - osim_body_name - [char]
            %   * Name of a body in the OpenSim model, on which the force
            %   is applied.
            %
            %   - location in body - [1x3 double] (optional) Default: [0, 0, 0]
            %   * Location of the point in the body where the force is applied.
            %
            %   - reference frame - [char] (optional) Default: osim_body_name
            %   * Name of a body in the OpenSim model, in whose reference
            %   frame the force vector is expressed.
            %

            % value should be a 3xNmesh matrix
            if size(value,1)~=3 || size(value,2)~=self.Nmesh
                error('Expected "%s" to have size %ix%i, but the size is %ix%i.',...
                    inputname(2),3,self.Nmesh,size(value,1),size(value,2));
            end

            F_name = ['BodyForce_',force_name];

            idx = find(cellfun(@(x)strcmp(x,F_name), self.names_res));
            if isempty(idx)
                self.res{end+1} = value;
                self.names_res{end+1} = F_name;

                self.BodyForces(end+1).body = osim_body_name;
                self.BodyForces(end).point_in_body = location_in_body;
                self.BodyForces(end).name = force_name;
                self.BodyForces(end).reference_frame = reference_frame;

                self.osimBodiesUsed{end+1} = osim_body_name;
                self.osimBodiesUsed{end+1} = reference_frame;

                self.meta_res(end+1).name = force_name;
                self.meta_res(end).type = 'Forces';
                self.meta_res(end).subtype = 'toExtFun';

            else
                self.res{idx} = self.res{idx} + value;
            end
        end % end of addBodyForce

        function [] = addBodyMoment(self,value,moment_name,osim_body_name,...
                reference_frame)
            arguments
                self Orthosis
                value (3,:)
                moment_name char
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                reference_frame char {inputExistsInOsimModel(self,reference_frame,'frame')} = osim_body_name;
            end
            % Add a moment vector acting on a body.
            %
            % EXAMPLE
            %   Add a moment on the tibia, and an equal but opposite moment
            %   on the hindfoot:
            %   orthosis_obj.addBodyMoment(T_ankle, 'T_exo_shank_r','tibia_r');
            %   orthosis_obj.addBodyMoment(-T_ankle, 'T_exo_foot_r','calcn_r','tibia_r']);
            %
            % INPUT:
            %   - value - [3x1 variable or double]
            %   * Value of the moment (in Nm)
            %
            %   - moment_name - [char]
            %   * Name of the moment.
            %
            %   - osim_body_name - [char]
            %   * Name of a body in the OpenSim model, on which the moment
            %   is applied.
            %
            %   - reference frame - [char] (optional) Default: osim_body_name
            %   * Name of a body in the OpenSim model, in whose reference
            %   frame the moment vector is expressed.
            %


            % value should be a 3xNmesh matrix
            if size(value,1)~=3 || size(value,2)~=self.Nmesh
                error('Expected "%s" to have size %ix%i, but the size is %ix%i.',...
                    inputname(2),3,self.Nmesh,size(value,1),size(value,2));
            end

            F_name = ['BodyMoment_',moment_name];

            idx = find(cellfun(@(x)strcmp(x,F_name), self.names_res));
            if isempty(idx)
                self.res{end+1} = value;
                self.names_res{end+1} = F_name;

                self.BodyMoments(end+1).body = osim_body_name;
                self.BodyMoments(end).name = moment_name;
                self.BodyMoments(end).reference_frame = reference_frame;

                self.osimBodiesUsed{end+1} = osim_body_name;
                self.osimBodiesUsed{end+1} = reference_frame;

                self.meta_res(end+1).name = moment_name;
                self.meta_res(end).type = 'Moments';
                self.meta_res(end).subtype = 'toExtFun';
            else
                self.res{idx} = self.res{idx} + value;
            end
        end % end of addBodyMoment

    %% analysis
        function [] = addVarToPostProcessing(self,var,label)
            arguments
                self Orthosis
                var
                label char
            end
            % Add internal variable to be evaluated in post-processing
            %   The values at each mesh point will be stored in 
            %   R.orthosis.separate{i}.(label)
            %
            % INPUT:
            %   - var - [variable]
            %   * Variable, or result of calculation based on variables.
            %
            %   - label - [char]
            %   * Label used to store the values of this variable.
            %

            self.res_pp{end+1} = var;
            self.names_res_pp{end+1} = label;
        end % end of addVarToPostProcessing


    %% OpenSim
        function [] = updatePropertiesFromOsimModel(self)
            % Read properties from OpenSim model and update the stored info

            if isempty(self.osimPath)
                error('Please use %s.setOsimPath to set the path to an OpenSim model file.',inputname(1));
            end

            import org.opensim.modeling.*;
            model = Model(self.osimPath);

            % coordinates
            n_coord = model.getCoordinateSet().getSize();
            osim_coords = cell(1,n_coord);
            for i=1:n_coord
                osim_coords{1,i} = char(model.getCoordinateSet().get(i-1).getName());
            end
            self.osimCoordsAll = osim_coords;

            % bodies
            n_bodies = model.getBodySet().getSize();
            osim_bodies = cell(1,n_bodies);
            for i=1:n_bodies
                osim_bodies{1,i} = char(model.getBodySet().get(i-1).getName());
            end
            self.osimBodiesAll = osim_bodies;

            % muscles
            n_mus = model.getMuscles().getSize();
            osim_muscles = cell(1,n_mus);
            for i=1:n_mus
                osim_muscles{1,i} = char(model.getMuscles().get(i-1).getName());
            end
            self.osimMusclesAll = osim_muscles;

            % contact spheres
            n_contact = model.getContactGeometrySet().getSize();
            osim_contacts = {};
            for i=1:n_contact
                contact_geom_i = model.getContactGeometrySet().get(i-1);
                if strcmp(contact_geom_i.getConcreteClassName(),'ContactSphere')
                    osim_contacts{1,end+1} = char(contact_geom_i.getName());
                end
            end
            osim_contacts{1,end+1} = 'left_total';
            osim_contacts{1,end+1} = 'right_total';
            self.osimContactsAll = osim_contacts;

        end % end of updatePropertiesFromOsimModel

    %% CasADi
        function [] = createCasadiFunction(self)
            % Create a CasADi Function for the added forces in function of the defined variables
            self.fun = casadi.Function(['f_orthosis_',self.name],...
                self.arg, [self.res, self.res_pp],...
                self.names_arg, [self.names_res, self.names_res_pp]);
        end

        function [wrap_fun, wrap_fun_pp] = wrapCasadiFunction(self,ExtFunIO,muscleNames)
            arguments
                self
                ExtFunIO
                muscleNames
            end
            % Wrap the CasADi Function of the orthosis for easy integration with PredSim
            %
            % INPUTS:
            %   - ExtFunIO - [struct]
            %   * model_info.ExtFunIO
            %
            %   - muscleNames - [cell array of chars]
            %   * model_info.muscle_info.muscle_names
            %
            % OUTPUTS:
            %   - wrap_fun - [casadi.Function]
            %   * CasADi Function for use in OCP formulation
            %
            %   - wrap_fun_pp - [casadi.Function]
            %   * CasADi Function for use in post-processing
            %

            import casadi.*

            if isempty(self.fun)
                self.createCasadiFunction();
            end

            % all inputs for wrapper function
            arg_SX.q = SX.sym('q',ExtFunIO.jointi.nq.all,self.Nmesh);
            arg_SX.qdot = SX.sym('qdot',ExtFunIO.jointi.nq.all,self.Nmesh);
            arg_SX.qddot = SX.sym('qddot',ExtFunIO.jointi.nq.all,self.Nmesh);
            arg_SX.act = SX.sym('a',length(muscleNames),self.Nmesh);
            arg_SX.fromExtFun = SX.sym('fromExtFun',ExtFunIO.nOutputs,self.Nmesh); % GRFs, point kinematics

            % all outputs for wrapper function
            res_SX.Mcoord = SX(ExtFunIO.jointi.nq.all,self.Nmesh);
            res_SX.toExtFun = SX(ExtFunIO.input.nInputs,self.Nmesh); % bodyforces, bodymoments

            % create inputs for function from inputs of wrapper function
            arg_fun = cell(1,self.fun.n_in);
            res_fun = cell(1,self.fun.n_out);
            for j=1:length(arg_fun)
                if strcmp(self.meta_arg(j).type,'muscle')
                    idx = find(strcmp(muscleNames,self.meta_arg(j).name));
                else
                    idx = ExtFunIO.(self.meta_arg(j).type).(self.meta_arg(j).name);
                end
                arg_fun{j} = arg_SX.(self.meta_arg(j).subtype)(idx,:);
            end

            % evaluate function
            if isempty(arg_fun)
                [res_0] = self.fun();
                for j=1:length(res_fun)
                    res_fun{j} = res_0.(self.fun.name_out(j-1));
                end
            else
                [res_fun{:}] = self.fun(arg_fun{:});
            end

            % assign outputs from function to outputs of wrapper function
            for j=1:length(self.meta_res)
                if strcmp(self.meta_res(j).type,'coordi')
                    idx = ExtFunIO.(self.meta_res(j).type).(self.meta_res(j).name);
                else
                    idx = ExtFunIO.input.(self.meta_res(j).type).(self.meta_res(j).name);
                end
                res_SX.(self.meta_res(j).subtype)(idx,:) = res_fun{j};
            end

            % collect outputs for post-processing function
            wrap_res_pp = res_fun(length(self.meta_res)+1:end);

            % input arguments for wrapper function
            arg_names = fieldnames(arg_SX);
            for i=1:length(arg_names)
                wrap_arg{i} = arg_SX.(arg_names{i});
            end

            % output arguments for wrapper function
            res_names = fieldnames(res_SX);
            for i=1:length(res_names)
                wrap_res{i} = res_SX.(res_names{i});
            end

            % create wrapper function
            wrap_fun = Function(['f_orthosis_',self.name,'_wrapped'],...
                wrap_arg, wrap_res, arg_names, res_names);

            wrap_fun_pp = Function(['f_orthosis_',self.name,'_wrapped_postprocessing'],...
                wrap_arg, wrap_res_pp, arg_names, self.names_res_pp);
            
        end % end of wrapCasadiFunction

    end % end of methods

%%
    methods(Access=private) % private — Access by class methods only (not from subclasses)

        function indentation = calcContactIndentation(self,contact_name)
            % get a variable representing the indentation of a contact sphere

            % open model
            if isempty(self.osimPath)
                error('Please use %s.setOsimPath to set the path to an OpenSim model file.',self.name);
            end
            import org.opensim.modeling.*;
            model = Model(self.osimPath);
            model.finalizeConnections();
            model.initSystem();
            
            % get contact sphere info
            csp1 = model.getContactGeometrySet().get(contact_name);
            csp1 = ContactSphere.safeDownCast(csp1);

            csp1_loc = csp1.getLocation().getAsMat;
            csp1_body = char(csp1.getBody.getName());
            
            csp1_r = csp1.getRadius();

            % get transformation matrix between ground and contact plane
            fl = model.getContactGeometrySet().get('floor');
            fl_tr = fl.getTransform();
            fl_T = sparsify(casadi.DM(fl_tr.T().getAsMat));
            fl_R = fl_tr.R();
            
            fl_rot = casadi.DM(4,4);
            for ii=1:3
                for jj=1:3
                    R_ij = fl_R.get(ii-1,jj-1);
                    if abs(R_ij) > eps
                        fl_rot(ii,jj) = R_ij;
                    end
                end
            end 
            fl_rot(4,1:3) = fl_T;
            fl_rot = sparsify(fl_rot);

            % add input variable for position of contact sphere centre
            pointPos = var_point(self,contact_name,csp1_body,csp1_loc,'pos');

            % express contact sphere centre relative to contact plane
            point = [pointPos; 1];
            point_inFl = fl_rot*point;
            
            % calculate indentation (negative means no contact)
            indentation = -(point_inFl(1,:) - csp1_r);
            
        end % end of getContactIndentation

    end % end of private methods

    methods(Access=protected, Hidden=true) 
        % protected — Access from methods in class or subclasses
        % hidden - do not show in documentation
        
        function inputExistsInOsimModel(self,input,type)
            % Assert that the input argument is valid for the OpenSim model

            switch type
                case 'coord'
                    osimComponents = self.osimCoordsAll;
                case 'body'
                    osimComponents = self.osimBodiesAll;
                case 'contact'
                    osimComponents = self.osimContactsAll;
                case 'muscle'
                    osimComponents = self.osimMusclesAll;
                case 'frame'
                    osimComponents = self.osimBodiesAll;
                    osimComponents{1,end+1} = 'ground';
                otherwise
                    osimComponents = {};
            end

            if ~isempty(osimComponents)
                try
                    mustBeMember(input,osimComponents)
                catch ex
                    throwAsCaller(ex);
                end
            elseif ~self.warningOsimPathNotSet
                self.warningOsimPathNotSet = true;
                warning(['Please use %s.setOsimPath to set ',...
                    'the path to an OpenSim model file to enable full ',...
                    'input validation'],self.name);
            end
        end % end of inputExistsInOsimModel

    end % end of protected methods


    methods(Hidden=true)
        % hidden - do not show in documentation

        % Hide methods inherited from handle so they are not in the
        % documentation for Orthosis
        function argout = eq(varargin)
            argout = eq@handle(varargin{:});
        end
        function argout = ge(varargin)
            argout = ge@handle(varargin{:});
        end
        function argout = gt(varargin)
            argout = gt@handle(varargin{:});
        end
        function argout = le(varargin)
            argout = le@handle(varargin{:});
        end
        function argout = lt(varargin)
            argout = lt@handle(varargin{:});
        end
        function argout = ne(varargin)
            argout = ne@handle(varargin{:});
        end
        function argout = listener(varargin)
            argout = listener@handle(varargin{:});
        end
        function argout = addlistener(varargin)
            argout = addlistener@handle(varargin{:});
        end
        function argout = notify(varargin)
            argout = notify@handle(varargin{:});
        end
        function argout = findobj(varargin)
            argout = findobj@handle(varargin{:});
        end
        function argout = findprop(varargin)
            argout = findprop@handle(varargin{:});
        end
    end

end % end of classdef