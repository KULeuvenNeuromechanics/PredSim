classdef Orthosis < handle
% --------------------------------------------------------------------------
% Orthosis
%   Interface to define custom orthosis or exoskeleton devices to be used
%   in predictive simulations. For examples on how to use this, see
%   parametricAFO, ankleExoZhang2017.
%
%   See also parametricAFO, ankleExoZhang2017 
% 
% Original author: Lars D'Hondt
% Original date: January 2024
% --------------------------------------------------------------------------

    properties(Access = private)
        % general properties
        name = []; % name of the orthosis
        Nmesh = 1; % number of meshpoints used to describe the time-varying 
            % behaviour of the orthosis. Set to 1 if time-independent.

        % properties of CasADi Function describing orthosis mechanics
        arg = {}; % input arguments
        res = {}; % output arguments
        names_arg = {}; % names of input arguments
        names_res = {}; % names of output arguments
        fun = []; % handle of CasADi Function

        % properties of interface with OpenSimAD. 
            % These will be passed to generateExternalFunction
        BodyForces = {}; % input3DBodyForces
        BodyMoments = {}; % input3DBodyMoments
        PointPositions = {}; % export3DPositions
        PointVelocities = {}; % export3DVelocities
        
        % properties of OpenSim model that is to be used with orthosis
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

    % constructor
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
            %
            % OUTPUT:
            %   - self - [Orthosis]
            %   * brief description of output_1
            
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

    % create CasADi Function
        function [] = createCasadiFunction(self)
            self.fun = casadi.Function(['f_ortosis_',self.name],...
                self.arg, self.res, self.names_arg, self.names_res);
        end

    % getters and setters
        function name = getName(self)
            name = self.name;
        end

        function N_mesh = getNmesh(self)
            N_mesh = self.Nmesh;
        end

        function fun = getFunction(self)
            if isempty(self.fun)
                createCasadiFunction(self);
            end
            fun = self.fun;
        end

        function BodyForces = getBodyForces(self)
            BodyForces = self.BodyForces;
        end

        function BodyMoments = getBodyMoments(self)
            BodyMoments = self.BodyMoments;
        end

        function PointPositions = getPointPositions(self)
            PointPositions = self.PointPositions;
        end

        function PointVelocities = getPointVelocities(self)
            PointVelocities = self.PointPositions;
        end

        function [] = setOsimPath(self,osimPath)
            self.osimPath = osimPath;
            updatePropertiesFromOsimModel(self);
        end

    % create a variable 
        % coordinate position, velocity or acceleration
        function coord = var_coord(self,osim_coord_name,pos_vel_acc)
            arguments
                self Orthosis
                osim_coord_name char {inputExistsInOsimModel(self,osim_coord_name,'osim_coord_name')}
                pos_vel_acc char {mustBeMember(pos_vel_acc,{'pos','vel','acc'})} = 'pos';
            end
            var_name = ['coord_',osim_coord_name,'_',pos_vel_acc];
            coord = casadi.SX.sym(var_name,1,self.Nmesh);
            self.arg{end+1} = coord;
            self.names_arg{end+1} = var_name;
            self.osimCoordsUsed{end+1} = osim_coord_name;
        end

        % point position or velocity
        function point = var_point(self,point_name,osim_body_name,location_in_body,pos_vel)
            arguments
                self Orthosis
                point_name char
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                location_in_body (1,3) double = [0, 0, 0];
                pos_vel char {mustBeMember(pos_vel,{'pos','vel'})} = 'pos';
            end
            var_name = ['point_',point_name,'_',pos_vel];
            point = casadi.SX.sym(var_name,3,self.Nmesh);
            self.arg{end+1} = point;
            self.names_arg{end+1} = var_name;
            self.osimBodiesUsed{end+1} = osim_body_name;
            switch pos_vel
                case 'pos'
                    self.PointPositions(end+1).body = osim_body_name;
                    self.PointPositions(end).point_in_body = location_in_body;
                    self.PointPositions(end).name = point_name;
                case 'vel'
                    self.PointVelocities(end+1).body = osim_body_name;
                    self.PointVelocities(end).point_in_body = location_in_body;
                    self.PointVelocities(end).name = point_name;
            end
        end

        % ground reaction force
        function GRF = var_GRF(self,osim_contact_name,F_d)
            % GRF_name can be "left_total", "right_total", or the name of a
            % contact sphere
            arguments
                self Orthosis
                osim_contact_name char {inputExistsInOsimModel(self,osim_contact_name,'contact')}
                F_d char {mustBeMember(F_d,{'Force','indentation'})} = 'Force';
            end
            var_name = ['GRF_',osim_contact_name,'_',F_d];
            switch F_d
                case 'Force'
                    GRF = casadi.SX.sym(var_name,3,self.Nmesh);
                    self.arg{end+1} = GRF;
                    self.names_arg{end+1} = var_name;
                    self.osimContactsUsed{end+1} = osim_contact_name;
                case 'indentation'
                    GRF = getContactIndentation(self,osim_contact_name);
            end
        end

        % muscle activity
        function act = var_act(self,osim_muscle_name)
            arguments
                self Orthosis
                osim_muscle_name char {inputExistsInOsimModel(self,osim_muscle_name,'muscle')}
            end
            var_name = ['muscle_',osim_muscle_name,'_act'];
            act = casadi.SX.sym(var_name,1,self.Nmesh);
            self.arg{end+1} = GRF;
            self.names_arg{end+1} = var_name;
            self.osimMusclesUsed{end+1} = osim_muscle_name;
        end

        % any variable
        function var = var(self,varargin)

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
        end


    % add a force/moment
        % coordinate force or moment
        function [] = addCoordForce(self,value,osim_coord_name)
            arguments
                self Orthosis
                value (1,:)
                osim_coord_name char {inputExistsInOsimModel(self,osim_coord_name,'coord')}
            end
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
            else
                self.res{idx} = self.res{idx} + value;
            end
        end

        % body force
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
            else
                self.res{idx} = self.res{idx} + value;
            end
        end

        % body moment
        function [] = addBodyMoment(self,value,force_name,osim_body_name,...
                reference_frame)
            arguments
                self Orthosis
                value (3,:)
                force_name char
                osim_body_name char {inputExistsInOsimModel(self,osim_body_name,'body')}
                reference_frame char {inputExistsInOsimModel(self,reference_frame,'frame')} = osim_body_name;
            end
            % value should be a 3xNmesh matrix
            if size(value,1)~=3 || size(value,2)~=self.Nmesh
                error('Expected "%s" to have size %ix%i, but the size is %ix%i.',...
                    inputname(2),3,self.Nmesh,size(value,1),size(value,2));
            end

            F_name = ['BodyMoment_',force_name];

            idx = find(cellfun(@(x)strcmp(x,F_name), self.names_res));
            if isempty(idx)
                self.res{end+1} = value;
                self.names_res{end+1} = F_name;

                self.BodyMoments(end+1).body = osim_body_name;
                self.BodyMoments(end).name = force_name;
                self.BodyMoments(end).reference_frame = reference_frame;

                self.osimBodiesUsed{end+1} = osim_body_name;
            else
                self.res{idx} = self.res{idx} + value;
            end
        end

        % any force
        function [] = addForce(self,varargin)

            n_in = length(varargin);

            idx = find(cellfun(@(x)(ischar(x)||isstring(x))&&contains(x,'CoordForce'), varargin));
            if ~isempty(idx)
                addCoordForce(self,varargin{setdiff(1:n_in,idx)});
                return
            end

            idx = find(cellfun(@(x)(ischar(x)||isstring(x))&&contains(x,'BodyForce'), varargin));
            if ~isempty(idx)
                addBodyForce(self,varargin{setdiff(1:n_in,idx)});
                return
            end

            idx = find(cellfun(@(x)(ischar(x)||isstring(x))&&contains(x,'BodyMoment'), varargin));
            if ~isempty(idx)
                addBodyMoment(self,varargin{setdiff(1:n_in,idx)});
                return
            end

            switch n_in
                case 3
                    addCoordForce(self,varargin{:});
                case 6
                    addBodyForce(self,varargin{:});
                case 5
                    addBodyMoment(self,varargin{:});
                otherwise
                    error('Unable to add an orthosis force based on these input arguments');
            end
        end

        % test that all references to components of the OpenSim model are valid
        function [] = testOsimModel(self)
            if isempty(self.osimPath)
                error('Please use %s.setOsimPath to set the path to an OpenSim model file.',inputname(1));
            end

            import org.opensim.modeling.*;
            model = Model(self.osimPath);

            errstr = [];

            % coordinates
            osim_coords = getUnique(self,self.osimCoordsUsed);
            errstr_coords = [];
            for name_osim = string(osim_coords)
                try
                    model.getCoordinateSet().get(name_osim);
                catch
                    errstr_coords = [errstr_coords,', ',char(name_osim)];
                end
            end
            if ~isempty(errstr_coords)
                errstr = [errstr, sprintf('\t Could not resolve coordinate names: %s\n',errstr_coords)];
            end

            % bodies
            osim_bodies = getUnique(self,self.osimBodiesUsed);
            errstr_bodies = [];
            for name_osim = string(osim_bodies)
                if strcmpi(name_osim,'ground')
                    continue
                end
                try
                    model.getBodySet().get(name_osim);
                catch
                    errstr_bodies = [errstr_bodies,', ',char(name_osim)];
                end
            end
            if ~isempty(errstr_bodies)
                errstr = [errstr, sprintf('\t Could not resolve body names: %s\n',errstr_bodies)];
            end

            % muscles
            osim_muscles = getUnique(self,self.osimMusclesUsed);
            errstr_muscles = [];
            for name_osim = string(osim_muscles)
                try
                    model.getMuscleSet().get(name_osim);
                catch
                    errstr_muscles = [errstr_muscles,', ',char(name_osim)];
                end
            end
            if ~isempty(errstr_muscles)
                errstr = [errstr, sprintf('\t Could not resolve muscle names: %s\n',errstr_muscles)];
            end

            % contact spheres
            osim_contacts = getUnique(self,self.osimContactsUsed);
            errstr_contacts = [];
            for name_osim = string(osim_contacts)
                if strcmpi(name_osim,'left_total') || strcmpi(name_osim,'right_total')
                    continue
                end
                try
                    model.getContactGeometrySet().get(name_osim);
                catch
                    errstr_contacts = [errstr_contacts,', ',char(name_osim)];
                end
            end
            if ~isempty(errstr_contacts)
                errstr = [errstr, sprintf('\t Could not resolve contact sphere names: %s\n',errstr_contacts)];
            end

            if ~isempty(errstr)
                error([sprintf(['Unable to resolve all OpenSim references for Orthosis %s',...
                ' based on the OpenSim model file "%s"\n'],inputname(1),self.osimPath),errstr]);
            end
        end

        function [] = updatePropertiesFromOsimModel(self)
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

        end

    end % end of methods

    methods(Access=private)

        % get only unique entrise of a cell array
        function cell_array = getUnique(self,cell_array)

            isUnique = true(size(cell_array));

            for ii = 1:length(cell_array)-1
                for jj = ii+1:length(cell_array)
                    if isequal(cell_array(ii),cell_array(jj))
                        isUnique(jj) = false;
                        break;
                    end
                end
            end
    
            cell_array(~isUnique) = [];
        end

        % get a variable representing the indentation of a contact sphere
        function indentation = getContactIndentation(self,contact_name)

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
            
        end

        function inputExistsInOsimModel(self,input,type)

%             disp(inputname(1))

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
                    'the path to an OpenSim model file to enable full',...
                    'input validation'],self.name);
            end
        end

    end % end of private methods

end % end of classdef