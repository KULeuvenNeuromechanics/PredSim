function [symmetry, jointi] = identify_kinematic_chains(S,osim_path,model_info)
% --------------------------------------------------------------------------
% getCoordinateSymmetry
%   This function returns arrays of coordinate indices needed to formulate
%   the periodicity constraint in case of a half gait cycle. 
%   With coordinate positions "Qs" and velocities "Qdots", the periodicity
%   is implemented as:
%       opti.subject_to(Qs(QsInvA,end) - Qs(QsInvB,1) == 0);
%       opti.subject_to(Qdots(QdotsInvA,end) - Qdots(QdotsInvB,1) == 0);
%       opti.subject_to(Qs(QsOpp,end) + Qs(QsOpp,1) == 0);
%       opti.subject_to(Qdots(QsOpp,end) + Qdots(QsOpp,1) == 0);
%   
%   Coordinates of left arm and leg are mapped to their corresponding
%   coordinate on the right side, and vice versa. The coordinates that do
%   not belong to a limb are assigned to Inv or Opp by evaluating their
%   effect on displacement. 
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - symmetry -
%   * structure with fields QsInvA, QsInvB, QdotsInvA, QdotsInvB, and QsOpp
%   each containing the array of coordinate indices as indicated above.
%       
%   - jointi -
%   * update structure with fields containing the kinematic chain for each limb
% 
% Original author: Lars D'Hondt
% Original date: 17/April/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import org.opensim.modeling.*;


model = Model(osim_path);
state = model.initSystem;
jointSet = model.getJointSet();
n_joints = jointSet.getSize();

jointi = model_info.ExtFunIO.jointi;

%% Create table
joint_struct = [];
for i=1:n_joints
    joint_i = jointSet.get(i-1);
    if ~isfield(jointi,char(joint_i.getName()))
        % skip joints that are not in jointi
        continue
    end
    joint_struct(end+1).joint = char(joint_i.getName());
    joint_struct(end).parent = char(joint_i.getParentFrame().findBaseFrame().getName());
    joint_struct(end).child = char(joint_i.getChildFrame().findBaseFrame().getName());

    crd = [];
    for j=1:joint_i.numCoordinates()
        crd{j} = char(joint_i.get_coordinates(j-1).getName());
    end
    joint_struct(end).coordinates = crd;

end

joint_table = struct2table(joint_struct);

%% Find limbs
base_legs = S.subject.base_joints_legs;
base_arms = S.subject.base_joints_arms;
if ~iscell(base_legs) && ~isempty(base_legs)
    base_legs = {base_legs};
end
if ~iscell(base_arms) && ~isempty(base_arms)
    base_arms = {base_arms};
end

limbs.all.joints = [];
limbs.all.coords = [];
leg_l = [];
leg_r = [];
arm_l = [];
arm_r = [];

for i=1:length(base_legs)
    base_joint_idx = find(strcmp(joint_table.joint,base_legs{i}));
    if isempty(base_joint_idx)
        base_joint_idx = find(strcmpi(joint_table.joint,[base_legs{i} '_r'])|strcmpi(joint_table.joint,['r_' base_legs{i}]));
        if isempty(base_joint_idx)
            error([base_legs{i} ' is not recognised as a joint name.']);
        end
    end
    [base_joint_l, base_joint_r] = mirrorName(joint_table.joint{base_joint_idx});

    [joints_limb_i, coords_limb_i] = getKinematicChain(joint_table,base_joint_r);

    coordi = nan(length(coords_limb_i),1);
    for j=1:length(coords_limb_i)
        coordi(j) = model_info.ExtFunIO.coordi.(coords_limb_i{j});
    end
    
    jointi.(['leg' num2str(i) '_r']) = coordi';
    leg_r = [leg_r coordi];

    limbs.all.joints = [limbs.all.joints; joints_limb_i];
    limbs.all.coords = [limbs.all.coords; coords_limb_i];
    
    [joints_limb_i, coords_limb_i] = getKinematicChain(joint_table,base_joint_l);

    coordi = nan(length(coords_limb_i),1);
    for j=1:length(coords_limb_i)
        coordi(j) = model_info.ExtFunIO.coordi.(coords_limb_i{j});
    end

    jointi.(['leg' num2str(i) '_l']) = coordi';
    leg_l = [leg_l coordi];

    limbs.all.joints = [limbs.all.joints; joints_limb_i];
    limbs.all.coords = [limbs.all.coords; coords_limb_i];
end

for i=1:length(base_arms)
    base_joint_idx = find(strcmp(joint_table.joint,base_arms{i}));
    if isempty(base_joint_idx)
        base_joint_idx = find(strcmpi(joint_table.joint,[base_arms{i} '_r'])|strcmpi(joint_table.joint,['r_' base_arms{i}]));
        if isempty(base_joint_idx)
            error(['"' base_arms{i} '" is not recognised as a joint name on the right side.']);
        end
    end
    [base_joint_l, base_joint_r] = mirrorName(joint_table.joint{base_joint_idx});

    [joints_limb_i, coords_limb_i] = getKinematicChain(joint_table,base_joint_r);

    coordi = nan(length(coords_limb_i),1);
    for j=1:length(coords_limb_i)
        coordi(j) = model_info.ExtFunIO.coordi.(coords_limb_i{j});
    end
    
    jointi.(['arm' num2str(i) '_r']) = coordi';
    arm_r = [arm_r coordi];

    limbs.all.joints = [limbs.all.joints; joints_limb_i];
    limbs.all.coords = [limbs.all.coords; coords_limb_i];
    
    [joints_limb_i, coords_limb_i] = getKinematicChain(joint_table,base_joint_l);

    coordi = nan(length(coords_limb_i),1);
    for j=1:length(coords_limb_i)
        coordi(j) = model_info.ExtFunIO.coordi.(coords_limb_i{j});
    end
    
    jointi.(['arm' num2str(i) '_l']) = coordi';
    arm_l = [arm_l coordi];

    limbs.all.joints = [limbs.all.joints; joints_limb_i];
    limbs.all.coords = [limbs.all.coords; coords_limb_i];
end

jointi.leg_r = leg_r';
jointi.leg_l = leg_l';
jointi.arm_r = arm_r';
jointi.arm_l = arm_l';

%% Find floating base
jointi.floating_base = [];
jointi.base_forward = [];
jointi.base_vertical = [];
jointi.base_lateral = [];

for i=1:size(joint_table,1)

    idx = find(contains(joint_table.child,joint_table.parent(i)));
    if isempty(idx)
        floating_base = jointSet.get(joint_table.joint(i));
        continue
    end
end

fl_b = floating_base.getConcreteClassName();
if strcmp(fl_b,'CustomJoint')
    floating_base = CustomJoint.safeDownCast(floating_base);

    sptr = floating_base.getSpatialTransform();

    for j=1:6
    
        tr1 = sptr.getTransformAxis(j-1);
        ax_j = tr1.get_axis().getAsMat;
        CoordNames = tr1.getCoordinateNames;
        if ~strcmp(char(CoordNames),'()')
            crd_j = char(tr1.get_coordinates(0));
            istr = strcmp(model.getCoordinateSet.get(crd_j).getMotionType(),'Translational');
    
            jointi.floating_base(end+1) = model_info.ExtFunIO.coordi.(crd_j);
            if istr
                if ax_j(1) == 1
                    jointi.base_forward = model_info.ExtFunIO.coordi.(crd_j);
                elseif ax_j(2) == 1
                    jointi.base_vertical = model_info.ExtFunIO.coordi.(crd_j);
                elseif ax_j(3) == 1
                    jointi.base_lateral = model_info.ExtFunIO.coordi.(crd_j);
                end
            end
        end
    
    end

elseif strcmp(fl_b,'PlanarJoint')
    floating_base = PlanarJoint.safeDownCast(floating_base);

    % rotation around z
    crd_rz = char(floating_base.get_coordinates(0).getName());
    jointi.floating_base(end+1) = model_info.ExtFunIO.coordi.(crd_rz);

    % translation along x
    crd_tx = char(floating_base.get_coordinates(1).getName());
    jointi.base_forward = model_info.ExtFunIO.coordi.(crd_tx);
    jointi.floating_base(end+1) = model_info.ExtFunIO.coordi.(crd_tx);

    % translation along y
    crd_ty = char(floating_base.get_coordinates(2).getName());
    jointi.base_vertical = model_info.ExtFunIO.coordi.(crd_ty);
    jointi.floating_base(end+1) = model_info.ExtFunIO.coordi.(crd_ty);

end





jointi.torso = setdiff(1:numel(fields(model_info.ExtFunIO.coordi)),...
    [jointi.floating_base,jointi.leg_r,jointi.leg_l,jointi.arm_r,jointi.arm_l]);

%% Find symmetry
joints_not_limbs = setdiff(joint_table.joint,limbs.all.joints);
coords_inverse = [];
coords_opposite = [];

for i=1:length(joints_not_limbs)

    % Set all states to 0
    state_vars = model.getStateVariableValues(state);
    state_vars.setToZero();
    model.setStateVariableValues(state,state_vars);
    model.realizePosition(state);

    body_frame = jointSet.get(joints_not_limbs{i}).getChildFrame().findBaseFrame();

    % Define 2 stations, symmetric w.r.t. global reference frame
    body_orig_in_gnd = body_frame.getPositionInGround(state).getAsMat;
    
    station_l_in_gnd = Vec3.createFromMat(body_orig_in_gnd + [0.2; 0.2; -0.2]);
    station_r_in_gnd = Vec3.createFromMat(body_orig_in_gnd + [0.2; 0.2; 0.2]);
    station_l = model.getGround.findStationLocationInAnotherFrame(state,station_l_in_gnd,body_frame);
    station_r = model.getGround.findStationLocationInAnotherFrame(state,station_r_in_gnd,body_frame);

    coords_i = joint_table.coordinates(strcmp(joint_table.joint,joints_not_limbs{i}));
    coords_i = coords_i{1};

    for j=1:length(coords_i)

        % Set all states to 0
        state_vars = model.getStateVariableValues(state);
        state_vars.setToZero();
        model.setStateVariableValues(state,state_vars);
        model.realizePosition(state);

        
        % Find the location of all 3 stations for the 0 state
        loc_l_0 = body_frame.findStationLocationInGround(state,station_l).getAsMat;
        loc_r_0 = body_frame.findStationLocationInGround(state,station_r).getAsMat;

        model.getCoordinateSet.get(coords_i{j}).setValue(state,0.3);
        model.realizePosition(state);
        loc_l_p = body_frame.findStationLocationInGround(state,station_l).getAsMat - loc_l_0;
        loc_r_p = body_frame.findStationLocationInGround(state,station_r).getAsMat - loc_r_0;

        model.getCoordinateSet.get(coords_i{j}).setValue(state,-0.3);
        model.realizePosition(state);
        loc_l_n = body_frame.findStationLocationInGround(state,station_l).getAsMat - loc_l_0;
        loc_r_n = body_frame.findStationLocationInGround(state,station_r).getAsMat - loc_r_0;

        loc_m_n = loc_l_n;
        loc_m_n(3) = -loc_m_n(3);
        loc_m_p = loc_l_p;
        loc_m_p(3) = -loc_m_p(3);

        t_inv = [loc_r_n(3), loc_r_p(3), loc_l_n(3), loc_l_p(3)]; % inverse: no effect on z
%         t_inv = [loc_l_n-loc_r_n; loc_l_p-loc_r_p]; % inverse: left and right loc are the same for a motion
        t_opp = [loc_r_p-loc_m_n; loc_r_n-loc_m_p]; % opposite: 


        if max(abs(t_inv)) < 1e-6
            coords_inverse{end+1,1} = coords_i{j};

        elseif max(abs(t_opp)) < 1e-6
            coords_opposite{end+1,1} = coords_i{j};
        else
            warning(['Could not identify symmetry for coordinate "' coords_i{j} '".'])

        end
    end

end

QsInv = [];
QsOpp = [];
for j=1:length(coords_inverse)
    QsInv(j,1) = model_info.ExtFunIO.coordi.(coords_inverse{j});
end
for j=1:length(coords_opposite)
    QsOpp(j,1) = model_info.ExtFunIO.coordi.(coords_opposite{j});
end
QsInv = setdiff(QsInv,jointi.base_forward);

%% Make lists of coordinates
% Arm and leg coordinates are mapped to their symmetric counterpart.
% Remaining Inverse coordinates are mapped to themselves.
QsInvA = [QsInv; leg_l; leg_r; arm_l; arm_r];
QsInvB = [QsInv; leg_r; leg_l; arm_r; arm_l];

% The forward motion of the floating base is periodic in velocity, but not
% in position.
QdotsInvA = [jointi.base_forward';QsInvA];
QdotsInvB = [jointi.base_forward';QsInvB];

% Assemble struct
symQs.QsInvA = QsInvA;
symQs.QsInvB = QsInvB;
symQs.QdotsInvA = QdotsInvA;
symQs.QdotsInvB = QdotsInvB;
symQs.QsOpp = sort(QsOpp)';

symmetry = symQs;

end
%%
function [joints_chain, coords_chain] = getKinematicChain(joint_table,starting_joint)
% --------------------------------------------------------------------------
% getKinematicChain
%   Get the names of the joints and coordinates that belong to the
%   kinematic chain starting at a given joint.
% 
%
%   See also identify_kinematic_chains
%
% INPUT:
%   - joint_table -
%   * Table with a row for each joint and columns:
%       - joint: name of the joint
%       - child: name of the body that has the child frame of the joint
%       - parent: name of the body that has the paret frame of the joint
%       - coordinates: cell array with names of all coordinates of the joint
% 
%   - starting_joint -
%   * Name of the joint that acts as the start of the kinematic chain.
%
%
% OUTPUT:
%   - joints_chain -
%   * Cell array with names of the joints of the kinematic chain.
%
%   - coords_chain -
%   * Cell array with names of the coordinates of the kinematic chain.
% 
% Original author: Lars D'Hondt
% Original date: May/2024
% --------------------------------------------------------------------------

    joints_chain = [];
    coords_chain = [];
    has_next_joint = true;
    next_joints_array = {starting_joint};
    if ~contains(joint_table.joint,next_joints_array)
        return
    end
    while has_next_joint
        % take joint to analyse from FIFO array
        next_joint = next_joints_array{1};
        % add joint to kinematic chain
        joints_chain{end+1,1} = next_joint;
        % add coordinates to kinematic chain
        coords_i = joint_table.coordinates(strcmp(joint_table.joint,next_joint));
        coords_chain = [coords_chain(:); coords_i{:}];
        % get the child frame of this joint
        next_frame = joint_table.child(strcmp(joint_table.joint,next_joint));
        % get joints that have this frame as their parent frame and add
        % them to FIFO array
        next_joint_i = joint_table.joint(strcmp(joint_table.parent,next_frame));
        next_joints_array = [next_joints_array(:); next_joint_i(:)];
        % remove 1st joint from array, since that one is done
        next_joints_array = next_joints_array(2:end);
        % set stop criterion if there ar no more joints to analyse
        if isempty(next_joints_array)
            has_next_joint = false;
        end
    end
end