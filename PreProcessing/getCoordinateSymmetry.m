function [symQs] = getCoordinateSymmetry(S,osim_path,model_info)
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
%   - symQs -
%   * structure with fields QsInvA, QsInvB, QdotsInvA, QdotsInvB, and QsOpp
%   each containing the array of coordinate indices as indicated above.
%       
% 
% Original author: Lars D'Hondt
% Original date: 11/April/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% Exclude coordinate indices that correspond to limbs
coords_wo_limbs = [model_info.ExtFunIO.jointi.floating_base, model_info.ExtFunIO.jointi.torso];
% Names of selected coordinates
coord_names = model_info.ExtFunIO.coord_names.all(coords_wo_limbs);

% Helper variables to store indices while sorting
QsInv = [];
QsOpp = [];
Qs_forward = [];

% Initialise model
import org.opensim.modeling.*;
model = Model(osim_path);
state = model.initSystem;

% Get state variables vector
state_vars = model.getStateVariableValues(state);

% Get the reference frame of the rigid body that is presumably the trunk
idx_torso = model_info.ExtFunIO.jointi.torso(1);
coord_tmp = model.getCoordinateSet().get(model_info.ExtFunIO.coord_names.all{idx_torso});
joint_tmp = coord_tmp.getJoint();
body_frame = joint_tmp.getChildFrame().findBaseFrame();

% Define a left, right, and origin station in trunk reference frame
station_l = Vec3.createFromMat([0,0.5,-0.2]);
station_r = Vec3.createFromMat([0,0.5,0.2]);
station_o = Vec3.createFromMat([0,0,0]);

% Set all states to 0
state_vars.setToZero();
model.setStateVariableValues(state,state_vars);
model.realizePosition(state);

% Find the location of all 3 stations for the 0 state
loc_l_0 = body_frame.findStationLocationInGround(state,station_l).getAsMat;
loc_r_0 = body_frame.findStationLocationInGround(state,station_r).getAsMat;
loc_o_0 = body_frame.findStationLocationInGround(state,station_o).getAsMat;


%% Loop through selected coordinates
for i=1:length(coord_names)
    % Set all states to 0
    state_vars.setToZero();

    % Set position of current coordinate to 0.3 (rad or m)
    state_vars.set(model_info.ExtFunIO.coordi_OpenSimAPIstate.(coord_names{i}),0.3);
    model.setStateVariableValues(state,state_vars);
    model.realizePosition(state);
    % Find displacement of all 3 stations for the positive position value
    loc_l_pos = body_frame.findStationLocationInGround(state,station_l).getAsMat - loc_l_0;
    loc_r_pos = body_frame.findStationLocationInGround(state,station_r).getAsMat - loc_r_0;
    loc_o_pos = body_frame.findStationLocationInGround(state,station_o).getAsMat - loc_o_0;

    % Set position of current coordinate to 0.3 (rad or m)
    state_vars.set(model_info.ExtFunIO.coordi_OpenSimAPIstate.(coord_names{i}),-0.3);
    model.setStateVariableValues(state,state_vars);
    model.realizePosition(state);
    % Find displacement of all 3 stations for the negative position value
    loc_l_neg = body_frame.findStationLocationInGround(state,station_l).getAsMat - loc_l_0;
    loc_r_neg = body_frame.findStationLocationInGround(state,station_r).getAsMat - loc_r_0;
    loc_o_neg = body_frame.findStationLocationInGround(state,station_o).getAsMat - loc_o_0;

    % Evaluate the station displacement as a result of setting the current coordinate value
    if max([abs(loc_l_pos-loc_o_pos); abs(loc_r_pos-loc_o_pos); abs(loc_l_neg-loc_o_neg);...
            abs(loc_r_neg-loc_o_neg)])<eps && loc_l_pos(1) >=0.1
        % If the resulting displacement of all 3 stations is identical, and
        % nonzero along x, this coordinate defines forward motion.
        Qs_forward(end+1) = coords_wo_limbs(i);

    elseif max(abs([loc_l_pos(3);loc_l_neg(3);loc_r_pos(3);loc_r_neg(3)]))<eps
        % If there is no displacement normal to the sagittal plane, the
        % current coordinate value should be inverted.
        QsInv(end+1) = coords_wo_limbs(i);

    else
        % Remaining coordinates should be opposite to ensure continuity.
        QsOpp(end+1) = coords_wo_limbs(i);
    end
    
end

%% Make lists of coordinates
% Arm and leg coordinates are mapped to their symmetric counterpart.
% Remaining Inverse coordinates are mapped to themselves.
QsInvA = [QsInv,...
    model_info.ExtFunIO.jointi.leg_l,model_info.ExtFunIO.jointi.leg_r,...
    model_info.ExtFunIO.jointi.arm_l,model_info.ExtFunIO.jointi.arm_r]';
QsInvB = [QsInv,...
    model_info.ExtFunIO.jointi.leg_r,model_info.ExtFunIO.jointi.leg_l,...
    model_info.ExtFunIO.jointi.arm_r,model_info.ExtFunIO.jointi.arm_l]';

% The forward motion of the floating base is periodic in velocity, but not
% in position.
QdotsInvA = [Qs_forward';QsInvA];
QdotsInvB = [Qs_forward';QsInvB];

% There is no distinction between positions and velocities for the opposite
% coordinates
orderQsOpp = QsOpp';

% Assemble struct
symQs.QsInvA = QsInvA;
symQs.QsInvB = QsInvB;
symQs.QdotsInvA = QdotsInvA;
symQs.QdotsInvB = QdotsInvB;
symQs.orderQsOpp = orderQsOpp;


end