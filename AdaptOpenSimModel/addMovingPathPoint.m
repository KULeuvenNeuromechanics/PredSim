function [model] = addMovivingPathPoint(model,muscle,name,parent,joint,poly)
% --------------------------------------------------------------------------
%addMovingPathPoint 
%   This function adds a path point to a osim model muscle.
% 
% INPUT:
%   -model-
%   * org.opensim.modeling.Model
%
%   -muscle-
%   * char name of muscle in which you want to add a path point
%
%   -name-
%   * char name of the path point to add
%
%   -parent-
%   * parent to which path point will be connected
%
%   -joint-
%   * joint to which the point has a relation
%
%   -poly-
%   * struct with 3 fields (x, y and z), each containing the factors for
%   the polynomial
% 
% OUTPUT:
%   -model-
%   * updated org.opensim.modeling.Model
% 
% Original author: Bram Van Den Bosch
% Original date: 20/05/2022
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 20/05/2022
% --------------------------------------------------------------------------
% clear all

import org.opensim.modeling.*
% 
% %Create the Original OpenSim model from a .osim file
% rootdir_model = 'C:\Users\u0138016\OneDrive - KU Leuven\SimCP_2\Subjects\CP3\T0\Model';
% gen_model = 'SimCP_2_to_scale.osim'; % base model
% filename = fullfile(rootdir_model,gen_model);
% model = Model(filename);
% state = model.initSystem;
% 
% muscle = 'rect_fem_r';
% parent = 'tibia_r';
% name = 'test';
% joint = 'knee_r';
% poly.x = [1 1 1];
% poly.y = [0.2 0.5 0.4];
% poly.z = [1 1 1];

%%
currentMuscle = model.getMuscles.get(muscle);
currentParent = model.getBodySet.get(parent);
currentJoint  = model.getJointSet.get(joint);

% get absolute path to parent body and joint coordinate
parent_path = currentParent.getAbsolutePathString;
coordinate_path = currentJoint.getCoordinate.getAbsolutePathString;

% create MovingPathPoint
pathpoint = MovingPathPoint;
pathpoint.setName(name);

% set parent and joint coordinate names
pathpoint.getSocket('parent_frame').setConnecteePath(parent_path);
pathpoint.getSocket('x_coordinate').setConnecteePath(coordinate_path);
pathpoint.getSocket('y_coordinate').setConnecteePath(coordinate_path);
pathpoint.getSocket('z_coordinate').setConnecteePath(coordinate_path);

% write polynomial values to locations
polyfunc = PolynomialFunction;
polyfunc.setName('function');
coeff = Vector.createFromMat(poly.x);
polyfunc.setCoefficients(coeff);
poly_dc = PolynomialFunction.safeDownCast(polyfunc);
pathpoint.set_x_location(poly_dc);

polyfunc = PolynomialFunction;
polyfunc.setName('function');
coeff = Vector.createFromMat(poly.y);
polyfunc.setCoefficients(coeff);
poly_dc = PolynomialFunction.safeDownCast(polyfunc);
pathpoint.set_y_location(poly_dc);

polyfunc = PolynomialFunction;
polyfunc.setName('function');
coeff = Vector.createFromMat(poly.z);
polyfunc.setCoefficients(coeff);
poly_dc = PolynomialFunction.safeDownCast(polyfunc);
pathpoint.set_z_location(poly_dc);

% add pathpoint to muscle
pp_dc = MovingPathPoint.safeDownCast(pathpoint);
currentMuscle.getGeometryPath.updPathPointSet.adoptAndAppend(pp_dc);

% finalize
model.updForceSet;
% model.finalizeConnections;

end