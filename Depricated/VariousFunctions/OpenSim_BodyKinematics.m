function [BodyKin] = OpenSim_BodyKinematics(filename_model,results_dir,kinematics_file,varargin)
%OpenSim_BodyKinematics Uses the opensim API to run the bodykinematics
%analysis
%   Input arguments:
%       1) filename_model: full path to .osim model used for bodykinematics
%       analysis
%       2) results_dir: directory to save the results
%       3) kinematics_file: full path to .mot file with kinematics
%       4) variable input arguments:
%           4.1) array with start and end time for the analysis (by default
%           the start and end time of the IK file be used).
%           4.2) delete .mot file files from bodykinematics
%   Output arguments:
%       BodyKin.Pos: segment position information
%       BodyKin.Vel: segment velocity information
%       BodyKin.Acc: segment acceleration information
%       BodyKin.header: colheaders
%   Author:
%       Maarten Afschrift

% variable input arguments
if ~isempty(varargin) && ~isempty(varargin{1})
    event = varargin{1};
else
    IK = ReadMotFile(kinematics_file);
    event = [IK.data(1,1) IK.data(end,1)];
end

% Delete .mot files
BoolDeleteFile = false;
if length(varargin)>1
    BoolDeleteFile = varargin{2};
end

% test if the results directory is a folder
if ~isfolder(results_dir)
    mkdir(results_dir);
end

% import the opensim API
import org.opensim.modeling.*

% import the opensim model
osimModel=Model(filename_model);

% construct the bodykinematics analysis tool
bk = BodyKinematics();
bk.setStartTime(event(1));
bk.setEndTime(event(2));
bk.setOn(true);
bk.setStepInterval(1);
bk.setInDegrees(true);

% add the analysis tool to the model (this is needed if we want to add this
% to the analyzetool)
osimModel.addAnalysis(bk);
osimModel.initSystem();

% create the analysis tool
tool=AnalyzeTool(osimModel);
tool.setLoadModelAndInput(true);    % needed to updata analysis set from model   
tool.setResultsDir(results_dir);
tool.setInitialTime(event(1));
tool.setFinalTime(event(2));
[path,name,ext]=fileparts(kinematics_file);
tool.setName(name);

% run the analysis
tool.setCoordinatesFileName(kinematics_file);
tool.run();

% % import the output files
Acc = ReadMotFile(fullfile(results_dir,[name '_BodyKinematics_acc_global.sto']));
Vel = ReadMotFile(fullfile(results_dir,[name '_BodyKinematics_vel_global.sto']));
Pos = ReadMotFile(fullfile(results_dir,[name '_BodyKinematics_pos_global.sto']));
BodyKin.Pos = Pos.data;
BodyKin.Vel = Vel.data;
BodyKin.Acc = Acc.data;
BodyKin.header = Pos.names;

if BoolDeleteFile
    delete(fullfile(results_dir,[name '_BodyKinematics_acc_global.sto']));
    delete(fullfile(results_dir,[name '_BodyKinematics_vel_global.sto']));
    delete(fullfile(results_dir,[name '_BodyKinematics_pos_global.sto']));
end
end

