function [] = OpenSim_Muscle_Analysis(motion_file,model_sel,output_path,event,varargin)
% --------------------------------------------------------------------------
% OpenSim_Muscle_Analysis 
%   Executes a muscle analsysis from the command line
%   Note: check compatibility with Opensim 4.0
%   
% INPUT:
%   - motion_file -
%   * 
%
%   - model_sel -
%   * 
% 
%   - output_path -
%   * 
%
%   - event -
%   * 
% 
%   - varargin -
%   * 
%
% OUTPUT:
%   * no outputs
% 
% Original author: 
% Original date: 
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


% run only muscle analysis for the selected DOFs
bool_Limit = 0;
if  ~isempty(varargin)
    CoordNames = varargin{1};
    bool_Limit = 1;
end

import org.opensim.modeling.*

% get the analysis tool and change variables
[FunctionPath,~]=fileparts(mfilename('fullpath'));
path_generic_file=fullfile(FunctionPath,'settings_Muscle_analysis.xml');
tool=AnalyzeTool(path_generic_file,false);
tool.setLoadModelAndInput(true)
osimModel=Model(model_sel);
tool.setModel(osimModel);
tool.setResultsDir(output_path);
tool.setInitialTime(event(1));
tool.setFinalTime(event(2));
[~, name, ~]=fileparts(motion_file);
tool.setName(name);
if bool_Limit
    mAnalysis =tool.getAnalysisSet().get(0);
    mA = MuscleAnalysis.safeDownCast(mAnalysis);
    as = ArrayStr();
    for i=1:length(CoordNames)
        as.set(i-1,CoordNames{i});
    end
    mA.setCoordinates(as);
    mA.setStartTime(event(1));
    mA.setEndTime(event(2));
end

% export the analysis
tool.setModelFilename(model_sel);
tool.setCoordinatesFileName(motion_file);
out_path_xml=fullfile(output_path,['muscle_analysis_' name '.xml']);
tool.print(out_path_xml);

% run the tool
tool2 = AnalyzeTool(out_path_xml);
tool2.print([out_path_xml(1:end-4) '_final.xml']);
tool2.run();


end

