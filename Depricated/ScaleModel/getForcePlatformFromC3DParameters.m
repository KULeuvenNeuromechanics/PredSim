function f = getForcePlatformFromC3DParameters(ParameterGroup)
% Purpose:  Retrieve FORCE_PLATFORM information from a C3D file            
%
% Usage:    f = getForcePlatformFromC3DParameters(C3DParameterGroup) 
%
% Input:    C3DParamterGroup is the read in C3D file structure containing 
%           data collection parameters 
%
% Output:   f returns a structure with the following format: 
%             f(FP#).corners(corner#, 3) - stores xyz coordinates of the 
%                   4 corners of each FP, in the lab CS
%             f(FP#).origin - stores xyz elements of the vector from the
%                   center of the FP surface to the FP origin, in the
%                   FP coordinate system (Type 2, AMTI plates).
%               
% ASeth 09-07 adapted from ASA 9-05

% Retrieve number of FPs from FORCE_PLATFORM:USED.
% Get details about analog data
I = 1; J = 1;
while ~strcmp(ParameterGroup(I).name, 'FORCE_PLATFORM');
        I = I+1;
end
   
for J = 1:length(ParameterGroup(I).Parameter),
    category = upper(ParameterGroup(I).Parameter(J).name);
    switch category{1}
        case 'USED'
            nFP = ParameterGroup(I).Parameter(J).data;
        case 'CORNERS'
            corners = ParameterGroup(I).Parameter(J).data;
        case 'ORIGIN'
            origin = ParameterGroup(I).Parameter(J).data;  
        case 'CHANNEL'
            channel= ParameterGroup(I).Parameter(J).data; 
    end            
end

for fpIndex = 1:nFP
    f(fpIndex).corners = corners(:,:,fpIndex)';
    f(fpIndex).origin = origin(:,fpIndex)';
    f(fpIndex).channel = channel(:,fpIndex);
end
