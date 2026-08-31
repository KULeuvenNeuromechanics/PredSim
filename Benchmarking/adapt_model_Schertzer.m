function [schertzer2014] = adapt_model_Schertzer(S,osim_path)
%adapt_model_Browning Adds mass to the opensim model based on the study of
%Browning 2008 (addref)
%   input arguments:
%       (1) S: default settings structure
%       (2) S: osim_path
%   output arguments:
%       (1) output structure with path information about the exported
%       models

import org.opensim.modeling.*;
Addedmass = [0.5; 1; 2;...
    0.5; 1; 2;...
    2; 7; 10; 16; 22;...
    0];

lTibia = 0.39; % this is hard coded for now. Get this from the opensim model in the future
lFemur = 0.34;% this is hard coded for now. Get this from the opensim model in the future

SegmentAdded = {'tibia_l','tibia_r';...
    'tibia_l','tibia_r';...
    'tibia_l','tibia_r';...
    'femur_l','femur_r';...
    'femur_l','femur_r';...
    'femur_l','femur_r';...
    'torso',[];...
    'torso',[];...
    'torso',[];...
    'torso',[];...
    'torso',[];...
    'torso',[]};
COMlocation = {[0 -lTibia 0], [0 -lTibia 0]; ...
    [0 -lTibia 0], [0 -lTibia 0];...
    [0 -lTibia 0], [0 -lTibia 0];...
    [0 -lFemur 0], [0 -lFemur 0]; ...
    [0 -lFemur 0], [0 -lFemur 0];...
    [0 -lFemur 0], [0 -lFemur 0];
    [0 0 0], [];...
    [0 0 0], [];...
    [0 0 0], [];...
    [0 0 0], [];...
    [0 0 0], [];...
    [0 0 0], []}; % deviation in COM location from current COM location

    ModelOut = {'Ankle_05';'Ankle_1';'Ankle_2';...
        'Knee_05';'Knee_1';'Knee_2';
        'Torso_2';'Torso_7';'Torso_10'; 'Torso_16'; 'Torso_22';...
        'Fall22_def'};
    location_addedmass_str = {'ankle','ankle','ankle',...
        'knee','knee','knee',...
        'torso','torso','torso','torso','torso',''};
    
ct = 1;
for i=1:length(Addedmass)
    % open model
    mSel = Model(osim_path);
    for j =1:2 % left and right leg
        segSel = SegmentAdded{i,j};
        if ~isempty(segSel)
            % added mass, COM lcation
            madd = Addedmass(i);
            madd_COM = COMlocation{i,j};
            madd_Inertia = [0, 0, 0];

            % adapt the segment mass
            BodySel = mSel.getBodySet().get(segSel);
            mBodySel = BodySel.getMass();
            BodySel.setMass(mBodySel + madd);
            % get COM location and inertia
            I = BodySel.getInertia;
            I_diag = I.getMoments;
            I_diag_mat = [I_diag.get(0) I_diag.get(1) I_diag.get(2)];
            COM = BodySel.getMassCenter();
            COM_mat = [COM.get(0) COM.get(1) COM.get(2)];
            % compute new COM location and inertia
            [COM_new,I_new] = AdaptSegmentCOMAndInetia(COM_mat,I_diag_mat,mBodySel,...
                COM_mat+madd_COM,madd_Inertia, madd);
            % update the model
            for ii=1:3
                I_diag.set(ii-1,I_new(ii));
                COM.set(ii-1,COM_new(ii))
            end
            BodySel.setInertia(I);
            BodySel.setMassCenter(COM);
        end
    end
    % init the model (not sure if this is needed)
    mSel.initSystem();
    
    % save model
    model_name = [S.subject.name '_' ModelOut{i}];
    out_folder =  fullfile(S.misc.main_path,'Subjects',model_name);
    out_modelname = fullfile(out_folder,[model_name '.osim']);
    if ~isfolder(out_folder)
        mkdir(out_folder)
    end
    mSel.print(out_modelname);

    % update output structure
    schertzer2014.modelnames{ct} = model_name;
    schertzer2014.osim_path{ct} = out_modelname;
    schertzer2014.added_mass{ct} = Addedmass(i);
    schertzer2014.location_added_mass{ct} =location_addedmass_str{i};
    ct = ct+1;

    % save 
    clear mSel;
end


end