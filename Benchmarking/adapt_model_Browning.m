function [browning2008] = adapt_model_Browning(S,osim_path,metadata_only)
%adapt_model_Browning Adds mass to the opensim model based on the study of
%Browning 2008 (addref)
%   input arguments:
%       (1) S: default settings structure
%       (2) S: osim_path
%   metadata_only (optional): true returns model names/conditions without
%   loading OpenSim or writing models. Used by the ID debugging tool.
%   output arguments:
%       (1) output structure with path information about the exported
%       models

if nargin < 3
    metadata_only = false;
end
Addedmass = [2; 4; 8;...
    2; 4;...
    4; 8; 12; 16;...
    2; 4];
SegmentAdded = {'femur_l','femur_r';...
    'femur_l','femur_r';...
    'femur_l','femur_r';...
    'calcn_l','calcn_r';...
    'calcn_l','calcn_r';...
    'pelvis',[];...
    'pelvis',[];...
    'pelvis',[];...
    'pelvis',[];
    'tibia_l','tibia_r';...
    'tibia_l','tibia_r'};
COMlocation = {[0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]; ...
    [0 0 0], [0 0 0]}; % deviation in COM location from current COM location
ModelOut = {'_femur2kg'; '_femur4kg'; '_femur8kg';...
    '_foot2kg';'_foot4kg';'_pelvis_4kg';'_pelvis_8kg';...
    'pelvis_12kg';'_pelvis_16kg';'_tibia_2kg';'_tibia_4kg'};
ct = 1;
for i=1:length(Addedmass)
    browning2008.modelnames{ct} = [S.subject.name ModelOut{i}];
    model_name = browning2008.modelnames{ct};
    browning2008.osim_path{ct} = fullfile(S.misc.main_path,'Subjects',model_name,[model_name '.osim']);
    browning2008.added_mass{ct} = Addedmass(i) * sum(~cellfun(@isempty,SegmentAdded(i,:)));
    browning2008.location_added_mass{ct} = regexprep(SegmentAdded{i,1},'_l$','');
    if metadata_only
        ct = ct + 1;
        continue
    end

    % open model
    mSel = org.opensim.modeling.Model(osim_path);
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
    model_name = [S.subject.name ModelOut{i}];
    out_folder =  fullfile(S.misc.main_path,'Subjects',model_name);
    out_modelname = fullfile(out_folder,[model_name '.osim']);
    changed = save_benchmark_model(mSel,out_modelname);

    % update output structure
    browning2008.modelnames{ct} = model_name;
    browning2008.osim_path{ct} = out_modelname;
    browning2008.model_changed{ct} = changed;
    ct = ct+1;

    % save 
    clear mSel;
end


end
