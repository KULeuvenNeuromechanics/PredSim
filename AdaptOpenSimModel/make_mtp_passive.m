% --------------------------------------------------------------------------
% make_mtp_passive
%   This script adapts an OpenSim model to use for simulations with a
%   passive mtp joint. The muscle pathpoints of long toe flexors and
%   extensors are fixed in calcaneus frame to remove their interaction with
%   the mtp joint.
% 
% 
% Original author: Lars D'Hondt
% Original date: 09/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


clear
clc

import org.opensim.modeling.*;
[pathHere,~,~] = fileparts(mfilename('fullpath'));

osim_filename = 'Fal_s1';

model = Model([pathHere '\' osim_filename '.osim']);
s = model.initSystem;




calcn = model.getBodySet().get('calcn_r').findBaseFrame().getPositionInGround(s).getAsMat;
toes = model.getBodySet().get('toes_r').findBaseFrame().getPositionInGround(s).getAsMat;


muscleNames = {'flex_dig_r','flex_hal_r','ext_dig_r','ext_hal_r'};
N_mus = length(muscleNames);
for i=1:N_mus
    muscleNames{N_mus+i} = [muscleNames{i}(1:end-1) 'l'];
end
N_mus = length(muscleNames);

Muscles = model.getMuscles();

shift_toes_2_calcn_r = calcn - toes;
shift_toes_2_calcn_l = shift_toes_2_calcn_r;
shift_toes_2_calcn_l(3) = -shift_toes_2_calcn_l(3);

calcn_body_r = model.getBodySet().get('calcn_r');
calcn_body_l = model.getBodySet().get('calcn_l');


path_points = {};
for i=1:N_mus
    path_points_i = Muscles.get(muscleNames{i}).getGeometryPath().getPathPointSet();
    N_pp = path_points_i.getSize();
    remove_pp = zeros(N_pp);
    for j=1:N_pp
        N_pp_new = N_pp;
        path_point_j = path_points_i.get(j-1);
        parent_socket = char(path_point_j.getBodyName());

        name_j = char(path_point_j.getName());
        loc_old = path_point_j.getLocation(s).getAsMat();
        ppt = PathPoint.safeDownCast(path_point_j);
        if strcmp(parent_socket,'toes_r')
            loc_new = loc_old - shift_toes_2_calcn_r;
            loc_osim = Vec3.createFromMat(loc_new);
            ppt.setLocation(loc_osim);
            ppt.setBody(calcn_body_r);
            ppt.connectSocket_parent_frame(calcn_body_r);
            upd=1;
        elseif strcmp(parent_socket,'toes_l')
            loc_new = loc_old - shift_toes_2_calcn_l;
            loc_osim = Vec3.createFromMat(loc_new);
            ppt.setLocation(loc_osim);
            ppt.setBody(calcn_body_l);
            upd=1;
        else
            upd=0;
        end
        if upd
            path_points{end+1,1} = name_j;
            path_points{end,2} = parent_socket;
            path_points{end,3} = loc_new;
        end

    end

end


model.finalizeConnections();
model.initSystem;
model.print([pathHere '\' osim_filename '_passive_mtp.osim']);


















