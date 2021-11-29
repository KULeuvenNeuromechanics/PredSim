% This function returns the ID given as inputs the path to the
% file containing the output of the ID and the joints of
% interest
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function ID = getID(pathID,joints)

load(pathID,'IDall')
ID.time = IDall.data(:,strcmp(IDall.colheaders,{'time'}));
ID.all(:,1) = ID.time;
count = 1;
for i = 1:length(joints)
    count = count + 1;
    if strcmp(joints{i},'pelvis_tx') || strcmp(joints{i},'pelvis_ty') || strcmp(joints{i},'pelvis_tz')
        ID.(joints{i}) = IDall.data(:,strcmp(IDall.colheaders,[joints{i},'_force']));
    else
        ID.(joints{i}) = IDall.data(:,strcmp(IDall.colheaders,[joints{i},'_moment']));
    end
    ID.all(:,count) = ID.(joints{i});
end