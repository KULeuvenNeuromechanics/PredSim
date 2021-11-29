%% Muscle analysis through API
%------------------------------

import org.opensim.modeling.*
modelpath = fullfile(pwd,'s1_Poggensee_Pinjoint_Exo.osim');
IKfile = fullfile(pwd,'s1Pog_dis4.mot');

IK = ReadMotFile(IKfile);
IK.data(:,[2:4 8:end]) =  IK.data(:,[2:4 8:end])*pi/180;



M = Model(modelpath);
x = M.initSystem();
muscles = M.getMuscles();
coords = M.getCoordinateSet();

m=1;
msel = muscles.get(m-1);        
LMT = msel.getLength(x);



iRow = 1;
IndexSel = nan(2,coords.getSize());
for i=1:coords.getSize()
    DOF = coords.get(i-1).getName();
    iSel = find(strcmp(char(DOF),IK.names));
    if ~isempty(iSel)
        IndexSel(1,i) = iSel;
        IndexSel(2,i) = i;
    end
end
iDel = find(isnan(IndexSel(2,:)));
IndexSel(:,iDel) = [];

nfr = 100;
LMT = nan(nfr,muscles.getSize());
dM = nan(nfr,muscles.getSize());
for iRow=1:nfr
    for i=IndexSel(2,:)
        Qval = IK.data(iRow,i+1);
        M.updCoordinateSet().get(i).setValue(x, Qval);
    end
    for m = 1:muscles.getSize()
        msel = muscles.get(m-1);        
        LMT(iRow,m) = msel.getLength(x);
        qs = M.updCoordinateSet().get(6);
        dM(iRow,m) =msel.computeMomentArm(x,qs);
    end
end




