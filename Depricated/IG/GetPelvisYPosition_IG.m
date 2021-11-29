
%% Debug pelvis Y position
import casadi.*;

% external function
cd('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\ExternalFunctions');
S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp 
F  = external('F',S.ExternalFunc2);
cd('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\IG');
% indexes
GRFi = [32:37];
jointi.pelvis.ty    = 5;

PyVect = 0.88:0.001:0.92;
OutVect = nan(43,length(PyVect));
for i=1:length(PyVect)
    Input = [zeros(31,1); zeros(31,1) ;zeros(31,1)];
    Input(jointi.pelvis.ty*2-1) = PyVect(i);
    Out = F(Input);
    OutVect(:,i) = full(Out);
end

figure();
FyPelvis = OutVect(jointi.pelvis.ty,:);
plot(PyVect,FyPelvis);
xlabel('Pelvis position');
ylabel('Residual Fy');
PelvisPos = PyVect(FyPelvis>=0);
disp(['Optimal pelvis position for this subjects is ' num2str(PelvisPos(1))]);

