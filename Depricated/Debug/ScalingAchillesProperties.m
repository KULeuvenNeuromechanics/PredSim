

%% Scale achilles tendon properties
%-----------------------------------

clear all; clc;
import casadi.*;

CasadiFunc = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\CasADiFunctions\Casadi_s1Pog_mtp';
OrigPath = pwd;
cd(CasadiFunc);
f_FiberLength_TendonForce_tendon = Function.load('f_FiberLength_TendonForce_tendon');
f_FiberVelocity_TendonForce_tendon = Function.load('f_FiberVelocity_TendonForce_tendon');
f_forceEquilibrium_FtildeState_all_tendon = Function.load('f_forceEquilibrium_FtildeState_all_tendon');
f_lMT_vMT_dM = Function.load('f_lMT_vMT_dM');
cd(OrigPath);

% mt properties
muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
    'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
    'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
    'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
    'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
    'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
    'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
    'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
    'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
    'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
    'intobl_l','extobl_l'};
musi = MuscleIndices(muscleNames(1:end-3));
NMuscle = length(muscleNames(1:end-3))*2;
ExtPoly = '_mtp';
subject = 's1_Poggensee';
pathmusclemodel = fullfile('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim','MuscleModel',subject);
load([pathmusclemodel,'/MTparameters_',subject, ExtPoly, '.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];

load Fvparam
load Fpparam
load Faparam

% calf muscles
IndexCalf = [32 33 34];
params = MTparameters_m(:,IndexCalf);

% anatomical posture
q   = zeros(1,10);
qd  =  zeros(1,10);
[lMT, vMT, ~] = f_lMT_vMT_dM(q,qd);
lMT     = full(lMT(IndexCalf))';
vMT     = full(vMT(IndexCalf))';

% muscle tension
tension = getSpecificTensions(muscleNames(IndexCalf))';

% muscle state
fse =  0.2*ones(1,3);
dfse = [0 0 0];

% defaul tendon info
aTendon = [35 35 35];
shift = [0 0 0];

% optimization variables
opti =casadi.Opti();
a = opti.variable(1,3);
error = MX(1,3);
for m=1:3
    error(m) = ForceEquilibrium_FtildeState_all_tendon(a(m),fse(m),dfse(m),lMT(m),vMT(m),params(:,m),Fvparam,Fpparam,Faparam,tension(m),aTendon(m),shift(m));
end
opti.subject_to(error == 0);


% options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.max_iter              = 10000;
options.ipopt.linear_solver         = 'mumps';
options.ipopt.tol                   = 1*10^(-6);
opti.solver('ipopt', options);

sol = opti.solve();
aOpt = full(sol.value(a));
disp(aOpt);
[lM,lMtilde ] = FiberLength_TendonForce_tendon(fse,params,lMT,aTendon,shift);

%% Find scale values for lTs and lMo

opti2 =casadi.Opti();

scale   = opti2.variable(1,3);
aTendon = opti2.parameter(1,3);
shift   = getShift(aTendon);
FMo     = params(1,:);
lMo     = params(2,:).*scale;
lTs     = params(3,:).*scale;
alphao  = params(4,:);
vMmax   = params(5,:);
Ferr    = MX(1,3);
for m=1:3
    Ferr(m) = ForceEquilibrium_FtildeState_all_tendonOptParam(aOpt(m),fse(m),dfse(m),...
        lMT(m),vMT(m),FMo(m),lMo(m),lTs(m),alphao(m),vMmax(m),Fvparam,Fpparam,Faparam,tension(m),aTendon(m),shift(m));
end
% opti2.subject_to(Ferr == 0);
opti2.solver('ipopt', options);

% initial guess
opti2.set_initial(scale,1);
opti2.subject_to(0.7 < scale <1.3);
% opti2.minimize(sumsqr(Ferr));
opti2.subject_to(Ferr == 0);

% change tendon stiffness and find scale parameters
ct = 1;
kVect = 15:1:40;
for k = kVect
    opti2.set_value(aTendon,[k k k]);
    sol = opti2.solve();
    ScaleSol(ct,:) = value(sol,scale);
    ct = ct+1;
end

save('lMo_lTsScale.mat','ScaleSol','kVect');



% figure();
hold on;
plot(kVect,ScaleSol(:,1));




% 
% % fiber length should be equal for a given force in the anatomical position
% FT = 1000;
% dfse = [0 0 0];
% fse = [
%     lMT
% opti =casadi.Opti();
% a = opti.variable(3);
% 
% 
% ForceEquilibrium_FtildeState_all(a,fse,dfse,lMT,vMT,params,Fvparam,...
%     Fpparam,Faparam,tension)
% 
% [err, FT, Fce, Fpass, Fiso, vMmax, massM] = ...
%     ForceEquilibrium_FtildeState_all_tendon(a,fse,dfse,lMT,vMT,params,...
%     Fvparam,Fpparam,Faparam,tension,aTendon,shift)


% now find one value for slck length and optimal fiber length such that
% lMtilde is constant for different values of ATendon







