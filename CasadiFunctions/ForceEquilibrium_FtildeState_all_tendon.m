function [err, FT, Fce, Fpass, Fiso, vMmax, massM] = ...
    ForceEquilibrium_FtildeState_all_tendon(a,fse,dfse,lMT,vMT,FMo_in,lMo_in,...
    lTs_in,alphao_in,vMmax_in,Fvparam,Fpparam,Faparam,tension,aTendon,shift,...
    MuscMoAsmp,d,stiffness,strength)

% This function derives the Hill-equilibrium.
% More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 

FMo = ones(size(a,1),1)*FMo_in;
lMo = ones(size(a,1),1)*lMo_in;
lTs = ones(size(a,1),1)*lTs_in;
alphao = ones(size(a,1),1)*alphao_in;
vMmax = ones(size(a,1),1)*vMmax_in;
Atendonsc = aTendon;
Atendon = ones(size(a,1),1)*Atendonsc;
volM = FMo.*lMo;
massM = volM.*(1059.7)./(tension*1e6);

% Inverse tendon force-length characteristic
lTtilde = log(5*(fse + 0.25 - shift))./Atendon + 0.995;

% Hill-type muscle model: geometric relationships
if(MuscMoAsmp == 0) % b = cst
    lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
else    % alpha = cst = alphao
   lM = (lMT-lTs.*lTtilde)./cos(alphao);
end
lMtilde = lM./lMo;

% Active muscle force-length characteristic
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);
b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);
num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);
num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);
FMltilde = FMtilde1+FMtilde2+FMtilde3;
Fiso = FMltilde;
% Active muscle force-velocity characteristic
vT = lTs.*dfse./(0.2*Atendonsc*exp(Atendonsc*(lTtilde-0.995)));
if(MuscMoAsmp == 0) % b = cst
    cos_alpha = (lMT-lTs.*lTtilde)./lM;
else    % alpha = cst = alphao
    cos_alpha = cos(alphao);
end
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;
e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);
FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;
% Active muscle force
Fcetilde = strength.*a.*FMltilde.*FMvtilde +d*vMtilde;
Fce = FMo.*Fcetilde;

% Passive muscle force-length characteristic
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilde - stiffness) / e0);
% Passive muscle force
Fpetilde = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);
Fpass = FMo.*Fpetilde;

% Muscle force (non-normalized)
% FM = FMo.*(Fcetilde+Fpetilde);

% Tendon force
FT = fse .* FMo;

% Equilibrium between muscle and tendon forces
% err =  FM.*cos_alpha-FT;
err =  (Fcetilde+Fpetilde).*cos_alpha-fse;

end
