% --------------------------------------------------------------------------
% getHelperMatrices
% Code to create helper matricies used in making CasADi function that
% determine the muscle tendon length, velocity and moment arm using
% polynomial fit.
% 
% Original author: Dhruv Gupta
% Original date: 30/November/2021
%
% Last update:
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

clear
close all
clc

nOrder = 10;
nDOF = 6;
nCoeffMat = nan(nOrder,nDOF);
expoVal_all = cell(nOrder,nDOF);
for order=1:nOrder
    for n_dof=1:nDOF
        [nCoeffMat(order,n_dof),expoVal_all{order,n_dof}]=calcCoeff(order,n_dof);
    end
end
save('nCoeffMat.mat','nCoeffMat')
save('expoVal_all.mat','expoVal_all')
        
        