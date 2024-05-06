clear
clc


import casadi.*

addpath('./Debugging/')

%% reference function
% example: normalised sum of 8 squared values
etemp8 = SX.sym('etemp8',8);
Jtemp8 = 0;
for i=1:length(etemp8)
    Jtemp8 = Jtemp8 + etemp8(i).^2;
end
Jtemp8 = Jtemp8/8;
f_J8_ref = Function('f_J8',{etemp8},{Jtemp8});

clear etemp8 Jtemp8

%% functions that needs testing
% same as reference (but MX function instead of SX)
etemp8 = MX.sym('etemp8',8);
Jtemp8 = 0;
for i=1:length(etemp8)
    Jtemp8 = Jtemp8 + etemp8(i).^2;
end
Jtemp8 = Jtemp8/8;
f_J8_test = Function('f_J8',{etemp8},{Jtemp8});

clear etemp8 Jtemp8

% small difference on the output
etemp8 = SX.sym('etemp8',8);
Jtemp8 = 0;
for i=1:length(etemp8)
    Jtemp8 = Jtemp8 + etemp8(i).^2;
end
Jtemp8 = Jtemp8/8 + 1e-10;
f_J8_wrong = Function('f_J8',{etemp8},{Jtemp8});

clear etemp8 Jtemp8

%% matlab function
x = sym('x',[8,1]);
f_J8_matlab = @(x) sum(x.^2)/8;

clear x
%% compare test functions to reference
% minimal input: both function handles
compareCasADiFunctions(f_J8_ref,f_J8_test)

% repeat the test 10 times to avoid rng-based false positive
diff = compareCasADiFunctions(f_J8_ref,f_J8_wrong,[],10);

% use a custom tolerance on the element-wise relative error
compareCasADiFunctions(f_J8_ref,f_J8_wrong,1e-7)

% test a matlab function vs a CasADi-function
compareCasADiFunctions(f_J8_ref,f_J8_matlab)



%% external functions

% path_ext = 'C:\Users\u0150099\Documents\master_thesis\3dpredictsim\ExternalFunctions';
% name1 = 'old\F_CP3_T0_scaled_MRI_v7_scaledMT_left.dll';
% name2 = 'old\F_CP3_T0_scaled_sf_v2_right.dll';
% 
% f_ext1 = external('F',fullfile(path_ext,name1));
% f_ext2 = external('F',fullfile(path_ext,name2));
% 
% compareCasADiFunctions(f_ext1,f_ext1)
% 
% compareCasADiFunctions(f_ext1,f_ext2)
% 
% compareCasADiFunctions(f_ext1,f_ext2,[],10)
% 
% compareCasADiFunctions(f_ext1,f_ext2,[],1,[0,1])


%% lookup table

% pathpolynomial = 'C:\Users\u0150099\Documents\master_thesis\3dpredictsim\Polynomials\Fal_s1_mtj_sc';
% f_table1 = Function.load((fullfile(pathpolynomial,'f_getMtjLigamentMoment_exp')));
% f_table2 = Function.load((fullfile(pathpolynomial,'f_getMtjLigamentMoment_exp_v2')));
% 
% compareCasADiFunctions(f_table1,f_table1)
% 
% compareCasADiFunctions(f_table1,f_table2)
% 
% compareCasADiFunctions(f_table1,f_table2,[],1,0)










