clear
close all
clc
import casadi.*

addpath('../VariousFunctions/')

%% 3D ellipsoid
a = 0.5;
b = 0.2;
c = 0.3;

fun =@(x1,x2) c*sqrt( 1 - (x1/a).^2 - (x2/b).^2 );


x1_SX = SX.sym('x1');
x2_SX = SX.sym('x2');
f_jac = Function('jac_fun',{x1_SX,x2_SX},{jacobian(fun(x1_SX,x2_SX),[x1_SX;x2_SX])'});

jac_fun =@(x1,x2) full(f_jac(x1(:)',x2(:)'))';


%% sample from function

% random samples
x12 = [a;b]/2.*(2*lhsdesign(100,2)'-1);
x1 = x12(1,:);
x2 = x12(2,:);

% grid samples
% [x1, x2] = meshgrid( linspace(-1,1,10)*a/2, linspace(-1,1,10)*b/2 );


x_s = [x1(:), x2(:)];
y_s = fun(x1(:),x2(:))';
ydx_s = jac_fun(x1(:),x2(:));

%% plot reference

[x_1, x_2] = meshgrid( linspace(-1,1,100)*a/2, linspace(-1,1,100)*b/2 );
y = reshape(fun(x_1(:),x_2(:)),100,100);
ydx = jac_fun(x_1(:),x_2(:));

figure
surf(x_1,x_2,y,"EdgeColor","none","FaceColor",[1,1,1]*0.8)
hold on
plot3(x_s(:,1),x_s(:,2),y_s,'.k')


%% fit 2nd order

[coeff1] = mvpolyfit(x_s, y_s, 2);

[y_fit1,jac_y_fit1] = mvpolyval(coeff1,[x_1(:),x_2(:)]);
y_fit1 = reshape(y_fit1,100,100);

rmse1 = rms(y-y_fit1,'all');
rmse_jac1 = rms(ydx-jac_y_fit1,'all');

surf(x_1,x_2,y_fit1,"EdgeColor","none","FaceColor",[0 0.4470 0.7410])

%% fit order in range 2-5, use Jacobian information

[coeff2, stats2] = mvpolyfit(x_s, y_s, [2,5], ydx_s);

[y_fit2,jac_y_fit2] = mvpolyval(coeff2,[x_1(:),x_2(:)]);
y_fit2 = reshape(y_fit2,100,100);

rmse2 = rms(y-y_fit2,'all');
rmse_jac2 = rms(ydx-jac_y_fit2,'all');

surf(x_1,x_2,y_fit2,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980])

%% centering and scaling of x

[coeff3, stats3, mu3] = mvpolyfit(x_s, y_s, [2,5], ydx_s);

[y_fit3,jac_y_fit3] = mvpolyval(coeff3,[x_1(:),x_2(:)],mu3);
y_fit3 = reshape(y_fit3,100,100);

rmse3 = rms(y-y_fit3,'all');
rmse_jac3 = rms(ydx-jac_y_fit3,'all');

surf(x_1,x_2,y_fit3,"EdgeColor","none","FaceColor",[0.9290 0.6940 0.1250])

%% custom mldivide

% mldivide_impl = 'mldivide';
mldivide_impl =@(A,b) A\b; 

[coeff4] = mvpolyfit(x_s, y_s, 5, ydx_s, "mldivide_impl",mldivide_impl);

[y_fit4,jac_y_fit4] = mvpolyval(coeff4,[x_1(:),x_2(:)]);
y_fit4 = reshape(y_fit4,100,100);

rmse4 = rms(y-y_fit4,'all');
rmse_jac4 = rms(ydx-jac_y_fit4,'all');

surf(x_1,x_2,y_fit4,"EdgeColor","none","FaceColor",[0.4940 0.1840 0.5560])







