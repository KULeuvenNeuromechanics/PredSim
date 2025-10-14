clear
close all
clc
import casadi.*



a = 0.5;
b = 0.2;
c = 0.3;

fun =@(x1,x2) c*sqrt( 1 - (x1/a).^2 - (x2/b).^2 );

n_points = 100;

% [x1, x2] = meshgrid( linspace(-1,1,10)*a/2, linspace(-1,1,10)*b/2 );

x12 = [a;b]/2.*(2*lhsdesign(n_points,2)'-1);
x1 = x12(1,:);
x2 = x12(2,:);


x1_SX = SX.sym('x1');
x2_SX = SX.sym('x2');
jac_fun = Function('jac_fun',{x1_SX,x2_SX},...
    {jacobian(fun(x1_SX,x2_SX),x1_SX), jacobian(fun(x1_SX,x2_SX),x2_SX)});

x3 = fun(x1,x2);
[ydx1,ydx2] = jac_fun(x1(:),x2(:));



x = [x1(:), x2(:)];
y = x3(:);
ydx = [full(ydx1), full(ydx2)];

order = [3];

[coeff] = mvpolyfit(x, y, order, ydx);

[coeff2, stats2, mu2] = mvpolyfit(x, y, order);



%%

[x_1, x_2] = meshgrid( linspace(-1,1,100)*a/2, linspace(-1,1,100)*b/2 );
x_3 = reshape(fun(x_1(:),x_2(:)),100,100);
[ydx1,ydx2] = jac_fun(x_1(:),x_2(:));
ydx = [full(ydx1), full(ydx2)];

[x_3_fit1,jac_x_3_fit1] = mvpolyval(coeff,[x_1(:),x_2(:)]);
x_3_fit1 = reshape(x_3_fit1,100,100);

[x_3_fit2,jac_x_3_fit2] = mvpolyval(coeff2,[x_1(:),x_2(:)],mu2);
x_3_fit2 = reshape(x_3_fit2,100,100);

rmse1 = rms(x_3-x_3_fit1,'all');
rmse2 = rms(x_3-x_3_fit2,'all');
rmse_jac1 = rms(ydx-jac_x_3_fit1,'all');
rmse_jac2 = rms(ydx-jac_x_3_fit2,'all');

figure
surf(x_1,x_2,x_3,"EdgeColor","none","FaceColor",[0 0.4470 0.7410])
hold on
surf(x_1,x_2,x_3_fit1,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980])
surf(x_1,x_2,x_3_fit2,"EdgeColor","none","FaceColor",[0.9290 0.6940 0.1250])






