%% DummyModel Example
startpath = pwd;
cd('C:\OpenSimGit\opensim-core-build\RelWithDebInfo');
import casadi.*
F = external('F','DummyModel.dll',struct('enable_fd',true,...
        'enable_forward',false,'enable_reverse',false,...
        'enable_jacobian',false,'fd_method','backward')); 
    
% test .dll file
x = zeros(16,1);
u = zeros(8,1);


x(5) = 0;           % base orientation z ?
x(7) = 2;          % base x position
x(9) = 0.85;        % base position y
x(11) = 0;         % base z position
x(13) = 0; % q1
x(15) = 0;  % q2
    

y = full(F(x,u))';
disp(y);
cd(startpath);

Force = y(12:14);
Moment = y(15:17);
COPx = Moment(3)./Force(2);
COPz = -Moment(1)./Force(2);

disp(['COPx ' num2str(COPx) '  COPz ' num2str(COPz)]);
% clear F;

Force2 = y(18:20);
Moment2 =  y(21:23);
MomentCross = y(24:26);

% figure();
% subplot(1,2,1)
% bar([Force; Force2]);
% subplot(1,2,2)
% bar([Moment; Moment2]);
disp(['Force ' num2str(Force)]);
disp(['Force 2 ' num2str(Force2)]);
disp(['Moment ' num2str(Moment)]);
disp(['Moment 2 ' num2str(Moment2)]);
disp(['Moment Cross ' num2str(MomentCross)]);