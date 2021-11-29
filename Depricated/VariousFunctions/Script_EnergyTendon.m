

% get equations using symbolic toolbox
% note lTr is the deviation of the tendon length w.r.t. the slack length 
syms FMo lTs Atendon shift lTr 'real'
lT = lTr+lTs;           % tendon length
lTtilde = lT./lTs;      % norm tendon length
fse = (exp(Atendon.*(lTtilde - 0.995)))/5-0.25+shift;   % tendon F/L properties
FT = fse*FMo;           % tendon force
energy = int(FT,lTr);   % tendon energie
matlabFunction(energy,'file',fullfile(pwd,'energy_tendon.m'));
clear FMo lTs Atendon shift lTr


% test for a simple case
lT = 0.2:0.001:0.22;
lTs = 0.201;
lTr = lT-lTs;
Atendon = 35;
shift = 0;
FMo = 1000;
energy = energy_tendon(Atendon,FMo,lTr,lTs,shift);
lTtilde = lT./lTs;
fse = (exp(Atendon.*(lTtilde - 0.995)))/5-0.25+shift;
FT = fse*FMo;
figure(); 
subplot(1,2,1);
plot(lTr,FT);
subplot(1,2,2)
plot(lTr,energy);
