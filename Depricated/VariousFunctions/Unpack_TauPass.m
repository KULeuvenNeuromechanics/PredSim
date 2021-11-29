function [Tau_passj,Tau_passj_J] = Unpack_TauPass(Tau_passj_all)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

Tau_passj.hip.flex.l = Tau_passj_all(1);
Tau_passj.hip.flex.r = Tau_passj_all(2);
Tau_passj.hip.add.l = Tau_passj_all(3);
Tau_passj.hip.add.r = Tau_passj_all(4);
Tau_passj.hip.rot.l = Tau_passj_all(5);
Tau_passj.hip.rot.r = Tau_passj_all(6);
Tau_passj.knee.l = Tau_passj_all(7);
Tau_passj.knee.r = Tau_passj_all(8);
Tau_passj.ankle.l = Tau_passj_all(9);
Tau_passj.ankle.r = Tau_passj_all(10);
Tau_passj.subt.l = Tau_passj_all(11);
Tau_passj.subt.r = Tau_passj_all(12);
Tau_passj.mtp.all = Tau_passj_all(13:14);
Tau_passj.trunk.ext = Tau_passj_all(15);
Tau_passj.trunk.ben = Tau_passj_all(16);
Tau_passj.trunk.rot = Tau_passj_all(17);
Tau_passj.arm = Tau_passj_all(18:25);
Tau_passj_J = Tau_passj_all([1:12 15:end]);
    
    
end

