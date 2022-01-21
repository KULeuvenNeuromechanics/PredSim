% test

% f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
%     Qdotsj,a_ak,a_aj,a_mtpk,a_mtpj,vAk,e_ak,e_mtpk,dFTtildej,Aj},...
%     {eq_constrV,ineq_constr,J});

% tfk_r       = rand(1);
tfk_r       = 1;
ak_r        = rand(NMuscle,1);
aj_r        = rand(NMuscle,d);
FTtildek_r  = rand(NMuscle,1);
FTtildej_r  = rand(NMuscle,d);
Qsk_r       = rand(nq.all,1);
Qsj_r       = rand(nq.all,d);
Qdotsk_r    = rand(nq.all,1);
Qdotsj_r    = rand(nq.all,d);
a_ak_r      = rand(nq.arms,1);
a_aj_r      = rand(nq.arms,d);
a_mtpk_r    = rand(nq.mtp,1);
a_mtpj_r    = rand(nq.mtp,d);
vAk_r       = rand(NMuscle,1);
e_ak_r      = rand(nq.arms,1);
e_mtpk_r    = rand(nq.mtp,1);
dFTtildej_r = rand(NMuscle,d);
Aj_r        = rand(nq.all,d);

[eq_constrV,ineq_constr,J] = f_coll(tfk_r,ak_r,aj_r,FTtildek_r,FTtildej_r,Qsk_r,Qsj_r,Qdotsk_r,...
    Qdotsj_r,a_ak_r,a_aj_r,a_mtpk_r,a_mtpj_r,vAk_r,e_ak_r,e_mtpk_r,dFTtildej_r,Aj_r)
 