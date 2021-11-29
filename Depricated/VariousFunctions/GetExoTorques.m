function [ExoVect] = GetExoTorques(S,pathRepo,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

body_mass = S.mass;
if S.ExoBool
    if strcmp(S.DataSet,'Zhang2017')
        % load the data from Zhang 2017
        Zhang = load([pathRepo,'\Data\Zhang_2017\opt_tau.mat']);
        Tankle = nanmean(Zhang.opt_tau)*-1.*body_mass; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle Tankle]);
    elseif strcmp(S.DataSet,'PoggenSee2020_AFO')
        Poggensee = load([pathRepo,'\Data\Poggensee_2020\torque_profile.mat']);
        Tankle = Poggensee.torque*-1*body_mass; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle' Tankle']);
    elseif strcmp(S.DataSet,'PoggenSee2020_Exp')
        Poggensee = load([pathRepo,'\Data\Poggensee_2020_Exp\torque_profile.mat']);
        Tankle = Poggensee.torque*-1; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle' Tankle']);
     elseif strcmp(S.DataSet,'PoggenSee2020_ExpPass')
        Poggensee = load([pathRepo,'\Data\Poggensee_2020_ExpPass\torque_profile.mat']);
        Tankle = Poggensee.torque*-1; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle' Tankle']);
    elseif strcmp(S.DataSet,'Poggensee_2020_Subj2_T1')
        % second subject send by Katie pogensee
        DPog = load([pathRepo,'\Data\Poggensee_2020_Subj2\torque_profile.mat']);
        Tankle = DPog.torque1*-1; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle' Tankle']);
    elseif strcmp(S.DataSet,'Poggensee_2020_Subj2_T2')
        % second subject send by Katie pogensee
        DPog = load([pathRepo,'\Data\Poggensee_2020_Subj2\torque_profile.mat']);
        Tankle = DPog.torque2*-1; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle' Tankle']);
    else
        error(['Could not find the dataset ' S.DataSet ' to prescribe the exoskeleton torques']);
    end
    
    % adapt timing of exoskeleton torque to "align with natural" timing of gait
    % phases in the simulation
    if S.PercStance.bool
        % adapt ExoVect to express exoskeleton assistance as a percentage of
        % the stance phase instead of stride
        nfr = length(Tankle);
        Tankle_stance = Tankle(1:nfr*S.PercStance.xStanceOr);
        Tankle_swing = Tankle(nfr*S.PercStance.xStanceOr:nfr);
        SplineStance = spline(linspace(0,1,length(Tankle_stance)),Tankle_stance');
        SplineSwing = spline(linspace(0,1,length(Tankle_swing)),Tankle_swing');
        xStance = linspace(0,1,nfr*S.PercStance.xStanceNew);
        xSwing = linspace(0,1,nfr*(1-S.PercStance.xStanceNew));
        TStance = ppval(SplineStance,xStance);
        TSwing = ppval(SplineSwing,xSwing);
        Tankle = [TStance TSwing(2:end)];
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle Tankle]);
    end
    
    % express exoskeleton on discretisation simulation
    ExoControl.Tankle_r = ppval(ExoSpline.Tankle,linspace(0,0.5,N));
    ExoControl.Tankle_l = ppval(ExoSpline.Tankle,linspace(0.5,1,N));
    
    % scale exoskeleton support
    if isfield(S,'ExoScale')
        ExoControl.Tankle_r = ExoControl.Tankle_r*S.ExoScale;
        ExoControl.Tankle_l = ExoControl.Tankle_l*S.ExoScale;
    end
    ExoVect = [ExoControl.Tankle_l; ExoControl.Tankle_r];
else
    % zero exoskeleton assitance
    ExoVect = zeros(2,N);
end



end

