function [] = PlotResults_3DSim(ResultsFile,Cs,LegName,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Notes: We have to update the computations of the biological joint moments
% because this is currently based on subtracting the exo torque from the
% total ankle moment. Since we actuate more than 1 joint, we will have to
% do this computation based on R.Exodiff_id



set(0,'defaultTextInterpreter','none');

if ~exist('LegName','var')
    [~,filename,~]= fileparts(ResultsFile);
    LegName = filename;
end
if exist(ResultsFile,'file')
    load(ResultsFile);
    
    if length(varargin)<2
        xParam = R.Sopt.ExoScale;
        xParamLab = 'Level assistance';
    else
        xParam = varargin{2};
        xParamLab = varargin{3};
    end
    
    boolFirst = 1;
    % create figure handles
    if ~isempty(varargin) && ~isempty(varargin{1})
        if ~isempty(varargin{1}.Name)
            h = varargin{1};
            tab3 = h.Parent.Children(1).Children(1).Children(1);
            tab6 = h.Parent.Children(1).Children(1).Children(2);
            tab7 = h.Parent.Children(1).Children(1).Children(3);
            boolFirst = 0;
        else
            h = varargin{1};
            hTabGroup = uitabgroup;
            
            tab3 = uitab(hTabGroup, 'Title', 'COT');
            tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force');
            tab7 = uitab(hTabGroup, 'Title', 'Objective Function');
            h.Name = 'Sim3D_Results';
            set(h,'Color','w');
        end
    else
        h = figure();
        h.Name = 'Sim3D_Results';
        hTabGroup = uitabgroup;
        tab3 = uitab(hTabGroup, 'Title', 'COT');
        tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force');
        tab7 = uitab(hTabGroup, 'Title', 'Objective Function');
        set(h,'Color','w');
    end
    
    %% Helper names
    joints_ref = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
        'hip_flexion','hip_adduction','hip_rotation',...
        'knee_angle','ankle_angle','subtalar_angle','mtp_angle',...
        'lumbar_extension','lumbar_bending','lumbar_rotation',...
        'arm_flex','arm_add','arm_rot','elbow_flex'};
    joints_tit = {'Pelvis tilt','Pelvis list','Pelvis rotation','Pelvis tx',...
        'Pelvis ty','Pelvis tz','Hip flexion L','Hip adduction L',...
        'Hip rotation L','Hip flexion R','Hip adduction R','Hip rotation R',...
        'Knee L','Knee R','Ankle L','Ankle R',...
        'Subtalar L','Subtalar R','MTP L','MTP R',...
        'Lumbar extension','Lumbar bending','Lumbar rotation',...
        'Arm flexion L','Arm adduction L','Arm rotation L',...
        'Arm flexion R','Arm adduction R','Arm rotation R',...
        'Elbow flexion L','Elbow flexion R'};
    
    
    
    %% Plot general information
    axes('parent', tab3);
    
    % Plot COT as a function of exo assitance
    subplot(2,2,1); hold on;
    plot(xParam,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('COT');  xlabel(xParamLab);
    title('Cost of Transport');
    
    % Plot stride frequency
    subplot(2,2,2); hold on;
    dt = R.t(end);
    plot(xParam,1./dt,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('Stride Frequency'); xlabel(xParamLab);
    title('Stride Frequency');
    
    subplot(2,2,3);  hold on;
    plot(xParam,R.dt_exoShift,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('Time shift exo torque[s]');  xlabel(xParamLab);
    
    
    subplot(2,2,4); hold on;
    plot(R.T_exo(:,2),'-','Color',Cs,'DisplayName',LegName);
    ylabel('Exo Moment [Nm]');  xlabel('% stride');
    title('Right');
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        %         lhPos = lh.Position;
        %         lhPos(1) = lhPos(1)+0.2;
        %         set(lh,'position',lhPos);
    end
    
    
    %% Ground reaction force
    axes('parent', tab6);
    for i=1:6
        subplot(2,3,i);
        if i== 6
            plot(R.GRFs(:,i),'-','Color',Cs,'DisplayName',LegName); hold on;
        else
            plot(R.GRFs(:,i),'-','Color',Cs); hold on;
        end
        title(R.colheaders.GRF{i});
        xlabel('% stride');
    end
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Objective function
    axes('parent', tab7);
    if isfield(R,'Obj')
        Fields = fieldnames(R.Obj);
        nf = length(Fields);
        for i= 1:nf
            subplot(3,4,i)
            if i== nf
                plot(xParam,R.Obj.(Fields{i}),'o','Color',Cs,'MarkerFaceColor',Cs,'DisplayName',LegName); hold on;
            else
                plot(xParam,R.Obj.(Fields{i}),'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
            end
            xlabel(xParamLab);
            title(Fields{i});
        end
        if boolFirst
            lh=legend('-DynamicLegend','location','east');
            lh.Interpreter = 'none';
        end
    end
    
else
    warning(['File not found: ' ResultsFile]);
end

end

