function [] = PlotResults_3DSim(ResultsFile,Cs,LegName,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Notes: We have to update the computations of the biological joint moments
% because this is currently based on subtracting the exo torque from the
% total ankle moment. Since we actuate more than 1 joint, we will have to
% do this computation based on R.Exodiff_id

% varargin:
%   (1) handle figure
%   (2) value x-axis (COT graph)
%   (3) label x-axis (COT graph)
%   (4) Structure with Reference data (ExperimentalData) you want to plot
%   on this figure

% Read variable input arguments
if length(varargin) > 3
    RefSNames = varargin{4}; % default subject from Falisse 2019
else
    RefSNames = {'subject1'};
end

% set colors for the experimentalData
CRef = parula(length(RefSNames)+3);
CRef = CRef(2:end,:);


lw = 2;

set(0,'defaultTextInterpreter','none');

if ~exist('LegName','var')
    [~,filename,~]= fileparts(ResultsFile);
    LegName = filename;
end
if exist(ResultsFile,'file')
    load(ResultsFile,'-mat','R');
    
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
            tab1 = h.Parent.Children(1).Children(1).Children(1);
            tab2 = h.Parent.Children(1).Children(1).Children(2);
            tab3 = h.Parent.Children(1).Children(1).Children(3);
            tab4 = h.Parent.Children(1).Children(1).Children(4);
            tab5 = h.Parent.Children(1).Children(1).Children(5);
            tab6 = h.Parent.Children(1).Children(1).Children(6);
            tab7 = h.Parent.Children(1).Children(1).Children(7);
            tab8 = h.Parent.Children(1).Children(1).Children(8);
            tab9 = h.Parent.Children(1).Children(1).Children(9);
            tab10 = h.Parent.Children(1).Children(1).Children(10);
            tab11 = h.Parent.Children(1).Children(1).Children(11);
            tab12 = h.Parent.Children(1).Children(1).Children(12);
            tab13 = h.Parent.Children(1).Children(1).Children(13);
            
            boolFirst = 0;
        else
            h = varargin{1};
            hTabGroup = uitabgroup;
            tab1 = uitab(hTabGroup, 'Title', 'Kinematics','BackgroundColor',[1 1 1]);
            tab2 = uitab(hTabGroup, 'Title', 'Kinetics','BackgroundColor',[1 1 1]);
            tab3 = uitab(hTabGroup, 'Title', 'COT','BackgroundColor',[1 1 1]);
            tab4 = uitab(hTabGroup, 'Title', 'ExoInfo','BackgroundColor',[1 1 1]);
            tab5 = uitab(hTabGroup, 'Title', 'CalfM','BackgroundColor',[1 1 1]);
            tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force','BackgroundColor',[1 1 1]);
            tab7 = uitab(hTabGroup, 'Title', 'Objective Function','BackgroundColor',[1 1 1]);
            tab8 = uitab(hTabGroup, 'Title', 'Ankle detailed','BackgroundColor',[1 1 1]);
            tab9 = uitab(hTabGroup, 'Title', 'SpatioTemporal','BackgroundColor',[1 1 1]);
            tab10 = uitab(hTabGroup, 'Title', 'Mech. Energy','BackgroundColor',[1 1 1]);
            tab11 = uitab(hTabGroup, 'Title', 'SolverInfo','BackgroundColor',[1 1 1]);
            tab12 = uitab(hTabGroup, 'Title', 'Body Kinematics','BackgroundColor',[1 1 1]);
            tab13 = uitab(hTabGroup, 'Title', 'DesignAssistance','BackgroundColor',[1 1 1]);
            h.Name = 'Sim3D_Results';
            set(h,'Color','w');
        end
    else
        h = figure('Color',[1 1 1]);
        h.Name = 'Sim3D_Results';
        hTabGroup = uitabgroup;
        tab1 = uitab(hTabGroup, 'Title', 'Kinematics','BackgroundColor',[1 1 1]);
        tab2 = uitab(hTabGroup, 'Title', 'Kinetics','BackgroundColor',[1 1 1]);
        tab3 = uitab(hTabGroup, 'Title', 'COT','BackgroundColor',[1 1 1]);
        tab4 = uitab(hTabGroup, 'Title', 'ExoInfo','BackgroundColor',[1 1 1]);
        tab5 = uitab(hTabGroup, 'Title', 'CalfM','BackgroundColor',[1 1 1]);
        tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force','BackgroundColor',[1 1 1]);
        tab7 = uitab(hTabGroup, 'Title', 'Objective Function','BackgroundColor',[1 1 1]);
        tab8 = uitab(hTabGroup, 'Title', 'Ankle detailed','BackgroundColor',[1 1 1]);
        tab9 = uitab(hTabGroup, 'Title', 'SpatioTemporal','BackgroundColor',[1 1 1]);
        tab10 = uitab(hTabGroup, 'Title', 'Mech. Energy','BackgroundColor',[1 1 1]);
        tab11 = uitab(hTabGroup, 'Title', 'SolverInfo','BackgroundColor',[1 1 1]);
        tab12 = uitab(hTabGroup, 'Title', 'Body Kinematics','BackgroundColor',[1 1 1]);
        tab13 = uitab(hTabGroup, 'Title', 'DesignAssistance','BackgroundColor',[1 1 1]);
        set(h,'Color','w');
        Pos = [331         221        1452         884];
        ScreenSize = get(0,'ScreenSize');
        if Pos(3)<ScreenSize(3) && Pos(4)<ScreenSize(4)
            set(gcf,'Position',Pos);
        else
            set(gcf,'Position',ScreenSize);
        end
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
    
    %% Plot default figure kinematics (Antoine his code)
    axes('parent', tab1);
    
    if boolFirst
        pathF = which('PlotResults_3DSim.m');
        pathrepo = pathF(1:end-26);
        pathReferenceData = [pathrepo,'/ExperimentalData'];
        load([pathReferenceData,'/ExperimentalData.mat'],'ExperimentalData');
        Qref = ExperimentalData.Q;
    end
    idx_Qs = [1,2,3,10,11,12,14,16,18,20,21,22,23,27,28,29,31];
    NumTicks = 2;
    ww = 10;
    
    label_fontsize  = 12;
    line_linewidth  = 2;
    for i = 1:length(idx_Qs)
        subplot(3,6,i)
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        % Experimental data
        if ~(strcmp(joints_ref{i},'mtp_angle')) && boolFirst == 1
           for is = 1:length(RefSNames)
              subject = RefSNames{is};
              if isfield(Qref,subject) && isfield(Qref.(subject),'Qs')
                  idx_jref = strcmp(Qref.(subject).Qs.colheaders,joints_ref{i});
                  meanPlusSTD = Qref.(subject).Qs.mean(:,idx_jref) + 2*Qref.(subject).Qs.std(:,idx_jref);
                  meanMinusSTD = Qref.(subject).Qs.mean(:,idx_jref) - 2*Qref.(subject).Qs.std(:,idx_jref);
                  stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                  intervalQ = 1:stepQ:size(R.Qs,1);
                  sampleQ = 1:size(R.Qs,1);
                  meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
                  meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
                  hold on
                  fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],CRef(is,:),'DisplayName',subject);
                  alpha(.25);
              end
           end            
        end
        
        % Simulation results
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        hold on;
        if i == length(idx_Qs)
            plot(x,R.Qs(:,idx_Qs(i)),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
        else
            plot(x,R.Qs(:,idx_Qs(i)),'color',Cs,'linewidth',line_linewidth);
        end
        
        % Plot settings
        if boolFirst
            set(gca,'Fontsize',label_fontsize);
            title(joints_tit{idx_Qs(i)},'Fontsize',label_fontsize);
            % Y-axis
            if i == 1 || i == 5 || i == 9 ||i == 13
                ylabel('Angle (°)','Fontsize',label_fontsize);
            end
            % X-axis
            %             L = get(gca,'XLim');
            if i > 12
                %                 set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            else
                set(gca,'XTick',[]);
            end
        end
    end
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos);
    end
    
    %% Plot joint kinetics
    axes('parent', tab2);
    label_fontsize  = 12;
    line_linewidth  = 2;
    for i = 1:length(idx_Qs)
        subplot(3,6,i)
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        % Experimental data
        % Experimental data
        if ~(strcmp(joints_ref{i},'mtp_angle')) && boolFirst == 1
            for is = 1:length(RefSNames)
                subject = RefSNames{is};
                IDref = ExperimentalData.Torques;
                if isfield(IDref,subject) && ~isempty(IDref.(subject).mean)                    
                    idx_jref = strcmp(IDref.(subject).colheaders,joints_ref{i});
                    meanPlusSTD = IDref.(subject).mean(:,idx_jref) + 2*IDref.(subject).std(:,idx_jref);
                    meanMinusSTD = IDref.(subject).mean(:,idx_jref) - 2*IDref.(subject).std(:,idx_jref);
                    stepID = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                    intervalID = 1:stepID:size(R.Qs,1);
                    sampleID = 1:size(R.Qs,1);
                    meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
                    meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID);
                    hold on
                    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],CRef(is,:),'DisplayName','MoCap');
                    alpha(.25);
                end
            end
        end
        
        % Simulation results
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        hold on;
        if i == length(idx_Qs)
            plot(x,R.Tid(:,idx_Qs(i)),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
        else
            plot(x,R.Tid(:,idx_Qs(i)),'color',Cs,'linewidth',line_linewidth);
        end
        
        % Plot settings
        if boolFirst
            % Plot settings
            set(gca,'Fontsize',label_fontsize);
            title(joints_tit{idx_Qs(i)},'Fontsize',label_fontsize);
            % Y-axis
            if i == 1 || i == 5 || i == 9 ||i == 13
                ylabel('Torque (Nm)','Fontsize',label_fontsize);
            end
            % X-axis
            L = get(gca,'XLim');
            NumTicks = 2;
            if i > 9
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            else
                set(gca,'XTick',[]);
            end
        end
    end
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos);
    end
    
    
    %% Plot general information
    axes('parent', tab3);
    
    % Plot COT as a function of exo assitance
    subplot(2,2,1); hold on;
    plot(xParam,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs,'DisplayName',LegName);
    ylabel('COT');  xlabel(xParamLab);
    title('Cost of Transport');
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    set(gca,'box','off');
    
    % Plot stride frequency
    subplot(2,2,2); hold on;
    plot(xParam,1./(R.tf_step*2),'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('Stride Frequency'); xlabel(xParamLab);
    title('Stride Frequency');
    set(gca,'box','off');
    if boolFirst
        if isfield(ExperimentalData,'StrideFreq')
            SF_mean = ExperimentalData.StrideFreq.Pilot1_AddedMass.mean;
            SF_std = ExperimentalData.StrideFreq.Pilot1_AddedMass.std;
            errorbar(xParam,SF_mean,SF_std,'o','Color',[0.5 0 0],'MarkerFaceColor',[0.5 0 0]);
        end
    end
    
    subplot(2,2,3);  hold on;
    plot(xParam,R.dt_exoShift,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('Time shift exo torque[s]');  xlabel(xParamLab);
    set(gca,'box','off');
    
    subplot(2,2,4); hold on;
    plot(x,R.T_exo(:,2),'-','Color',Cs);
    ylabel('Exo Moment [Nm]');  xlabel('% stride');
    title('Right');
    set(gca,'box','off');
    
    %% Plot Torque information
    % update this here
    boolActuation = 0;
    if isfield(R,'Exodiff_id')
        boolActuation = 1;
    end
    axes('parent', tab4);
    iAnkle = strcmp(R.colheaders.joints,'ankle_angle_r');
    iSubtalar = strcmp(R.colheaders.joints,'subtalar_angle_r');
    subplot(3,2,1)
    if boolActuation
        plot(x,R.Exodiff_id(:,iAnkle),'-','Color',Cs,'LineWidth',lw); hold on;
    else
        plot(x,R.T_exo(:,2),'-','Color',Cs); hold on;
    end
    ylabel('Exo Moment - Ankle [Nm]');  xlabel('% stride');
    set(gca,'box','off');
    
    subplot(3,2,2);
    if boolActuation
        plot(x,R.Exodiff_id(:,iSubtalar),'-','Color',Cs,'LineWidth',lw); hold on;
    else
    end
    ylabel('Exo Moment- Subtalar [Nm]'); xlabel('% stride');
    set(gca,'box','off');
    
    subplot(3,2,3)
    plot(x,R.Tid(:,iAnkle),'-','Color',Cs,'LineWidth',lw); hold on;
    ylabel('Ankle moment [Nm]'); xlabel('% stride');
    set(gca,'box','off');
    
    subplot(3,2,4)
    plot(x,R.Tid(:,iSubtalar),'-','Color',Cs,'LineWidth',lw); hold on;
    ylabel('Subtalar moment [Nm]'); xlabel('% stride');
    set(gca,'box','off');
    
    subplot(3,2,5)
    if boolActuation
        plot(x,R.Tid(:,iAnkle)-R.Exodiff_id(:,iAnkle),'-','Color',Cs,'LineWidth',lw); hold on;
    else
        plot(x,R.Tid(:,iAnkle)-R.T_exo(:,2),'-','Color',Cs,'LineWidth',lw); hold on;
    end
    ylabel('Biological ankle moment [Nm]'); xlabel('% stride');
    title('Left');
    set(gca,'box','off');
    
    subplot(3,2,6)
    if boolActuation
        plot(x,R.Tid(:,iSubtalar)-R.Exodiff_id(:,iSubtalar),'-','Color',Cs,'DisplayName',LegName,'LineWidth',lw); hold on;
    else
        plot(x,R.Tid(:,iSubtalar),'-','Color',Cs,'DisplayName',LegName,'LineWidth',lw); hold on;
    end
    ylabel('Biological subtalar moment [Nm]'); xlabel('% stride');
    title('Left');
    set(gca,'box','off');
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Plot ankle muscle energetics
    axes('parent', tab5);
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
    
    subplot(5,2,1)
    plot(x,R.a(:,iSol),'-','Color',Cs,'LineWidth',lw); hold on; title('Soleus');
    xlabel('% stride'); ylabel('activity');
    set(gca,'box','off');
    
    subplot(5,2,2)
    plot(x,R.a(:,iGas),'-','Color',Cs,'LineWidth',lw); hold on; title('Gastrocnemius');
    xlabel('% stride'); ylabel('activity');
    set(gca,'box','off');
    
    subplot(5,2,3)
    plot(x,R.MetabB.Etot(:,iSol),'-','Color',Cs,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Muscle metab power');
    set(gca,'box','off');
    
    subplot(5,2,4)
    plot(x,R.MetabB.Etot(:,iGas),'-','Color',Cs,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Muscle metab power (W)');
    set(gca,'box','off');
    
    subplot(5,2,5)
    plot(x,R.lMtilde(:,iSol),'-','Color',Cs,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Norm fiber length');
    set(gca,'box','off');
    
    subplot(5,2,6)
    plot(x,R.lMtilde(:,iGas),'-','Color',Cs,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Norm fiber length');
    set(gca,'box','off');
    
    subplot(5,2,7)
    plot(x,R.MetabB.Wdot(:,iSol),'-','Color',Cs,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Wdot');
    set(gca,'box','off');
    
    subplot(5,2,8)
    plot(x,R.MetabB.Wdot(:,iGas),'-','Color',Cs,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Wdot');
    set(gca,'box','off');
    
    subplot(5,2,9)
    plot(x,R.FT(:,iSol),'-','Color',Cs,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Norm muscle force');
    set(gca,'box','off');
    
    subplot(5,2,10)
    plot(x,R.FT(:,iGas),'-','Color',Cs,'DisplayName',LegName,'LineWidth',lw); hold on;
    xlabel('% stride'); ylabel('Norm muscle force');
    set(gca,'box','off');
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Ground reaction force
    axes('parent', tab6);
    if boolFirst == 1
        for iP = 1:3
            subplot(2,3,iP)
            for is = 1:length(RefSNames)
                subject = RefSNames{is};
                if isfield(ExperimentalData.GRFs,subject) && ~isempty(ExperimentalData.GRFs.(subject).mean)
                    GRFref = ExperimentalData.GRFs.(subject);
                    meanPlusSTD = GRFref.mean(:,iP) + 2*GRFref.std(:,iP);
                    meanMinusSTD = GRFref.mean(:,iP) - 2*GRFref.std(:,iP);
                    stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                    intervalQ = 1:stepQ:size(R.Qs,1);
                    sampleQ = 1:size(R.Qs,1);
                    meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
                    meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
                    hold on
                    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],CRef(is,:),'DisplayName','MoCap - P2');
                    alpha(.25);
                end
            end
        end
    end
    
    
    for i=1:6
        subplot(2,3,i);
        l = plot(x,R.GRFs(:,i),'-','Color',Cs,'LineWidth',lw); hold on;
        title(R.colheaders.GRF{i});
        xlabel('% stride');
        set(gca,'box','off');
        if i==1
            ylabel('Force [%BodyWeight]');
        end
    end
    l.DisplayName = LegName;
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
            l = plot(xParam,R.Obj.(Fields{i}),'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
            xlabel(xParamLab);
            title(Fields{i});
            set(gca,'box','off');
        end
        l.DisplayName = LegName;
        if boolFirst
            lh=legend('-DynamicLegend','location','east');
            lh.Interpreter = 'none';
        end
    end
    
    
    %% Plot ankle muscle energetics
    axes('parent', tab8);
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
    iGas2 = find(strcmp(R.colheaders.muscles,'med_gas_r'));
    iTib = find(strcmp(R.colheaders.muscles,'tib_ant_r'));
    if isempty(iGas)
        iGas = find(strcmp(R.colheaders.muscles,'gaslat_r'));
    end
    if isempty(iGas2)
        iGas2 = find(strcmp(R.colheaders.muscles,'gasmed_r'));
    end
    if isempty(iTib)
        iTib = find(strcmp(R.colheaders.muscles,'tibant_r'));
    end
    mVect = {'Soleus','Gas-lat','Gas-med','Tib-ant'};
    
    iM = [iSol iGas iGas2 iTib];
    
    for i=1:4
        subplot(5,4,i)
        plot(x,R.a(:,iM(i)),'-','Color',Cs,'LineWidth',lw); hold on; title(mVect{i});
        xlabel('% stride'); ylabel('activity');
        set(gca,'box','off');
        
        subplot(5,4,i+4)
        plot(x,R.MetabB.Etot(:,iM(i)),'-','Color',Cs,'LineWidth',lw); hold on;
        xlabel('% stride'); ylabel('Muscle metab power');
        set(gca,'box','off');
        
        subplot(5,4,i+8)
        plot(x,R.lMtilde(:,iM(i)),'-','Color',Cs,'LineWidth',lw); hold on;
        xlabel('% stride'); ylabel('Norm fiber length');
        set(gca,'box','off');
        
        subplot(5,4,i+12)
        plot(x,R.FT(:,iM(i)),'-','Color',Cs,'LineWidth',lw); hold on;
        xlabel('% stride'); ylabel('Norm muscle force');
        set(gca,'box','off');
        
    end
    % plot (biological) joint moments
    subplot(5,4,17)
    plot(x,R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r')),'-','Color',Cs,'LineWidth',lw); hold on;
    ylabel('Ankle moment [Nm]'); xlabel('% stride');
    set(gca,'box','off');
    
    subplot(5,4,18)
    plot(x,R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r'))-R.T_exo(:,2),'-','Color',Cs,'LineWidth',lw); hold on;
    ylabel('Muscle ankle [Nm]'); xlabel('% stride');
    set(gca,'box','off');
    
    subplot(5,4,19)
    plot(x,R.Tid(:,strcmp(R.colheaders.joints,'subtalar_angle_r')),'-','Color',Cs,'LineWidth',lw); hold on;
    ylabel('subtalar moment [Nm]'); xlabel('% stride');
    set(gca,'box','off');
    
    subplot(5,4,20)
    l = plot(x,R.Tid(:,strcmp(R.colheaders.joints,'mtp_angle_r')),'-','Color',Cs,'LineWidth',lw); hold on;
    ylabel('mtp moment [Nm]'); xlabel('% stride');
    l.DisplayName = LegName;
    set(gca,'box','off');
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Spatiotemporal results
    
    axes('parent', tab9);
    
    subplot(2,3,1);
    plot(xParam,R.StrideLength,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
    ylabel('Stride length [m]');  xlabel(xParamLab);
    set(gca,'box','off');
    
    if boolFirst
        for is = 1:length(RefSNames)
            subject = RefSNames{is};
            if isfield(ExperimentalData,'StrideLength') && isfield(ExperimentalData.StrideLength,subject)
                SL_mean = ExperimentalData.StrideLength.(subject).mean;
                SL_std =   ExperimentalData.StrideLength.(subject).std;
                errorbar(xParam,SL_mean,SL_std,'o','Color',CRef(is,:),'MarkerFaceColor',CRef(is,:));
            end
        end
    end
    
    if isfield(R,'StepWidth_COP')
        subplot(2,3,2);
        plot(xParam,R.StepWidth_COP,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        ylabel('Stride width [m]');  xlabel(xParamLab);
        set(gca,'box','off');
        if boolFirst
            for is = 1:length(RefSNames)
                subject = RefSNames{is};
                if isfield(ExperimentalData,'StepWidth') && isfield(ExperimentalData.StepWidth,subject)
                    SW_mean = ExperimentalData.StepWidth.(subject).mean;
                    SW_std = ExperimentalData.StepWidth.(subject).std;
                    errorbar(xParam,SW_mean,SW_std,'o','Color',CRef(is,:),'MarkerFaceColor',CRef(is,:));
                end
            end
        end
    end
    
    subplot(2,3,3);
    plot(xParam,1./(R.tf_step*2),'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
    ylabel('stride frequency [Hz]');  xlabel(xParamLab);
    set(gca,'box','off');
    
    if boolFirst
        for is = 1:length(RefSNames)
            subject = RefSNames{is};
            if isfield(ExperimentalData,'StrideFreq') && isfield(ExperimentalData.StrideFreq,subject)
                SF_mean = ExperimentalData.StrideFreq.(subject).mean;
                SF_std = ExperimentalData.StrideFreq.(subject).std;
                errorbar(xParam,SF_mean,SF_std,'o','Color',CRef(is,:),'MarkerFaceColor',CRef(is,:));
            end
        end
    end
    
    subplot(2,3,4);
    plot(xParam,R.Event.Stance,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
    ylabel('% stance');  xlabel(xParamLab);
    set(gca,'box','off');

    if boolFirst
        for is = 1:length(RefSNames)
            subject = RefSNames{is};
            if isfield(ExperimentalData,'PercStance') && isfield(ExperimentalData.PercStance,subject)
                Smean = ExperimentalData.PercStance.(subject).mean;
                Sstd = ExperimentalData.PercStance.(subject).std;
                errorbar(xParam,Smean,Sstd,'o','Color',CRef(is,:),'MarkerFaceColor',CRef(is,:));
            end
        end
    end
    
    subplot(2,3,5);
    plot(xParam,R.Event.Swing,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
    ylabel('% swing');  xlabel(xParamLab);
    set(gca,'box','off');
    
    if boolFirst
        for is = 1:length(RefSNames)
            subject = RefSNames{is};
            if isfield(ExperimentalData,'PercSwing') && isfield(ExperimentalData.PercSwing,subject)
                Smean = ExperimentalData.PercSwing.(subject).mean;
                Sstd = ExperimentalData.PercSwing.(subject).std;
                errorbar(xParam,Smean,Sstd,'o','Color',CRef(is,:),'MarkerFaceColor',CRef(is,:));
            end
        end
    end
    
    subplot(2,3,6);
    l = plot(xParam,R.Event.DS,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
    ylabel('% double support');  xlabel(xParamLab);
    set(gca,'box','off');
    l.DisplayName = LegName;
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    
    %% Plot mechanical energy results
    axes('parent', tab10);
    if isfield(R,'MechE')
        nc = 3;
        nr = 2;
        subplot(nr,nc,1);
        plot(xParam,R.MechE.MuscleWork,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Net. Fiber Work gc [J]');
        xlabel(xParamLab);
        
        subplot(nr,nc,2);
        plot(xParam,R.MechE.MusclePosWorkTotal,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Positive Fiber Work gc [J]');
        xlabel(xParamLab);
        
        subplot(nr,nc,3);
        plot(xParam,R.MechE.MuscleNegWorkTotal,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Negative Fiber Work gc [J]');
        xlabel(xParamLab);
        
        subplot(nr,nc,4);
        plot(xParam,R.MechE.TotalWork,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Total mechanical work [J]');
        xlabel(xParamLab);
        
        subplot(nr,nc,5:6);
        l = plot(x,R.MechE.TotalPower,'Color',Cs,'MarkerFaceColor',Cs,'LineWidth',lw); hold on;
        set(gca,'box','off');
        ylabel('Total Power [W]');
        xlabel('% gait cycle');
        l.DisplayName = LegName;
        if boolFirst
            lh=legend('-DynamicLegend','location','east');
            lh.Interpreter = 'none';
        end
    end
    
    %% Plot solver info
    axes('parent', tab11);
    if isfield(R,'stats')
        nc = 2;
        nr = 2;
        
        subplot(nr,nc,1);
        plot(xParam,R.stats.t_proc_total,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Total time [s]');
        xlabel(xParamLab);
        
        subplot(nr,nc,2);
        plot(xParam,R.stats.iter_count,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('N iterations []');
        xlabel(xParamLab);
        
        dt_NLP = R.stats.t_proc_nlp_f +R.stats.t_proc_nlp_g +R.stats.t_proc_nlp_grad_f +R.stats.t_proc_nlp_jac_g;
        subplot(nr,nc,3);
        plot(xParam,dt_NLP,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('time Func Eval [s]');
        xlabel(xParamLab);
        
        subplot(nr,nc,4);
        plot(xParam,R.stats.success,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Boolean succes [s]');
        xlabel(xParamLab);
        
    end
    
    %% Plot Body kinematics
    
    axes('parent', tab12);
    if isfield(R,'BodyKin')
        nc = 3;
        nr = 2;
        
        % plot the experimental data
        PlotHeaders = {'pelvis_Oz','femur_r_Oz','tibia_r_Oz','talus_r_Oz','calcn_r_Oz','toes_r_Oz'};
        if boolFirst == 1
            for iP = 1:length(PlotHeaders)                
                subplot(nr,nc,iP)
                for is = 1:length(RefSNames)
                subject = RefSNames{is};
                    if isfield(ExperimentalData,'BK') && isfield(ExperimentalData.BK,subject)
                        BK = ExperimentalData.BK;                        
                        idx_jref = strcmp(BK.(subject).Pos.colheaders,PlotHeaders{iP});
                        meanPlusSTD = BK.(subject).Pos.mean(:,idx_jref) + 2*BK.(subject).Pos.std(:,idx_jref);
                        meanMinusSTD = BK.(subject).Pos.mean(:,idx_jref) - 2*BK.(subject).Pos.std(:,idx_jref);
                        stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                        intervalQ = 1:stepQ:size(R.Qs,1);
                        sampleQ = 1:size(R.Qs,1);
                        meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
                        meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
                        hold on
                        fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],CRef(is,:),'DisplayName','MoCap - P2');
                        alpha(.25);
                        title(PlotHeaders{iP});
                    end
                end
            end            
            subplot(nr,nc,1)
            ylabel('orientation [deg]');
            subplot(nr,nc,4)
            ylabel('orientation [deg]');
            xlabel('%gait cycle');
            subplot(nr,nc,5)
            xlabel('%gait cycle');
            subplot(nr,nc,6)
            xlabel('%gait cycle');
        end
        
        % plot the simulation data
        if  ~isempty(R.BodyKin)
            for iP = 1:length(PlotHeaders)
                subplot(nr,nc,iP); hold on;
                iSel = strcmp(R.BodyKin.header,PlotHeaders{iP});
                nBK = length(R.BodyKin.Pos(:,iSel)); % rounding problems with bodykinematics
                if nBK >length(x)
                    nBK = length(x);
                end
                l = plot(x(1:nBK),R.BodyKin.Pos(:,iSel),'-','Color',Cs,'LineWidth',lw); hold on;
            end
        end
        
        l.DisplayName = LegName;
        if boolFirst
            lh=legend('-DynamicLegend','location','east');
            lh.Interpreter = 'none';
        end
    end
    if isfield(R,'MechE')
        plot(x,R.MechE.TotalPower,'Color',Cs,'MarkerFaceColor',Cs,'LineWidth',lw); hold on;
        set(gca,'box','off');
        ylabel('Total Power [W]');
        xlabel('% gait cycle');
    end
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Plot solver info
    axes('parent', tab11);
    if isfield(R,'stats')
        nc = 2;
        nr = 2;
        
        subplot(nr,nc,1);
        plot(xParam,R.stats.t_proc_total,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Total time [s]');
        xlabel(xParamLab);
        
        subplot(nr,nc,2);
        plot(xParam,R.stats.iter_count,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('N iterations []');
        xlabel(xParamLab);
        
        dt_NLP = R.stats.t_proc_nlp_f +R.stats.t_proc_nlp_g +R.stats.t_proc_nlp_grad_f +R.stats.t_proc_nlp_jac_g;
        subplot(nr,nc,3);
        plot(xParam,dt_NLP,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('time Func Eval [s]');
        xlabel(xParamLab);
        
        subplot(nr,nc,4);
        plot(xParam,R.stats.success,'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
        set(gca,'box','off');
        ylabel('Boolean succes [s]');
        xlabel(xParamLab);
        
    end
    
    %% Design assistance
    axes('parent', tab13);
    if isfield(R.S,'ExoDesignBool') && R.S.ExoDesignBool
        StrVect = {'ankle_angle_r','knee_angle_r','hip_flexion_r'};
        YLab = {'Ankle','Knee','Hip'};
        for ij = 1:length(StrVect)            
            iJoint = strcmp(R.colheaders.joints,StrVect{ij});
            T = R.Exodiff_id(:,iJoint);
            qd = R.Qdots(:,iJoint)*pi/180;
            x = linspace(1,100,length(T));
            ExoPower = T.*qd;    

            subplot(3,2,ij*2-1)
            plot(x,T,'Color',Cs,'LineWidth',lw); hold on;
            ylabel(YLab{ij});
            if ij ==1
                title('Moment');
            elseif ij == 3
                xlabel('% gait cycle');
            end
            set(gca,'box','off');

            subplot(3,2,ij*2);        
            ls = plot(x,ExoPower,'Color',Cs,'LineWidth',lw); hold on;            
            if ij ==1
                title('Power');
            elseif ij == 3
                xlabel('% gait cycle');
            end
            set(gca,'box','off');
        end
    end

else
    warning(['File not found: ' ResultsFile]);
end

end

