function [] = PlotResults_3DSim(ResultsFile,Cs,LegName,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('LegName','var')
    LegName = '';
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
        h = varargin{1};
        figure(h);
        tab0 = h.Parent.Children(1).Children(1).Children(1);
        tab1 = h.Parent.Children(1).Children(1).Children(2);
        tab2 = h.Parent.Children(1).Children(1).Children(3);
        tab3 = h.Parent.Children(1).Children(1).Children(4);
        tab4 = h.Parent.Children(1).Children(1).Children(5);
        tab5 = h.Parent.Children(1).Children(1).Children(6);
        tab6 = h.Parent.Children(1).Children(1).Children(7);
        tab7 = h.Parent.Children(1).Children(1).Children(8);
        tab8 = h.Parent.Children(1).Children(1).Children(9);
        tab9 = h.Parent.Children(1).Children(1).Children(10);
        tab10 = h.Parent.Children(1).Children(1).Children(11);
        boolFirst = 0;
    else
        h = figure();
%         set(h,'Position',get(0,'ScreenSize'));
        hTabGroup = uitabgroup;
        tab0 = uitab(hTabGroup, 'Title', 'Kinematics');
        tab1 = uitab(hTabGroup, 'Title', 'COT');
        tab2 = uitab(hTabGroup, 'Title', 'ExoInfo');
        tab3 = uitab(hTabGroup, 'Title', 'CalfM');
        tab4 = uitab(hTabGroup, 'Title', 'Kinematics - Sag');
        tab5 = uitab(hTabGroup, 'Title', 'Kinetics - Sag');
        tab6 = uitab(hTabGroup, 'Title', 'Kinematics - Front');
        tab7 = uitab(hTabGroup, 'Title', 'Kinetics - Front');
        tab8 = uitab(hTabGroup, 'Title', 'Ground reaction force');
        tab9 = uitab(hTabGroup, 'Title', 'Objective Function');
        tab10 = uitab(hTabGroup, 'Title', 'Ankle detailed');
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
    axes('parent', tab0); 
    
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
            subject = 'subject1';
            idx_jref = strcmp(Qref.(subject).Qs.colheaders,joints_ref{i});
            meanPlusSTD = Qref.(subject).Qs.mean(:,idx_jref) + 2*Qref.(subject).Qs.std(:,idx_jref);
            meanMinusSTD = Qref.(subject).Qs.mean(:,idx_jref) - 2*Qref.(subject).Qs.std(:,idx_jref);
            stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
            intervalQ = 1:stepQ:size(R.Qs,1);
            sampleQ = 1:size(R.Qs,1);
            meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
            meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
            hold on
            fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName','MoCap');
            alpha(.25);
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
        set(gca,'Fontsize',label_fontsize);
        title(joints_tit{idx_Qs(i)},'Fontsize',label_fontsize);
        % Y-axis
        if i == 1 || i == 5 || i == 9 ||i == 13
            ylabel('Angle (°)','Fontsize',label_fontsize);
        end
        % X-axis
        L = get(gca,'XLim');
        if i > 12
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        else
            set(gca,'XTick',[]);
        end
    end
%     subplot(3,6,18)
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos)
    else
%         leg = findobj(tab0,'Type','axes','Tag','legend');
%         legappend(LegName);
    end
    


    %% Plot general information
    axes('parent', tab1);
    
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
    
    %
%     subplot(2,2,3);  hold on;
%     plot(R.T_exo(:,1),'-','Color',Cs);
%     ylabel('Exo Moment [Nm]');  xlabel('% stride');
%     title('Left');
     subplot(2,2,3);  hold on;
     plot(xParam,R.dt_exoShift,'o','Color',Cs,'MarkerFaceColor',Cs);
     ylabel('Time shift exo torque[s]');  xlabel(xParamLab);
    % text(20,-10,['Time shift: ' num2str(round(R.dt_exoShift,3)) ' s']);
    
    
    subplot(2,2,4); hold on;
    plot(R.T_exo(:,2),'-','Color',Cs);
    ylabel('Exo Moment [Nm]');  xlabel('% stride');
    title('Right');
    
    
    %% Plot Torque information
    axes('parent', tab2);
    
    subplot(3,2,1)
    plot(R.T_exo(:,1),'-','Color',Cs); hold on;
    ylabel('Exo Moment [Nm]');  xlabel('% stride');
    title('Left');
    
    subplot(3,2,2);
    plot(R.T_exo(:,2),'-','Color',Cs); hold on;
    ylabel('Exo Moment [Nm]'); xlabel('% stride');
    title('Right');
    
    subplot(3,2,3)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_l')),'-','Color',Cs); hold on;
    ylabel('Ankle moment [Nm]'); xlabel('% stride');
    title('Left');
    
    subplot(3,2,4)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r')),'-','Color',Cs); hold on;
    ylabel('Ankle moment [Nm]'); xlabel('% stride');
    title('Right');
    
    subplot(3,2,5)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_l'))-R.T_exo(:,1),'-','Color',Cs); hold on;
    ylabel('Biological ankle moment [Nm]'); xlabel('% stride');
    title('Left');
    
    subplot(3,2,6)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r'))-R.T_exo(:,2),'-','Color',Cs); hold on;
    ylabel('Biological ankle moment [Nm]'); xlabel('% stride');
    title('Left');
    
    ax =[];
    ax2 = [];
    for i=1:3
        ax(i) = subplot(3,2,i*2-1);
        ax2(i) = subplot(3,2,i*2);
    end
    linkaxes(ax,'x');
    linkaxes(ax2,'x');
    
    %% Plot ankle muscle energetics
    axes('parent', tab3);
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
    
    subplot(5,2,1)
    plot(R.a(:,iSol),'-','Color',Cs); hold on; title('Soleus');
    xlabel('% stride'); ylabel('activity');
    
    subplot(5,2,2)
    plot(R.a(:,iGas),'-','Color',Cs); hold on; title('Gastrocnemius');
    xlabel('% stride'); ylabel('activity');
    
    subplot(5,2,3)
    plot(R.MetabB.Etot(:,iSol),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Muscle metab power');
    
    subplot(5,2,4)
    plot(R.MetabB.Etot(:,iGas),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Muscle metab power (W)');
    
    subplot(5,2,5)
    plot(R.lMtilde(:,iSol),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Norm fiber length');
    
    subplot(5,2,6)
    plot(R.lMtilde(:,iGas),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Norm fiber length');
    
    subplot(5,2,7)
    plot(R.MetabB.Wdot(:,iSol),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Wdot');
    
    subplot(5,2,8)
    plot(R.MetabB.Wdot(:,iGas),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Wdot');
    
    subplot(5,2,9)
    plot(R.FT(:,iSol),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Norm muscle force');
    
    subplot(5,2,10)
    plot(R.FT(:,iGas),'-','Color',Cs); hold on;
    xlabel('% stride'); ylabel('Norm muscle force');
    
    ax =[];
    ax2 = [];
    for i=1:5
        ax(i) = subplot(5,2,i*2-1);
        ax2(i) = subplot(5,2,i*2);
    end
    linkaxes(ax,'x');
    linkaxes(ax2,'x');
    
    %% Plot joint kinematics sagital plane
    
    axes('parent', tab4);
    
    % headers
    header =  R.colheaders.joints;
    
    % Plot pelvis kinematics
    nr = 4; nc= 3;
    for i=1:6
        subplot(nr,nc,i);
        plot(R.Qs(:,i),'-','Color',Cs); hold on;
        title(header{i});
    end
    
    % plot sagital plane kinematics
    subplot(nr,nc,7);
    plot(R.Qs(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,8);
    plot(R.Qs(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee flexion');
    subplot(nr,nc,9);
    plot(R.Qs(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle flexion');
    subplot(nr,nc,10);
    plot(R.Qs(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    subplot(nr,nc,11);
    plot(R.Qs(:,strcmp(header,'lumbar_extension')),'-','Color',Cs); hold on;
    title('lumbar extension');
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('Joint kinematics [deg or m]');
    end
    
    
    %% Plot joint kinetics sagital plane
    
    axes('parent', tab5);
    
    % headers
    header =  R.colheaders.joints;
    
    % Plot pelvis kinematics
    nr = 4; nc= 3;
    for i=1:6
        subplot(nr,nc,i);
        plot(R.Tid(:,i),'-','Color',Cs); hold on;
        title(header{i});
    end
    
    % plot sagital plane kinematics
    subplot(nr,nc,7);
    plot(R.Tid(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,8);
    plot(R.Tid(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee flexion');
    subplot(nr,nc,9);
    plot(R.Tid(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle flexion');
    subplot(nr,nc,10);
    plot(R.Tid(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    subplot(nr,nc,11);
    plot(R.Tid(:,strcmp(header,'lumbar_extension')),'-','Color',Cs); hold on;
    title('lumbar extension');
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('Joint kinetics [Nm or N]');
    end
    
    
    %% Plot other kinematics
    
    axes('parent', tab6);
    
    nr = 3; nc = 3;
    subplot(nr,nc,1);
    plot(R.Qs(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,2);
    plot(R.Qs(:,strcmp(header,'hip_adduction_r')),'-','Color',Cs); hold on;
    title('hip adduction');
    subplot(nr,nc,3);
    plot(R.Qs(:,strcmp(header,'hip_rotation_r')),'-','Color',Cs); hold on;
    title('hip rotation');
    
    
    subplot(nr,nc,5);
    plot(R.Qs(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee angle');
    
    subplot(nr,nc,7);
    plot(R.Qs(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle angle');
    subplot(nr,nc,8);
    plot(R.Qs(:,strcmp(header,'subtalar_angle_r')),'-','Color',Cs); hold on;
    title('subtalar angle');
    subplot(nr,nc,9);
    plot(R.Qs(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('angle [deg]');
    end
    
    
    %% Plot other kinetics
    
    axes('parent', tab7);
    
    nr = 3; nc = 3;
    subplot(nr,nc,1);
    plot(R.Tid(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,2);
    plot(R.Tid(:,strcmp(header,'hip_adduction_r')),'-','Color',Cs); hold on;
    title('hip adduction');
    subplot(nr,nc,3);
    plot(R.Tid(:,strcmp(header,'hip_rotation_r')),'-','Color',Cs); hold on;
    title('hip rotation');
    
    
    subplot(nr,nc,5);
    plot(R.Tid(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee angle');
    
    subplot(nr,nc,7);
    plot(R.Tid(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle angle');
    subplot(nr,nc,8);
    plot(R.Tid(:,strcmp(header,'subtalar_angle_r')),'-','Color',Cs); hold on;
    title('subtalar angle');
    subplot(nr,nc,9);
    plot(R.Tid(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('Joint moment [Nm]');
    end
    
    %% Ground reaction force  
     axes('parent', tab8);
    for i=1:6
        subplot(2,3,i);
        plot(R.GRFs(:,i),'-','Color',Cs); hold on;
        title(R.colheaders.GRF{i});
        xlabel('% stride');
    end
    
    %% Objective function
    axes('parent', tab9);
    if isfield(R,'Obj')
        Fields = fieldnames(R.Obj);
        nf = length(Fields);
        for i= 1:nf
            subplot(3,4,i)
            plot(xParam,R.Obj.(Fields{i}),'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
            xlabel(xParamLab);
            title(Fields{i});
        end        
    end    


    %% detailed ankle energetics
    %% Plot ankle muscle energetics
    axes('parent', tab10);
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
    iGas2 = find(strcmp(R.colheaders.muscles,'med_gas_r'));
    iTib = find(strcmp(R.colheaders.muscles,'tib_ant_r'));

    mVect = {'Soleus','Gas-lat','Gas-med','Tib-ant'};

    iM = [iSol iGas iGas2 iTib];

    for i=1:4
        subplot(5,4,i)
        plot(R.a(:,iM(i)),'-','Color',Cs); hold on; title(mVect{i});
        xlabel('% stride'); ylabel('activity');

        subplot(5,4,i+4)
        plot(R.MetabB.Etot(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Muscle metab power');

        subplot(5,4,i+8)
        plot(R.lMtilde(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Norm fiber length');

        subplot(5,4,i+12)
        plot(R.FT(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Norm muscle force');

    end
    % plot (biological) joint moments
    subplot(5,4,17)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r')),'-','Color',Cs); hold on;
    ylabel('Ankle moment [Nm]'); xlabel('% stride');
    
    subplot(5,4,18)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r'))-R.T_exo(:,2),'-','Color',Cs); hold on;
    ylabel('Muscle ankle [Nm]'); xlabel('% stride');

    subplot(5,4,19)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'subtalar_angle_r')),'-','Color',Cs); hold on;
    ylabel('subtalar moment [Nm]'); xlabel('% stride');

    subplot(5,4,20)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'mtp_angle_r')),'-','Color',Cs); hold on;
    ylabel('mtp moment [Nm]'); xlabel('% stride');
   
   
   
    
    ax =[];
    ax2 = [];
    ax3 =[];
    ax4 = [];

    for i=1:5
        ax(i) = subplot(5,4,i*4-3);
        ax2(i) = subplot(5,4,i*4-2);
        ax3(i) = subplot(5,4,i*4-1);
        ax4(i) = subplot(5,4,i*4);
    end
    linkaxes(ax,'x');
    linkaxes(ax2,'x');
    linkaxes(ax3,'x');
    linkaxes(ax4,'x');




else
    warning(['File not found: ' ResultsFile]);
end

end

