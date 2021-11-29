function [] = FitPolynomials(MainPath,ModelName,Modelpath,PolyFolder,Bool_RunMA,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


max_order = 9;
if ~isempty(varargin)
    max_order = varargin{1};
end

threshold = 0.003; % 3mm
if length(varargin)>1
    threshold = varargin{2};
end


%% Settings

% Import the opensim API
import org.opensim.modeling.*

% Output folder for this subject
OutFolder = PolyFolder;
SubjFolder = fullfile(MainPath,'Polynomials',OutFolder);
MA_path=fullfile(SubjFolder,'MuscleAnalysis');
if ~isfolder(MA_path)
    mkdir(MA_path);
end
%% Create dummy motion and run muscle analysis
if Bool_RunMA
    
    % bounds on dofs to create training dataset
    Bound_hipflex = [-50 50];
    Bound_hipadd = [-30 30];
    Bound_hiprot = [-30 30];    
    Bound_ankle = [-30 30];
    Bound_subt = [-30 30];
    Bound_mtp = [-30 10];
    if strcmp(ModelName,'Gait92')
        Bound_knee = [-90 0];
    elseif strcmp(ModelName,'Rajagopal')
        Bound_knee = [0 90];
    else
        Bound_knee = [-90 0]; % default settings
    end
    
    % get the model coordinates
    m = Model(Modelpath);
    CoordSet = m.getCoordinateSet();
    nc = CoordSet.getSize();
    NamesCoordinates = cell(1,nc);
    for i = 1:nc
        NamesCoordinates{i} = char(CoordSet.get(i-1).getName());
    end
    % construct a file with generalized coordinates
    headers=[{'time'} NamesCoordinates];
    
    if strcmp(ModelName,'Rajagopal')
        % create the dummy motion
        n=5000; p=7;
        X = lhsdesign(n,p);
        X_scale=[diff(Bound_hipflex) diff(Bound_hipadd ) diff(Bound_hiprot),...
            diff(Bound_knee) diff(Bound_ankle) diff(Bound_subt) diff(Bound_mtp)];
        X_min=[Bound_hipflex(1) Bound_hipadd(1) Bound_hiprot(1) Bound_knee(1),...
            Bound_ankle(1) Bound_subt(1) Bound_mtp(1)];
        Angles=X.*(ones(n,1)*X_scale)+(ones(n,1)*X_min);
        IndexAngles = [7 8 9 10 12 13 14]+1; % +1 because of time vector
        
        % path with dummy motion
        time=(1:n)./100;
        data=zeros(length(time),length(headers));
        data(:,1)=time;
        data(:,IndexAngles) = Angles;   % right leg
        data(:,IndexAngles+8) = Angles; % left leg
        data(:,12) = Angles(:,4)*pi./180;   % Needed for the patella in the model
        data(:,20) = Angles(:,4)*pi./180;   % left leg
        pathDummyMotion = fullfile(SubjFolder,'dummy_motion.mot');
        generateMotFile(data,headers,pathDummyMotion);
        
        % select coordinate names to run muscle analysis
        CoordNames = headers(IndexAngles);
        
    elseif strcmp(ModelName,'Gait92')
        % additional bounds for the lumbar joint
        Bound_LumbExt = [-30 30];
        Bound_LumbBend = [-30 30];
        Bound_LumbRot = [-30 30];
        
        % create the dummy motion
        n=5000; p=10;
        X = lhsdesign(n,p);
        X_scale=[diff(Bound_hipflex) diff(Bound_hipadd ) diff(Bound_hiprot),...
            diff(Bound_knee) diff(Bound_ankle) diff(Bound_subt) diff(Bound_mtp),...
            diff(Bound_LumbExt) diff(Bound_LumbBend) diff(Bound_LumbRot)];
        X_min=[Bound_hipflex(1) Bound_hipadd(1) Bound_hiprot(1) Bound_knee(1),...
            Bound_ankle(1) Bound_subt(1) Bound_mtp(1), Bound_LumbExt(1), ...
            Bound_LumbBend(1), Bound_LumbRot(1)];
        Angles=X.*(ones(n,1)*X_scale)+(ones(n,1)*X_min);
        DofNames = {'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
            'knee_angle_r','ankle_angle_r','subtalar_angle_r','mtp_angle_r',...
            'lumbar_extension','lumbar_bending','lumbar_rotation'};
        IndexAngles = nan(1,length(DofNames));
        for iq = 1:length(DofNames)
            IndexAngles(iq) = find(strcmp(headers,DofNames{iq}));
        end                
        % path with dummy motion
        time=(1:n)./100;
        data=zeros(length(time),length(headers));
        data(:,1)=time;
        data(:,IndexAngles) = Angles;   % right leg
        data(:,IndexAngles+8) = Angles; % left leg
        pathDummyMotion = fullfile(SubjFolder,'dummy_motion.mot');
        generateMotFile(data,headers,pathDummyMotion);
        
        % select coordinate names to run muscle analysis
        CoordNames = headers(IndexAngles);
    end    
    %Run a muscle analysis on the dummy motion
    disp('Muscle analysis running....');
    OpenSim_Muscle_Analysis(pathDummyMotion,Modelpath,MA_path,[time(1) time(end)],CoordNames);    
end

%% import the muscle analysis data
% subject pre-fix
SubjPre = 'dummy_motion';
lMT = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_Length.sto']);
MA.hip.flex = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_hip_flexion_r.sto']);
MA.hip.add = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_hip_adduction_r.sto']);
MA.hip.rot = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_hip_rotation_r.sto']);
MA.knee.flex = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_knee_angle_r.sto']);
MA.ankle.flex = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_ankle_angle_r.sto']);
MA.sub = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_subtalar_angle_r.sto']);
MA.mtp = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_mtp_angle_r.sto']);
if strcmp(ModelName,'Gait92')
    % lumbar joint
    MA.trunk.ext = ReadMotFile([MA_path,'/',SubjPre,'_MuscleAnalysis_MomentArm_lumbar_extension.sto']);
    MA.trunk.ben = ReadMotFile([MA_path,'/',SubjPre,'_MuscleAnalysis_MomentArm_lumbar_bending.sto']);
    MA.trunk.rot = ReadMotFile([MA_path,'/',SubjPre,'_MuscleAnalysis_MomentArm_lumbar_rotation.sto']);

%% Fit polynomial functions
if strcmp(ModelName,'Rajagopal')
    %% Load the dummy motion
    pathDummyMotion = fullfile(SubjFolder,'dummy_motion.mot');
    dummy_motion = ReadMotFile(pathDummyMotion);

    % 15 dofs (mtp locked)
    % Order of dofs: hip flex r, hip add r, hip rot r, knee flex r, ankle flex
    % r, hip flex l, hip add l, hip rot l, knee flex l, ankle flex l, lumbar
    % ext, lumbar bend, lumbar rot, subtalar r, subtalar l, mtp_r, mtp_l
    order_Qs = [7 8 9 10 12 13 14]+1;
    q = dummy_motion.data(:,order_Qs).*(pi/180);
    % adapt the angle the knee such that it's similar to the definition in
    % gait92 model
    q(:,4) = -q(:,4);
    % changes sign moment arms knee joint
    MA.knee.flex.data(:,2:end) = -MA.knee.flex.data(:,2:end);
    
    %% Organize MuscleData   
    MuscleData.dof_names = dummy_motion.names(order_Qs);

    %% Organize MuscleData
    if ~isfield(dummy_motion,'colheaders')
        dummy_motion.colheaders = strsplit(dummy_motion.textdata{end});
    end
    MuscleData.dof_names = dummy_motion.colheaders(order_Qs);
    muscleNames = {'addbrev_r','addlong_r','addmagDist_r','addmagIsch_r','addmagMid_r','addmagProx_r',...
        'bflh_r','bfsh_r','edl_r','ehl_r','fdl_r','fhl_r','gaslat_r','gasmed_r','glmax1_r','glmax2_r',...
        'glmax3_r','glmed1_r','glmed2_r','glmed3_r','glmin1_r','glmin2_r','glmin3_r','grac_r','iliacus_r',...
        'perbrev_r','perlong_r','piri_r','psoas_r','recfem_r','sart_r','semimem_r','semiten_r','soleus_r',...
        'tfl_r','tibant_r','tibpost_r','vasint_r','vaslat_r','vasmed_r'};
    MuscleData.muscle_names = muscleNames;
    for m = 1:length(muscleNames)
        MuscleData.lMT(:,m)     = lMT.data(:,strcmp(lMT.names,muscleNames{m}));            % lMT
        MuscleData.dM(:,m,1)    = MA.hip.flex.data(:,strcmp(lMT.names,muscleNames{m}));    % hip_flex
        MuscleData.dM(:,m,2)    = MA.hip.add.data(:,strcmp(lMT.names,muscleNames{m}));     % hip_add
        MuscleData.dM(:,m,3)    = MA.hip.rot.data(:,strcmp(lMT.names,muscleNames{m}));     % hip_rot
        MuscleData.dM(:,m,4)    = MA.knee.flex.data(:,strcmp(lMT.names,muscleNames{m}));   % knee
        MuscleData.dM(:,m,5)    = MA.ankle.flex.data(:,strcmp(lMT.names,muscleNames{m}));  % ankle
        MuscleData.dM(:,m,6)    = MA.sub.data(:,strcmp(lMT.names,muscleNames{m}));         % sub
        MuscleData.dM(:,m,7)    = MA.mtp.data(:,strcmp(lMT.names,muscleNames{m}));         % mtp
    end
    MuscleData.q = q;
    
    %% Call PolynomialFit
    [muscle_spanning_joint_INFO,MuscleInfo] = PolynomialFit_mtp(MuscleData);
    save(fullfile(SubjFolder,'MuscleData.mat'),'MuscleData')
    save(fullfile(SubjFolder,'muscle_spanning_joint_INFO.mat'),'muscle_spanning_joint_INFO')
    save(fullfile(SubjFolder,'MuscleInfo.mat'),'MuscleInfo');
    
elseif strcmp(ModelName,'Gait92')
    pathDummyMotion = fullfile(SubjFolder,'dummy_motion.mot');
    dummy_motion = ReadMotFile(pathDummyMotion);
    % Order of dofs: hip flex r, hip add r, hip rot r, knee flex r, ankle flex
    % r, hip flex l, hip add l, hip rot l, knee flex l, ankle flex l, lumbar
    % ext, lumbar bend, lumbar rot, subtalar r, subtalar l, mtp_r, mtp_l
    DofNames = {'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
        'knee_angle_r','ankle_angle_r','subtalar_angle_r','mtp_angle_r',...
        'lumbar_extension','lumbar_bending','lumbar_rotation'};
    IndexAngles = nan(1,length(DofNames));
    for iq = 1:length(DofNames)
        IndexAngles(iq) = find(strcmp(dummy_motion.names,DofNames{iq}));
    end    
    q = dummy_motion.data(:,IndexAngles).*(pi/180);    
    MuscleData.dof_names = dummy_motion.names(IndexAngles);

    muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
        'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
        'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
        'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
        'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
        'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
        'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
        'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
        'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
        'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
        'intobl_l','extobl_l'};
    MuscleData.muscle_names = muscleNames;
    for m = 1:length(muscleNames)
        MuscleData.lMT(:,m)     = lMT.data(:,strcmp(lMT.names,muscleNames{m}));            % lMT
        MuscleData.dM(:,m,1)    = MA.hip.flex.data(:,strcmp(lMT.names,muscleNames{m}));    % hip_flex
        MuscleData.dM(:,m,2)    = MA.hip.add.data(:,strcmp(lMT.names,muscleNames{m}));     % hip_add
        MuscleData.dM(:,m,3)    = MA.hip.rot.data(:,strcmp(lMT.names,muscleNames{m}));     % hip_rot
        MuscleData.dM(:,m,4)    = MA.knee.flex.data(:,strcmp(lMT.names,muscleNames{m}));   % knee
        MuscleData.dM(:,m,5)    = MA.ankle.flex.data(:,strcmp(lMT.names,muscleNames{m}));  % ankle
        MuscleData.dM(:,m,6)    = MA.sub.data(:,strcmp(lMT.names,muscleNames{m}));         % sub
        MuscleData.dM(:,m,7)    = MA.mtp.data(:,strcmp(lMT.names,muscleNames{m}));         % mtp
        MuscleData.dM(:,m,8)    = MA.trunk.ext.data(:,strcmp(lMT.names,muscleNames{m}));   % lumbar ext
        MuscleData.dM(:,m,9)    = MA.trunk.ben.data(:,strcmp(lMT.names,muscleNames{m}));   % lumbar bend
        MuscleData.dM(:,m,10)    = MA.trunk.rot.data(:,strcmp(lMT.names,muscleNames{m}));  % lumbar rot
    end
    MuscleData.q = q;
    MuscleData.qdot = zeros(size(q));
    
    %% Call PolynomialFit
    [muscle_spanning_joint_INFO,MuscleInfo] = ...
        PolynomialFit_mtp(MuscleData,max_order,threshold);
    save(fullfile(SubjFolder,'MuscleData.mat'),'MuscleData')
    save(fullfile(SubjFolder,'muscle_spanning_joint_INFO.mat'),'muscle_spanning_joint_INFO')
    save(fullfile(SubjFolder,'MuscleInfo.mat'),'MuscleInfo');
end






end

