%% Convert MoCap data for scaling musculoskeletal model
%------------------------------------------------------

% settings
datapath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\Data\Poggensee_2020_Subj2';
c3dFile = fullfile(datapath,'static_exo_full.c3d');
outPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_3D\OpenSimModel\s2_Poggensee';
OutTRC = fullfile(outPath,'static_Exo.trc');

% read the c3d file
[Markers,MLabels,VideoFrameRate,AnalogSignals,ALabels,AUnits,AnalogFrameRate,Event,ParameterGroup,CameraInfo]=...
    readC3D(c3dFile);


% Convert markerNames
c3dNames    = {'STER','LWRI','RWRI','LTHI','LTH2','LTH3','RTHI','RTH2','RTH3',...
    'RTIB','RSH2','RSH3','LTIB','LSH2','LSH3',...
    'RMML','LMML','RMEP','LMEP','LKNM','RKNM',...
    'LMED','RMED'};    
OsimNames   = {'CLAV','LW1','RW1','LTHI3','LTH1','LTH2','RTHI3','RTH1','RTH2','RTIB3',...
    'RTIB1','RTIB2','LTIB3','LTIB1','LTIB2',...
    'RANK_med','LANK_med','RKNE_med','LKNE_med','LKNE_med','RKNE_med',...
    'LANK_med','RANK_med'};

disp('');
for i=1:length(c3dNames)
   Ind = find(strcmp(MLabels,c3dNames{i}));
   MLabels(Ind) = OsimNames(i);  
   if isempty(Ind)
       disp(c3dNames{i})
   end
end
disp('');

% export to a trc file
R = rotx(90)*roty(90); R=R';
Markers = rot3DVectors(R, Markers); % rotate markers
[nvF, nm] = size(Markers);           % get number of frames
vFrms = (1:nvF)';                   % frame index
vTime = (1/VideoFrameRate*(vFrms)); % Time vector
% set zeros to NaNs
for i=1:nm
    msel = Markers(:,i);
    Markers(msel == 0,i)= NaN;
    if any(msel==0)
%         disp(MLabels{ceil(i/3)});
    end
end
Markers=Markers * 0.001;          % convert from mm to m
writeMarkersToTRC(OutTRC, Markers, MLabels, VideoFrameRate, vFrms, vTime, 'm');  % write to TRC file
