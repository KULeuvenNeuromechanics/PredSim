function [Markers,MLabels,VideoFrameRate,AnalogSignals,ALabels,AUnits,AnalogFrameRate,Event,ParameterGroup,CameraInfo]=...
    readC3D(FullFileName, varargin)
% ReadC3D:	Read 3D video coordinate and analog (FP) data from a C3D file
%
% Input:             FullFileName - file (including path) to be read
%       (optional)   force number of markers
%       (optional)   subtractAnalogOffset - indices of columns (in order of ALabels)
% Output:
% Markers            3D-marker data [NvideoFrames x Nmarkers*Ndim(=3)]
% MLabels            marker column identifiers
% VideoFrameRate     Frames/sec
% AnalogSignals      Analog signals [Nsignals x NanalogSamples ]
% ALabels            anolog column (signal) idenitifiers
% AnalogFrameRate    Samples/sec
% Event              Event(Nevents).time ..value  ..name
% ParameterGroup     ParameterGroup(Ngroups).Parameters(Nparameters).data ..etc.
% CameraInfo         MarkerRelated CameraInfo [Nmarkers x NvideoFrames]
% ResidualError      MarkerRelated ErrorInfo  [Nmarkers x NvideoFrames]
%
% AUTHOR(S) AND VERSION-HISTORY
% Ver. 1.0 Creation (Alan Morris, Toronto, October 1998) [originally named "getc3d.m"]
% Ver. 2.0 Revision (Jaap Harlaar, Amsterdam, april 2002)
% Ver. 3.0 Revision (Ajay Seth, Stanford University, August 2007)
% Ver. 3.1 Revision (Sjoerd Bruijn, KU leuven, May 2012. Speedup.)
% Note: removed ability to limit number of markers
% note: assumes scale <0
%==========================================================================

Markers=[];
MLabels=[];
VideoFrameRate=0;
AnalogSignals=[];
ALabels=[];
AnalogFrameRate=0;
Event=[];
ParameterGroup=[];
CameraInfo=[];
ResidualError=[];

% ###############################################
% ##                                           ##
% ##    open the file                          ##
% ##                                           ##
% ###############################################

ind=findstr(FullFileName,'\');
if ind>0, FileName=FullFileName(ind(length(ind))+1:length(FullFileName)); else FileName=FullFileName; end

fid=fopen(FullFileName,'r','n'); % native format (PC-intel)

if fid==-1,
    h=errordlg(['File: ',FileName,' could not be opened'],'application error');
    uiwait(h)
    return
end

NrecordFirstParameterblock=fread(fid,1,'int8');     % Reading record number of parameter section
key=fread(fid,1,'int8');                           % key = 80;

if key~=80,
    h=errordlg(['File: ',FileName,' does not comply to the C3D format'],'application error');
    uiwait(h)
    fclose(fid)
    return
end


fseek(fid,512*(NrecordFirstParameterblock-1)+3,'bof'); % jump to processortype - field
proctype=fread(fid,1,'int8')-83;                       % proctype: 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)

if proctype==2,
    fclose(fid);
    fid=fopen(FullFileName,'r','d'); % DEC VAX D floating point and VAX ordering
end

% ###############################################
% ##                                           ##
% ##    read header                            ##
% ##                                           ##
% ###############################################

%NrecordFirstParameterblock=fread(fid,1,'int8');     % Reading record number of parameter section
%key1=fread(fid,1,'int8');                           % key = 80;

fseek(fid,2,'bof');


Nmarkers=fread(fid,1,'int16');			        %number of markers

if nargin > 1,
    if ~isempty(varargin{1}),
        Nmarkersinput = varargin{1};
    end
end

NanalogSamplesPerVideoFrame=fread(fid,1,'int16');			%number of analog channels x #analog frames per video frame
StartFrame=fread(fid,1,'int16');		        %# of first video frame

EndFrame=fread(fid,1,'int16');			        %# of last video frame

MaxInterpolationGap=fread(fid,1,'int16');		%maximum interpolation gap allowed (in frame)

Scale=fread(fid,1,'float32');			        %floating-point scale factor to convert 3D-integers to ref system units
if Scale>0
    warndlg('Program was written to use scale < 0. Apparently, this is not the case for this file. Reverting back to old readc3d code, will be slow!')
end
NrecordDataBlock=fread(fid,1,'int16');			%starting record number for 3D point and analog data

NanalogFramesPerVideoFrame=fread(fid,1,'int16');
if NanalogFramesPerVideoFrame > 0,
    NanalogChannels=NanalogSamplesPerVideoFrame/NanalogFramesPerVideoFrame;
else
    NanalogChannels=0;
end


VideoFrameRate=fread(fid,1,'float32');
AnalogFrameRate=VideoFrameRate*NanalogFramesPerVideoFrame;

% ###############################################
% ##                                           ##
% ##    read events                            ##
% ##                                           ##
% ###############################################

fseek(fid,298,'bof');
EventIndicator=fread(fid,1,'int16');
if EventIndicator==12345,
    Nevents=fread(fid,1,'int16');
    fseek(fid,2,'cof'); % skip one position/2 bytes
    if Nevents>0,
        for i=1:Nevents,
            Event(i).time=fread(fid,1,'float');
        end
        fseek(fid,188*2,'bof');
        for i=1:Nevents,
            Event(i).value=fread(fid,1,'int8');
        end
        fseek(fid,198*2,'bof');
        for i=1:Nevents,
            Event(i).name=cellstr(char(fread(fid,4,'char')'));
        end
    end
end


% ###############################################
% ##                                           ##
% ##    read 1st parameter block               ##
% ##                                           ##
% ###############################################

fseek(fid,512*(NrecordFirstParameterblock-1),'bof');

dat1=fread(fid,1,'int8');
key2=fread(fid,1,'int8');                   % key = 80;
NparameterRecords=fread(fid,1,'int8');
proctype=fread(fid,1,'int8')-83;            % proctype: 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)


Ncharacters=fread(fid,1,'int8');   			% characters in group/parameter name
GroupNumber=fread(fid,1,'int8');				% id number -ve=group / +ve=parameter


while Ncharacters > 0 % The end of the parameter record is indicated by <0 characters for group/parameter name
    
    if GroupNumber<0 % Group data
        GroupNumber=abs(GroupNumber);
        GroupName=fread(fid,[1,Ncharacters],'char');
        ParameterGroup(GroupNumber).name=cellstr(char(GroupName));	%group name
        offset=fread(fid,1,'int16');							%offset in bytes
        deschars=fread(fid,1,'int8');							%description characters
        GroupDescription=fread(fid,[1,deschars],'char');
        ParameterGroup(GroupNumber).description=cellstr(char(GroupDescription)); %group description
        
        ParameterNumberIndex(GroupNumber)=0;
        fseek(fid,offset-3-deschars,'cof');
        
        
    else % parameter data
        clear dimension;
        ParameterNumberIndex(GroupNumber)=ParameterNumberIndex(GroupNumber)+1;
        ParameterNumber=ParameterNumberIndex(GroupNumber);              % index all parameters within a group
        
        ParameterName=fread(fid,[1,Ncharacters],'char');				% name of parameter
        
        % read parameter name
        if size(ParameterName)>0
            ParameterGroup(GroupNumber).Parameter(ParameterNumber).name=cellstr(char(ParameterName));	%save parameter name
        end
        
        % read offset
        offset=fread(fid,1,'int16');							%offset of parameters in bytes
        filepos=ftell(fid);										%present file position
        nextrec=filepos+offset(1)-2;							%position of beginning of next record
        
        
        % read type
        type=fread(fid,1,'int8');     % type of data: -1=char/1=byte/2=integer*2/4=real*4
        ParameterGroup(GroupNumber).Parameter(ParameterNumber).datatype=type;
        
        
        % read number of dimensions
        dimnum=fread(fid,1,'int8');
        if dimnum==0
            datalength=abs(type);								%length of data record
        else
            mult=1;
            for j=1:dimnum
                dimension(j)=fread(fid,1,'uint8');
                mult=mult*dimension(j);
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).dim(j)=dimension(j);  %save parameter dimension data
            end
            datalength=abs(type)*mult;							%length of data record for multi-dimensional array
        end
        
        
        if type==-1 %datatype=='char'
            
            wordlength=dimension(1);	%length of character word
            if dimnum==2 & datalength>0 %& parameter(idnumber,index,2).dim>0
                for j=1:dimension(2)
                    data=fread(fid,[1,wordlength],'char');	%character word data record for 2-D array
                    ParameterGroup(GroupNumber).Parameter(ParameterNumber).data(j)=cellstr(char(data));
                end
                
            elseif dimnum==1 & datalength>0
                data=fread(fid,[1,wordlength],'char');		%numerical data record of 1-D array
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data=cellstr(char(data));
            end
            
        elseif type==1    %1-byte for boolean
            
            Nparameters=datalength/abs(type);
            data=fread(fid,Nparameters,'int8');
            ParameterGroup(GroupNumber).Parameter(ParameterNumber).data=data;
            
        elseif type==2 & datalength>0			%integer
            
            Nparameters=datalength/abs(type);
            data=fread(fid,Nparameters,'int16');
            if dimnum>1
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data=reshape(data,dimension);
            else
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data=data;
            end
            
        elseif type==4 & datalength>0
            
            Nparameters=datalength/abs(type);
            data=fread(fid,Nparameters,'float');
            if dimnum>1
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data=reshape(data,dimension);
            else
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data=data;
            end
        else
            % error
        end
        
        deschars=fread(fid,1,'int8');							%description characters
        if deschars>0
            description=fread(fid,[1,deschars],'char');
            ParameterGroup(GroupNumber).Parameter(ParameterNumber).description=cellstr(char(description));
        end
        %moving ahead to next record
        fseek(fid,nextrec,'bof');
    end
    
    % check group/parameter characters and idnumber to see if more records present
    Ncharacters=fread(fid,1,'int8');   			% characters in next group/parameter name
    GroupNumber=fread(fid,1,'int8');				% id number -ve=group / +ve=parameter
end


% ###############################################
% ##                                           ##
% ##    read data block                        ##
% ##                                           ##
% ###############################################
%  Get the coordinate and analog data

fseek(fid,(NrecordDataBlock-1)*512,'bof');


NvideoFrames=EndFrame - StartFrame + 1;
% h = waitbar(0,[FileName,' is loading...']);
% Original marker structure is annoying to work with and to output to file
% so changing here so that it is listed three coordinates at a time
% in order of the marker labels. -aseth
if Scale<0
    reps = [int2str(4*Nmarkers),'*float32'];
    if Nmarkers>0
        tmpMKR = fread(fid,inf,reps,4*NanalogFramesPerVideoFrame*NanalogChannels);   % note four bytes per analog channel data skipped
        
        % Reshape/resize/reorder data
        % First do the markers
        tmpMKR = reshape(tmpMKR,Nmarkers*4,length(tmpMKR)./(4*Nmarkers)) ;           % go from one long column to a 4xNMxNVF matrix
        
        Markers = tmpMKR';
        
        Markers(:,4:4:end)=[];% trim off the residual/camera contribution stuff%
    else
        Markers=ones(0,3)*nan;
    end
%     waitbar(0.3,h)
    % Get camera contribution/residual
    if Nmarkers>0
        a = fix(tmpMKR(4:4:end,:));
        highbyte=fix(a/256);
        lowbyte=a-highbyte*256;
        CameraInfo=highbyte;
        ResidualError=lowbyte*abs(Scale);
    else
        CameraInfo=nan;
        ResidualError=nan;
    end
%      waitbar(0.6,h)
    if NanalogChannels > 0
        % Read in analog data as repeated blocks
        fseek(fid,(NrecordDataBlock-1)*512+4*4*Nmarkers,'bof');
        reps = [int2str(NanalogFramesPerVideoFrame*NanalogChannels),'*float32'];
        tmpANL = fread(fid,inf,reps,4*4*Nmarkers);   % note: four bytes per 3d data skipped
        % Reshape the analog data
        AnalogSignals = reshape(tmpANL,NanalogChannels,length(tmpANL)/NanalogChannels)';
    else
        AnalogSignals = 0;
    end
%      waitbar(1,h)
else
    for i=1:NvideoFrames
        for j=1:Nmarkers
            Markers(i,3*j-2:3*j)=fread(fid,3,'int16')'.*Scale;
            ResidualError(i,j)=fread(fid,1,'int8');
            CameraInfo(i,j)=fread(fid,1,'int8');
        end
%         waitbar(i/NvideoFrames,h)
        for j=1:NanalogFramesPerVideoFrame,
            AnalogSignals(j+NanalogFramesPerVideoFrame*(i-1),1:NanalogChannels)=...
                fread(fid,NanalogChannels,'int16')';
        end
    end
end
% close(h) % waitbar
fclose(fid);



I = 1; J = 1;
while ~strcmp(ParameterGroup(I).name, 'POINT');
    I = I+1;
end

while ~strcmp(ParameterGroup(I).Parameter(J).name, 'LABELS');
    J = J+1;
end


MLabels = ParameterGroup(I).Parameter(J).data;


% Get details about analog data
I = 1; J = 1;
while ~strcmp(ParameterGroup(I).name, 'ANALOG');
    I = I+1;
end

for J = 1:length(ParameterGroup(I).Parameter),
    category = upper(ParameterGroup(I).Parameter(J).name);
    switch category{1}
        case 'GEN_SCALE'
            gen_scale = ParameterGroup(I).Parameter(J).data;
        case 'SCALE'
            scales = ParameterGroup(I).Parameter(J).data;
        case 'OFFSET'
            offsets = ParameterGroup(I).Parameter(J).data;
        case 'LABELS'
            ALabels = ParameterGroup(I).Parameter(J).data;
        case 'UNITS'
            AUnits = ParameterGroup(I).Parameter(J).data;
    end
end
% number of analog frames and individual signals
[nAF, nAS] = size(AnalogSignals);
% Apply the appropriate scaling and offsets to convert analog voltage to
% appropriate units

if nAS,
    offsets = offsets(1:nAS);
    scales = scales(1:nAS);
    AnalogSignals = (AnalogSignals-ones(nAF, 1)*offsets').*(ones(nAF, 1)*scales')*gen_scale;
end

if nargin > 2,
    if ~isempty(varargin{2}),
        offInds = varargin{2};
        offSignals = AnalogSignals(:,offInds);
        % super smooth and rectify the data
        offSignals = abs(smooth(offSignals, 10, AnalogFrameRate));
        dOffSignals = deriv(offSignals,1/AnalogFrameRate);
        % Use a 2% threshold
        thresh = 0.02*max(dOffSignals);
        for I = 1:length(offInds),
            % find the starting index where data begins
            startInd = min(find(dOffSignals(:,I) >= thresh(I)));
            background = mean(AnalogSignals(1:startInd,offInds(I)));
            AnalogSignals(:,offInds(I)) = AnalogSignals(:,offInds(I))-background;
        end
    end
end

% ======================
% end readc3D.m

