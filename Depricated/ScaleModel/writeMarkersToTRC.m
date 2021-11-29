function err = writeMarkersToTRC(trcfile, Markers, MLabels, Rate, Frames, Time, Units)
% Write 3D Markers trajectories (real or virtual) to a .trc file                        
% USAGE: error = writeMarkersToTRC(trcFile, Markers, MLabels, Rate, Frames, Time, Units)
%
% ASeth 10-07 
% Stanford University

err = 0;
[nvF, nc] = size(Markers);
nM = length(MLabels);

if (nM >= nc/3),
    % Maybe more labels than we need
    nM = nc/3;
    MLabels = MLabels(1:nM);
else
    % number of labels does not correspond to the number of Markers
    error('number of labels does not correspond to the number of Markers');
    err = 1;
end

if isempty(Frames),
    vFrms = [1:nvF]';
    vTime = 1/Rate*(vFrms);
else
    if (length(Frames) ~= nvF),
        error('number of frames does not correspond to the length of Markers');
        err = 1;
    end
    vFrms = Frames;
    vTime = Time;
end
    
% Assemble Marker data for writing out to .trc file
data = [vFrms vTime Markers];


% Generate the header for the .trc file
fid = fopen(trcfile, 'wt');


fprintf(fid, 'PathFileType\t4\t(X/Y/Z)\t%s\n', trcfile);
fprintf(fid, 'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid, '%f\t%f\t%d\t%d\t%s\t%f\t%d\t%d\n', ...
Rate, Rate, nvF, nM, Units, Rate, vFrms(1), vFrms(end));
fprintf(fid, 'Frame#\tTime');

fprintf(fid,'\t%s',MLabels{1});

for I = 2:nM,
    fprintf(fid,'\t\t\t%s', MLabels{I});
end
fprintf(fid, '\n\t');

for I = 1:nM,
    fprintf(fid,'\tX%i\tY%i\tZ%i', I,I,I);
end
fprintf(fid, '\n\n');

fclose(fid);

% Now append the data to the file now that header has been written out.
dlmwrite(trcfile, data, '-append', 'delimiter', '\t');