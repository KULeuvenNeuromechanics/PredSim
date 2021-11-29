function [] = Plot3D_pwd(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dpath = pwd;
if ~isempty(varargin)
    dpath = varargin{1};
end


% Get the names of the results files
MatFiles = dir(fullfile(dpath,'*.mat'));
nFil = length(MatFiles);
ct = 1;
for i = 1:nFil
    filename = MatFiles(i).name;
    FileEnd = filename(end-6:end);
    OutName = fullfile(dpath,[filename(1:end-4) '_pp.mat']);
    if ~strcmp(FileEnd,'_pp.mat')
        Names{ct} = OutName;
        ct= ct+1;
    end
end

% plot influence Obj A
disp('Plotting: ');
DispHeader(Names);
disp('....');
CsV = hsv(ct);

C_NoExo     = [66 143 155]./255;
C_Passive   = [0.4167    0.2604    0.1658];
C_Active    = [0.6944    0.4340    0.2764];

CsV(2,:) = C_NoExo;
CsV(3,:) = C_Passive;
CsV(1,:) = C_Active;

h = figure();
set(h,'Position',[82         151        1497         827]);
for i = 1:ct-1
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim(Names{i},CsV(i,:),name,h);
end

filename = fullfile(dpath,'FigureResults.fig');
% if ~exist(filename,'file')
saveas(h,fullfile(dpath,'FigureResults.fig'));
% end


end

