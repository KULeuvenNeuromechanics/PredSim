% cutting the output data from predsim in gait cycles.
% assuming you have already run convertResultsForOpenSim.m



clear ; clc ; close



GRFdir = 'C:/Users/u0130218/Documents/MATLAB/PredSimResults/trackPred_new/Falisse_et_al_2022_20240821T120216_grf.mot';
fn = strsplit(GRFdir ,'/');
fn = char(fn{end});
fn = strrep(fn, '_grf.mot' , '.mot'); 

kindir = strrep(GRFdir, '_grf.mot' , '.mot'); 

% load the specific trial GRF
grfFile = read_motionFile_v40(GRFdir);

% load the motion file
kinFile = read_motionFile_v40(kindir);

% in theory there is only 1 gait cycle per foot - so lets ujsut cut the
% kinematics 

grfTags = {'groundr_force_vy', 'groundl_force_vy'};

for gt = 1:length(grfTags)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define the points based on vGRF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
   % find the correct column
   ci= 1;
   while strcmp(char(grfFile.labels{ci}), char(grfTags{gt})) == 0 
       ci = ci +1 ;
   end

    tempGRF = grfFile.data(:,ci);
    tempGRFtime = grfFile.data(:,1);
    
    thr = ones(length(tempGRFtime)) * 10;
    
    h = figure();
    plot(tempGRFtime,tempGRF, 'b', 'LineWidth', 2);
    hold on
    plot(tempGRFtime,thr, 'r', 'LineWidth', 1 , 'LineStyle', '--')
    [xi,yi] = getpts(h);
    waitfor(h)
    
    % lets jsut assume I dont mess up
    % start index
    diftimes = abs(tempGRFtime - xi(1)) ;
    % now find the minimum
    [ ~,indStart] = min(diftimes);
   
    % start index
    diftimes = abs(tempGRFtime - xi(2)) ;
    % now find the minimum
    [ ~ ,indEnd] = min(diftimes);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now cut
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    cutKinFile = kinFile ;
    cutKinFile.data = kinFile.data(indStart:indEnd,:) ;
    cutKinFile.nr = (indEnd - indStart) + 1;
    
    newKinFileName = strrep(kindir, '.mot', strcat('_' , char(grfTags(gt)), '_kinematics.mot'));
    write_motionFile_v40(cutKinFile, newKinFileName);   
end