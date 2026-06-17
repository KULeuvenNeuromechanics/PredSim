% take the outputs from pred sim and pull them out for use in opensim
% specifically I want to pull out muscle activations, forces, and ground
% reactions for re-running kineatics, and verification of joint contact
% loading reductions


% typcial
clear all
close all
clc
%% what do I need to do .

reps = {'A','B','C','D'};
stackedData = struct();
stackedData.stack = struct();

for x = 1:length(reps)

% create an entry in the struct
stackedData.(reps{x}) =struct();

kinMotDir = fullfile('C:/Users/u0130218/Documents/Data/HappyJoints_PredSimTesting/S204/COMAK/Gait/Right', char(reps{x}) ,'comak-inverse-kinematics/ik_generic/S204_ik_forPredSim.mot');
% start by laoding the kinemaitcs file to get teh time.
kinMotData = read_motionFile_v40(kinMotDir);

    for li = 1:length(kinMotData.labels)
        stackedData.(reps{x}).(kinMotData.labels{li}) =  kinMotData.data(:,li);
        if li == 1 % this means it is the time
            stackedData.stack.time(1,x) = kinMotData.data(end,li)-kinMotData.data(1,li);
            
        elseif strcmp((kinMotData.labels{li}), 'pelvis_tx')
            td = kinMotData.data(:,li)+10;
            stackedData.stack.pelvis_tx(1,x) = max(td) - min(td);
            
        elseif strcmp((kinMotData.labels{li}), 'pelvis_rotation')
            % NEED TO CHANGE TO BE BETWEN -90 and 90 
            if  kinMotData.data(1,li) < -90
                td = kinMotData.data(:,li) + 180;
                stackedData.stack.(kinMotData.labels{li})(x,:) = interpolation_100(td',100); 
            elseif kinMotData.data(1,li) > 90
                td = kinMotData.data(:,li) - 180;
                stackedData.stack.(kinMotData.labels{li})(x,:) = interpolation_100(td',100); 
            else
                stackedData.stack.(kinMotData.labels{li})(x,:) = interpolation_100(kinMotData.data(:,li)',100);  
            end

            
        else
            stackedData.stack.(kinMotData.labels{li})(x,:) = interpolation_100(kinMotData.data(:,li)',100);
        end
    end
end

% now loop throug the stacks and find the average.
stackedData.average=struct();

cols = fieldnames(stackedData.stack);
avgMotionData=[];

for slab = 1:length(cols)

    if slab == 1
        avgtime = mean(stackedData.stack.time);
        stackedData.average.time = linspace(0,avgtime,100)';
        
    elseif strcmp(cols{slab}, 'pelvis_tx')
        avglen = mean(stackedData.stack.pelvis_tx);
        stackedData.average.pelvis_tx = linspace(0,avglen,100)';
  
    else
        stackedData.average.(cols{slab}) = mean(stackedData.stack.(cols{slab}),1)';
    end
    avgMotionData(:,slab) = stackedData.average.(cols{slab});

end

outMotionData = struct();
outMotionData.nr = 100;
outMotionData.nc = kinMotData.nc;
outMotionData.inDeg = kinMotData.inDeg;
outMotionData.labels = cols';


outMotionData.data = avgMotionData;

% write out
write_motionFile_v40(outMotionData, 'AverageKineamtics.mot');
%

disp('done cuh')
