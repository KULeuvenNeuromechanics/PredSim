
clear
close all
clc


%% reference result
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathTests,~,~] = fileparts(pathHere);
[pathRepo,~,~] = fileparts(pathTests);

ref_file = fullfile(pathRepo,'Tests','ReferenceResults','Falisse_et_al_2022',...
    'Falisse_et_al_2022_paper.mat');
load(ref_file,'R','model_info');
R_ref = R;
model_info_ref = model_info;

x = 1:(100-1)/(size(R_ref.kinematics.Qs,1)-1):100;


test_file = fullfile(); % change this to path of result you want to test
load(test_file,'R','model_info')

figure
for i=1:31
    subplot(4,8,i)
    hold on
    plot(x,R_ref.kinematics.Qs(:,i),'DisplayName','ref')
    title(model_info_ref.ExtFunIO.coord_names.all{i},'Interpreter','none')
end

for i=1:31
    idx = model_info.ExtFunIO.coordi.(model_info_ref.ExtFunIO.coord_names.all{i});
    subplot(4,8,idx)
    hold on
    plot(x,R.Qs(:,i),'--','DisplayName','test')
end
legend

%%
figure
idx_r = [];
for i=1:length(model_info_ref.muscle_info.muscle_names)
    if strcmp(model_info_ref.muscle_info.muscle_names{i}(end),'r')
        idx_r(end+1) = i;
    end
end
for i=1:length(idx_r)
    subplot(7,7,i)
    hold on
    plot(x,R_ref.muscles.a(:,idx_r(i)),'DisplayName','ref')
    title(model_info_ref.muscle_info.muscle_names{idx_r(i)},'Interpreter','none')
end

idx_r2 = [];
for i=1:length(model_info.muscle_info.muscle_names)
    if strcmp(model_info.muscle_info.muscle_names{i}(end),'r')
        idx_r2(end+1) = i;
    end
end
for i=1:length(idx_r2)
    idx0 = find(strcmp(model_info.muscle_info.muscle_names(idx_r2),...
        model_info_ref.muscle_info.muscle_names(idx_r(i))));
%     idx = find(idx_r(:)==idx0);
    subplot(7,7,i)
    hold on
    plot(x,R.muscles.a(:,idx_r2(idx0)),'--','DisplayName','test')
end

legend


