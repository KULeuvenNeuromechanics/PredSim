
clear
clc
[pathTests,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathTests);
addpath([pathRepo '\VariousFunctions'])

%
names_all = {'name1','name2','name3','name4','name5'};

%% proper example call 1
test_input = {{'name1','name4'},[5,3],0.5,'name2',[-3;200],1.3};
[values1,values2] = unpack_name_value_combinations(test_input,names_all,[2,1]);

%% error
% test_input = {{'name1','name4'},[5,3],0.5,'nam2',[-3;200],1.3};
% [values1,values2] = unpack_name_value_combinations(test_input,names_all,[2,1]);

%% error
% test_input = {{'name1','name4'},[5,3],0.5,'name1',[-3;200],1.3};
% [values1,values2] = unpack_name_value_combinations(test_input,names_all,[2,1]);

%% error
% test_input = {{'name1','name4'},[5,3],'name1',[-3;200],1.3};
% [values1,values2] = unpack_name_value_combinations(test_input,names_all,[2,1]);

%% proper example call 2
test_input = {{'name1','name4'},[],0.5,'name2',1.3,[]};
[values1,values2] = unpack_name_value_combinations(test_input,names_all,[1,1]);





