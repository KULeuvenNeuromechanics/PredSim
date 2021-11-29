%% CreateFoldersFor .dll file
%-----------------------------

% path information
PathOsimSource = 'C:\opensim-ad-core-source';
SoftwareOCP_predsim = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim';
PathExtFunc = 'C:\opensim-ExternalFunc';
OpenSimBuild = 'C:\opensim-ad-core-build2';
VisualStudioInstall = 'C:\Program Files (x86)\Microsoft Visual Studio 14.0';
nInputDll = 93;

% 
% name of the new project
Name = 'Example_AutomatedCreation';

% path info
PathStart = pwd;

%% OpenSim source - external function 
% create folder in the opensim source folder
ExtFuncPath = fullfile(PathOsimSource,'OpenSim\External_Functions',Name);
mkdir(ExtFuncPath);

% create the CMakeList file
fid = fopen(fullfile(ExtFuncPath,'CMakeLists.txt'),'wt');
fprintf( fid, '%s\n', ['set(TEST_TARGET ' Name ')']);
fprintf( fid, '%s\n', 'add_executable(${TEST_TARGET} ${TEST_TARGET}.cpp)');
fprintf( fid, '%s\n', 'target_link_libraries(${TEST_TARGET} osimSimulation)');
fprintf( fid, '%s\n', 'set_target_properties(${TEST_TARGET} PROPERTIES');
fprintf( fid, '%s\n', '    FOLDER "External_Functions")');
fclose(fid);
% copy the cpp files
PathCpp = fullfile(SoftwareOCP_predsim,'ExternalFunctions\CppFiles');
copyfile(fullfile(PathCpp,[Name '.cpp']),fullfile(ExtFuncPath,[Name '.cpp']));
% 
% % create folder in
% ExtFuncPath = fullfile(PathOsimSource,'OpenSim\External_Functions',[Name '_pp']);
% mkdir(ExtFuncPath);
% % create the CMakeList file
% fid = fopen(fullfile(ExtFuncPath,'CMakeLists.txt'),'wt');
% fprintf( fid, '%s\n', ['set(TEST_TARGET ' Name '_pp']);
% fprintf( fid, '%s\n', 'add_executable(${TEST_TARGET} ${TEST_TARGET}.cpp)');
% fprintf( fid, '%s\n', 'target_link_libraries(${TEST_TARGET} osimSimulation)');
% fprintf( fid, '%s\n', 'set_target_properties(${TEST_TARGET} PROPERTIES');
% fprintf( fid, '%s\n', 'FOLDER "External_Functions"');
% fclose(fid);
% % copy the cpp files
% PathCpp = fullfile(PathSoftware,'ExternalFunctions\CppFiles');
% copyfile(fullfile(PathCpp,[Name '_pp.cpp']),fullfile(ExtFuncPath,[Name '_pp.cpp']));

% add the project to the current cmakelist of the the external projects
CmakeFile = fullfile(PathOsimSource,'OpenSim\External_Functions','CMakeLists.txt');
ListExtFunc = importdata(CmakeFile);
delete(CmakeFile);
fid = fopen(fullfile(PathOsimSource,'OpenSim\External_Functions','CMakeLists.txt'),'wt');
for i =1:length(ListExtFunc)-1
    fprintf( fid, '%s\n', ListExtFunc{i});
end
% check if external function is not already in the list
if ~any(strcmp(Name,ListExtFunc))
    fprintf( fid, '%s\n', ['	add_subdirectory(' Name ')']);
%     fprintf( fid, '%s\n', ['	add_subdirectory(' Name '_pp)']);
else
    disp('External function was already included in the CmakeList');
end
fprintf( fid, '%s\n', ListExtFunc{end});
fclose(fid);

%% Add new project to opensim build
% update visual studio project with cmake
% go to folder where you want to build the projects
% then run for example: "cmake C:\sourceTest -G "Visual Studio 14 2015
% Win64"
disp('... running cmake to update opensim build');
cd(OpenSimBuild);
TextCommand = ['cmake ' PathOsimSource ' -G "Visual Studio 14 2015 Win64"'];
system(TextCommand);

% Now Run visual studio from the command line 
% TextCommand = [OpenSimBuild '\OpenSim.sln'];
% system(TextCommand);
% disp(['Select RelWidthDebInfo, right click on ' Name 'project and build project']);
% wait('Hit enter when building is finished');

% or run it automatically using devnenv
% information: devenv mysln.sln /build Debug /project proj1
disp('... building project in visual studio');
ExeDevPath = fullfile(VisualStudioInstall,'Common7\IDE');
cd(ExeDevPath);
TextCommand = ['devenv.exe ' OpenSimBuild '\OpenSim.sln ' '/Build ' 'RelWithDebInfo ' '/Project ' Name ];
% TextCommand2 = ['devenv.exe ' OpenSimBuild '\OpenSim.sln ' '/Build ' 'RelWithDebInfo ' '/Project ' Name '_pp'];
system(TextCommand);
% system(TextCommand2);

%% Create executable, foo.m and foo_jac
% - first: check if build is finished
ExePath = fullfile(OpenSimBuild,'RelWithDebInfo');
File1 = fullfile(ExePath,[Name '.exe']);
% File2 = fullfile(ExePath,[Name '_pp.exe']);
% BoolExist = exist(File1,'file') && exist(File2,'file');
% BoolExist = exist(File1,'file');
% while ~BoolExist
%     %    BoolExist = exist(File1,'file') && exist(File2,'file');
%     BoolExist = exist(File1,'file');    
% end

% - second: create path with cgeneration info:
CgenDir1 = fullfile(PathOsimSource,'cgeneration',Name);
% CgenDir2 = fullfile(PathOsimSource,'cgeneration',Name);
mkdir(CgenDir1);
% mkdir(CgenDir2);

% - third: run the exectuables and copy foo.m
disp('... Running exectuble');

cd(ExePath);
if exist(fullfile(ExePath,'foo.m'),'file')
    delete(fullfile(ExePath,'foo.m'));
end
system([Name '.exe']);
% while ~exist(fullfile(ExePath,'foo.m'),'file')
% end
copyfile(fullfile(ExePath,'foo.m'),fullfile(CgenDir1,'foo.m'));
delete(fullfile(ExePath,'foo.m'));
% 
% system([Name '_pp.exe']);
% while ~exist(fullfile(ExePath,'foo.m'),'file')
% end
% copyfile(fullfile(ExePath,'foo.m'),fullfile(CgenDir2,'foo.m'));
% delete(fullfile(ExePath,'foo.m'));

% - fourth: create the foo_jac.m using a matlab function
disp('... Creating foo_jac.c');
% cd(CgenDir1)
generate_foo_jac(nInputDll,CgenDir1);
% cd(CgenDir2)
% disp('Creating foo_jac.m 2/2');
% generate_foo_jac(nInputDll);

%% Run cmake to build external function project

% First: create the cmake file

fid = fopen(fullfile(CgenDir1,'CMakeLists.txt'),'wt');
fprintf( fid, '%s\n', ['set(TEST_TARGET ' Name ')']);
fprintf( fid, '%s\n', 'project(${TEST_TARGET})');
fprintf( fid, '%s\n', 'add_library(${TEST_TARGET} SHARED foo_jac.c)');
fprintf( fid, '%s\n', 'install(TARGETS ${TEST_TARGET}');
fprintf( fid, '%s\n', '	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE');
fprintf( fid, '%s\n', '				GROUP_READ GROUP_WRITE GROUP_EXECUTE');
fprintf( fid, '%s\n', '				WORLD_READ WORLD_EXECUTE');
fprintf( fid, '%s\n', '	LIBRARY DESTINATION lib');
fprintf( fid, '%s\n', '	ARCHIVE DESTINATION lib');
fprintf( fid, '%s\n', '	RUNTIME DESTINATION bin');
fprintf( fid, '%s\n', '	)');
fclose(fid);
% 
% fid = fopen(fullfile(CgenDir2,'CMakeLists.txt'),'wt');
% fprintf( fid, '%s\n', ['set(TEST_TARGET ' Name '_pp']);
% fprintf( fid, '%s\n', 'project(${TEST_TARGET})');
% fprintf( fid, '%s\n', 'add_library(${TEST_TARGET} SHARED foo_jac.c)');
% fprintf( fid, '%s\n', 'install(TARGETS ${TEST_TARGET}');
% fprintf( fid, '%s\n', '	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE');
% fprintf( fid, '%s\n', '				GROUP_READ GROUP_WRITE GROUP_EXECUTE');
% fprintf( fid, '%s\n', '				WORLD_READ WORLD_EXECUTE');
% fprintf( fid, '%s\n', '	LIBRARY DESTINATION lib');
% fprintf( fid, '%s\n', '	ARCHIVE DESTINATION lib');
% fprintf( fid, '%s\n', '	RUNTIME DESTINATION bin');
% fprintf( fid, '%s\n', '	)');
% fclose(fid);

% second: create folder for external functions
mkdir(fullfile(PathExtFunc,Name,'build'));
mkdir(fullfile(PathExtFunc,Name,'install'));
% mkdir(fullfile(PathExtFunc,[Name '_pp'],'build'));
% mkdir(fullfile(PathExtFunc,[Name '_pp'],'install'));

% Run cmake again
disp('... running cmake to create project external function');
cd(fullfile(PathExtFunc,Name,'build'));
TextCommand = ['cmake ' CgenDir1 ' -G "Visual Studio 14 2015 Win64"'];
system(TextCommand);

% cd(fullfile(PathExtFunc,[Name '_pp'],'build'));
% TextCommand = ['cmake ' CgenDir2 ' -G "Visual Studio 14 2015 Win64"'];
% system(TextCommand);


%% create .dll file using the visual studio

% go to path with executable to run visual studio from command line
disp('... Build and install external func project in visual studio ');
cd(ExeDevPath);
% run for firs tproject
ExtFunc_Build = fullfile(PathExtFunc,Name,'build');
TextCommand = ['devenv.exe ' ExtFunc_Build '\' Name '.sln /Build ' 'RelWithDebInfo ' '/Project ALL_BUILD' ];
system(TextCommand);
TextCommand = ['devenv.exe ' ExtFunc_Build '\' Name '.sln /Build ' 'RelWithDebInfo ' '/Project INSTALL' ];
system(TextCommand);
% ExtFunc_Build = fullfile(PathExtFunc,[Name '_pp'],'build');
% TextCommand = ['devenv.exe ' ExtFunc_Build '\' Name '_pp.sln /Build ' 'RelWithDebInfo ' '/Project ALL_BUILD' ];
% system(TextCommand);

%% Copy the .dll files to your predictive simulation folder

if ~exist(fullfile(SoftwareOCP_predsim,'ExternalFunctions',[Name '.dll']),'file')
    copyfile(fullfile(PathExtFunc,Name,'build','RelWithDebInfo',[Name '.dll']),...
        fullfile(SoftwareOCP_predsim,'ExternalFunctions'));
else
    disp(['cannot copy file ' Name '.dll because this is already in your external function folder ', ...
        fullfile(SoftwareOCP_predsim,'ExternalFunctions')]);
end
cd(PathStart);









