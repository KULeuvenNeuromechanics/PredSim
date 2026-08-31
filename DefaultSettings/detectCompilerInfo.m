function info = detectCompilerInfo()

info = struct();
info.ConfiguredMexCompiler = [];
info.InstalledVisualStudios = {};
info.InstalledMSVCVersions = {};
info.DetectedCLOnPath = [];
info.CMakeGenerators = {};
info.LatestCMakeGenerator = [];

%% ==========================================================
% 1) Detect configured MEX compiler
%% ==========================================================
try
    cfg = mex.getCompilerConfigurations('C++','Selected');
    if ~isempty(cfg)
        info.ConfiguredMexCompiler = struct( ...
            'Name', cfg.Name, ...
            'Version', cfg.Version, ...
            'Location', cfg.Location);
    else
        info.ConfiguredMexCompiler = 'No MEX C++ compiler configured';
    end
catch ME
    info.ConfiguredMexCompiler = ['Error: ' ME.message];
end


%% ==========================================================
% 2) Detect installed Visual Studio instances
%% ==========================================================
if ispc
    
    vswherePath = fullfile( ...
        getenv('ProgramFiles(x86)'), ...
        'Microsoft Visual Studio', ...
        'Installer', ...
        'vswhere.exe');
    
    if exist(vswherePath,'file')
        cmd = sprintf('"%s" -products * -format json', vswherePath);
        [status, cmdout] = system(cmd);
        
        if status == 0 && ~isempty(cmdout)
            vsData = jsondecode(cmdout);
            
            for k = 1:numel(vsData)
                vs = struct();
                vs.DisplayName = vsData(k).displayName;
                vs.Version = vsData(k).installationVersion;
                vs.Path = vsData(k).installationPath;
                
                % Extract major version
                major = str2double(extractBefore(vs.Version,'.'));
                vs.MajorVersion = major;
                
                % Map to CMake generator
                vs.CMakeGenerator = mapToCMakeGenerator(major);
                
                info.InstalledVisualStudios{end+1} = vs;
                
                if ~isempty(vs.CMakeGenerator)
                    info.CMakeGenerators{end+1} = vs.CMakeGenerator;
                end
            end
        end
    end
end


%% ==========================================================
% 3) Select newest available CMake generator
%% ==========================================================
if ~isempty(info.InstalledVisualStudios)
    
    majors = cellfun(@(x) x.MajorVersion, info.InstalledVisualStudios);
    [~, idx] = max(majors);
    
    info.LatestCMakeGenerator = ...
        info.InstalledVisualStudios{idx}.CMakeGenerator;
end


%% ==========================================================
% 4) Pretty print
%% ==========================================================
if nargout == 0
    
    fprintf('\n========== COMPILER DETECTION ==========\n\n');
    
    fprintf('Configured MEX Compiler:\n');
    disp(info.ConfiguredMexCompiler)
    
    fprintf('\nInstalled Visual Studios:\n');
    disp(info.InstalledVisualStudios)
    
    fprintf('\nAvailable CMake Generators:\n');
    disp(info.CMakeGenerators)
    
    fprintf('\nNewest Recommended CMake Generator:\n');
    disp(info.LatestCMakeGenerator)
    
    fprintf('\n========================================\n\n');
end

end


%% ==========================================================
% Helper: Map VS major version to CMake generator string
%% ==========================================================
function gen = mapToCMakeGenerator(major)

switch major
    
    case 14
        gen = 'Visual Studio 14 2015 Win64';
        
    case 15
        gen = 'Visual Studio 15 2017 Win64';
        
    case 16
        % VS 2019 — modern CMake style preferred
        gen = 'Visual Studio 16 2019';
        
    case 17
        % VS 2022
        gen = 'Visual Studio 17 2022';
     
    case 18
        % VS 2026
        gen = 'Visual Studio 18 2026';
        
    otherwise
        gen = [];
end

end