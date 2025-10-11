function [init_path_savefolder,init_path_geom,init_path_casadi] = getPaths()
    name = getenv('COMPUTERNAME');

    % fill these in based on your computer
    if strcmp(name,'GBW-L-W2122')
        init_path_savefolder = 'C:\Users\u0150099\OneDrive - KU Leuven\Results_PredSim_App\TGCS';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.3\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
    elseif strcmp(name, 'GBW-L-W2075')
        init_path_savefolder = 'C:\Users\u0138016\OneDrive - KU Leuven\Outreach\Kinderuniversiteit\2022\Resultaten';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.3\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
    elseif strcmp(name,'GBW-L-W2195')
        init_path_savefolder = 'C:\GBW_MyPrograms\KinderuniversiteitApp\Resultaten';
        init_path_geom = 'C:\OpenSim 4.4\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
    elseif strcmp(name,'GBW-L-W2109')
        init_path_savefolder = 'C:\GBW_MyPrograms\PredSim\Subjects';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim\OpenSim 4.3\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\Casadi\casadi-windows-matlabR2016a-v3.5.5';
    elseif strcmp(name,'GBW-L-W4394')
        init_path_savefolder = 'C:\Users\u0125183\OneDrive - KU Leuven\Kinderuniversiteit';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.5\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\Casadi';
 elseif strcmp(name, 'GBW-L-W2099')
        init_path_savefolder = 'C:\Users\u0167448\OneDrive - KU Leuven\Kinderuniversiteit2025';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.5\Geometry';
        init_path_casadi = 'C:\Users\u0167448\Documents\GitHub\casadi-windows-matlabR2016a-v3.5.5';
     elseif strcmp(name, 'GBW-L-W4275')
        init_path_savefolder = 'C:\GBW_MyPrograms\Kinderuniversiteit2025';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.4\Geometry';
        init_path_casadi = 'C:\Users\u0046458\Documents\Code\casadi-3.6.3-windows64-matlab2018b';
           elseif strcmp(name, 'GBW-L-W2239')
        init_path_savefolder = 'C:\Users\u0167448\OneDrive - KU Leuven\Kinderuniversiteit2025';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.5\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\casadi-3.7.1-windows64-matlab2018b';

            elseif strcmp(name, 'LAPTOP-E7OQMBKV')
               init_path_savefolder = 'C:\Users\timvd\OneDrive - KU Leuven\Kinderuniversiteit2025';
        init_path_geom = 'C:\OpenSim 4.5\Geometry';
        init_path_casadi = 'C:\Users\timvd\Documents\casadi-3.7.1-windows64-matlab2018b';
    % elseif strcmp(name,'MSI')
    %     init_path_savefolder = 'D:\OneDrive - KU Leuven\Results_PredSim_App\Testing';
    %     init_path_geom = 'D:\software\OpenSim 4.3\Geometry';
    %     init_path_casadi = 'D:\software\casadi 3.5.5';
    % else
        % init_path_savefolder = 'C:\Users\u0138016\OneDrive - KU Leuven\Outreach\Kinderuniversiteit\2022\Resultaten';
        % init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.3\Geometry';
        % init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
    end


end