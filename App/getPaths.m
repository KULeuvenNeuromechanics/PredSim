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
        init_path_savefolder = 'O:\DemoVitruvianMan\ProfielStudent2024\Results';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.4\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
    elseif strcmp(name,'MSI')
        init_path_savefolder = 'D:\OneDrive - KU Leuven\Results_PredSim_App\Testing';
        init_path_geom = 'D:\software\OpenSim 4.3\Geometry';
        init_path_casadi = 'D:\software\casadi 3.5.5';
    else
        init_path_savefolder = 'C:\Users\u0138016\OneDrive - KU Leuven\Outreach\Kinderuniversiteit\2022\Resultaten';
        init_path_geom = 'C:\GBW_MyPrograms\OpenSim 4.3\Geometry';
        init_path_casadi = 'C:\GBW_MyPrograms\casadi_3_5_5';
    end


end