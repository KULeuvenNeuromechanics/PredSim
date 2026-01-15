% determine specific tension for the muscles in the 2D model

muscles_2D = {'hamstring', 'glut_max', 'iliopsoas', 'vasti', 'gastroc'};
muscles_3D = struct();
muscles_3D(1).names = {'semimem', 'semiten', 'bifemlh'};
muscles_3D(2).names = {'glut_max1', 'glut_max2', 'glut_max3'};
muscles_3D(3).names = {'iliacus', 'psoas'};
muscles_3D(4).names = {'vas_med', 'vas_int','vas_lat'};
muscles_3D(5).names = {'med_gas', 'lat_gas'};

ST_wa = nan(1, length(muscles_2D));

for k = 1:length(muscles_2D)
    mnames = muscles_3D(k).names;

    Fmo = nan(1, length(mnames));
    ST = nan(1, length(mnames));

    for i = 1:length(muscle_info.parameters)
        for j = 1:length(mnames)
            if strcmp(muscle_info.parameters(i).muscle_name, [mnames{j}, '_r'])
                Fmo(j) = muscle_info.parameters(i).FMo;
                ST(j) = specific_tension(i);
            end
        end
    end

    ST_wa(k) = sum(Fmo .* ST) / sum(Fmo);
end


