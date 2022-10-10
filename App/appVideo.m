classdef appVideo

properties (Access = private)

    % path to .osim model file
    path_osim
    % path to .mot motion file
    path_mot
    % path to folder with geometries to display model
    path_geom

    % cell array with coordinate names in api order
    CoordNamesAPI

    % opensim ModelVisualizer
    vis
    % opensim state
    model_state

end % end of private properties

methods (Access = public)

    % Constructor
    function obj = appVideo(osim_file,geom_folder,coord_names)
        obj.path_osim = osim_file;
        obj.path_geom = geom_folder;
        obj.CoordNamesAPI = coord_names;
    end

    % Destructor
    function delete(obj)
        obj.vis.delete();
        obj.model_state.delete();
        clear obj
    end

    % Initialise window with model in default pose
    function obj = initWindow(obj)
        import org.opensim.modeling.*
        osimModel = Model(obj.path_osim);
        osimModel.setUseVisualizer(true);
        obj.model_state = osimModel.initSystem();
        obj.vis = osimModel.getVisualizer();
%         h_vis = objectHandle(osimModel.getVisualizer());
%         obj.vis = h_vis;
        obj.vis.addDirToGeometrySearchPaths(obj.path_geom);
        pause(0.1);
        obj.vis.show(obj.model_state);
        simbodyVis = obj.vis.updSimbodyVisualizer();
        simbodyVis.setShowSimTime(true);
    end

    % Play motion once
    function playVideo(obj,mot_file)
        import org.opensim.modeling.*
        obj.path_mot = mot_file;

        dat = read_motionFile_v40(obj.path_mot);

        dt = dat.data(2,1)-dat.data(1,1);
        N = length(dat.data(:,1));                
        
        % column index in datastrcture for every coordName
        IndexCoord = nan(length(obj.CoordNamesAPI),1);
        for i=1:length(obj.CoordNamesAPI)
            IndexCoord(i) = find(strcmp(obj.CoordNamesAPI{i},dat.labels))-1;
        end
        dat.data(:,IndexCoord(2)+1) = dat.data(:,IndexCoord(2)+1) - mean(dat.data(:,IndexCoord(2)+1));
        
%         simbodyVis = obj.vis.updSimbodyVisualizer();
%         simbodyVis.setShowSimTime(true);
        Qvect = obj.model_state.getQ();
        for i=[1:N 1]
            dSel = dat.data(i,2:end);
            for j=1:length(obj.CoordNamesAPI)
                if j==2 || j==3
                    Qvect.set(j-1,dSel(IndexCoord(j)));
                else
                    Qvect.set(j-1,dSel(IndexCoord(j))*pi/180);
                end
            end
            obj.model_state.setQ(Qvect);
            obj.model_state.setTime(dat.data(i,1));
            pause(dt);
            obj.vis.show(obj.model_state);
        end
        
        Qvect.delete();
        clear Qvect
    end

end % end of public methods


end % end of classdef