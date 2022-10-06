

 % Button pushed function: VideoButton
        function VideoButtonPushed(app, event)
            global S
            global vis
            global model_state;
            try
                motFile = S.MotFile;
                if ~exist(motFile,'file')
                    app.SimulationResults.Value = [app.SimulationResults.Value; {'Please run simulation first, using default simulation as visualisation'}];
                    motFile = fullfile(S.pathRepo,'Results',S.ResultsFolder,'DefaultSim_q.mot'); % file with outputs
                end
                dat = ReadMotFile(motFile);
                dt = dat.data(2,1)-dat.data(1,1);
                N = length(dat.data(:,1));                
                CoordNamesAPI = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_r',...
                    'hip_flexion_l','lumbar_extension','knee_angle_r','knee_angle_l','ankle_angle_r','ankle_angle_l','mtp_angle_r','mtp_angle_l'};
                
                % column index in datastrcture for every coordName
                IndexCoord = nan(length(CoordNamesAPI),1);
                for i=1:length(CoordNamesAPI)
                    IndexCoord(i) = find(strcmp(CoordNamesAPI{i},dat.names))-1;
                end
                
                simbodyVis = vis.updSimbodyVisualizer();
                simbodyVis.setShowSimTime(true);
                Qvect = model_state.getQ();
                for i=[1:N 1]
                    dSel = dat.data(i,2:end);
                    for j=1:length(CoordNamesAPI)
                        if j==2 || j==3
                            Qvect.set(j-1,dSel(IndexCoord(j)));
                        else
                            Qvect.set(j-1,dSel(IndexCoord(j))*pi/180);
                        end
                    end
                    model_state.setQ(Qvect);
                    model_state.setTime(dat.data(i,1));
                    pause(dt);
                    vis.show(model_state);
                end
            catch
                app.SimulationResults.Value = [app.SimulationResults.Value; {'Unknown error in visualisation'}];
            end           
        end

        