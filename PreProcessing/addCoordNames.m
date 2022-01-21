function model_info = addCoordNames(model_info,group)
coordinates=fields(model_info.ExtFunIO.coordi);
for i = 1:length(model_info.ExtFunIO.jointi.(group))
    for c = 1:length(coordinates)
        if model_info.ExtFunIO.jointi.(group)(i) == model_info.ExtFunIO.coordi.(coordinates{c})
            model_info.ExtFunIO.jointi.names.(group){i} = coordinates{c};
        end
    end
end
