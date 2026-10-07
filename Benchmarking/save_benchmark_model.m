function changed = save_benchmark_model(model,filename)
% Preserve unchanged model files and report when dynamics must be regenerated.
folder = fileparts(filename);
if ~isfolder(folder)
    mkdir(folder);
end
temporary = [tempname(folder) '.osim'];
cleanup = onCleanup(@() removeTemporary(temporary));
model.print(temporary);
changed = ~isfile(filename) || ~strcmp(fileread(filename),fileread(temporary));
if changed
    movefile(temporary,filename,'f');
end
end

function removeTemporary(filename)
if isfile(filename)
    delete(filename);
end
end
