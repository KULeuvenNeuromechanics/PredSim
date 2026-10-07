function idx = match_benchmark_id(id, identifiers)
% Match exact IDs, or a unique Schertzer condition within 0.005 m/s.
% Digitized Schertzer speeds vary slightly around the prescribed 4/5/6 km/h.
% Never relax the load location or mass, or choose among ambiguous matches.
idx = find(strcmp(id,identifiers));
if ~isempty(idx)
    return
end
pattern = '^schertzer2014_([0-9]+p[0-9]+)ms_(ankle|knee|torso)_([0-9]+)kg$';
condition = regexp(char(id),pattern,'tokens','once');
if isempty(condition)
    return
end
speed = str2double(strrep(condition{1},'p','.'));
for k = 1:numel(identifiers)
    candidate = regexp(char(identifiers{k}),pattern,'tokens','once');
    if ~isempty(candidate) && strcmp(condition{2},candidate{2}) && ...
            strcmp(condition{3},candidate{3}) && ...
            abs(speed-str2double(strrep(candidate{1},'p','.'))) <= 0.005
        idx(end+1) = k;
    end
end
end
