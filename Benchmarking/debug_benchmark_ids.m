function [matches,experiments,simulations] = debug_benchmark_ids(S,osim_path,S_benchmark,varargin)
%DEBUG_BENCHMARK_IDS Audit production ID assignments without running PredSim.
% [matches,experiments,simulations] = debug_benchmark_ids(S,osim_path,S_benchmark)
% Optional 'DataFolder': local benchmark JSON folder (no download when given).
% Optional 'ReportFolder': write two CSV tables; default '' writes nothing.
% The default data download uses download_benchmarkdata(false).
% No model conversion, CasADi, OpenSim, batch jobs or result files are needed.
% matches contains one row per simulation and includes ALL candidate IDs for
% ambiguous matches. experiments gives reverse coverage, including unused data.
% Slope columns are percentages: the five supported datasets store e.g. 8,
% whereas model metadata uses 0.08. AddedMass is the bilateral model total;
% ExperimentalAddedMass is the JSON value, without silently rescaling it.
% Gomenuka loads are compared as fractions of unloaded subject/model mass.
% Load locations absent from JSON are checked against the experimental ID
% where its format permits it, and labelled accordingly in LocationSource.

p = inputParser;
addParameter(p,'DataFolder','',@(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p,'ReportFolder','',@(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p,varargin{:});
if ~isfield(S_benchmark,'studies')
    S_benchmark.studies = {};
end
datafolder = char(p.Results.DataFolder);
if isempty(datafolder)
    datafolder = download_benchmarkdata(false);
end
assert(isfolder(datafolder),'PredSim:MissingBenchmarkData', ...
    'Experimental data folder does not exist: %s',datafolder);
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'Benchmarking','data','functions'));
[data,studies,identifiers] = get_all_benchmarkdata(datafolder);
assert(~isempty(data),'PredSim:MissingBenchmarkData', ...
    'No experimental JSON files found in %s.',datafolder);
files = dir(fullfile(datafolder,'**','*.json'));

S_benchmark.dry_run = true;
simulations = benchmark_predsim(S,osim_path,S_benchmark);
matches = struct2table(simulations);
n = height(matches);
matches.SlopePercent = 100*matches.Slope;
matches.MatchCount = zeros(n,1);
matches.MatchType = strings(n,1);
matches.ExperimentalID = strings(n,1);
matches.ExperimentalFile = strings(n,1);
matches.ExperimentalSpeed = nan(n,1);
matches.ExperimentalSlopePercent = nan(n,1);
matches.ExperimentalAddedMass = nan(n,1);
matches.ExperimentalMassFraction = nan(n,1);
matches.ExperimentalLoadLocation = strings(n,1);
matches.LocationSource = strings(n,1);
matches.Issues = strings(n,1);
matches.Unchecked = strings(n,1);
matches.Status = strings(n,1);
usage = zeros(numel(data),1);

for k = 1:n
    sim = simulations(k);
    duplicate_issues = strings(0,1);
    if sum(strcmp(sim.SimulationID,{simulations.SimulationID})) > 1
        duplicate_issues(end+1) = "duplicate_simulation_id";
    end
    if sum(strcmp(sim.SaveFolder,{simulations.SaveFolder})) > 1
        duplicate_issues(end+1) = "duplicate_output_folder";
    end
    matches.Issues(k) = strjoin(duplicate_issues,', ');
    idx = match_benchmark_id(sim.SimulationID,identifiers);
    matches.MatchCount(k) = numel(idx);
    matches.ExperimentalID(k) = strjoin(string(identifiers(idx)),' | ');
    matches.ExperimentalFile(k) = strjoin(string(arrayfun(@(j) ...
        fullfile(files(j).folder,files(j).name),idx,'UniformOutput',false)),' | ');
    if isempty(idx)
        matches.MatchType(k) = "none";
        if strcmp(sim.Study,'gait_speeds')
            matches.Status(k) = "no_experiment_expected";
        else
            matches.Status(k) = "missing";
        end
        continue
    elseif numel(idx) > 1
        matches.MatchType(k) = "ambiguous";
        matches.Status(k) = "ambiguous";
        continue % Production also refuses to attach ambiguous candidates.
    end
    usage(idx) = usage(idx) + 1;
    matches.MatchType(k) = "exact";
    if ~strcmp(sim.SimulationID,identifiers{idx})
        matches.MatchType(k) = "speed_tolerance";
    end
    exp = data{idx};
    matches.ExperimentalSpeed(k) = number(exp,'speed');
    matches.ExperimentalSlopePercent(k) = number(exp,'slope');
    matches.ExperimentalAddedMass(k) = number(exp,'added_mass');
    matches.ExperimentalMassFraction(k) = number(exp,'added_mass')/number(exp,'subject_mass');
    issues = duplicate_issues;
    unchecked = strings(0,1);
    if ~strcmpi(sim.Study,exp.study)
        issues(end+1) = "study";
    end
    speed_tolerance = 1e-4;
    if strcmp(sim.Study,'schertzer2014')
        speed_tolerance = 0.005; % Same tolerance as match_benchmark_id.
    end
    [issues,unchecked] = compare(issues,unchecked,"speed", ...
        sim.Speed,matches.ExperimentalSpeed(k),speed_tolerance);
    [issues,unchecked] = compare(issues,unchecked,"slope", ...
        matches.SlopePercent(k),matches.ExperimentalSlopePercent(k),1e-6);
    if strcmp(sim.Study,'gomenuka2014')
        [issues,unchecked] = compare(issues,unchecked,"mass_fraction", ...
            sim.MassFraction,matches.ExperimentalMassFraction(k),1e-6);
    else
        [issues,unchecked] = compare(issues,unchecked,"added_mass", ...
            sim.AddedMass,matches.ExperimentalAddedMass(k),1e-6);
        [issues,unchecked] = compare(issues,unchecked,"settings_added_mass", ...
            sim.AddedMass,sim.SettingsAddedMass,1e-6);
    end
    [location,source] = experimental_location(exp);
    matches.ExperimentalLoadLocation(k) = location;
    matches.LocationSource(k) = source;
    loaded = sim.AddedMass > 0 || sim.MassFraction > 0;
    if loaded
        if strlength(location) == 0
            unchecked(end+1) = "load_location";
        elseif ~strcmp(normalize_location(sim.LoadLocation),normalize_location(location))
            issues(end+1) = "load_location";
        end
    end
    matches.Issues(k) = strjoin(issues,', ');
    matches.Unchecked(k) = strjoin(unchecked,', ');
    matches.Status(k) = "matched";
    if ~isempty(issues)
        matches.Status(k) = "condition_mismatch";
    elseif ~isempty(unchecked)
        matches.Status(k) = "matched_unchecked";
    end
end

% Reverse audit: distinguish unique attachments from merely being a candidate.
candidate_count = zeros(numel(data),1);
for k = 1:n
    idx = match_benchmark_id(simulations(k).SimulationID,identifiers);
    candidate_count(idx) = candidate_count(idx) + 1;
end
selected = ismember(lower(string(studies)),lower(string(S_benchmark.studies)));
experiments = table(string(studies(:)),string(identifiers(:)), ...
    string(arrayfun(@(f) fullfile(f.folder,f.name),files,'UniformOutput',false)), ...
    selected(:),candidate_count,usage, ...
    'VariableNames',{'Study','ExperimentalID','ExperimentalFile', ...
    'SelectedStudy','CandidateCount','UniqueMatchCount'});
experiments.DuplicateID = false(height(experiments),1);
for k = 1:height(experiments)
    experiments.DuplicateID(k) = sum(strcmp(identifiers{k},identifiers)) > 1;
end
experiments.Status = repmat("unused",height(experiments),1);
experiments.Status(candidate_count > 0 & usage == 0) = "ambiguous_candidate";
experiments.Status(usage == 1) = "matched";
experiments.Status(usage > 1) = "reused";
experiments.Status(~selected) = "study_not_selected";

disp(matches(:,{'Study','SimulationID','ExperimentalID','MatchType','Status','Issues','Unchecked'}));
fprintf('ID audit: %d planned simulations; %d experimental records.\n',n,numel(data));
fprintf('Missing: %d; ambiguous: %d; condition mismatches: %d; unchecked: %d.\n', ...
    sum(matches.Status == "missing"),sum(matches.Status == "ambiguous"), ...
    sum(matches.Status == "condition_mismatch"),sum(strlength(matches.Unchecked) > 0));
fprintf('Selected experimental records without a unique assignment: %d.\n',sum(selected(:) & usage == 0));
reportfolder = char(p.Results.ReportFolder);
if ~isempty(reportfolder)
    if ~isfolder(reportfolder)
        mkdir(reportfolder);
    end
    writetable(matches,fullfile(reportfolder,'simulation_id_matches.csv'));
    writetable(experiments,fullfile(reportfolder,'experimental_id_coverage.csv'));
end
end

function value = number(data,field)
value = NaN;
if isfield(data,field) && isnumeric(data.(field)) && isscalar(data.(field))
    value = data.(field);
end
end

function [issues,unchecked] = compare(issues,unchecked,label,sim,exp,tolerance)
if ~isfinite(sim) || ~isfinite(exp)
    unchecked(end+1) = label;
elseif abs(sim-exp) > tolerance
    issues(end+1) = label;
end
end

function [location,source] = experimental_location(exp)
location = "";
source = "unavailable";
if isfield(exp,'location_added_mass') && ~isempty(exp.location_added_mass)
    location = string(exp.location_added_mass);
    source = "metadata";
    return
end
token = regexp(exp.identifier,'^browning2008_(femur|foot|pelvis|tibia)[0-9]+kg$', ...
    'tokens','once');
if isempty(token)
    token = regexp(exp.identifier,'^schertzer2014_[0-9]+p[0-9]+ms_(ankle|knee|torso)_[0-9]+kg$', ...
        'tokens','once');
end
if ~isempty(token)
    location = string(token{end});
    source = "experimental_id";
end
end

function location = normalize_location(location)
location = regexprep(lower(string(location)),'_[lr]$','');
if location == "calcn"
    location = "foot";
end
end
