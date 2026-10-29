% Hypervolume of the E2 fronts for several reference points (writes E2_HV_RefPoints.csv).
clc; clear;

Experiment = 'E2';
OutRoot = 'Results_TRC';
RefPoints = [1.05 1.1 1.2 1.5];

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);
inDir = fullfile(thisDir, OutRoot, Experiment);
files = dir(fullfile(inDir, [Experiment '_*.mat']));
if isempty(files), error('No result files in %s', inDir); end

rows = {}; fronts = {};
for q = 1:numel(files)
    S = load(fullfile(inDir, files(q).name));
    rows = [rows; S.rows]; fronts = [fronts; S.fronts];
end
key = cellfun(@(a, b, c, d, e, g) sprintf('%s|%d|%d|%s|%g|%d', a, b, c, d, e, g), ...
    rows(:, 2), rows(:, 3), rows(:, 4), rows(:, 5), rows(:, 6), rows(:, 7), 'UniformOutput', false);
fullKey = strcat(key, '|', rows(:, 8));
[~, ia] = unique(fullKey, 'stable');
rows = rows(ia, :); fronts = fronts(ia); key = key(ia);

hv = nan(size(rows, 1), numel(RefPoints));
groups = unique(key);
for g = 1:numel(groups)
    idx = find(strcmp(key, groups{g}));
    allF = cell2mat(fronts(idx));
    lo = min(allF, [], 1); hi = max(allF, [], 1); rg = max(hi - lo, 1e-12);
    for r = idx(:)'
        Fn = bsxfun(@rdivide, bsxfun(@minus, fronts{r}, lo), rg);
        for k = 1:numel(RefPoints)
            hv(r, k) = hv_mc(Fn, RefPoints(k), 50000, 777);
        end
    end
end

fid = fopen(fullfile(inDir, 'E2_HV_RefPoints.csv'), 'w');
fprintf(fid, 'Workload,Tasks,Nodes,Factor,Level,Run,Method');
fprintf(fid, ',HV_%g', RefPoints); fprintf(fid, '\n');
for r = 1:size(rows, 1)
    fprintf(fid, '%s,%d,%d,%s,%g,%d,%s', rows{r, 2:8});
    fprintf(fid, ',%.6g', hv(r, :)); fprintf(fid, '\n');
end
fclose(fid);
fprintf('Written %s\n', fullfile(inDir, 'E2_HV_RefPoints.csv'));
