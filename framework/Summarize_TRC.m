clc; clear;

Experiment = 'E2';
OutRoot = 'Results_TRC';
KeyMetrics = {'makespan', 'slaViolation', 'deadlineMiss', 'energy', 'critLoss', ...
              'unrecovered', 'maxTemp', 'activeNodes', 'planTime', 'aging', 'meanAF'};

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);
inDir = fullfile(thisDir, OutRoot, Experiment);
files = dir(fullfile(inDir, [Experiment '_*.mat']));
if isempty(files), error('No result files in %s', inDir); end

rows = {}; data = []; fronts = {};
for q = 1:numel(files)
    S = load(fullfile(inDir, files(q).name));
    rows = [rows; S.rows]; data = [data; S.data]; fronts = [fronts; S.fronts];
    Metrics = S.Metrics; cfg = S.cfg;
    fprintf('loaded %s (%d rows)\n', files(q).name, size(S.rows, 1));
end
cfgFull = experiment_config(Experiment);
Methods = cfgFull.Methods(ismember(cfgFull.Methods, unique(rows(:, 8))));
RefMethod = cfgFull.RefMethod;

key = cellfun(@(a, b, c, d, e, g) sprintf('%s|%d|%d|%s|%g|%d', a, b, c, d, e, g), ...
    rows(:, 2), rows(:, 3), rows(:, 4), rows(:, 5), rows(:, 6), rows(:, 7), 'UniformOutput', false);
fullKey = strcat(key, '|', rows(:, 8));
[~, ia] = unique(fullKey, 'stable');
if numel(ia) < numel(fullKey)
    fprintf('Removed %d duplicated rows\n', numel(fullKey) - numel(ia));
end
rows = rows(ia, :); data = data(ia, :); fronts = fronts(ia); key = key(ia);

if cfgFull.useHV
    hv = nan(size(rows, 1), 1);
    groups = unique(key);
    for g = 1:numel(groups)
        idx = find(strcmp(key, groups{g}));
        allF = cell2mat(fronts(idx));
        lo = min(allF, [], 1); hi = max(allF, [], 1); rg = max(hi - lo, 1e-12);
        for r = idx(:)'
            Fn = bsxfun(@rdivide, bsxfun(@minus, fronts{r}, lo), rg);
            hv(r) = hv_mc(Fn, 1.1, 50000, 777);
        end
    end
    data = [data, hv]; Metrics = [Metrics, {'HV'}];
    KeyMetrics = [{'HV'}, KeyMetrics];
end
if any(strcmp(Metrics, 'alarms')) && any(data(:, strcmp(Metrics, 'alarms')) > 0)
    KeyMetrics = [KeyMetrics, {'alarms', 'moves'}];
end
col = @(name) find(strcmp(Metrics, name));

fid = fopen(fullfile(inDir, [Experiment '_AllRuns.csv']), 'w');
fprintf(fid, 'Experiment,Workload,Tasks,Nodes,Factor,Level,Run,Method'); fprintf(fid, ',%s', Metrics{:}); fprintf(fid, '\n');
for r = 1:size(rows, 1)
    fprintf(fid, '%s,%s,%d,%d,%s,%g,%d,%s', rows{r, :}); fprintf(fid, ',%.6g', data(r, :)); fprintf(fid, '\n');
end
fclose(fid);

facCol = rows(:, 5); levCol = cell2mat(rows(:, 6)); methCol = rows(:, 8);
facs = unique(facCol, 'stable');
fid = fopen(fullfile(inDir, [Experiment '_Summary.csv']), 'w');
fprintf(fid, 'Factor,Level,Method'); fprintf(fid, ',%s', KeyMetrics{:}); fprintf(fid, '\n');
fprintf('\n================ %s summary (mean over workloads, scales, runs) ================\n', Experiment);
for fa = 1:numel(facs)
    levs = unique(levCol(strcmp(facCol, facs{fa})))';
    for lv = [levs, NaN]
        if isnan(lv), sel0 = strcmp(facCol, facs{fa}); lvTxt = 'all';
        else, sel0 = strcmp(facCol, facs{fa}) & levCol == lv; lvTxt = num2str(lv); end
        fprintf('\n%s = %s\n%-18s', facs{fa}, lvTxt, 'Method');
        fprintf(' %12s', KeyMetrics{:}); fprintf('\n');
        for k = 1:numel(Methods)
            s = sel0 & strcmp(methCol, Methods{k});
            fprintf(fid, '%s,%s,%s', facs{fa}, lvTxt, Methods{k});
            fprintf('%-18s', Methods{k});
            for mm = 1:numel(KeyMetrics)
                v = data(s, col(KeyMetrics{mm}));
                fprintf(fid, ',%.5g (%.3g)', mean(v), std(v));
                fprintf(' %12.4g', mean(v));
            end
            fprintf(fid, '\n'); fprintf('\n');
        end
    end
end
fclose(fid);

better = struct('HV', 'max');
fid = fopen(fullfile(inDir, [Experiment '_Stats.csv']), 'w');
fprintf(fid, 'Reference,Other,Factor,Level,Pairs,Metric,MeanRef,MeanOther,p_value,RefWinRate,Significant\n');
fprintf('\n================ paired Wilcoxon signed-rank: %s vs others ================\n', RefMethod);
refRows = find(strcmp(methCol, RefMethod));
for k = 1:numel(Methods)
    if strcmp(Methods{k}, RefMethod), continue; end
    pr = zeros(0, 2);
    for r = refRows(:)'
        o = find(strcmp(key, key{r}) & strcmp(methCol, Methods{k}));
        if ~isempty(o), pr(end+1, :) = [r, o(1)]; end
    end
    if isempty(pr), continue; end
    for fa = 1:numel(facs)
        inFac = strcmp(facCol(pr(:, 1)), facs{fa});
        levs = unique(levCol(pr(inFac, 1)))';
        for lv = [NaN, levs]
            if isnan(lv), s = inFac; lvTxt = 'all';
            else, s = inFac & levCol(pr(:, 1)) == lv; lvTxt = num2str(lv); end
            if ~any(s), continue; end
            a = data(pr(s, 1), :); b = data(pr(s, 2), :);
            fprintf('%s vs %s | %s = %s (%d pairs)\n', RefMethod, Methods{k}, facs{fa}, lvTxt, size(a, 1));
            for mm = 1:numel(KeyMetrics)
                c = col(KeyMetrics{mm});
                p = trc_signrank(a(:, c), b(:, c));
                if isfield(better, KeyMetrics{mm}), win = mean(a(:, c) > b(:, c)); else, win = mean(a(:, c) < b(:, c)); end
                sig = p < 0.05;
                fprintf(fid, '%s,%s,%s,%s,%d,%s,%.6g,%.6g,%.3g,%.3f,%d\n', RefMethod, Methods{k}, facs{fa}, lvTxt, ...
                    size(a, 1), KeyMetrics{mm}, mean(a(:, c)), mean(b(:, c)), p, win, sig);
                fprintf('   %-13s ref %.4g  other %.4g  p=%.3g  ref wins %.0f%%\n', KeyMetrics{mm}, ...
                    mean(a(:, c)), mean(b(:, c)), p, 100 * win);
            end
        end
    end
end
fclose(fid);

plotM = KeyMetrics(1:min(6, numel(KeyMetrics)));
for fa = 1:numel(facs)
    try
        levs = unique(levCol(strcmp(facCol, facs{fa})))';
        fig = figure('Visible', 'off', 'Position', [100 100 1300 800]);
        for pm = 1:numel(plotM)
            subplot(2, 3, pm); c = col(plotM{pm});
            Mtx = zeros(numel(levs), numel(Methods));
            for a2 = 1:numel(levs)
                for k = 1:numel(Methods)
                    Mtx(a2, k) = mean(data(strcmp(facCol, facs{fa}) & levCol == levs(a2) & strcmp(methCol, Methods{k}), c));
                end
            end
            if numel(levs) > 2
                plot(levs, Mtx, '-o', 'LineWidth', 1.5); xlabel(facs{fa});
            else
                bar(Mtx);
                set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%s=%g', facs{fa}, x), levs, 'UniformOutput', false));
            end
            grid on; title(plotM{pm});
            if pm == 1, legend(Methods, 'Location', 'best'); end
        end
        print(fig, fullfile(inDir, sprintf('%s_Figure_%s.png', Experiment, facs{fa})), '-dpng', '-r150');
        close(fig);
    catch ME
        fprintf('Figure for %s skipped: %s\n', facs{fa}, ME.message);
    end
end
fprintf('\nTables written to %s\n', inDir);
