% Experiments of the TRC framework: set Experiment and run, then run Summarize_TRC.m.
clc; clear; close all;

Experiment     = 'E2';
RunTag         = 'A';
WorkloadsToRun = {};
MethodsToRun   = {};
numRunsOverride = [];
PopSize        = 100;
MaxIter        = 100;
QuickTest      = false;

P.Ea = 0.7;
P.Tref = 60;
P.p0 = 0.02;
P.throttleSlope = 0.03;
P.minSpeed = 0.40;
P.thermalIters = 3;
P.sleepFrac = 0.10;
P.standbyFrac = 0.05;
P.replicaCostFrac = 0.50;
P.ambientShift = 0;

P.epsMin = 0.05;
P.epsMax = 0.30;
P.epsFixed = 0.10;
P.pscGain = 0.05;
P.pscRiskTol = 1.0;

P.numMC = 200;
P.bgNoise = 0.05;
P.etaError = 0.10;
P.EaError = 0.15;
P.epsSLA = 0.10;
P.detectDelay = 2; P.replicaDelay = 1; P.wakeDelay = 10;
DeadlineFactor = 1.25;

P.carsW = [0.25 0.25 0.25 0.25];
P.carsTop = [0.50 0.20];

P.tauLow = 0.30; P.pThermal = 0.30; P.pFTRS = 0.30; P.pDRFA = 0.20;

HistoryLength = 20;
OutRoot = 'Results_TRC';

cfg = experiment_config(Experiment);
if ~isempty(WorkloadsToRun), cfg.Workloads = WorkloadsToRun; end
if ~isempty(MethodsToRun),   cfg.Methods   = MethodsToRun;   end
if ~isempty(numRunsOverride), cfg.numRuns = numRunsOverride; end
if QuickTest
    cfg.Workloads = cfg.Workloads(1); cfg.numRuns = 1;
    PopSize = 10; MaxIter = 4; P.numMC = 20;
    if strcmp(Experiment, 'E6')
        for f = 1:numel(cfg.Factors), cfg.Factors(f).levels = cfg.Factors(f).levels(1:2); end
    else
        cfg.Scales = cfg.Scales(1, :); cfg.Factors(1).levels = cfg.Factors(1).levels(end);
    end
end

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);
addpath(fullfile(thisDir, '..', 'optimizers'));
outDir = fullfile(thisDir, OutRoot, Experiment);
if ~exist(outDir, 'dir'), mkdir(outDir); end

allW = {'Montage', 'Epigenomics', 'CyberShake', 'Inspiral', 'SIPHT'};
Metrics = {'makespan','makespanP95','deadlineMiss','slaViolation','slowdown','energy','critLoss', ...
           'migrated','replicaUsed','unrecovered','nodeFailures','maxTemp','throttledNodes','activeNodes', ...
           'planTime','promised','nReplicas','forecastErr','nCritical','alarms','moves','aging','meanAF'};
rows = {}; data = zeros(0, numel(Metrics)); fronts = {};

fprintf('==== %s | methods: %s | runs %d | N=%d T=%d ====\n', Experiment, ...
    strjoin(cfg.Methods, ', '), cfg.numRuns, PopSize, MaxIter);

for w = 1:numel(cfg.Workloads)
    wl = cfg.Workloads{w};
    wIdx = find(strcmp(allW, wl));
    for sc = 1:size(cfg.Scales, 1)
        for f = 1:numel(cfg.Factors)
            fac = cfg.Factors(f);
            if strcmp(fac.name, 'deadline')
                levelSets = {fac.levels};
            else
                levelSets = num2cell(fac.levels);
            end
            for ls = 1:numel(levelSets)
                lev = levelSets{ls};
                Pl = P; dlFactors = DeadlineFactor; scale = cfg.Scales(sc, :);
                switch fac.name
                    case 'ambient',  Pl.ambientShift = lev;
                    case 'p0',       Pl.p0 = lev;
                    case 'Ea',       Pl.Ea = lev;
                    case 'kappa',    Pl.throttleSlope = lev;
                    case 'deadline', dlFactors = lev;
                    case 'scale',    scale = cfg.ScaleLevels(lev, :);
                    case 'repCost',  Pl.replicaCostFrac = lev;
                    case 'weibullAmb'
                        Pl.weibullK2 = floor(lev / 100) + mod(floor(lev / 10), 10) / 10;
                        Pl.ambientShift = mod(lev, 10);
                    case 'carsAmb'
                        Pl.ambientShift = mod(lev, 10);
                end
                nT = scale(1); nN = scale(2);

                wf     = generate_workflow_DAG(wl, nT, 42 + nT);
                fogEnv = generate_fog_nodes(nN, 2024 + nN);
                histAll = simulate_resource_history(nN, HistoryLength + 1, 99 + nN + wIdx);
                hist   = histAll(:, 1:HistoryLength);
                bgTrue = histAll(:, end)';
                rs = rng(); rng(500000 + 1000 * wIdx + 10 * nT + nN, 'twister');
                bgARS  = ARSENAL_forecast(hist, 5);
                rng(rs);
                bgLast = hist(:, end)';
                Pl.Href = sum(wf.taskFLOPs) / sum(fogEnv.nodeCapacity .* max(0.15, 1 - bgARS));
                [crit, psi] = CARS_forecast(wf, fogEnv, bgARS, Pl);
                critPlan = crit; psiPlan = psi;
                if strcmp(fac.name, 'carsAmb')
                    Pw = Pl; Pw.carsW = cfg.carsWeights(floor(lev / 10), :);
                    [critPlan, psiPlan] = CARS_forecast(wf, fogEnv, bgARS, Pw);
                end
                ctx = struct('wf', wf, 'fogEnv', fogEnv, 'nT', nT, 'psi', psiPlan, 'crit', critPlan, ...
                             'bgARSENAL', bgARS, 'bgLast', bgLast, 'P', Pl);
                aEFT = eft_schedule(wf, fogEnv, bgARS);
                [~, DE] = trc_evaluate_schedule(aEFT, zeros(1, nT), wf, fogEnv, bgARS, psi, Pl, 'coupled');
                deadlines = dlFactors * DE.makespan;
                ctx.deadline = DeadlineFactor * DE.makespan;
                fprintf('\n--- %s | %d tasks %d nodes | %s = %s | critical tasks %d | ARSENAL MAE %.3f (last value %.3f)\n', ...
                    wl, nT, nN, fac.name, mat2str(lev), sum(crit == 3), mean(abs(bgARS - bgTrue)), mean(abs(bgLast - bgTrue)));

                for run = 1:cfg.numRuns
                    seed = 100000 * wIdx + 1000 * nT + 10 * run + round(100 * ls) + 7 * f;
                    simSeed = 7000 + seed;
                    planKeys = {}; planOuts = {};
                    for k = 1:numel(cfg.Methods)
                        meth = cfg.Methods{k};
                        mdl = trc_models(meth);
                        q = find(strcmp(planKeys, mdl.planner), 1);
                        if isempty(q)
                            rng(seed, 'twister');
                            out = run_method(mdl.planner, ctx, PopSize, MaxIter);
                            planKeys{end+1} = mdl.planner; planOuts{end+1} = out;
                        else
                            out = planOuts{q};
                        end
                        rt = trc_runtime_plan(out, ctx, mdl);
                        for d = 1:numel(deadlines)
                            if isfield(cfg, 'truth') && strcmp(cfg.truth, 'transient')
                                R = simulate_TRC_truth2(out.assign, out.rep, wf, fogEnv, bgTrue, psi, Pl, ...
                                                        deadlines(d), simSeed, out.promised);
                            else
                                R = simulate_TRC_truth(out.assign, out.rep, wf, fogEnv, bgTrue, psi, Pl, ...
                                                       deadlines(d), simSeed, out.promised, rt);
                            end
                            R.planTime = out.time; R.promised = out.promised;
                            R.nReplicas = sum(out.rep > 0);
                            if strcmp(trc_models(meth).forecast, 'last'), bgUsed = bgLast; else, bgUsed = bgARS; end
                            R.forecastErr = mean(abs(bgUsed - bgTrue));
                            R.nCritical = sum(crit == 3);
                            v = zeros(1, numel(Metrics));
                            for mm = 1:numel(Metrics), v(mm) = R.(Metrics{mm}); end
                            if strcmp(fac.name, 'deadline'), levOut = dlFactors(d); else, levOut = lev; end
                            rows(end+1, :) = {Experiment, wl, nT, nN, fac.name, levOut, run, meth};
                            data(end+1, :) = v;
                            fronts{end+1, 1} = out.F;
                        end
                        fprintf('   run %d  %-17s makespan %.4f  SLA-v %.2f  energy %.4f  critLoss %.4f  (%.1f s)\n', ...
                            run, meth, R.makespan, R.slaViolation, R.energy, R.critLoss, out.time);
                    end
                end
            end
        end
    end
end

file = fullfile(outDir, sprintf('%s_%s.mat', Experiment, RunTag));
save(file, 'rows', 'data', 'fronts', 'Metrics', 'P', 'cfg', 'PopSize', 'MaxIter');
fprintf('\nSaved %s\nRun Summarize_TRC.m (Experiment = ''%s'').\n', file, Experiment);
