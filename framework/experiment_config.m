function cfg = experiment_config(Experiment)
% Design of the experiments.
allW = {'Montage', 'Epigenomics', 'CyberShake', 'Inspiral', 'SIPHT'};
cfg.Experiment = Experiment;
cfg.Workloads = allW;
cfg.Scales = [30 6; 60 10];
cfg.Factors = struct('name', {'ambient'}, 'levels', {[0 8]});
cfg.numRuns = 10;
cfg.useHV = false;

switch Experiment
    case 'E2'
        cfg.Methods = {'IMOWO', 'MOWO', 'NSGA-II', 'MO-GWO', 'MO-PSO', 'SPEA2'};
        cfg.RefMethod = 'IMOWO';
        cfg.useHV = true;
    case 'E3'
        cfg.Methods = {'TRC', 'EFT', 'Thermal-only', 'Failure-only', 'PSB', 'TDRM', 'TRC-EqE'};
        cfg.RefMethod = 'TRC';
    case 'E3T'
        cfg.Methods = {'TRC', 'EFT', 'Thermal-only', 'Failure-only', 'PSB', 'TDRM', 'TRC-EqE'};
        cfg.RefMethod = 'TRC';
        cfg.truth = 'transient';
    case 'E4'
        cfg.Methods = {'TRC', 'TRC-noThrottle', 'TRC-noArrhenius', 'TRC-noReplica', ...
                       'TRC-noSleep', 'TRC-noForecast', 'Sequential'};
        cfg.RefMethod = 'TRC';
    case 'E5'
        cfg.Methods = {'IMOWO', 'IMOWO-noOps', 'IMOWO-noLHS', 'IMOWO-noOBL', 'IMOWO-noCDLeader', ...
                       'IMOWO-noKNN', 'IMOWO-noLevy', 'MOWO-core'};
        cfg.RefMethod = 'IMOWO';
        cfg.useHV = true;
    case 'E6'
        cfg.Methods = {'TRC', 'EFT', 'Thermal-only', 'Failure-only'};
        cfg.RefMethod = 'TRC';
        cfg.Workloads = {'Montage', 'CyberShake', 'Epigenomics'};
        cfg.Scales = [30 6];
        cfg.numRuns = 5;
        cfg.Factors = struct( ...
            'name',   {'p0', 'Ea', 'ambient', 'kappa', 'deadline', 'scale', 'repCost'}, ...
            'levels', {[0.005 0.01 0.02 0.05 0.10], [0.5 0.7 0.9 1.1], [0 4 8 12], ...
                       [0.02 0.03 0.05], [1.10 1.25 1.50 2.00], [1 2 3], [0.2 0.5 1.0]});
        cfg.ScaleLevels = [30 6; 60 10; 100 15];
    case 'E6p'
        cfg.Methods = {'TRC', 'EFT', 'Thermal-only', 'Failure-only'};
        cfg.RefMethod = 'TRC';
        cfg.Scales = [30 6];
        cfg.numRuns = 5;
        cfg.Factors = struct('name', {'p0'}, 'levels', {[0.0001 0.001]});
    case 'E3W'
        cfg.Methods = {'TRC', 'EFT', 'Thermal-only', 'Failure-only', 'PSB', 'TDRM', 'TRC-EqE'};
        cfg.RefMethod = 'TRC';
        cfg.truth = 'transient';
        cfg.Factors = struct('name', {'weibullAmb'}, 'levels', {[100 108 300 308]});
    case 'E8'
        cfg.Methods = {'TRC'};
        cfg.RefMethod = 'TRC';
        cfg.carsWeights = [0.25 0.25 0.25 0.25; 0.55 0.15 0.15 0.15; 0.15 0.55 0.15 0.15; ...
                           0.15 0.15 0.55 0.15; 0.15 0.15 0.15 0.55];
        cfg.Factors = struct('name', {'carsAmb'}, 'levels', {[10 18 20 28 30 38 40 48 50 58]});
    case 'E7'
        cfg.Methods = {'TRC', 'TRC+Local', 'TRC+PSC-U', 'TRC+PSC', 'EFT', 'EFT+PSC'};
        cfg.RefMethod = 'TRC+PSC';
    case 'E7b'
        cfg.Methods = {'TRC', 'TRC+PSC', 'TRC+PSC-R', 'TRC+PSC@T', 'TRC+PSC-R@T', ...
                       'TRC+Local@T', 'EFT', 'EFT+PSC@T'};
        cfg.RefMethod = 'TRC+PSC-R@T';
    otherwise
        error('Unknown experiment %s (use E2..E7, E7b)', Experiment);
end
end
