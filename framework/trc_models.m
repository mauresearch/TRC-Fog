function model = trc_models(name)
% Planning model of every method.

k = strfind(name, '+');
if ~isempty(k)
    model = trc_models(name(1:k(1)-1));
    model.name = name;
    model.runtime = name(k(1)+1:end);
    return;
end

model = struct('name', name, 'objective', 'TRC', 'mode', 'coupled', 'throttle', true, ...
    'useRep', true, 'sleepPlan', true, 'forecast', 'ARSENAL', ...
    'opThermal', true, 'opFTRS', true, 'opDRFA', true, 'drfaRule', 'analytic', ...
    'optimizer', 'IMOWO', 'imowoOpts', struct(), 'special', '');
model.planner = name;
model.runtime = 'none';

switch name
    case {'IMOWO', 'TRC', 'TRC (proposed)'}
    case 'MOWO'
        model.optimizer = 'MOWO';
        model.opThermal = false; model.opFTRS = false; model.opDRFA = false;
    case {'NSGA-II', 'MO-GWO', 'MO-PSO', 'SPEA2'}
        model.optimizer = strrep(name, '-', '_');
        if strcmp(name, 'NSGA-II'), model.optimizer = 'NSGA2'; end
        model.opThermal = false; model.opFTRS = false; model.opDRFA = false;

    case 'EFT'
        model.special = 'EFT';
        model.objective = 'TRC'; model.mode = 'uncoupled'; model.throttle = false; model.useRep = false;
    case 'Thermal-only'
        model.objective = 'Thermal'; model.mode = 'uncoupled'; model.throttle = false;
        model.useRep = false; model.opFTRS = false; model.drfaRule = 'thermal';
    case 'TRC-EqE'
        model.special = 'EqE';
    case 'PSB'
        model.special = 'PSB';
        model.mode = 'uncoupled'; model.throttle = false; model.useRep = true;
    case 'TDRM'
        model.special = 'TDRM';
        model.mode = 'uncoupled'; model.throttle = false; model.useRep = false;
    case 'Failure-only'
        model.mode = 'uncoupled'; model.throttle = false;
        model.opThermal = false; model.drfaRule = 'failure';

    case 'TRC-noThrottle'
        model.throttle = false;
    case 'TRC-noArrhenius'
        model.mode = 'uncoupled'; model.drfaRule = 'failure';
    case 'TRC-noReplica'
        model.useRep = false; model.opFTRS = false;
    case 'TRC-noSleep'
        model.sleepPlan = false;
    case 'TRC-noForecast'
        model.forecast = 'last';
    case 'Sequential'
        model.special = 'Sequential';
        model.objective = 'Thermal'; model.mode = 'uncoupled'; model.throttle = false;
        model.useRep = false; model.opFTRS = false; model.opDRFA = false;

    case 'IMOWO-noOps'
        model.opThermal = false; model.opFTRS = false; model.opDRFA = false;
    case 'IMOWO-noLHS'
        model.imowoOpts.useLHS = false;
    case 'IMOWO-noOBL'
        model.imowoOpts.useOBL = false;
    case 'IMOWO-noKNN'
        model.imowoOpts.archive = 'cd';
    case 'IMOWO-noCDLeader'
        model.imowoOpts.leader = 'first';
    case 'IMOWO-noLevy'
        model.imowoOpts.useLevy = false;
    case 'MOWO-core'
        model.opThermal = false; model.opFTRS = false; model.opDRFA = false;
        model.imowoOpts.useLHS = false; model.imowoOpts.useOBL = false;
        model.imowoOpts.leader = 'first'; model.imowoOpts.archive = 'cd';

    otherwise
        error('trc_models: unknown method "%s"', name);
end
end
